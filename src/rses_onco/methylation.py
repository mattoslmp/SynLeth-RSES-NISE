"""Utilities for TCGA/GDC gene-associated DNA methylation evidence."""
from __future__ import annotations

import math
import re
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .utils import canonical_gene_name

CANCER_PROJECTS = {
  "colon": ("TCGA-COAD", "TCGA-READ"),
  "stomach": ("TCGA-STAD",),
  "lung": ("TCGA-LUAD", "TCGA-LUSC"),
}
METHYLATION_SUBWEIGHTS = {
  "tumor_profile_divergence": 0.70,
  "tumor_normal_delta_divergence": 0.30,
}


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if math.isfinite(result) else None


def tcga_sample_type(sample: str) -> str:
  parts = str(sample).split("-")
  if len(parts) < 4:
    return ""
  match = re.match(r"^(\d{2})", parts[3])
  return match.group(1) if match else ""


def tcga_case_id(sample: str) -> str:
  parts = str(sample).split("-")
  return "-".join(parts[:3]) if len(parts) >= 3 else str(sample)


def position_rows(position: object) -> list[dict[str, Any]]:
  if isinstance(position, list):
    return [dict(item) for item in position if isinstance(item, dict)]
  if not isinstance(position, dict):
    return []
  lengths = [
    len(value) for value in position.values() if isinstance(value, list)
  ]
  if not lengths:
    return []
  rows = []
  for index in range(max(lengths)):
    row = {
      key: (
        value[index]
        if isinstance(value, list) and index < len(value)
        else value
      )
      for key, value in position.items()
    }
    if isinstance(row.get("position"), dict):
      row.update(row["position"])
    rows.append(row)
  return rows


def probe_identifier(row: dict[str, Any]) -> str:
  return str(
    row.get("name")
    or row.get("probe")
    or row.get("id")
    or row.get("probe_id")
    or ""
  )


def query_gene_frame(
  position: object,
  values: object,
  samples: list[str],
  gene: str,
) -> pd.DataFrame:
  rows = position_rows(position)
  matrices = values if isinstance(values, list) else []
  records = []
  for index, row in enumerate(rows):
    if index >= len(matrices):
      continue
    probe = probe_identifier(row)
    if not probe:
      continue
    for sample, value in zip(samples, matrices[index], strict=False):
      beta = numeric(value)
      if beta is None:
        continue
      records.append({
        "gene": canonical_gene_name(gene),
        "probe_id": probe,
        "sample_id": str(sample),
        "case_id": tcga_case_id(str(sample)),
        "sample_type_code": tcga_sample_type(str(sample)),
        "beta_value": float(np.clip(beta, 0.0, 1.0)),
      })
  return pd.DataFrame(records)


def choose_dataset(
  xena: Any,
  host: str,
  project: str,
  explicit: str | None = None,
) -> tuple[str | None, str, int]:
  if explicit:
    return explicit, "explicit_dataset_override", 1
  cohorts = xena.all_cohorts(host, [])
  matching = [
    str(item.get("cohort") if isinstance(item, dict) else item)
    for item in cohorts
    if project.casefold()
    in str(
      item.get("cohort") if isinstance(item, dict) else item
    ).casefold()
  ]
  candidates = []
  for cohort in matching:
    for record in xena.dataset_list(host, [cohort]):
      text = " ".join(
        str(record.get(key) or "")
        for key in (
          "name", "longtitle", "type", "datasubtype", "text"
        )
      ).casefold()
      if "methyl" in text:
        candidates.append((record, text))
  if not candidates:
    return None, "no_methylation_dataset_discovered", 0

  def rank(item: tuple[dict[str, Any], str]) -> tuple[int, int]:
    _, text = item
    platform = (
      4 if "450" in text
      else 3 if "epic v2" in text
      else 2 if "epic" in text
      else 1 if "27" in text
      else 0
    )
    return platform, int("sesame" in text)

  selected = max(candidates, key=rank)[0]
  return (
    str(selected.get("name") or ""),
    "discovered_gdc_xena_dataset",
    len(candidates),
  )


def coverage_subscore(
  components: dict[str, float | None],
) -> dict[str, float | int]:
  numerator = 0.0
  observed_weight = 0.0
  observed = 0
  for name, weight in METHYLATION_SUBWEIGHTS.items():
    value = numeric(components.get(name))
    if value is None:
      continue
    numerator += weight * float(np.clip(value, 0, 1))
    observed_weight += weight
    observed += 1
  raw = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / sum(METHYLATION_SUBWEIGHTS.values())
  adjusted = raw * coverage if math.isfinite(raw) else float("nan")
  return {
    "raw": raw,
    "coverage": coverage,
    "adjusted": adjusted,
    "observed": observed,
  }


def summarize_genes(long: pd.DataFrame) -> pd.DataFrame:
  if long.empty:
    return pd.DataFrame()
  rows = []
  for keys, group in long.groupby(["cancer", "project", "gene"]):
    cancer, project, gene = keys
    tumor = group.loc[
      group["sample_type_code"].eq("01"), "beta_value"
    ]
    normal = group.loc[
      group["sample_type_code"].eq("11"), "beta_value"
    ]
    rows.append({
      "cancer": cancer,
      "project": project,
      "gene": gene,
      "dataset": group["dataset"].iloc[0],
      "probe_count": int(group["probe_id"].nunique()),
      "primary_tumor_sample_count": int(
        group.loc[
          group["sample_type_code"].eq("01"), "case_id"
        ].nunique()
      ),
      "solid_tissue_normal_sample_count": int(
        group.loc[
          group["sample_type_code"].eq("11"), "case_id"
        ].nunique()
      ),
      "primary_tumor_median_beta": (
        float(tumor.median()) if len(tumor) else np.nan
      ),
      "solid_tissue_normal_median_beta": (
        float(normal.median()) if len(normal) else np.nan
      ),
      "tumor_normal_delta_beta": (
        float(tumor.median() - normal.median())
        if len(tumor) and len(normal)
        else np.nan
      ),
    })
  return pd.DataFrame(rows)


def build_pair_metrics(
  long: pd.DataFrame,
  summaries: pd.DataFrame,
  pairs: pd.DataFrame,
  cancer_status: dict[str, dict[str, Any]],
  min_samples: int,
) -> pd.DataFrame:
  summary_index = {
    (str(row["cancer"]), str(row["gene"])): row
    for row in summaries.to_dict("records")
  }
  rows = []
  for cancer in CANCER_PROJECTS:
    status = cancer_status.get(cancer, {})
    source_eligible = status.get("status") == "available"
    subset = (
      long.loc[
        long["cancer"].eq(cancer)
        & long["sample_type_code"].eq("01")
      ]
      if not long.empty
      else pd.DataFrame()
    )
    gene_case = (
      subset.groupby(["gene", "case_id"], as_index=False)[
        "beta_value"
      ].median()
      if not subset.empty
      else pd.DataFrame(columns=["gene", "case_id", "beta_value"])
    )
    observed_genes = set(gene_case["gene"].astype(str))
    for pair in pairs.to_dict("records"):
      lost = canonical_gene_name(pair.get("lost_gene"))
      target = canonical_gene_name(pair.get("target_gene"))
      lost_values = gene_case.loc[
        gene_case["gene"].eq(lost), ["case_id", "beta_value"]
      ].rename(columns={"beta_value": "lost_beta"})
      target_values = gene_case.loc[
        gene_case["gene"].eq(target), ["case_id", "beta_value"]
      ].rename(columns={"beta_value": "target_beta"})
      overlap = lost_values.merge(target_values, on="case_id")
      rho = median_difference = profile = None
      if len(overlap) >= min_samples:
        rho = numeric(
          spearmanr(
            overlap["lost_beta"],
            overlap["target_beta"],
            nan_policy="omit",
          ).statistic
        )
        median_difference = float(
          np.median(
            np.abs(overlap["lost_beta"] - overlap["target_beta"])
          )
        )
        pieces = [
          float(np.clip(median_difference / 0.5, 0, 1))
        ]
        if rho is not None:
          pieces.append((1 - rho) / 2)
        profile = float(np.mean(pieces))
      lost_delta = numeric(
        summary_index.get((cancer, lost), {}).get(
          "tumor_normal_delta_beta"
        )
      )
      target_delta = numeric(
        summary_index.get((cancer, target), {}).get(
          "tumor_normal_delta_beta"
        )
      )
      delta_divergence = (
        float(
          np.clip(abs(lost_delta - target_delta) / 0.25, 0, 1)
        )
        if lost_delta is not None and target_delta is not None
        else None
      )
      score = coverage_subscore({
        "tumor_profile_divergence": profile,
        "tumor_normal_delta_divergence": delta_divergence,
      })
      reason = ""
      if not source_eligible:
        reason = str(
          status.get("reason") or "methylation_source_unavailable"
        )
      elif lost not in observed_genes and target not in observed_genes:
        reason = (
          "both_genes_unmapped_or_without_gene_associated_probes"
        )
      elif len(overlap) < min_samples:
        reason = "insufficient_overlapping_primary_tumor_samples"
      rows.append({
        "pair_id": pair.get("pair_id"),
        "cancer": cancer,
        "lost_gene": lost,
        "target_gene": target,
        "methylation_source_eligible": source_eligible,
        "methylation_source_status": status.get(
          "status", "unavailable"
        ),
        "methylation_absence_reason": reason,
        "methylation_primary_tumor_overlap_n": len(overlap),
        "methylation_pair_spearman_rho": rho,
        "methylation_pair_median_absolute_beta_difference": (
          median_difference
        ),
        "methylation_tumor_profile_divergence": profile,
        "methylation_lost_tumor_normal_delta_beta": lost_delta,
        "methylation_target_tumor_normal_delta_beta": target_delta,
        "methylation_tumor_normal_delta_divergence": (
          delta_divergence
        ),
        "methylation_context_raw": score["raw"],
        "methylation_context_coverage": score["coverage"],
        "component_promoter_methylation_context": score["adjusted"],
        "methylation_observed_subcomponents": score["observed"],
        "methylation_source": "NCI_GDC_via_UCSC_Xena_GDC_Hub",
        "methylation_interpretation": (
          "gene-associated CpG methylation context; not direct proof "
          "of promoter silencing, compensation or causality"
        ),
      })
  return pd.DataFrame(rows)
