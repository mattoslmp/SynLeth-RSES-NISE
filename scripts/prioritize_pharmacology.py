#!/usr/bin/env python3
"""Prioritize target-drug hypotheses for experimental cancer validation.

The script merges RSES-Onco vulnerabilities with Open Targets, ChEMBL, DGIdb,
MyChem, Pharos, CIViC and optional PRISM/GDSC/CTRP evidence. Outputs are research
priorities, not treatment recommendations or claims of clinical efficacy.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

from rses_onco.pharmacology import (
  PHARMACOLOGY_WEIGHTS,
  coverage_aware_pharmacology_score,
  evidence_rows_to_components,
  normalize_phase,
  therapeutic_hypothesis_score,
)
from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
CANCER_TERMS = {
  "colon": {"colon", "colorectal", "rectal", "bowel"},
  "stomach": {"stomach", "gastric", "gastroesophageal"},
  "lung": {"lung", "pulmonary", "non-small cell", "small cell"},
}
SUPPRESSIVE_TERMS = {
  "inhibitor", "inhibition", "antagonist", "blocker", "degrader",
  "negative modulator", "suppressor", "antisense", "sirna", "rnai",
}
ACTIVATING_TERMS = {"agonist", "activator", "positive modulator", "stimulator"}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists():
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t")
  except EmptyDataError:
    return pd.DataFrame()


def drug_key(value: object) -> str:
  if value is None or pd.isna(value):
    return ""
  return re.sub(r"[^A-Z0-9]+", "", str(value).upper())


def suppression_compatibility(rows: pd.DataFrame) -> float | None:
  texts: list[str] = []
  for column in (
    "interaction_type", "action_type", "mechanism_of_action", "directionality"
  ):
    if column in rows:
      texts.extend(rows[column].dropna().astype(str).str.casefold().tolist())
  if not texts:
    return None
  joined = " ; ".join(texts)
  if any(term in joined for term in SUPPRESSIVE_TERMS):
    return 1.0
  if any(term in joined for term in ACTIVATING_TERMS):
    return 0.0
  return 0.50


def cancer_relevance(rows: pd.DataFrame, cancer: str) -> float | None:
  terms = CANCER_TERMS.get(cancer, set())
  diseases = rows.get("disease_name", pd.Series(index=rows.index, dtype=object))
  disease_text = " ; ".join(diseases.dropna().astype(str).str.casefold())
  if not disease_text:
    source = rows.get("source", pd.Series(index=rows.index, dtype=object)).astype(str)
    civic = rows.loc[source.eq("civic")]
    civic_values = civic.get(
      "civic_gene_record", pd.Series(index=civic.index, dtype=object)
    )
    if (
      not civic.empty
      and civic_values.astype(str).str.casefold().isin({"1", "true", "yes"}).any()
    ):
      return 0.45
    return None
  if any(term in disease_text for term in terms):
    return 1.0
  if any(term in disease_text for term in ("cancer", "carcinoma", "neoplasm", "tumor")):
    return 0.60
  return 0.20


def collect_drug_groups(evidence: pd.DataFrame) -> dict[str, list[dict[str, Any]]]:
  groups: dict[str, list[dict[str, Any]]] = {}
  if evidence.empty:
    return groups
  for record in evidence.to_dict("records"):
    target = canonical_gene_name(record.get("target_gene"))
    if not target:
      continue
    key = drug_key(record.get("drug_name")) or drug_key(record.get("drug_id"))
    if not key:
      key = "TARGET_ONLY"
    one = dict(record)
    one["drug_key"] = key
    groups.setdefault(target, []).append(one)
  return groups


def merge_mychem(evidence: pd.DataFrame) -> pd.DataFrame:
  if evidence.empty or "source" not in evidence:
    return evidence
  mychem = evidence.loc[evidence["source"].astype(str).eq("mychem")].copy()
  primary = evidence.loc[~evidence["source"].astype(str).eq("mychem")].copy()
  if mychem.empty or primary.empty:
    return evidence
  mychem["drug_key"] = mychem.get(
    "drug_id", pd.Series(index=mychem.index, dtype=object)
  ).map(drug_key)
  annotations = {
    key: record.get("raw_record")
    for key, record in (
      mychem.loc[mychem["drug_key"].ne("")]
        .drop_duplicates("drug_key")
        .set_index("drug_key")
        .to_dict("index")
        .items()
    )
  }
  primary["drug_key"] = primary.apply(
    lambda row: drug_key(row.get("drug_name")) or drug_key(row.get("drug_id")),
    axis=1,
  )
  primary["mychem_annotation"] = primary["drug_key"].map(annotations)
  return pd.concat([primary, mychem], ignore_index=True, sort=False)


def candidate_drug_rows(
  ranking: pd.DataFrame,
  evidence: pd.DataFrame,
  sensitivity: pd.DataFrame,
) -> pd.DataFrame:
  evidence = merge_mychem(evidence)
  drug_groups = collect_drug_groups(evidence)
  sensitivity_groups: dict[tuple[str, str, str], pd.DataFrame] = {}
  if not sensitivity.empty:
    sensitivity = sensitivity.copy()
    sensitivity["drug_key"] = sensitivity.apply(
      lambda row: drug_key(row.get("drug_name")) or drug_key(row.get("drug_id")),
      axis=1,
    )
    for group_key, group in sensitivity.groupby(
      ["pair_id", "cancer", "drug_key"], dropna=False
    ):
      sensitivity_groups[
        (str(group_key[0]), str(group_key[1]), str(group_key[2]))
      ] = group

  output_rows: list[dict[str, Any]] = []
  for candidate in ranking.to_dict("records"):
    target = canonical_gene_name(
      candidate.get("analysis_target_gene") or candidate.get("target_gene")
    )
    cancer = str(candidate.get("cancer") or "")
    pair_id = str(candidate.get("pair_id") or "")
    by_drug: dict[str, list[dict[str, Any]]] = {}
    for record in drug_groups.get(target, []):
      by_drug.setdefault(str(record.get("drug_key") or "TARGET_ONLY"), []).append(record)
    if not by_drug:
      by_drug = {"TARGET_ONLY": []}

    for key, records in by_drug.items():
      rows = pd.DataFrame(records)
      sensitivity_rows = sensitivity_groups.get((pair_id, cancer, key), pd.DataFrame())
      if not sensitivity_rows.empty:
        one = sensitivity_rows.copy()
        one["target_gene"] = target
        rows = pd.concat([rows, one], ignore_index=True, sort=False)
      components = evidence_rows_to_components(rows)
      relevance = cancer_relevance(rows, cancer)
      if relevance is not None:
        components["cancer_relevance"] = relevance
      compatibility = suppression_compatibility(rows)
      if compatibility is not None:
        direct = components.get("direct_interaction")
        components["direct_interaction"] = (
          compatibility
          if direct is None
          else float(np.mean([direct, compatibility]))
        )
      pharmacology = coverage_aware_pharmacology_score(components)
      vulnerability = candidate.get("coverage_adjusted_rses")
      combined = therapeutic_hypothesis_score(
        vulnerability, pharmacology.adjusted_score
      )
      drug_names = rows.get(
        "drug_name", pd.Series(index=rows.index, dtype=object)
      ).dropna().astype(str)
      drug_ids = rows.get(
        "drug_id", pd.Series(index=rows.index, dtype=object)
      ).dropna().astype(str)
      drug_name = drug_names.mode().iloc[0] if not drug_names.empty else None
      drug_id = drug_ids.mode().iloc[0] if not drug_ids.empty else None
      sources = ";".join(sorted(set(
        rows.get("source", pd.Series(index=rows.index, dtype=object))
          .dropna().astype(str)
      )))
      phase_values = [
        normalize_phase(value)
        for value in rows.get("max_phase", pd.Series(index=rows.index, dtype=object))
      ]
      phase_values = [value for value in phase_values if value is not None]
      output = dict(candidate)
      output.update({
        "drug_key": key,
        "drug_name": drug_name,
        "drug_id": drug_id,
        "pharmacology_sources": sources,
        "pharmacology_source_count": len(set(sources.split(";"))) if sources else 0,
        **{
          f"pharmacology_component_{domain}": value
          for domain, value in components.items()
        },
        "pharmacology_rses": pharmacology.observed_score,
        "pharmacology_coverage": pharmacology.coverage,
        "pharmacology_adjusted": pharmacology.adjusted_score,
        "pharmacology_n_domains": pharmacology.n_domains,
        "pharmacology_class": pharmacology.interpretation,
        "therapeutic_hypothesis_score": combined,
        "suppression_compatibility": compatibility,
        "maximum_clinical_phase_normalized": (
          max(phase_values) if phase_values else None
        ),
        "repurposing_candidate": bool(
          drug_name
          and max(phase_values, default=0.0) >= 0.5
          and (compatibility is None or compatibility >= 0.5)
        ),
        "research_only": True,
        "interpretation_boundary": (
          "Computational experimental-priority hypothesis; not evidence of clinical "
          "efficacy, patient benefit or cure. Requires biomarker-matched validation."
        ),
      })
      output_rows.append(output)
  return pd.DataFrame(output_rows)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--evidence",
    default="data/processed/pharmacology/pharmacology_evidence_long.tsv",
  )
  parser.add_argument(
    "--sensitivity",
    default="data/processed/pharmacology/drug_response_selectivity.tsv",
  )
  parser.add_argument(
    "--output-dir",
    default="results/expanded_26Q1/pharmacology",
  )
  args = parser.parse_args()

  ranking = pd.read_csv(resolve_path(args.ranking), sep="\t")
  evidence = read_optional(resolve_path(args.evidence))
  sensitivity = read_optional(resolve_path(args.sensitivity))
  output_dir = resolve_path(args.output_dir)
  output_dir.mkdir(parents=True, exist_ok=True)

  result = candidate_drug_rows(ranking, evidence, sensitivity)
  if result.empty:
    raise SystemExit("No pharmacology hypotheses could be constructed")
  result = result.sort_values(
    [
      "therapeutic_hypothesis_score",
      "coverage_adjusted_rses",
      "pharmacology_adjusted",
    ],
    ascending=[False, False, False],
    na_position="last",
  )
  result.to_csv(
    output_dir / "pharmacology_ranked_hypotheses.tsv", sep="\t", index=False
  )

  actionability_columns = [
    column for column in [
      "analysis_target_gene", "target_gene", "drug_key", "drug_name", "drug_id",
      "pharmacology_sources", "pharmacology_source_count",
      *[f"pharmacology_component_{domain}" for domain in PHARMACOLOGY_WEIGHTS],
      "pharmacology_rses", "pharmacology_coverage", "pharmacology_adjusted",
      "pharmacology_class", "repurposing_candidate",
    ] if column in result.columns
  ]
  dedup_columns = [
    column for column in ("analysis_target_gene", "target_gene", "drug_key")
    if column in actionability_columns
  ]
  target_actionability = (
    result[actionability_columns]
      .sort_values("pharmacology_adjusted", ascending=False)
      .drop_duplicates(dedup_columns or ["drug_key"])
  )
  target_actionability.to_csv(
    output_dir / "target_drug_actionability.tsv", sep="\t", index=False
  )

  repurposing = result.loc[result["repurposing_candidate"].fillna(False)].copy()
  repurposing.to_csv(
    output_dir / "repurposing_hypotheses.tsv", sep="\t", index=False
  )

  component_columns = [
    column for column in result.columns
    if column.startswith("pharmacology_component_")
  ]
  result[[
    column for column in [
      "cancer", "pair_id", "analysis_lost_gene", "lost_gene", "target_gene",
      "drug_key", "drug_name", "drug_id", "coverage_adjusted_rses",
      "pharmacology_adjusted", "therapeutic_hypothesis_score", *component_columns,
    ] if column in result.columns
  ]].to_csv(
    output_dir / "pharmacology_component_matrix.tsv", sep="\t", index=False
  )

  source_coverage = (
    evidence.groupby("source", dropna=False)
      .agg(
        evidence_rows=("source", "size"),
        targets=("target_gene", "nunique"),
        drugs=("drug_name", "nunique"),
      )
      .reset_index()
    if not evidence.empty and "source" in evidence
    else pd.DataFrame(columns=["source", "evidence_rows", "targets", "drugs"])
  )
  source_coverage.to_csv(
    output_dir / "pharmacology_source_coverage.tsv", sep="\t", index=False
  )

  summary = {
    "vulnerability_rows": len(ranking),
    "pharmacology_hypothesis_rows": len(result),
    "unique_targets": result.get(
      "target_gene", pd.Series(dtype=object)
    ).nunique(),
    "unique_drugs": result.get(
      "drug_name", pd.Series(dtype=object)
    ).nunique(),
    "repurposing_hypotheses": len(repurposing),
    "sensitivity_contrasts": len(sensitivity),
  }
  (output_dir / "pharmacology_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
  )
  print(json.dumps(summary, indent=2, sort_keys=True))
  print(f"Wrote pharmacology priorities to {output_dir}")


if __name__ == "__main__":
  main()
