#!/usr/bin/env python3
"""Aggregate GDC methylation beta files to promoter-level gene/sample matrices.

The workflow joins each GDC beta-value file to the official GDC GENCODE-v36
platform annotation (HM27, HM450 or EPIC), retains promoter-proximal probes, and
writes a long gene/sample table plus cancer-level gene summaries. Missing files,
unsupported platforms and unreadable inputs are reported in a status table.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROJECT_TO_CANCER = {
  "TCGA-COAD": "colon",
  "TCGA-READ": "colon",
  "TCGA-STAD": "stomach",
  "TCGA-LUAD": "lung",
  "TCGA-LUSC": "lung",
}


def resolve_path(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def first_column(frame: pd.DataFrame, aliases: Iterable[str]) -> str | None:
  lookup = {str(column).strip().casefold(): str(column) for column in frame.columns}
  for alias in aliases:
    if alias.casefold() in lookup:
      return lookup[alias.casefold()]
  return None


def platform_key(value: object) -> str | None:
  text = str(value or "").upper().replace("-", "").replace("_", "")
  if "EPIC" in text or "850" in text:
    return "EPIC"
  if "450" in text:
    return "HM450"
  if "27" in text:
    return "HM27"
  return None


def split_tokens(value: object) -> list[str]:
  if value is None or (isinstance(value, float) and np.isnan(value)):
    return []
  tokens = re.split(r"[;,|]", str(value))
  return [token.strip() for token in tokens if token.strip()]


def promoter_mask(frame: pd.DataFrame) -> pd.Series:
  mask = pd.Series(False, index=frame.index)
  distance_column = first_column(frame, (
    "distance_to_tss",
    "Distance_to_TSS",
    "pos_to_tss",
    "Position_to_TSS",
    "gene_TSS_distance",
  ))
  if distance_column:
    distances = pd.to_numeric(frame[distance_column], errors="coerce")
    mask |= distances.between(-2000, 500, inclusive="both")
  group_column = first_column(frame, (
    "gene_group",
    "UCSC_RefGene_Group",
    "regulatory_feature_group",
    "Relation_to_Island",
  ))
  if group_column:
    groups = frame[group_column].astype(str)
    mask |= groups.str.contains(
      r"TSS1500|TSS200|5['’]?UTR|1stExon|promoter",
      case=False,
      regex=True,
      na=False,
    )
  return mask


def load_annotation(path: Path) -> pd.DataFrame:
  frame = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
  probe_column = first_column(frame, (
    "probe_id",
    "Probe_ID",
    "IlmnID",
    "Composite Element",
    "Composite.Element",
    "ID",
  ))
  gene_column = first_column(frame, (
    "gene",
    "gene_name",
    "gene_symbol",
    "Gene Symbol",
    "Gene_Symbol",
    "UCSC_RefGene_Name",
    "gene_HGNC",
  ))
  if probe_column is None or gene_column is None:
    raise ValueError(
      f"Annotation {path} lacks a supported probe or gene column; "
      f"columns={list(frame.columns)}"
    )
  mask_column = first_column(frame, ("MASK_general", "mask_general"))
  keep = promoter_mask(frame)
  if not keep.any():
    raise ValueError(
      f"Annotation {path} lacks promoter/TSS information usable by the configured rule"
    )
  if mask_column:
    masked = frame[mask_column].astype(str).str.casefold().isin(
      {"true", "1", "yes", "y", "masked"}
    )
    keep &= ~masked
  selected = frame.loc[keep, [probe_column, gene_column]].copy()
  selected.columns = ["probe_id", "gene_symbols"]
  rows = []
  for record in selected.to_dict("records"):
    for gene in split_tokens(record["gene_symbols"]):
      symbol = gene.split(".")[0].strip().upper()
      if symbol:
        rows.append({"probe_id": str(record["probe_id"]), "gene": symbol})
  result = pd.DataFrame(rows).drop_duplicates()
  if result.empty:
    raise ValueError(f"No promoter probe-to-gene mappings remained for {path}")
  return result


def load_beta(path: Path) -> pd.DataFrame:
  frame = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
  probe_column = first_column(frame, (
    "Composite Element",
    "Composite.Element",
    "probe_id",
    "Probe_ID",
    "IlmnID",
    "ID",
  ))
  beta_column = first_column(frame, (
    "Beta Value",
    "Beta_Value",
    "beta_value",
    "beta",
    "value",
  ))
  if probe_column is None or beta_column is None:
    if frame.shape[1] >= 2:
      probe_column, beta_column = str(frame.columns[0]), str(frame.columns[1])
    else:
      raise ValueError(f"Beta file {path} lacks probe and beta columns")
  result = frame[[probe_column, beta_column]].copy()
  result.columns = ["probe_id", "beta"]
  result["probe_id"] = result["probe_id"].astype(str)
  result["beta"] = pd.to_numeric(result["beta"], errors="coerce")
  result = result.dropna(subset=["beta"])
  result = result.loc[result["beta"].between(0, 1, inclusive="both")]
  return result


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--manifest",
    default="data/raw/methylation/gdc_methylation_manifest.tsv",
  )
  parser.add_argument(
    "--annotation-dir",
    default="data/raw/methylation/annotations",
  )
  parser.add_argument(
    "--output",
    default="data/processed/methylation/gdc_promoter_methylation_gene_sample.tsv",
  )
  parser.add_argument(
    "--gene-summary",
    default="data/processed/methylation/gdc_promoter_methylation_gene_summary.tsv",
  )
  parser.add_argument(
    "--status-output",
    default="data/processed/methylation/gdc_promoter_methylation_aggregation_status.tsv",
  )
  args = parser.parse_args()

  manifest_path = resolve_path(args.manifest)
  annotation_dir = resolve_path(args.annotation_dir)
  output_path = resolve_path(args.output)
  summary_path = resolve_path(args.gene_summary)
  status_path = resolve_path(args.status_output)
  manifest = pd.read_csv(manifest_path, sep="\t", low_memory=False)
  annotation_paths = {
    "EPIC": annotation_dir / "EPIC.hg38.manifest.gencode.v36.tsv.gz",
    "HM27": annotation_dir / "HM27.hg38.manifest.gencode.v36.tsv.gz",
    "HM450": annotation_dir / "HM450.hg38.manifest.gencode.v36.tsv.gz",
  }
  annotations: dict[str, pd.DataFrame] = {}
  for platform, path in annotation_paths.items():
    if path.exists() and path.stat().st_size:
      annotations[platform] = load_annotation(path)

  rows: list[pd.DataFrame] = []
  status_rows = []
  for record in manifest.to_dict("records"):
    path = resolve_path(str(record.get("local_path", "")))
    platform = platform_key(record.get("platform"))
    status = "ok"
    reason = ""
    gene_count = 0
    try:
      if platform is None or platform not in annotations:
        raise ValueError("unsupported_or_missing_platform_annotation")
      if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError("methylation_beta_file_missing")
      beta = load_beta(path)
      merged = beta.merge(annotations[platform], on="probe_id", how="inner")
      if merged.empty:
        raise ValueError("no_beta_probes_matched_promoter_annotation")
      aggregated = (
        merged.groupby("gene", as_index=False)
        .agg(
          promoter_beta=("beta", "median"),
          n_promoter_probes=("probe_id", "nunique"),
        )
      )
      aggregated["cancer"] = PROJECT_TO_CANCER.get(str(record.get("project_id")))
      aggregated["project_id"] = record.get("project_id")
      aggregated["case_id"] = record.get("case_id")
      aggregated["sample_id"] = record.get("sample_id") or record.get("sample_submitter_id")
      aggregated["sample_submitter_id"] = record.get("sample_submitter_id")
      aggregated["platform"] = platform
      aggregated["file_id"] = record.get("file_id")
      aggregated["source"] = "GDC_SeSAMe_Methylation_Beta_Value"
      aggregated["promoter_definition"] = "GDC_GENCODE_v36 promoter/TSS annotation; -2000 to +500 bp or annotated promoter/TSS group"
      rows.append(aggregated)
      gene_count = len(aggregated)
    except Exception as exc:  # noqa: BLE001 - audit every unusable source file
      status = "skipped"
      reason = str(exc)
    status_rows.append({
      "file_id": record.get("file_id"),
      "project_id": record.get("project_id"),
      "sample_id": record.get("sample_id"),
      "platform": platform,
      "local_path": str(path),
      "status": status,
      "reason": reason,
      "aggregated_gene_count": gene_count,
    })

  if not rows:
    raise RuntimeError("No GDC methylation files produced promoter-level gene values")
  long = pd.concat(rows, ignore_index=True)
  long = long.dropna(subset=["cancer", "sample_id", "gene", "promoter_beta"])
  long = long.sort_values(["cancer", "sample_id", "gene", "platform"])
  output_path.parent.mkdir(parents=True, exist_ok=True)
  long.to_csv(output_path, sep="\t", index=False)

  summary = (
    long.groupby(["cancer", "gene"], as_index=False)
    .agg(
      median_promoter_beta=("promoter_beta", "median"),
      mean_promoter_beta=("promoter_beta", "mean"),
      q25_promoter_beta=("promoter_beta", lambda values: values.quantile(0.25)),
      q75_promoter_beta=("promoter_beta", lambda values: values.quantile(0.75)),
      n_samples=("sample_id", "nunique"),
      median_promoter_probe_count=("n_promoter_probes", "median"),
      platforms=("platform", lambda values: ";".join(sorted(set(values.astype(str))))),
    )
  )
  summary["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
  summary["source"] = "GDC_SeSAMe_Methylation_Beta_Value"
  summary_path.parent.mkdir(parents=True, exist_ok=True)
  summary.to_csv(summary_path, sep="\t", index=False)
  pd.DataFrame(status_rows).to_csv(status_path, sep="\t", index=False)
  metadata = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "gene_sample_rows": int(len(long)),
    "gene_summary_rows": int(len(summary)),
    "cancers": sorted(long["cancer"].astype(str).unique()),
    "interpretation_boundary": "Promoter methylation beta values are not direct proof of transcriptional silencing.",
  }
  output_path.with_suffix(".metadata.json").write_text(
    json.dumps(metadata, indent=2),
    encoding="utf-8",
  )
  print(f"Wrote promoter methylation gene/sample table: {output_path} ({len(long):,} rows)")
  print(f"Wrote promoter methylation gene summary: {summary_path} ({len(summary):,} rows)")


if __name__ == "__main__":
  main()
