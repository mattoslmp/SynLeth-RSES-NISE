#!/usr/bin/env python3
"""Annotate row-level eligibility of every scored extended multi-omics layer."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import yaml

from rses_onco.extended_multiomics import (
  canonical_gene,
  read_model_feature_matrix,
  read_table,
)


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--config", default="config/extended_multiomics_sources.yaml"
  )
  parser.add_argument("--input-dir", default="dmap_data")
  parser.add_argument("--models", default="data/raw/depmap/Model.csv")
  parser.add_argument(
    "--evidence",
    default=(
      "data/processed/extended_multiomics/"
      "extended_pair_evidence_by_cancer.tsv"
    ),
  )
  args = parser.parse_args()

  config = yaml.safe_load(
    resolve(args.config).read_text(encoding="utf-8")
  ) or {}
  sources = config.get("sources") or {}
  data_dir = resolve(args.input_dir)
  models = read_table(resolve(args.models))
  evidence_path = resolve(args.evidence)
  evidence = pd.read_csv(evidence_path, sep="\t", low_memory=False)

  def genes_for_source(key: str) -> set[str]:
    path = data_dir / str(sources[key]["filename"])
    if not path.exists() or path.stat().st_size == 0:
      return set()
    matrix = read_model_feature_matrix(path, models, gene_features=True)
    return {
      canonical_gene(value)
      for value in matrix.columns
      if canonical_gene(value)
    }

  dependency_genes = genes_for_source("crispr_dependency")
  rnai_genes = genes_for_source("rnai_demeter2")
  protein_genes = set()
  for key in (
    "proteomics_gygi",
    "proteomics_sanger",
    "proteomics_olink",
    "proteomics_rppa",
    "proteomics_rppa500",
  ):
    protein_genes |= genes_for_source(key)

  target = evidence["target_gene"].map(canonical_gene)
  evidence["eligible_integrated_functional_loss"] = (
    pd.to_numeric(
      evidence.get("integrated_event_source_count_median"),
      errors="coerce",
    ).fillna(0) > 0
  )
  evidence["eligible_dependency_probability"] = target.isin(
    dependency_genes
  )
  evidence["eligible_protein_compensation"] = target.isin(
    protein_genes
  )
  evidence["eligible_rnai_orthogonal_support"] = target.isin(
    rnai_genes
  )

  value_columns = {
    "integrated_functional_loss": "integrated_functional_loss_support",
    "dependency_probability": "dependency_probability_support",
    "protein_compensation": "protein_compensation_support",
    "rnai_orthogonal_support": "rnai_orthogonal_support",
  }
  for layer, value_column in value_columns.items():
    eligible = evidence[f"eligible_{layer}"]
    observed = pd.to_numeric(
      evidence.get(value_column), errors="coerce"
    ).notna()
    evidence[f"{layer}_evidence_status"] = (
      "not_eligible_feature_not_measured"
    )
    evidence.loc[
      eligible & ~observed,
      f"{layer}_evidence_status",
    ] = "eligible_missing_or_insufficient_group_size"
    evidence.loc[
      eligible & observed,
      f"{layer}_evidence_status",
    ] = "observed"

  temporary = evidence_path.with_suffix(evidence_path.suffix + ".tmp")
  evidence.to_csv(temporary, sep="\t", index=False)
  temporary.replace(evidence_path)
  print(f"Annotated eligibility for {len(evidence):,} pair-context rows")
  for column in (
    "eligible_integrated_functional_loss",
    "eligible_dependency_probability",
    "eligible_protein_compensation",
    "eligible_rnai_orthogonal_support",
  ):
    print(f"{column}: {int(evidence[column].sum()):,}")


if __name__ == "__main__":
  main()
