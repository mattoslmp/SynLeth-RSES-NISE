#!/usr/bin/env python3
"""Validate WGCNA, promoter and TF regulatory evidence and score integration."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SEMANTICS = "eligibility-aware-wgcna-regulatory-v2"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(f"Missing or empty WGCNA/regulatory file: {path}")
  return pd.read_csv(path, sep="\t", low_memory=False)


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
  missing = sorted(columns - set(frame.columns))
  if missing:
    raise ValueError(f"{label} missing columns: {missing}")


def bounded(frame: pd.DataFrame, columns: tuple[str, ...], label: str) -> None:
  for column in columns:
    if column not in frame:
      continue
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    if not values.between(0, 1).all():
      raise ValueError(f"{label} has values outside [0,1] in {column}")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--functional-evidence",
    default="data/processed/expanded_pair_functional_evidence.tsv",
  )
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()

  ranking = read(resolve_path(args.ranking))
  functional = read(resolve_path(args.functional_evidence))
  require_columns(
    ranking,
    {
      "scoring_semantics_version",
      "pairwise_expression_context",
      "wgcna_expression_network",
      "expression_context_subcoverage",
      "regulatory_tf_association_divergence",
      "regulatory_tf_expression_profile_divergence",
      "regulatory_promoter_motif_divergence",
      "regulatory_network_subcoverage",
      "direct_promoter_binding_claim",
    },
    "ranking",
  )
  versions = set(
    ranking["scoring_semantics_version"].dropna().astype(str)
  )
  if versions != {EXPECTED_SEMANTICS}:
    raise ValueError(
      f"Unexpected scoring semantics {sorted(versions)}; "
      f"expected {EXPECTED_SEMANTICS}"
    )
  if ranking["direct_promoter_binding_claim"].astype(str).str.casefold().isin(
    {"1", "true", "yes"}
  ).any():
    raise ValueError("Ranking incorrectly claims direct promoter binding")
  bounded(
    ranking,
    (
      "wgcna_expression_network",
      "expression_context_subcoverage",
      "regulatory_tf_association_divergence",
      "regulatory_tf_expression_profile_divergence",
      "regulatory_promoter_motif_divergence",
      "regulatory_network_subcoverage",
      "coverage_adjusted_rses",
      "evidence_coverage",
    ),
    "ranking",
  )

  require_columns(
    functional,
    {
      "pair_id",
      "component_wgcna_expression_network",
      "component_regulatory_network_composite",
      "regulatory_network_coverage",
      "promoter_evidence_type",
      "regulatory_layer_version",
    },
    "functional evidence",
  )
  bounded(
    functional,
    (
      "component_wgcna_expression_network",
      "wgcna_expression_network_coverage",
      "component_regulatory_network_composite",
      "regulatory_network_coverage",
    ),
    "functional evidence",
  )
  invalid_promoter = ~functional["promoter_evidence_type"].fillna("").astype(str).eq(
    "JASPAR_motif_prediction_not_direct_binding"
  )
  if invalid_promoter.any():
    raise ValueError("Promoter evidence is not explicitly labelled as motif prediction")

  root = resolve_path(args.article_root)
  manifest = read(
    root
    / "tables/supporting_evidence/expression_regulatory/"
      "wgcna_regulatory_supporting_evidence_manifest.tsv"
  )
  require_columns(
    manifest,
    {
      "evidence_family",
      "source_path",
      "output_path",
      "rows",
      "sha256",
      "interpretation_boundary",
    },
    "WGCNA/regulatory support manifest",
  )
  required_families = {
    "cancer_specific_wgcna_pair_metrics",
    "wgcna_input_preparation",
    "promoter_tf_regulatory_pair_metrics",
    "ensembl_canonical_promoters",
    "jaspar_promoter_motif_predictions",
  }
  missing_families = sorted(
    required_families - set(manifest["evidence_family"].astype(str))
  )
  if missing_families:
    raise ValueError(f"WGCNA/regulatory manifest missing: {missing_families}")
  for value in manifest["output_path"].astype(str):
    path = Path(value)
    if not path.is_absolute():
      path = ROOT / path
    read(path)

  print("WGCNA/promoter regulatory evidence validation passed.")
  print(f"Ranking rows: {len(ranking):,}")
  print(f"Functional pair rows: {len(functional):,}")


if __name__ == "__main__":
  main()
