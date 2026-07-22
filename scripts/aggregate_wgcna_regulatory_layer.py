#!/usr/bin/env python3
"""Collapse cancer-specific WGCNA/regulatory evidence to a consensus pair prior.

Cancer-specific measurements remain in the detailed source tables. The functional
microniche receives a consensus prior so the same expression matrix is not counted as
three independent evidence units in each cancer-specific RSES-Onco score.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def median_columns(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
  available = [column for column in columns if column in frame]
  aggregation = {column: "median" for column in available}
  aggregation.update({"cancer": "nunique"})
  result = frame.groupby("pair_id", as_index=False).agg(aggregation)
  return result.rename(columns={"cancer": "consensus_cancers_observed"})


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--base", required=True)
  parser.add_argument("--cancer-specific", required=True)
  parser.add_argument("--output", required=True)
  args = parser.parse_args()

  base = pd.read_csv(resolve_path(args.base), sep="\t", low_memory=False)
  detailed = pd.read_csv(
    resolve_path(args.cancer_specific),
    sep="\t",
    low_memory=False,
  )
  if "pair_id" not in base or "pair_id" not in detailed:
    raise ValueError("Both evidence tables require pair_id")

  columns = [
    "component_wgcna_expression_network",
    "wgcna_expression_network_raw",
    "wgcna_expression_network_coverage",
    "wgcna_tom_similarity",
    "wgcna_tom_divergence",
    "wgcna_module_divergence",
    "wgcna_kME_divergence",
    "component_regulatory_network_composite",
    "regulatory_network_raw",
    "regulatory_network_coverage",
    "regulatory_tf_association_divergence",
    "regulatory_tf_expression_profile_divergence",
    "regulatory_promoter_motif_divergence",
    "regulatory_lost_regulator_count",
    "regulatory_target_regulator_count",
    "regulatory_lost_promoter_tf_count",
    "regulatory_target_promoter_tf_count",
  ]
  consensus = median_columns(detailed, columns)
  if "component_regulatory_network" in base:
    base = base.rename(columns={
      "component_regulatory_network": (
        "component_regulatory_network_dorothea_pair_level"
      )
    })
  duplicated = [
    column
    for column in consensus.columns
    if column in base.columns and column != "pair_id"
  ]
  base = base.drop(columns=duplicated, errors="ignore")
  merged = base.merge(consensus, on="pair_id", how="left")
  merged["component_regulatory_network"] = merged[
    "component_regulatory_network_composite"
  ]
  merged["expression_network_method"] = (
    "pairwise_Spearman_plus_consensus_signed_WGCNA"
  )
  merged["regulatory_network_method"] = (
    "DoRothEA_TF_target_plus_TF_expression_consistency_plus_"
    "JASPAR_promoter_motif_prediction"
  )
  merged["promoter_evidence_type"] = (
    "JASPAR_motif_prediction_not_direct_binding"
  )
  merged["regulatory_layer_version"] = (
    "wgcna-promoter-regulatory-v1"
  )

  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  merged.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)
  print(
    f"Wrote consensus WGCNA/regulatory evidence: {output} "
    f"({len(merged):,} pair rows)"
  )


if __name__ == "__main__":
  main()
