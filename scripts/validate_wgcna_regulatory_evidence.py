#!/usr/bin/env python3
"""Validate WGCNA, promoter and TF regulatory evidence and score integration."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ELIGIBILITY_SEMANTICS = "eligibility-aware-v1"
EXPECTED_EXPRESSION_REGULATORY_SEMANTICS = (
  "eligibility-aware-wgcna-regulatory-methylation-v4"
)


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
      "expression_regulatory_semantics_version",
      "score_version",
      "pairwise_expression_context",
      "wgcna_expression_network",
      "expression_context_subcoverage",
      "regulatory_tf_association_divergence",
      "regulatory_tf_expression_profile_divergence",
      "regulatory_promoter_motif_divergence",
      "regulatory_promoter_methylation_context",
      "methylation_context_subcoverage",
      "regulatory_network_subcoverage",
      "direct_promoter_binding_claim",
    },
    "ranking",
  )
  eligibility_versions = set(
    ranking["scoring_semantics_version"].dropna().astype(str)
  )
  if eligibility_versions != {EXPECTED_ELIGIBILITY_SEMANTICS}:
    raise ValueError(
      f"Unexpected eligibility semantics {sorted(eligibility_versions)}; "
      f"expected {EXPECTED_ELIGIBILITY_SEMANTICS}"
    )
  expression_regulatory_versions = set(
    ranking["expression_regulatory_semantics_version"]
      .dropna()
      .astype(str)
  )
  if expression_regulatory_versions != {
    EXPECTED_EXPRESSION_REGULATORY_SEMANTICS
  }:
    raise ValueError(
      "Unexpected expression/regulatory semantics "
      f"{sorted(expression_regulatory_versions)}; expected "
      f"{EXPECTED_EXPRESSION_REGULATORY_SEMANTICS}"
    )
  score_versions = set(ranking["score_version"].dropna().astype(str))
  if score_versions != {"RSES-Onco-expanded-v0.11.1"}:
    raise ValueError(
      f"Unexpected score versions {sorted(score_versions)}; "
      "expected RSES-Onco-expanded-v0.11.1"
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
  invalid_promoter = ~functional[
    "promoter_evidence_type"
  ].fillna("").astype(str).eq(
    "JASPAR_motif_prediction_not_direct_binding"
  )
  if invalid_promoter.any():
    raise ValueError(
      "Promoter evidence is not explicitly labelled as motif prediction"
    )

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
    "wgcna_correlation_fallback_audit",
    "wgcna_run_diagnostics",
    "promoter_tf_regulatory_pair_metrics",
    "ensembl_canonical_promoters",
    "jaspar_promoter_motif_predictions",
    "tcga_gdc_methylation_pair_metrics",
    "tcga_gdc_methylation_source_status",
  }
  missing_families = sorted(
    required_families - set(manifest["evidence_family"].astype(str))
  )
  if missing_families:
    raise ValueError(f"WGCNA/regulatory manifest missing: {missing_families}")
  exported = {}
  for record in manifest.to_dict("records"):
    path = Path(str(record["output_path"]))
    if not path.is_absolute():
      path = ROOT / path
    exported[str(record["evidence_family"])] = read(path)

  diagnostics = exported["wgcna_run_diagnostics"]
  require_columns(
    diagnostics,
    {
      "cancer",
      "correlation",
      "correlation_policy",
      "max_p_outliers",
      "pearson_fallback",
      "signed_kme_correlation",
      "signed_kme_max_p_outliers",
      "signed_kme_pearson_fallback",
      "zero_mad_gene_count",
      "zero_mad_module_eigengene_count",
      "pearson_fallback_entity_count",
    },
    "WGCNA run diagnostics",
  )
  expected_policy = (
    "bicor_primary_with_individual_Pearson_fallback_for_zero_MAD_only"
  )
  if not diagnostics["correlation"].astype(str).eq("bicor").all():
    raise ValueError("WGCNA primary correlation must be bicor")
  if not diagnostics["correlation_policy"].astype(str).eq(
    expected_policy
  ).all():
    raise ValueError("Unexpected WGCNA correlation policy")
  if not diagnostics["pearson_fallback"].astype(str).eq(
    "individual"
  ).all():
    raise ValueError("WGCNA Pearson fallback must be individual")
  if not diagnostics["signed_kme_correlation"].astype(str).eq(
    "bicor"
  ).all():
    raise ValueError("signedKME primary correlation must be bicor")
  if not diagnostics["signed_kme_pearson_fallback"].astype(str).eq(
    "individual"
  ).all():
    raise ValueError("signedKME Pearson fallback must be individual")
  for column in ("max_p_outliers", "signed_kme_max_p_outliers"):
    values = pd.to_numeric(diagnostics[column], errors="coerce")
    if not values.eq(0.10).all():
      raise ValueError(f"{column} must equal 0.10")

  fallback = exported["wgcna_correlation_fallback_audit"]
  require_columns(
    fallback,
    {
      "cancer",
      "entity_type",
      "entity",
      "mad",
      "pearson_fallback_expected",
      "fallback_reason",
      "primary_correlation",
      "fallback_correlation",
      "pearson_fallback_policy",
      "max_p_outliers",
    },
    "WGCNA correlation fallback audit",
  )
  if not fallback["primary_correlation"].astype(str).eq("bicor").all():
    raise ValueError("Fallback audit primary correlation must be bicor")
  if not fallback["fallback_correlation"].astype(str).eq(
    "pearson"
  ).all():
    raise ValueError("Fallback audit fallback correlation must be Pearson")
  if not fallback["pearson_fallback_policy"].astype(str).eq(
    "individual"
  ).all():
    raise ValueError("Fallback audit policy must be individual")
  fallback_expected = (
    fallback["pearson_fallback_expected"]
      .astype(str)
      .str.casefold()
      .isin({"1", "true", "yes"})
  )
  invalid_reason = (
    fallback_expected
    & fallback["fallback_reason"]
      .fillna("")
      .astype(str)
      .str.strip()
      .eq("")
  )
  if invalid_reason.any():
    raise ValueError("Fallback-eligible entities lack a zero-MAD reason")

  print("WGCNA/promoter regulatory evidence validation passed.")
  print(f"Ranking rows: {len(ranking):,}")
  print(f"Functional pair rows: {len(functional):,}")


if __name__ == "__main__":
  main()
