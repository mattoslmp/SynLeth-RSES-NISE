#!/usr/bin/env python3
"""Validate WGCNA, promoter, TF and methylation evidence integration."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ELIGIBILITY_SEMANTICS = "eligibility-aware-v1"
EXPECTED_EXPRESSION_REGULATORY_SEMANTICS = (
  "eligibility-aware-wgcna-regulatory-v3"
)
EXPECTED_METHYLATION_SEMANTICS = "promoter-methylation-context-v1"
EXPECTED_SCORE_VERSION = "RSES-Onco-expanded-v0.10.10"
EXPECTED_REGULATORY_LAYER_VERSION = "wgcna-promoter-methylation-regulatory-v3"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(
      f"Missing or empty WGCNA/regulatory/methylation file: {path}"
    )
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


def false_claims(frame: pd.DataFrame, column: str, label: str) -> None:
  if frame[column].astype(str).str.casefold().isin({"1", "true", "yes"}).any():
    raise ValueError(label)


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
      "methylation_semantics_version",
      "pairwise_expression_context",
      "wgcna_expression_network",
      "expression_context_subcoverage",
      "regulatory_tf_association_divergence",
      "regulatory_tf_expression_profile_divergence",
      "regulatory_promoter_motif_divergence",
      "regulatory_methylation_context",
      "component_promoter_methylation_context",
      "methylation_pair_profile_divergence",
      "methylation_target_hypomethylation_support",
      "methylation_coverage",
      "methylation_absence_reason",
      "methylation_evidence_type",
      "regulatory_network_subcoverage",
      "regulatory_network_formula",
      "methylation_context_formula",
      "direct_promoter_binding_claim",
      "direct_methylation_causality_claim",
    },
    "ranking",
  )
  if set(ranking["scoring_semantics_version"].dropna().astype(str)) != {
    EXPECTED_ELIGIBILITY_SEMANTICS
  }:
    raise ValueError("Unexpected eligibility semantics")
  if set(
    ranking["expression_regulatory_semantics_version"].dropna().astype(str)
  ) != {EXPECTED_EXPRESSION_REGULATORY_SEMANTICS}:
    raise ValueError("Unexpected expression/regulatory semantics")
  if set(ranking["methylation_semantics_version"].dropna().astype(str)) != {
    EXPECTED_METHYLATION_SEMANTICS
  }:
    raise ValueError("Unexpected methylation semantics")
  if set(ranking["score_version"].dropna().astype(str)) != {
    EXPECTED_SCORE_VERSION
  }:
    raise ValueError("Unexpected score version")
  false_claims(
    ranking,
    "direct_promoter_binding_claim",
    "Ranking incorrectly claims direct promoter binding",
  )
  false_claims(
    ranking,
    "direct_methylation_causality_claim",
    "Ranking incorrectly claims causal methylation silencing",
  )
  if not ranking["regulatory_network_formula"].astype(str).str.contains(
    "0.20*CCLE_promoter_methylation_context",
    regex=False,
  ).all():
    raise ValueError("Ranking lacks methylation-aware regulatory formula")
  bounded(
    ranking,
    (
      "wgcna_expression_network",
      "expression_context_subcoverage",
      "regulatory_tf_association_divergence",
      "regulatory_tf_expression_profile_divergence",
      "regulatory_promoter_motif_divergence",
      "regulatory_methylation_context",
      "component_promoter_methylation_context",
      "methylation_pair_profile_divergence",
      "methylation_target_hypomethylation_support",
      "methylation_coverage",
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
      "component_promoter_methylation_context",
      "regulatory_methylation_context",
      "methylation_coverage",
      "methylation_absence_reason",
      "methylation_evidence_type",
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
      "component_promoter_methylation_context",
      "regulatory_methylation_context",
      "methylation_coverage",
      "regulatory_network_coverage",
    ),
    "functional evidence",
  )
  if (~functional["promoter_evidence_type"].fillna("").astype(str).eq(
    "JASPAR_motif_prediction_not_direct_binding"
  )).any():
    raise ValueError("Promoter evidence is not explicitly labelled as motif prediction")
  if (~functional["methylation_evidence_type"].fillna("").astype(str).eq(
    "CCLE_RRBS_weighted_1kb_upstream_TSS_promoter_methylation"
  )).any():
    raise ValueError("Unexpected methylation evidence type")
  if set(functional["regulatory_layer_version"].dropna().astype(str)) != {
    EXPECTED_REGULATORY_LAYER_VERSION
  }:
    raise ValueError("Unexpected regulatory layer version")
  missing_methylation = functional["component_promoter_methylation_context"].isna()
  missing_reason = functional["methylation_absence_reason"].fillna("").astype(str).str.strip().eq("")
  if (missing_methylation & missing_reason).any():
    raise ValueError("Missing methylation components lack an explicit absence reason")

  root = resolve_path(args.article_root)
  manifest = read(
    root
    / "tables/supporting_evidence/expression_regulatory/"
      "wgcna_regulatory_supporting_evidence_manifest.tsv"
  )
  require_columns(
    manifest,
    {
      "evidence_family", "source_path", "output_path", "rows", "sha256",
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
    "promoter_methylation_pair_metrics",
    "ensembl_canonical_promoters",
    "jaspar_promoter_motif_predictions",
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

  methylation = exported["promoter_methylation_pair_metrics"]
  require_columns(
    methylation,
    {
      "pair_id", "cancer", "methylation_pair_profile_divergence",
      "methylation_target_hypomethylation_support",
      "component_promoter_methylation_context", "methylation_coverage",
      "methylation_absence_reason", "methylation_evidence_type",
    },
    "promoter methylation metrics",
  )
  bounded(
    methylation,
    (
      "methylation_pair_profile_divergence",
      "methylation_target_hypomethylation_support",
      "component_promoter_methylation_context",
      "methylation_coverage",
    ),
    "promoter methylation metrics",
  )
  absent_metric = methylation["component_promoter_methylation_context"].isna()
  absent_reason = methylation["methylation_absence_reason"].fillna("").astype(str).str.strip().eq("")
  if (absent_metric & absent_reason).any():
    raise ValueError("Methylation metric missingness is not traceable")

  diagnostics = exported["wgcna_run_diagnostics"]
  require_columns(
    diagnostics,
    {
      "cancer", "correlation", "correlation_policy", "max_p_outliers",
      "pearson_fallback", "signed_kme_correlation",
      "signed_kme_max_p_outliers", "signed_kme_pearson_fallback",
      "zero_mad_gene_count", "zero_mad_module_eigengene_count",
      "pearson_fallback_entity_count",
    },
    "WGCNA run diagnostics",
  )
  expected_policy = (
    "bicor_primary_with_individual_Pearson_fallback_for_zero_MAD_only"
  )
  if not diagnostics["correlation"].astype(str).eq("bicor").all():
    raise ValueError("WGCNA primary correlation must be bicor")
  if not diagnostics["correlation_policy"].astype(str).eq(expected_policy).all():
    raise ValueError("Unexpected WGCNA correlation policy")
  if not diagnostics["pearson_fallback"].astype(str).eq("individual").all():
    raise ValueError("WGCNA Pearson fallback must be individual")
  if not diagnostics["signed_kme_correlation"].astype(str).eq("bicor").all():
    raise ValueError("signedKME primary correlation must be bicor")
  if not diagnostics["signed_kme_pearson_fallback"].astype(str).eq("individual").all():
    raise ValueError("signedKME Pearson fallback must be individual")
  for column in ("max_p_outliers", "signed_kme_max_p_outliers"):
    values = pd.to_numeric(diagnostics[column], errors="coerce")
    if not values.eq(0.10).all():
      raise ValueError(f"{column} must equal 0.10")

  fallback = exported["wgcna_correlation_fallback_audit"]
  require_columns(
    fallback,
    {
      "cancer", "entity_type", "entity", "mad", "pearson_fallback_expected",
      "fallback_reason", "primary_correlation", "fallback_correlation",
      "pearson_fallback_policy", "max_p_outliers",
    },
    "WGCNA correlation fallback audit",
  )
  if not fallback["primary_correlation"].astype(str).eq("bicor").all():
    raise ValueError("Fallback audit primary correlation must be bicor")
  if not fallback["fallback_correlation"].astype(str).eq("pearson").all():
    raise ValueError("Fallback audit fallback correlation must be Pearson")
  if not fallback["pearson_fallback_policy"].astype(str).eq("individual").all():
    raise ValueError("Fallback audit policy must be individual")
  fallback_expected = fallback["pearson_fallback_expected"].astype(str).str.casefold().isin({"1", "true", "yes"})
  invalid_reason = fallback_expected & fallback["fallback_reason"].fillna("").astype(str).str.strip().eq("")
  if invalid_reason.any():
    raise ValueError("Fallback-eligible entities lack a zero-MAD reason")

  print("WGCNA/promoter/methylation regulatory evidence validation passed.")
  print(f"Ranking rows: {len(ranking):,}")
  print(f"Functional pair rows: {len(functional):,}")
  print(f"Methylation pair-cancer rows: {len(methylation):,}")


if __name__ == "__main__":
  main()
