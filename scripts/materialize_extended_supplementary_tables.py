#!/usr/bin/env python3
"""Materialize source-backed supplementary Tables S26-S47.

Every table is copied or generated from a traceable pipeline source. Missing or
empty optional sources remain explicitly unavailable and are never replaced with
invented biological evidence.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil

import pandas as pd

from rses_onco.evidence_categories import evidence_category_definitions

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_tsv(frame: pd.DataFrame, destination: Path) -> None:
  destination.parent.mkdir(parents=True, exist_ok=True)
  temporary = destination.with_suffix(destination.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(destination)


def atomic_copy(source: Path, destination: Path) -> None:
  if not source.exists() or source.stat().st_size == 0:
    raise FileNotFoundError(
      f"Required supplementary source is missing or empty: {source}"
    )
  destination.parent.mkdir(parents=True, exist_ok=True)
  temporary = destination.with_suffix(destination.suffix + ".tmp")
  shutil.copy2(source, temporary)
  temporary.replace(destination)


def atomic_optional_copy(source: Path, destination: Path) -> None:
  if source.exists() and source.stat().st_size > 0:
    atomic_copy(source, destination)
    return
  atomic_tsv(
    pd.DataFrame([{
      "evidence_status": "unavailable_or_not_acquired",
      "source_path": str(source),
      "absence_reason": (
        "Optional methylation source was not available for this run."
      ),
      "interpretation_boundary": (
        "Unavailable promoter methylation is not negative biological "
        "evidence."
      ),
    }]),
    destination,
  )


def bootstrap_reproduction_registry(destination: Path) -> None:
  registry = destination / "Table_S44_asset_reproduction_registry.tsv"
  if registry.exists() and registry.stat().st_size > 0:
    return
  atomic_tsv(
    pd.DataFrame([{
      "asset_id": "registry_pending_finalization",
      "asset_type": "pipeline_status",
      "category": "supplementary",
      "scientific_title": "Asset reproduction registry",
      "script": "scripts/build_asset_reproduction_registry.py",
      "primary_inputs": "config/article_assets.yaml",
      "intermediate_data": "article_outputs/manifests",
      "reproduction_command": (
        "bash scripts/run_publication_pipeline.sh assets-only"
      ),
      "outputs": str(registry),
      "dependencies": "environment.yml; pyproject.toml",
      "document_location": "supplementary data",
      "status": (
        "bootstrap row; replaced by the complete registry later in the run"
      ),
    }]),
    registry,
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument(
    "--methylation-evidence",
    default=(
      "data/processed/methylation/"
      "pair_promoter_methylation_evidence.tsv"
    ),
  )
  parser.add_argument(
    "--methylation-gene-summary",
    default=(
      "data/processed/methylation/"
      "gdc_promoter_methylation_gene_summary.tsv"
    ),
  )
  parser.add_argument(
    "--methylation-source-status",
    default=(
      "data/processed/methylation/"
      "gdc_promoter_methylation_aggregation_status.tsv"
    ),
  )
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  destination = article_root / "tables/supplementary"

  mapping = {
    "Table_S26_missing_data_sensitivity.tsv": (
      article_root / "tables/robustness/missing_data_sensitivity.tsv"
    ),
    "Table_S27_evidence_category_assignments.tsv": (
      article_root / "tables/qc/evidence_category_assignments.tsv"
    ),
    "Table_S28_coexpression_by_event_group.tsv": (
      article_root
      / "tables/supporting_evidence/expression/coexpression_by_event_group.tsv"
    ),
    "Table_S29_compensation_and_dependency_contrasts.tsv": (
      article_root
      / "tables/supporting_evidence/expression/compensation_and_dependency_contrasts.tsv"
    ),
    "Table_S30_model_level_expression_crispr_copy_number.tsv": (
      article_root
      / "tables/supporting_evidence/model_level/model_level_expression_crispr_copy_number.tsv"
    ),
    "Table_S31_string_candidate_edges_all_channels.tsv": (
      article_root
      / "tables/supporting_evidence/networks/raw_sources/string_candidate_edges_all_channels.tsv"
    ),
    "Table_S32_dorothea_candidate_regulatory_edges.tsv": (
      article_root
      / "tables/supporting_evidence/networks/raw_sources/dorothea_candidate_regulatory_edges.tsv"
    ),
    "Table_S33_hpa_candidate_localization.tsv": (
      article_root
      / "tables/supporting_evidence/localization/hpa_candidate_localization.tsv"
    ),
    "Table_S34_uniprot_candidate_annotations.tsv": (
      article_root
      / "tables/supporting_evidence/structures/uniprot_candidate_annotations.tsv"
    ),
    "Table_S35_wgcna_pair_metrics_all_cancers.tsv": (
      ROOT
      / "data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv"
    ),
    "Table_S36_wgcna_correlation_fallback_audit.tsv": (
      ROOT
      / "data/processed/regulatory/wgcna/wgcna_correlation_fallback_all_cancers.tsv"
    ),
    "Table_S37_wgcna_run_diagnostics.tsv": (
      ROOT
      / "data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv"
    ),
    "Table_S38_jaspar_promoter_tf_summary.tsv": (
      ROOT / "data/processed/regulatory/jaspar_promoter_tf_summary.tsv"
    ),
    "Table_S39_conditional_dependency_contrasts.tsv": (
      article_root
      / "tables/supporting_evidence/phenotypes/conditional_dependency_contrasts.tsv"
    ),
    "Table_S40_tcga_gene_event_summary.tsv": (
      article_root
      / "tables/supporting_evidence/genomic_context/tcga_gene_event_summary.tsv"
    ),
    "Table_S41_pharmacology_evidence_long_support.tsv": (
      article_root
      / "tables/supporting_evidence/pharmacology/pharmacology_evidence_long.tsv"
    ),
    "Table_S42_pharmacology_source_coverage.tsv": (
      ROOT
      / "results/expanded_26Q1/pharmacology/pharmacology_source_coverage.tsv"
    ),
  }
  for name, source in mapping.items():
    atomic_copy(source, destination / name)

  atomic_tsv(
    evidence_category_definitions(),
    destination / "Table_S43_evidence_category_definitions.tsv",
  )
  bootstrap_reproduction_registry(destination)

  optional_mapping = {
    "Table_S45_pair_promoter_methylation_evidence.tsv": resolve_path(
      args.methylation_evidence
    ),
    "Table_S46_gene_promoter_methylation_summary.tsv": resolve_path(
      args.methylation_gene_summary
    ),
    "Table_S47_methylation_source_status.tsv": resolve_path(
      args.methylation_source_status
    ),
  }
  for name, source in optional_mapping.items():
    atomic_optional_copy(source, destination / name)

  print(
    f"Materialized {len(mapping) + 2 + len(optional_mapping)} extended "
    f"supplementary tables in {destination}"
  )


if __name__ == "__main__":
  main()
