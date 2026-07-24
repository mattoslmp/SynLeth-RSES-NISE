#!/usr/bin/env python3
"""Validate scientific and publication integrity of RSES-Onco v0.12.0."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import pandas as pd

from rses_onco.extended_multiomics import (
  DIRECT_SCORE_SOURCE_KEYS,
  sha256_file,
)


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def require(path: Path, errors: list[str]) -> None:
  if not path.exists() or path.stat().st_size == 0:
    errors.append(f"missing_or_empty:{path}")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--processed-dir",
    default="data/processed/extended_multiomics",
  )
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument("--strict-sources", action="store_true")
  parser.add_argument("--require-publication-assets", action="store_true")
  args = parser.parse_args()

  processed = resolve(args.processed_dir)
  ranking_path = resolve(args.ranking)
  article_root = resolve(args.article_root)
  errors: list[str] = []
  required = {
    "inventory": processed / "extended_multiomics_source_inventory.tsv",
    "status": processed / "extended_multiomics_source_status.tsv",
    "loss": processed / "functional_loss_states.tsv",
    "evidence": processed / "extended_pair_evidence_by_cancer.tsv",
    "protein": processed / "proteomics_pair_evidence_by_source.tsv",
    "covariates": processed / "extended_covariate_context.tsv",
    "drug": processed / "custom_drug_sensitivity_long.tsv",
    "combination": processed / "gdsc_combination_evidence_long.tsv",
    "provenance": processed / "extended_multiomics_source_provenance.tsv",
    "status_json": processed / "extended_multiomics_status.json",
    "ranking": ranking_path,
  }
  for path in required.values():
    require(path, errors)
  if errors:
    raise SystemExit("\n".join(errors))

  inventory = pd.read_csv(
    required["inventory"], sep="\t", low_memory=False
  )
  status = pd.read_csv(required["status"], sep="\t", low_memory=False)
  loss = pd.read_csv(required["loss"], sep="\t", low_memory=False)
  evidence = pd.read_csv(
    required["evidence"], sep="\t", low_memory=False
  )
  protein = pd.read_csv(
    required["protein"], sep="\t", low_memory=False
  )
  covariates = pd.read_csv(
    required["covariates"], sep="\t", low_memory=False
  )
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  payload = json.loads(
    required["status_json"].read_text(encoding="utf-8")
  )

  if inventory["source_key"].duplicated().any():
    errors.append("duplicate_source_keys")
  for record in inventory.loc[
    inventory["exists"].astype(bool)
  ].to_dict("records"):
    path = Path(str(record["path"]))
    path = path if path.is_absolute() else ROOT / path
    if (
      path.exists()
      and str(record.get("sha256") or "") != sha256_file(path)
    ):
      errors.append(f"source_sha256_mismatch:{record['source_key']}")
  if args.strict_sources:
    available = set(
      inventory.loc[
        inventory["exists"].astype(bool),
        "source_key",
      ].astype(str)
    )
    missing_direct = sorted(DIRECT_SCORE_SOURCE_KEYS - available)
    if missing_direct:
      errors.append(f"missing_direct_score_sources:{missing_direct}")
    standardized_direct = set(
      status.loc[
        status["status"].astype(str).eq("ok"),
        "source_key",
      ].astype(str)
    )
    failed_direct = sorted(
      DIRECT_SCORE_SOURCE_KEYS - standardized_direct
    )
    if failed_direct:
      errors.append(
        f"unstandardized_direct_score_sources:{failed_direct}"
      )

  required_loss_columns = {
    "ModelID",
    "gene",
    "functional_loss_score",
    "functional_loss_state",
    "event_source_count",
  }
  if not required_loss_columns.issubset(loss.columns):
    errors.append(
      "loss_table_missing_columns:"
      f"{sorted(required_loss_columns - set(loss.columns))}"
    )
  valid_states = {
    "biallelic_or_homdel",
    "probable_functional_loss",
    "partial_loss",
    "intact",
    "missing",
  }
  observed_states = set(loss["functional_loss_state"].astype(str))
  if not observed_states.issubset(valid_states):
    errors.append(
      "invalid_functional_loss_states:"
      f"{sorted(observed_states - valid_states)}"
    )
  if loss.duplicated(["ModelID", "gene"]).any():
    errors.append("duplicate_model_gene_loss_states")

  if evidence.duplicated(["pair_id", "cancer"]).any():
    errors.append("duplicate_extended_pair_evidence")
  scored_columns = {
    "integrated_functional_loss_support",
    "dependency_probability_support",
    "protein_compensation_support",
    "rnai_orthogonal_support",
  }
  if not scored_columns.issubset(evidence.columns):
    errors.append(
      "extended_evidence_missing_scored_columns:"
      f"{sorted(scored_columns - set(evidence.columns))}"
    )
  for column in scored_columns:
    values = pd.to_numeric(evidence[column], errors="coerce").dropna()
    if ((values < 0) | (values > 1)).any():
      errors.append(f"scored_component_outside_0_1:{column}")

  if not protein.empty and "source" in protein.columns:
    expected_protein_sources = {
      "proteomics_gygi",
      "proteomics_sanger",
      "proteomics_olink",
      "proteomics_rppa",
      "proteomics_rppa500",
    }
    observed = set(protein["source"].astype(str))
    if (
      args.strict_sources
      and not expected_protein_sources.issubset(observed)
    ):
      errors.append(
        "missing_protein_source_records:"
        f"{sorted(expected_protein_sources - observed)}"
      )

  if not covariates.empty:
    scored = covariates.get(
      "scored_in_primary_rses",
      pd.Series(False, index=covariates.index),
    )
    if scored.astype(bool).any():
      errors.append("exploratory_covariate_was_scored_in_primary_rses")

  required_ranking_columns = {
    "baseline_coverage_adjusted_rses",
    "coverage_adjusted_rses",
    "extended_score_delta",
    "baseline_rank_within_cancer",
    "extended_rank_within_cancer",
    "extended_rank_change",
    "extended_multiomics_semantics_version",
    "score_version",
    "ablation_without_integrated_functional_loss_coverage_adjusted_rses",
    "ablation_without_dependency_probability_coverage_adjusted_rses",
    "ablation_without_protein_compensation_coverage_adjusted_rses",
    "ablation_without_rnai_orthogonal_support_coverage_adjusted_rses",
  }
  if not required_ranking_columns.issubset(ranking.columns):
    errors.append(
      "ranking_missing_extended_columns:"
      f"{sorted(required_ranking_columns - set(ranking.columns))}"
    )
  elif not ranking["score_version"].astype(str).eq(
    "RSES-Onco-expanded-v0.12.0"
  ).all():
    errors.append("ranking_contains_non_v0120_score_rows")

  if not str(payload.get("version", "")).endswith("v0.12.0"):
    errors.append("extended_status_version_mismatch")
  if payload.get("pair_evidence_rows") != len(evidence):
    errors.append("extended_status_pair_count_mismatch")

  figure_manifest = article_root / "manifests/figure_manifest.tsv"
  table_manifest = article_root / "manifests/table_manifest.tsv"
  if args.require_publication_assets:
    require(figure_manifest, errors)
    require(table_manifest, errors)
    if figure_manifest.exists():
      figures = pd.read_csv(
        figure_manifest, sep="\t", low_memory=False
      )
      expected_figures = {
        f"Figure_S{index}" for index in range(71, 79)
      }
      observed = set(figures["figure_id"].astype(str))
      if not expected_figures.issubset(observed):
        errors.append(
          "missing_extended_figures:"
          f"{sorted(expected_figures - observed)}"
        )
    if table_manifest.exists():
      tables = pd.read_csv(
        table_manifest, sep="\t", low_memory=False
      )
      expected_tables = {
        f"Table_S{index}" for index in range(53, 65)
      }
      observed = set(
        tables["table_id"]
          .astype(str)
          .str.extract(r"(Table_S\d+)", expand=False)
          .dropna()
      )
      if not expected_tables.issubset(observed):
        errors.append(
          "missing_extended_tables:"
          f"{sorted(expected_tables - observed)}"
        )

  validation = {
    "status": "pass" if not errors else "fail",
    "errors": errors,
    "available_sources": int(inventory["exists"].sum()),
    "standardized_sources": int(
      status["status"].astype(str).eq("ok").sum()
    ),
    "functional_loss_rows": len(loss),
    "pair_evidence_rows": len(evidence),
    "protein_evidence_rows": len(protein),
    "covariate_rows": len(covariates),
    "ranking_rows": len(ranking),
  }
  validation_path = (
    processed / "extended_multiomics_integrity_validation.json"
  )
  validation_path.write_text(
    json.dumps(validation, indent=2), encoding="utf-8"
  )
  if errors:
    raise SystemExit(
      "Extended multi-omics validation failed:\n"
      + "\n".join(f"- {error}" for error in errors)
    )
  print("Extended multi-omics integrity validation passed.")
  print(json.dumps(validation, indent=2))


if __name__ == "__main__":
  main()
