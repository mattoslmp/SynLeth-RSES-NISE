#!/usr/bin/env python3
"""Validate model-level and raw-source supporting evidence without overclaiming."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(f"Missing or empty supporting-evidence file: {path}")
  return pd.read_csv(path, sep="\t", low_memory=False)


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
  missing = sorted(columns - set(frame.columns))
  if missing:
    raise ValueError(f"{label} missing columns: {missing}")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  root = resolve_path(args.article_root)
  support = root / "tables/supporting_evidence"

  model_manifest = read(
    support / "model_level_supporting_evidence_manifest.tsv"
  )
  raw_manifest = read(
    support / "networks/raw_sources/raw_functional_evidence_manifest.tsv"
  )
  require_columns(
    model_manifest,
    {"evidence_family", "path", "rows", "sha256", "status"},
    "model-level manifest",
  )
  require_columns(
    raw_manifest,
    {"evidence_family", "output_path", "rows", "sha256", "status"},
    "raw-source manifest",
  )
  for value in model_manifest["path"].astype(str):
    read(Path(value) if Path(value).is_absolute() else ROOT / value)
  for value in raw_manifest["output_path"].astype(str):
    read(Path(value) if Path(value).is_absolute() else ROOT / value)

  model_level = read(
    support
    / "model_level/model_level_expression_crispr_copy_number.tsv"
  )
  require_columns(
    model_level,
    {
      "cancer", "pair_id", "lost_gene", "target_gene", "ModelID",
      "lost_copy_number", "loss_group", "lost_expression",
      "target_expression", "lost_gene_effect", "target_gene_effect",
    },
    "model-level measurements",
  )
  invalid_groups = set(model_level["loss_group"].astype(str)) - {
    "loss", "intact", "copy_number_unavailable"
  }
  if invalid_groups:
    raise ValueError(f"Unexpected model-level loss groups: {sorted(invalid_groups)}")

  correlations = read(
    support / "expression/coexpression_by_event_group.tsv"
  )
  require_columns(
    correlations,
    {
      "measurement", "event_stratum", "n_models", "spearman_rho",
      "p_value", "q_value_bh", "interpretation_boundary",
    },
    "coexpression table",
  )
  if correlations["interpretation_boundary"].fillna("").str.contains(
    "automatic evidence",
    case=False,
    regex=False,
  ).all() is False:
    raise ValueError("Coexpression rows do not preserve the interpretation boundary")

  contrasts = read(
    support / "expression/compensation_and_dependency_contrasts.tsv"
  )
  require_columns(
    contrasts,
    {
      "analysis", "n_loss", "n_intact", "delta_loss_minus_intact",
      "p_value", "q_value_bh", "analysis_status",
    },
    "compensation/dependency contrasts",
  )
  available = contrasts["analysis_status"].eq("available")
  if (
    (pd.to_numeric(contrasts.loc[available, "n_loss"], errors="coerce") < 3).any()
    or (pd.to_numeric(contrasts.loc[available, "n_intact"], errors="coerce") < 3).any()
  ):
    raise ValueError("An available contrast has fewer than three models per group")

  string_edges = read(
    support
    / "networks/raw_sources/string_candidate_edges_all_channels.tsv"
  )
  if "combined_score_interpretation" in string_edges:
    if not string_edges["combined_score_interpretation"].fillna("").str.contains(
      "not direct experimental evidence",
      case=False,
      regex=False,
    ).all():
      raise ValueError("STRING combined-score interpretation is incomplete")

  dorothea = read(
    support
    / "networks/raw_sources/dorothea_candidate_regulatory_edges.tsv"
  )
  if "promoter_binding_evidence" in dorothea:
    claimed = ~dorothea["promoter_binding_evidence"].astype(str).eq(
      "not_available_from_this_table"
    )
    if claimed.any():
      raise ValueError("DoRothEA rows incorrectly claim promoter-binding evidence")

  promoter = read(
    support / "networks/raw_sources/promoter_evidence_status.tsv"
  )
  require_columns(
    promoter,
    {"status", "reason", "scientific_rule"},
    "promoter evidence status",
  )
  if not promoter["status"].eq(
    "not_available_in_current_pipeline_sources"
  ).all():
    raise ValueError("Promoter evidence availability is misrepresented")

  print("Extended supporting-evidence validation passed.")
  print(f"Model-level rows: {len(model_level):,}")
  print(f"Coexpression rows: {len(correlations):,}")
  print(f"Compensation/dependency contrast rows: {len(contrasts):,}")


if __name__ == "__main__":
  main()
