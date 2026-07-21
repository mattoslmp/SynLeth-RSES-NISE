#!/usr/bin/env python3
"""Create manuscript-ready empirical RSES-Onco summary tables."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]


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


def bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
  if column not in frame:
    return pd.Series(False, index=frame.index)
  values = frame[column]
  if values.dtype == bool:
    return values.fillna(False)
  return values.fillna(False).astype(str).str.casefold().isin({"1", "true", "yes"})


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--input",
    default="results/empirical_26Q1/full/empirical_rses_onco_by_cancer.tsv",
  )
  parser.add_argument(
    "--depmap-only",
    default="results/empirical_26Q1/depmap_only/empirical_rses_onco_by_cancer.tsv",
  )
  parser.add_argument("--top-n", type=int, default=20)
  parser.add_argument("--fdr", type=float, default=0.05)
  parser.add_argument(
    "--output-dir",
    default="results/empirical_26Q1/full/article_tables",
  )
  args = parser.parse_args()

  input_path = resolve_path(args.input)
  output_dir = resolve_path(args.output_dir)
  output_dir.mkdir(parents=True, exist_ok=True)
  if not input_path.exists():
    raise SystemExit(f"Input ranking not found: {input_path}")

  scores = pd.read_csv(input_path, sep="\t")
  if scores.empty:
    raise SystemExit(f"Input ranking is empty: {input_path}")
  required = {"cancer", "pair_id", "coverage_adjusted_rses"}
  missing = required - set(scores.columns)
  if missing:
    raise SystemExit(f"Ranking lacks columns: {sorted(missing)}")

  top = (
    scores.sort_values(
      ["cancer", "coverage_adjusted_rses"],
      ascending=[True, False],
    )
    .groupby("cancer", group_keys=False)
    .head(args.top_n)
  )
  top.to_csv(output_dir / "Table_main_top_rses_onco_by_cancer.tsv", sep="\t", index=False)

  empirical_mask = pd.Series(False, index=scores.index)
  for column in [
    "has_empirical_dependency",
    "has_empirical_expression",
    "has_empirical_tcga",
  ]:
    empirical_mask |= bool_series(scores, column)
  scores.loc[empirical_mask].sort_values(
    ["cancer", "coverage_adjusted_rses"],
    ascending=[True, False],
  ).to_csv(
    output_dir / "Table_main_candidates_with_empirical_evidence.tsv",
    sep="\t",
    index=False,
  )

  dependency = read_optional(input_path.with_name("dependency_contrasts.tsv"))
  if not dependency.empty:
    dependency = dependency.sort_values(
      [column for column in ["q_value_bh", "delta_effect"] if column in dependency],
      ascending=True,
    )
    dependency.to_csv(
      output_dir / "Table_S_DepMap_dependency_contrasts.tsv",
      sep="\t",
      index=False,
    )
    if {"q_value_bh", "delta_effect"}.issubset(dependency.columns):
      significant = dependency.loc[
        (dependency["q_value_bh"] < args.fdr)
        & (dependency["delta_effect"] < 0)
      ]
      significant.to_csv(
        output_dir / "Table_main_significant_synthetic_lethal_dependencies.tsv",
        sep="\t",
        index=False,
      )

  expression = read_optional(input_path.with_name("expression_contrasts.tsv"))
  if not expression.empty:
    expression = expression.sort_values(
      [column for column in ["q_value_bh", "delta_expression"] if column in expression],
      ascending=[True, False][:len([
        column for column in ["q_value_bh", "delta_expression"] if column in expression
      ])],
    )
    expression.to_csv(
      output_dir / "Table_S_DepMap_expression_contrasts.tsv",
      sep="\t",
      index=False,
    )
    if {"q_value_bh", "delta_expression"}.issubset(expression.columns):
      expression.loc[
        (expression["q_value_bh"] < args.fdr)
        & (expression["delta_expression"] > 0)
      ].to_csv(
        output_dir / "Table_main_significant_expression_compensation.tsv",
        sep="\t",
        index=False,
      )

  tcga_columns = [
    column for column in [
      "cancer", "pair_id", "lost_feature", "analysis_lost_gene", "target_gene",
      "tcga_homdel_n", "tcga_evaluable_n", "tcga_homdel_frequency",
      "component_tumor_event", "has_empirical_tcga",
    ] if column in scores.columns
  ]
  scores[tcga_columns].to_csv(
    output_dir / "Table_S_TCGA_homozygous_deletion_events.tsv",
    sep="\t",
    index=False,
  )

  skipped = read_optional(input_path.with_name("skipped_complex_biomarkers.tsv"))
  if not skipped.empty:
    skipped.to_csv(
      output_dir / "Table_S_complex_biomarkers_not_in_simple_CN_analysis.tsv",
      sep="\t",
      index=False,
    )

  depmap_path = resolve_path(args.depmap_only)
  if depmap_path.exists():
    depmap = pd.read_csv(depmap_path, sep="\t")
    comparison = depmap.merge(
      scores,
      on=["cancer", "pair_id"],
      suffixes=("_depmap", "_full"),
    )
    comparison["delta_adjusted_rses"] = (
      comparison["coverage_adjusted_rses_full"]
      - comparison["coverage_adjusted_rses_depmap"]
    )
    comparison.sort_values(
      ["cancer", "delta_adjusted_rses"],
      ascending=[True, False],
    ).to_csv(
      output_dir / "Table_S_DepMap_vs_TCGA_DepMap_score_comparison.tsv",
      sep="\t",
      index=False,
    )

  summary = pd.DataFrame({
    "metric": [
      "scored_rows",
      "unique_pairs",
      "cancers",
      "rows_with_dependency",
      "rows_with_expression",
      "rows_with_tcga",
    ],
    "value": [
      len(scores),
      scores["pair_id"].nunique(),
      scores["cancer"].nunique(),
      int(bool_series(scores, "has_empirical_dependency").sum()),
      int(bool_series(scores, "has_empirical_expression").sum()),
      int(bool_series(scores, "has_empirical_tcga").sum()),
    ],
  })
  summary.to_csv(output_dir / "analysis_summary.tsv", sep="\t", index=False)
  print(summary.to_string(index=False), flush=True)
  print(f"Wrote manuscript-ready tables to {output_dir}", flush=True)


if __name__ == "__main__":
  main()
