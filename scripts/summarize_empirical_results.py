#!/usr/bin/env python3
"""Create manuscript-ready benchmark or expanded RSES-Onco summary tables."""
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


def read_first(parent: Path, names: list[str]) -> pd.DataFrame:
  for name in names:
    frame = read_optional(parent / name)
    if not frame.empty:
      return frame
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
    scores.sort_values(["cancer", "coverage_adjusted_rses"], ascending=[True, False])
      .groupby("cancer", group_keys=False)
      .head(args.top_n)
  )
  top.to_csv(output_dir / "Table_main_top_rses_onco_by_cancer.tsv", sep="\t", index=False)

  empirical_mask = pd.Series(False, index=scores.index)
  for column in [
    "has_empirical_dependency",
    "has_empirical_expression",
    "has_empirical_expression_compensation",
    "has_empirical_expression_context",
    "has_empirical_phenotype_profile",
    "has_empirical_tcga",
  ]:
    empirical_mask |= bool_series(scores, column)
  scores.loc[empirical_mask].sort_values(
    ["cancer", "coverage_adjusted_rses"], ascending=[True, False]
  ).to_csv(
    output_dir / "Table_main_candidates_with_empirical_evidence.tsv",
    sep="\t",
    index=False,
  )

  parent = input_path.parent
  dependency = read_first(parent, [
    "dependency_contrasts.tsv",
    "expanded_dependency_contrasts.tsv",
  ])
  if not dependency.empty:
    sort_columns = [column for column in ["q_value_bh", "delta_effect"] if column in dependency]
    if sort_columns:
      dependency = dependency.sort_values(sort_columns, ascending=True)
    dependency.to_csv(
      output_dir / "Table_S_DepMap_dependency_contrasts.tsv", sep="\t", index=False
    )
    if {"q_value_bh", "delta_effect"}.issubset(dependency.columns):
      dependency.loc[
        (dependency["q_value_bh"] < args.fdr) & (dependency["delta_effect"] < 0)
      ].to_csv(
        output_dir / "Table_main_significant_synthetic_lethal_dependencies.tsv",
        sep="\t",
        index=False,
      )

  expression = read_first(parent, [
    "expression_contrasts.tsv",
    "expanded_expression_compensation.tsv",
  ])
  if not expression.empty:
    sort_columns = [column for column in ["q_value_bh", "delta_expression"] if column in expression]
    if sort_columns:
      ascending = [True if column == "q_value_bh" else False for column in sort_columns]
      expression = expression.sort_values(sort_columns, ascending=ascending)
    expression.to_csv(
      output_dir / "Table_S_DepMap_expression_compensation.tsv", sep="\t", index=False
    )
    if {"q_value_bh", "delta_expression"}.issubset(expression.columns):
      expression.loc[
        (expression["q_value_bh"] < args.fdr) & (expression["delta_expression"] > 0)
      ].to_csv(
        output_dir / "Table_main_significant_expression_compensation.tsv",
        sep="\t",
        index=False,
      )

  expression_profiles = read_optional(parent / "expanded_expression_context_profiles.tsv")
  if not expression_profiles.empty:
    expression_profiles.to_csv(
      output_dir / "Table_S_expression_context_microniches.tsv", sep="\t", index=False
    )
  phenotype_profiles = read_optional(parent / "expanded_crispr_phenotype_profiles.tsv")
  if not phenotype_profiles.empty:
    phenotype_profiles.to_csv(
      output_dir / "Table_S_CRISPR_mutant_phenotype_profiles.tsv", sep="\t", index=False
    )

  microniche_columns = [
    column for column in scores.columns
    if column.startswith("microniche_") or column.startswith("functional_microniche_")
  ]
  microniche_key_columns = [
    column for column in [
      "cancer", "pair_id", "source_class", "lost_gene", "analysis_lost_gene", "target_gene"
    ] if column in scores.columns
  ]
  if microniche_columns:
    scores[microniche_key_columns + microniche_columns].to_csv(
      output_dir / "Table_S_human_functional_microniche_components.tsv",
      sep="\t",
      index=False,
    )

  tcga_columns = [
    column for column in [
      "cancer", "pair_id", "source_class", "lost_feature", "analysis_lost_gene",
      "target_gene", "tcga_homdel_n", "tcga_evaluable_n", "tcga_homdel_frequency",
      "component_tumor_event", "has_empirical_tcga",
    ] if column in scores.columns
  ]
  scores[tcga_columns].to_csv(
    output_dir / "Table_S_TCGA_homozygous_deletion_events.tsv", sep="\t", index=False
  )

  skipped = read_first(parent, [
    "skipped_complex_biomarkers.tsv",
    "expanded_skipped_complex_biomarkers.tsv",
  ])
  if not skipped.empty:
    skipped.to_csv(
      output_dir / "Table_S_complex_biomarkers_not_in_simple_CN_analysis.tsv",
      sep="\t",
      index=False,
    )

  depmap_path = resolve_path(args.depmap_only)
  if depmap_path.exists():
    depmap = pd.read_csv(depmap_path, sep="\t")
    comparison = depmap.merge(scores, on=["cancer", "pair_id"], suffixes=("_depmap", "_full"))
    comparison["delta_adjusted_rses"] = (
      comparison["coverage_adjusted_rses_full"]
      - comparison["coverage_adjusted_rses_depmap"]
    )
    comparison.sort_values(
      ["cancer", "delta_adjusted_rses"], ascending=[True, False]
    ).to_csv(
      output_dir / "Table_S_DepMap_vs_TCGA_DepMap_score_comparison.tsv",
      sep="\t",
      index=False,
    )

  if "source_class" in scores:
    class_summary = (
      scores.groupby("source_class", dropna=False)
        .agg(
          scored_rows=("pair_id", "size"),
          unique_directions=("pair_id", "nunique"),
          unique_targets=("target_gene", "nunique"),
          median_adjusted_rses=("coverage_adjusted_rses", "median"),
        )
        .reset_index()
        .sort_values("unique_directions", ascending=False)
    )
    class_summary.to_csv(output_dir / "Table_S_candidate_class_summary.tsv", sep="\t", index=False)

  summary_metrics = {
    "scored_rows": len(scores),
    "unique_pairs": scores["pair_id"].nunique(),
    "cancers": scores["cancer"].nunique(),
    "rows_with_dependency": int(bool_series(scores, "has_empirical_dependency").sum()),
    "rows_with_expression_compensation": int(
      (bool_series(scores, "has_empirical_expression")
       | bool_series(scores, "has_empirical_expression_compensation")).sum()
    ),
    "rows_with_expression_context": int(bool_series(scores, "has_empirical_expression_context").sum()),
    "rows_with_crispr_phenotype_profile": int(bool_series(scores, "has_empirical_phenotype_profile").sum()),
    "rows_with_tcga": int(bool_series(scores, "has_empirical_tcga").sum()),
  }
  if "source_class" in scores:
    summary_metrics["unique_nise_directions"] = scores.loc[
      scores["source_class"].astype(str).eq("NISE"), "pair_id"
    ].nunique()
  summary = pd.DataFrame({"metric": list(summary_metrics), "value": list(summary_metrics.values())})
  summary.to_csv(output_dir / "analysis_summary.tsv", sep="\t", index=False)
  print(summary.to_string(index=False), flush=True)
  print(f"Wrote manuscript-ready tables to {output_dir}", flush=True)


if __name__ == "__main__":
  main()
