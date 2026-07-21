#!/usr/bin/env python3
"""Generate publication-ready benchmark or expanded RSES-Onco figures."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]
CANCER_LABELS = {
  "colon": "Colorectal",
  "stomach": "Gastric",
  "lung": "Lung",
}
FORMATS = ("png", "pdf", "svg")


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


def save_all(fig: plt.Figure, output_base: Path) -> None:
  output_base.parent.mkdir(parents=True, exist_ok=True)
  for extension in FORMATS:
    fig.savefig(
      output_base.with_suffix(f".{extension}"),
      dpi=600 if extension == "png" else None,
      bbox_inches="tight",
    )
  plt.close(fig)


def pair_label(frame: pd.DataFrame) -> pd.Series:
  lost = frame.get("analysis_lost_gene", frame.get("lost_gene", frame.get("lost_feature", "")))
  return lost.astype(str) + " → " + frame["target_gene"].astype(str)


def ranking_figure(scores: pd.DataFrame, output_dir: Path, top_n: int) -> None:
  cancers = [c for c in ["colon", "stomach", "lung"] if c in set(scores["cancer"])]
  if not cancers:
    return
  fig, axes = plt.subplots(1, len(cancers), figsize=(6.2 * len(cancers), 8), squeeze=False)
  for axis, cancer in zip(axes[0], cancers):
    subset = (
      scores.loc[scores["cancer"].eq(cancer)]
      .sort_values("coverage_adjusted_rses", ascending=False)
      .head(top_n)
      .copy()
      .sort_values("coverage_adjusted_rses", ascending=True)
    )
    labels = pair_label(subset)
    axis.barh(labels, subset["coverage_adjusted_rses"])
    axis.set_xlabel("Coverage-adjusted RSES-Onco")
    axis.set_title(CANCER_LABELS.get(cancer, cancer.title()))
    maximum = float(subset["coverage_adjusted_rses"].max()) if not subset.empty else 1.0
    axis.set_xlim(0, max(1.0, maximum * 1.10))
    axis.tick_params(axis="y", labelsize=9)
    axis.grid(axis="x", alpha=0.25)
  fig.suptitle("Cancer-specific RSES-Onco ranking", fontsize=16, fontweight="bold")
  fig.tight_layout()
  save_all(fig, output_dir / "Figure_empirical_1_cancer_specific_ranking")


def dependency_figure(dependency: pd.DataFrame, output_dir: Path) -> None:
  if dependency.empty or "delta_effect" not in dependency:
    return
  subset = dependency.copy()
  subset["label"] = pair_label(subset) + " | " + subset["cancer"].map(CANCER_LABELS).fillna(subset["cancer"])
  subset = subset.sort_values("delta_effect").head(30).sort_values("delta_effect", ascending=False)
  fig, axis = plt.subplots(figsize=(10, max(6, 0.34 * len(subset) + 2)))
  axis.barh(subset["label"], subset["delta_effect"])
  axis.axvline(0, linewidth=1)
  axis.set_xlabel("Δ CRISPR effect (loss − intact); more negative supports selectivity")
  axis.set_title("DepMap target-dependency contrasts")
  axis.tick_params(axis="y", labelsize=9)
  axis.grid(axis="x", alpha=0.25)
  fig.tight_layout()
  save_all(fig, output_dir / "Figure_empirical_2_dependency_contrasts")


def expression_figure(expression: pd.DataFrame, output_dir: Path) -> None:
  if expression.empty or "delta_expression" not in expression:
    return
  subset = expression.copy()
  subset["label"] = pair_label(subset) + " | " + subset["cancer"].map(CANCER_LABELS).fillna(subset["cancer"])
  subset = subset.reindex(subset["delta_expression"].abs().sort_values(ascending=False).index).head(30)
  subset = subset.sort_values("delta_expression")
  fig, axis = plt.subplots(figsize=(10, max(6, 0.34 * len(subset) + 2)))
  axis.barh(subset["label"], subset["delta_expression"])
  axis.axvline(0, linewidth=1)
  axis.set_xlabel("Δ target expression (loss − intact), log2(TPM+1)")
  axis.set_title("DepMap expression-compensation contrasts")
  axis.tick_params(axis="y", labelsize=9)
  axis.grid(axis="x", alpha=0.25)
  fig.tight_layout()
  save_all(fig, output_dir / "Figure_empirical_3_expression_compensation")


def event_dependency_figure(scores: pd.DataFrame, output_dir: Path) -> None:
  required = {"tcga_homdel_frequency", "component_selectivity", "coverage_adjusted_rses"}
  if not required.issubset(scores.columns):
    return
  subset = scores.dropna(subset=["tcga_homdel_frequency", "component_selectivity"]).copy()
  if subset.empty:
    return
  fig, axis = plt.subplots(figsize=(9, 7))
  for cancer, group in subset.groupby("cancer"):
    axis.scatter(
      group["tcga_homdel_frequency"],
      group["component_selectivity"],
      s=50 + 250 * group["coverage_adjusted_rses"].fillna(0),
      alpha=0.75,
      label=CANCER_LABELS.get(cancer, cancer.title()),
    )
  top = subset.sort_values("coverage_adjusted_rses", ascending=False).head(12)
  for _, row in top.iterrows():
    axis.annotate(
      f"{row['analysis_lost_gene']}→{row['target_gene']}",
      (row["tcga_homdel_frequency"], row["component_selectivity"]),
      xytext=(4, 4),
      textcoords="offset points",
      fontsize=8,
    )
  axis.set_xlabel("TCGA homozygous-deletion frequency")
  axis.set_ylabel("DepMap loss-selectivity component")
  axis.set_title("Tumor-event prevalence versus target selectivity")
  axis.grid(alpha=0.25)
  axis.legend(frameon=False)
  fig.tight_layout()
  save_all(fig, output_dir / "Figure_empirical_4_event_selectivity_integration")


def heatmap_figure(scores: pd.DataFrame, output_dir: Path) -> None:
  pivot = scores.pivot_table(
    index=["pair_id", "lost_feature", "target_gene"],
    columns="cancer",
    values="coverage_adjusted_rses",
    aggfunc="first",
  )
  if pivot.empty:
    return
  pivot["maximum"] = pivot.max(axis=1)
  pivot = pivot.sort_values("maximum", ascending=False).head(25).drop(columns="maximum")
  columns = [c for c in ["colon", "stomach", "lung"] if c in pivot.columns]
  pivot = pivot[columns]
  labels = [f"{lost} → {target}" for _, lost, target in pivot.index]
  fig, axis = plt.subplots(figsize=(8, max(7, 0.34 * len(pivot) + 2)))
  image = axis.imshow(pivot.to_numpy(), aspect="auto", vmin=0, vmax=1)
  axis.set_xticks(np.arange(len(columns)), [CANCER_LABELS.get(c, c.title()) for c in columns])
  axis.set_yticks(np.arange(len(labels)), labels)
  axis.tick_params(axis="y", labelsize=9)
  axis.set_title("Coverage-adjusted RSES-Onco across cancer types")
  colorbar = fig.colorbar(image, ax=axis)
  colorbar.set_label("Coverage-adjusted score")
  for row_index in range(pivot.shape[0]):
    for column_index in range(pivot.shape[1]):
      value = pivot.iloc[row_index, column_index]
      if pd.notna(value):
        axis.text(column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=8)
  fig.tight_layout()
  save_all(fig, output_dir / "Figure_empirical_5_cross_cancer_heatmap")


def microniche_figure(scores: pd.DataFrame, output_dir: Path, top_n: int = 30) -> None:
  domains = [
    "microniche_expression_context",
    "microniche_localization",
    "microniche_biochemical_structural",
    "microniche_genetic_phenotype",
    "microniche_interaction_network",
    "microniche_regulatory_network",
  ]
  available = [column for column in domains if column in scores.columns]
  if not available:
    return
  pair_level = (
    scores.sort_values("functional_microniche_adjusted", ascending=False)
      .drop_duplicates("pair_id")
      .head(top_n)
      .copy()
  )
  if pair_level.empty:
    return
  matrix = pair_level[available].astype(float)
  labels = pair_label(pair_level).tolist()
  domain_labels = [
    column.removeprefix("microniche_").replace("_", " ").title()
    for column in available
  ]
  fig, axis = plt.subplots(figsize=(11, max(8, 0.34 * len(pair_level) + 2)))
  image = axis.imshow(matrix.to_numpy(), aspect="auto", vmin=0, vmax=1)
  axis.set_xticks(np.arange(len(available)), domain_labels, rotation=35, ha="right")
  axis.set_yticks(np.arange(len(labels)), labels)
  axis.tick_params(axis="y", labelsize=8)
  axis.set_title("Human functional-microniche evidence domains")
  colorbar = fig.colorbar(image, ax=axis)
  colorbar.set_label("Normalized divergence / specialization evidence")
  for row_index in range(matrix.shape[0]):
    for column_index in range(matrix.shape[1]):
      value = matrix.iloc[row_index, column_index]
      if pd.notna(value):
        axis.text(column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=7)
  fig.tight_layout()
  save_all(fig, output_dir / "Figure_empirical_6_functional_microniche_domains")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--input",
    default="results/empirical_26Q1/full/empirical_rses_onco_by_cancer.tsv",
  )
  parser.add_argument("--output-dir", default="figures/empirical_26Q1")
  parser.add_argument("--top-n", type=int, default=15)
  args = parser.parse_args()

  input_path = resolve_path(args.input)
  output_dir = resolve_path(args.output_dir)
  scores = pd.read_csv(input_path, sep="\t")
  dependency = read_first(input_path.parent, [
    "dependency_contrasts.tsv",
    "expanded_dependency_contrasts.tsv",
  ])
  expression = read_first(input_path.parent, [
    "expression_contrasts.tsv",
    "expanded_expression_compensation.tsv",
  ])

  ranking_figure(scores, output_dir, args.top_n)
  dependency_figure(dependency, output_dir)
  expression_figure(expression, output_dir)
  event_dependency_figure(scores, output_dir)
  heatmap_figure(scores, output_dir)
  microniche_figure(scores, output_dir)
  print(f"Wrote empirical figures to {output_dir}", flush=True)


if __name__ == "__main__":
  main()
