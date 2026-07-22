#!/usr/bin/env python3
"""Source-backed coverage, missingness and robustness publication figures."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd

from rses_onco.publication import (
  FigureRecord,
  figure_record,
  save_figure_triplet,
  set_publication_style,
  wrap_label,
  write_source_data,
)

ROOT = Path(__file__).resolve().parents[1]

READABLE_STATUS = {
  "observed_evidence": "Observed evidence",
  "negative_evidence": "Observed negative evidence",
  "neutral_evidence": "Observed neutral evidence",
  "missing": "Evidence unavailable",
  "not_eligible": "Not eligible",
  "technical_failure": "Technical/source failure",
  "insufficient_sample": "Insufficient sample",
}


def read_tsv(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    return pd.DataFrame()
  return pd.read_csv(path, sep="\t", low_memory=False)


def panel_label(axis: plt.Axes, label: str) -> None:
  axis.text(
    0.01, 0.99, label, transform=axis.transAxes, ha="left", va="top",
    fontsize=13, fontweight="bold",
    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5}, zorder=20,
  )


def save_record(
  fig: plt.Figure,
  figure_id: str,
  file_name: str,
  title: str,
  caption: str,
  output_root: Path,
  source: pd.DataFrame,
  inputs: list[Path],
  script: str,
  strict: bool,
  category: str = "supplementary",
) -> FigureRecord:
  output_dir = output_root / "figures" / ("main" if category == "main" else "supplementary")
  source_dir = output_root / "source_data" / "figures" / ("main" if category == "main" else "supplementary")
  base = output_dir / file_name
  source_path = write_source_data(source, source_dir / f"{file_name}_source_data.tsv")
  audit = save_figure_triplet(fig, base, figure_id, strict_layout=strict)
  return figure_record(
    figure_id=figure_id,
    category=category,
    title=title,
    caption=caption,
    base_path=base,
    source_data_path=source_path,
    input_paths=inputs,
    audit=audit,
    script=script,
  )


def _bar_coverage(axis: plt.Axes, frame: pd.DataFrame, label_col: str, limit: int | None = None) -> None:
  work = frame.copy()
  if limit is not None:
    work = work.sort_values("coverage_fraction", ascending=False).head(limit)
  work = work.sort_values("coverage_fraction")
  positions = np.arange(len(work), dtype=float)
  labels = [wrap_label(value, 28) for value in work[label_col].astype(str)]
  axis.barh(positions, work["coverage_fraction"])
  axis.set_yticks(positions, labels)
  axis.set_xlim(0, 1.05)
  axis.set_xlabel("Eligible hypotheses with evidence")
  axis.grid(axis="x", alpha=0.25)
  for y, row in zip(positions, work.to_dict("records")):
    axis.text(
      min(float(row["coverage_fraction"]) + 0.015, 0.98), y,
      str(row["coverage_label"]), va="center", fontsize=8.5,
    )


def expanded_coverage_s1(
  module: Any,
  item: dict,
  ranking: pd.DataFrame,
  output_dir: Path,
  source_dir: Path,
  strict: bool,
  input_path: Path,
):
  """Replace S1 with domain/source/cancer/class/missingness coverage panels."""
  qc = ROOT / "article_outputs" / "tables" / "qc"
  domain_path = qc / "coverage_by_domain.tsv"
  source_path = qc / "coverage_by_source.tsv"
  cancer_path = qc / "coverage_by_cancer.tsv"
  class_path = qc / "coverage_by_mechanistic_class.tsv"
  reason_path = qc / "missingness_reasons.tsv"
  domain = read_tsv(domain_path)
  source = read_tsv(source_path)
  cancer = read_tsv(cancer_path)
  classes = read_tsv(class_path)
  reasons = read_tsv(reason_path)
  required_frames = {
    "domain": domain, "source": source, "cancer": cancer,
    "mechanistic class": classes, "missingness reason": reasons,
  }
  missing = [name for name, frame in required_frames.items() if frame.empty]
  if missing:
    raise RuntimeError(
      "Expanded Figure S1 requires evidence-audit tables; missing/empty: " + ", ".join(missing)
    )

  set_publication_style()
  fig = plt.figure(figsize=(20.0, 17.0), constrained_layout=True)
  grid = fig.add_gridspec(3, 2, height_ratios=[1.35, 1.0, 1.15])
  axes = [
    fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
    fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1]),
    fig.add_subplot(grid[2, :]),
  ]

  _bar_coverage(axes[0], domain, "domain_label")
  axes[0].set_title("Coverage by evidence domain")
  panel_label(axes[0], "A")

  _bar_coverage(axes[1], source, "evidence_source")
  axes[1].set_title("Coverage by data source")
  panel_label(axes[1], "B")

  cancer_display = cancer.copy()
  cancer_display["cancer_label"] = cancer_display["domain_family"].astype(str) + " | " + cancer_display["cancer"].astype(str).replace({"colon": "Colorectal", "stomach": "Gastric", "lung": "Lung"})
  _bar_coverage(axes[2], cancer_display, "cancer_label")
  axes[2].set_title("Coverage by cancer context")
  panel_label(axes[2], "C")

  class_display = classes.copy()
  class_display["class_label"] = class_display["domain_family"].astype(str) + " | " + class_display["mechanistic_class"].fillna("Unclassified").astype(str).str.replace("_", " ", regex=False)
  _bar_coverage(axes[3], class_display, "class_label", limit=20)
  axes[3].set_title("Coverage by mechanistic class")
  panel_label(axes[3], "D")

  missing_only = reasons.loc[~reasons["evidence_state"].astype(str).isin({"observed_evidence", "negative_evidence", "neutral_evidence"})].copy()
  missing_only["reason_label"] = missing_only["evidence_state"].map(READABLE_STATUS).fillna(missing_only["evidence_state"].astype(str)) + " | " + missing_only["absence_reason"].astype(str).str.replace("_", " ", regex=False)
  missing_only = missing_only.groupby("reason_label", as_index=False).agg(rows=("rows", "sum")).sort_values("rows").tail(20)
  positions = np.arange(len(missing_only), dtype=float)
  axes[4].barh(positions, missing_only["rows"])
  axes[4].set_yticks(positions, [wrap_label(value, 52) for value in missing_only["reason_label"]])
  axes[4].set_xlabel("Candidate–cancer–domain records")
  axes[4].set_title("Documented reasons for unavailable evidence")
  axes[4].grid(axis="x", alpha=0.25)
  for y, value in zip(positions, missing_only["rows"]):
    axes[4].text(float(value), y, f" {int(value)}", va="center", fontsize=8.5)
  panel_label(axes[4], "E")

  source_data = pd.concat([
    domain.assign(panel="A_domain"),
    source.assign(panel="B_source"),
    cancer_display.assign(panel="C_cancer"),
    class_display.assign(panel="D_class"),
    missing_only.assign(panel="E_missingness_reason"),
  ], ignore_index=True, sort=False)
  return module.save_record(
    fig, item, output_dir, source_dir, source_data,
    [input_path, domain_path, source_path, cancer_path, class_path, reason_path], strict,
  )


def figure_s33(output_root: Path, strict: bool) -> FigureRecord:
  path = output_root / "tables" / "qc" / "candidate_domain_evidence_audit.tsv"
  audit = read_tsv(path)
  if audit.empty:
    raise RuntimeError(f"Missing audit table: {path}")
  status_order = [
    "observed_evidence", "negative_evidence", "neutral_evidence",
    "missing", "insufficient_sample", "technical_failure", "not_eligible",
  ]
  matrix = audit.pivot_table(
    index="hypothesis_direction", columns="domain_label", values="evidence_state",
    aggfunc="first",
  )
  top_pairs = (
    audit.groupby("hypothesis_direction", as_index=False)
      .agg(observed=("evidence_present", "sum"), records=("domain", "size"))
      .sort_values(["observed", "records"], ascending=False)
      .head(80)["hypothesis_direction"]
  )
  matrix = matrix.reindex(top_pairs)
  mapping = {status: index for index, status in enumerate(status_order)}
  numeric = matrix.applymap(lambda value: mapping.get(str(value), np.nan)).to_numpy(dtype=float)
  set_publication_style()
  fig, axis = plt.subplots(figsize=(18.0, max(10.0, 0.32 * len(matrix))), constrained_layout=True)
  cmap = ListedColormap(["#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7", "#ef8a62", "#bdbdbd"])
  image = axis.imshow(numeric, aspect="auto", cmap=cmap, vmin=-0.5, vmax=len(status_order) - 0.5)
  axis.set_xticks(np.arange(len(matrix.columns)), [wrap_label(value, 20) for value in matrix.columns], rotation=28, ha="right")
  axis.set_yticks(np.arange(len(matrix)), [wrap_label(value, 38) for value in matrix.index])
  axis.set_xlabel("Evidence domain")
  axis.set_ylabel("Directed hypothesis")
  colorbar = fig.colorbar(image, ax=axis, fraction=0.025, pad=0.02, ticks=np.arange(len(status_order)))
  colorbar.ax.set_yticklabels([READABLE_STATUS[status] for status in status_order])
  colorbar.set_label("Evidence state")
  return save_record(
    fig, "Figure_S33", "Figure_S33_complete_missingness_matrix",
    "Complete evidence-state and missingness matrix",
    "Evidence states across candidate–cancer–domain combinations; unavailable, non-eligible and negative evidence are distinct.",
    output_root, audit, [path], "scripts/make_audit_supplementary_figures.py", strict,
  )


def figure_s34(output_root: Path, strict: bool) -> FigureRecord:
  path = output_root / "tables" / "robustness" / "raw_vs_coverage_adjusted.tsv"
  frame = read_tsv(path)
  if frame.empty:
    raise RuntimeError(f"Missing robustness table: {path}")
  set_publication_style()
  fig, axes = plt.subplots(1, 3, figsize=(18.0, 6.2), constrained_layout=True)
  cancers = ["colon", "stomach", "lung"]
  labels = {"colon": "Colorectal", "stomach": "Gastric", "lung": "Lung"}
  for index, (axis, cancer) in enumerate(zip(axes, cancers)):
    group = frame.loc[frame["cancer"].astype(str).eq(cancer)].copy()
    axis.scatter(group["recomputed_raw_score"], group["recomputed_adjusted_score"], alpha=0.45, s=18)
    axis.plot([0, 1], [0, 1], linestyle="--", linewidth=1, color="black")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("Observed-domain score")
    axis.set_ylabel("Coverage-adjusted score")
    axis.set_title(labels[cancer])
    axis.grid(alpha=0.2)
    panel_label(axis, chr(ord("A") + index))
  return save_record(
    fig, "Figure_S34", "Figure_S34_raw_vs_coverage_adjusted",
    "Observed-domain and coverage-adjusted RSES-Onco scores",
    "Coverage adjustment penalizes incomplete evidence without imputing missing components as zero.",
    output_root, frame, [path], "scripts/make_audit_supplementary_figures.py", strict,
  )


def figure_s35(output_root: Path, strict: bool) -> FigureRecord:
  path = output_root / "tables" / "robustness" / "ranking_correlations.tsv"
  frame = read_tsv(path)
  if frame.empty:
    raise RuntimeError(f"Missing rank-correlation table: {path}")
  selected = frame.loc[frame["scenario"].astype(str).str.startswith("leave_out_")].copy()
  matrix = selected.pivot_table(index="scenario", columns="cancer", values="spearman_rho", aggfunc="first")
  matrix = matrix[[column for column in ("colon", "stomach", "lung") if column in matrix]]
  set_publication_style()
  fig, axis = plt.subplots(figsize=(11.5, max(6.5, 0.55 * len(matrix))), constrained_layout=True)
  image = axis.imshow(matrix.to_numpy(dtype=float), aspect="auto", vmin=-1, vmax=1)
  axis.set_xticks(np.arange(len(matrix.columns)), [str(value).title() for value in matrix.columns])
  axis.set_yticks(np.arange(len(matrix)), [wrap_label(value.replace("leave_out_", "Without ").replace("_", " "), 32) for value in matrix.index])
  for y in range(len(matrix)):
    for x in range(len(matrix.columns)):
      value = matrix.iloc[y, x]
      axis.text(x, y, f"{value:.3f}" if pd.notna(value) else "NA", ha="center", va="center", fontsize=9)
  colorbar = fig.colorbar(image, ax=axis, fraction=0.035, pad=0.03)
  colorbar.set_label("Spearman rank correlation")
  return save_record(
    fig, "Figure_S35", "Figure_S35_leave_one_domain_out",
    "Leave-one-domain-out ranking robustness",
    "Correlation of the baseline ranking with rankings recomputed after excluding one domain at a time.",
    output_root, selected, [path], "scripts/make_audit_supplementary_figures.py", strict,
  )


def figure_s36(output_root: Path, strict: bool) -> FigureRecord:
  path = output_root / "tables" / "robustness" / "top_k_stability.tsv"
  frame = read_tsv(path)
  if frame.empty:
    raise RuntimeError(f"Missing stability table: {path}")
  summary = frame.groupby("scenario", as_index=False).agg(
    median_jaccard=("jaccard", "median"),
    minimum_jaccard=("jaccard", "min"),
    maximum_jaccard=("jaccard", "max"),
  ).sort_values("median_jaccard")
  set_publication_style()
  fig, axis = plt.subplots(figsize=(13.5, max(7.0, 0.35 * len(summary))), constrained_layout=True)
  positions = np.arange(len(summary), dtype=float)
  lower = summary["median_jaccard"] - summary["minimum_jaccard"]
  upper = summary["maximum_jaccard"] - summary["median_jaccard"]
  axis.errorbar(summary["median_jaccard"], positions, xerr=np.vstack([lower, upper]), fmt="o", capsize=3)
  axis.set_yticks(positions, [wrap_label(value.replace("_", " "), 42) for value in summary["scenario"]])
  axis.set_xlim(0, 1.02)
  axis.set_xlabel("Top-k Jaccard stability across cancers")
  axis.grid(axis="x", alpha=0.25)
  return save_record(
    fig, "Figure_S36", "Figure_S36_top_k_stability",
    "Stability of high-priority candidates",
    "Median and cancer-specific range of top-k overlap under leave-one-domain-out and controlled weight perturbations.",
    output_root, frame, [path], "scripts/make_audit_supplementary_figures.py", strict,
  )


def figure_s37(output_root: Path, strict: bool) -> FigureRecord:
  path = output_root / "tables" / "qc" / "evidence_overlap_summary.tsv"
  frame = read_tsv(path)
  if frame.empty:
    raise RuntimeError(f"Missing overlap table: {path}")
  summary = frame.groupby(["domain", "overlap_class"], as_index=False).agg(
    evidence_units=("deduplication_key", "nunique"),
    evidence_rows=("evidence_rows", "sum"),
  )
  pivot = summary.pivot_table(index="domain", columns="overlap_class", values="evidence_units", fill_value=0, aggfunc="sum")
  set_publication_style()
  fig, axis = plt.subplots(figsize=(13.5, max(6.5, 0.55 * len(pivot))), constrained_layout=True)
  left = np.zeros(len(pivot))
  positions = np.arange(len(pivot), dtype=float)
  for column in pivot.columns:
    values = pivot[column].to_numpy(dtype=float)
    axis.barh(positions, values, left=left, label=column.replace("_", " "))
    left += values
  axis.set_yticks(positions, [wrap_label(value.replace("_", " "), 30) for value in pivot.index])
  axis.set_xlabel("Independent or overlap-controlled evidence units")
  axis.legend(frameon=False, loc="lower right")
  axis.grid(axis="x", alpha=0.2)
  return save_record(
    fig, "Figure_S37", "Figure_S37_evidence_overlap_control",
    "Evidence overlap and deduplication control",
    "Evidence units are classified as independent or overlapping; repeated aggregations cannot receive more than one total evidence unit.",
    output_root, frame, [path], "scripts/make_audit_supplementary_figures.py", strict,
  )


def figure_s38(output_root: Path, strict: bool) -> FigureRecord:
  path = output_root / "tables" / "qc" / "evidence_category_assignments.tsv"
  frame = read_tsv(path)
  if frame.empty:
    raise RuntimeError(f"Missing category table: {path}")
  counts = frame.groupby(["cancer", "highest_evidence_category"], as_index=False).agg(
    hypotheses=("candidate_id", "nunique")
  )
  pivot = counts.pivot_table(index="highest_evidence_category", columns="cancer", values="hypotheses", fill_value=0, aggfunc="sum")
  set_publication_style()
  fig, axis = plt.subplots(figsize=(13.0, max(6.5, 0.65 * len(pivot))), constrained_layout=True)
  positions = np.arange(len(pivot), dtype=float)
  width = 0.24
  for index, cancer in enumerate([column for column in ("colon", "stomach", "lung") if column in pivot]):
    axis.barh(positions + (index - 1) * width, pivot[cancer], height=width, label=cancer.title())
  axis.set_yticks(positions, [wrap_label(value.replace("_", " "), 36) for value in pivot.index])
  axis.set_xlabel("Unique hypotheses")
  axis.legend(frameon=False)
  axis.grid(axis="x", alpha=0.2)
  return save_record(
    fig, "Figure_S38", "Figure_S38_evidence_categories",
    "Explicit evidence-category assignments",
    "Candidate-universe, prioritized, microniche-supported, conditional-dependency and FDR-supported categories are kept distinct.",
    output_root, frame, [path], "scripts/make_audit_supplementary_figures.py", strict,
  )


AUDIT_FIGURE_BUILDERS = [figure_s33, figure_s34, figure_s35, figure_s36, figure_s37, figure_s38]
