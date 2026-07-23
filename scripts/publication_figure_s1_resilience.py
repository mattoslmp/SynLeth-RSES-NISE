#!/usr/bin/env python3
"""Data-adaptive, strict-layout-safe implementation of Supplementary Figure S1.

The original fixed 3x2 grid allocated the same vertical space to panels with very
different numbers and lengths of labels. Real-data runs can therefore produce
valid content that fails the strict layout audit through y-tick collisions,
especially in the mechanistic-class panel. This module preserves the same source
tables and plotted records but derives row heights from the rendered label burden.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rses_onco.publication import set_publication_style, wrap_label
from scripts.publication_audit_figures import (
  READABLE_STATUS,
  panel_label,
  read_tsv,
)

ROOT = Path(__file__).resolve().parents[1]


def _coverage_work(
  frame: pd.DataFrame,
  label_col: str,
  limit: int | None = None,
) -> pd.DataFrame:
  work = frame.copy()
  work["coverage_fraction"] = pd.to_numeric(
    work["coverage_fraction"],
    errors="coerce",
  )
  work = work.dropna(subset=["coverage_fraction"])
  if limit is not None:
    work = (
      work.sort_values("coverage_fraction", ascending=False)
      .head(limit)
    )
  return work.sort_values("coverage_fraction").reset_index(drop=True)


def _wrapped(values: pd.Series, width: int) -> list[str]:
  return [wrap_label(value, width) for value in values.astype(str)]


def _label_burden_height(
  labels: list[str],
  *,
  minimum: float,
) -> float:
  """Estimate panel height from tick count and wrapped-line count.

  The constants deliberately leave room for 9-point multiline labels, the panel
  title and the x-axis label at Matplotlib's default canvas DPI. This keeps labels
  readable at 100% zoom rather than merely shrinking them until the audit passes.
  """
  if not labels:
    return minimum
  line_units = sum(max(1, label.count("\n") + 1) for label in labels)
  return max(minimum, 1.15 + 0.22 * line_units + 0.07 * len(labels))


def _bar_coverage(
  axis: plt.Axes,
  frame: pd.DataFrame,
  label_col: str,
  *,
  limit: int | None = None,
  wrap_width: int = 44,
) -> tuple[pd.DataFrame, list[str]]:
  work = _coverage_work(frame, label_col, limit)
  labels = _wrapped(work[label_col], wrap_width)
  positions = np.arange(len(work), dtype=float)
  axis.barh(positions, work["coverage_fraction"])
  axis.set_yticks(positions, labels)
  axis.tick_params(axis="y", labelsize=9.0, pad=5)
  for label in axis.get_yticklabels():
    label.set_linespacing(1.12)
  axis.set_xlim(0, 1.05)
  axis.set_xlabel("Eligible hypotheses with evidence")
  axis.grid(axis="x", alpha=0.25)
  if len(work):
    axis.set_ylim(-0.65, len(work) - 0.35)
  for y, row in zip(positions, work.to_dict("records")):
    axis.text(
      min(float(row["coverage_fraction"]) + 0.015, 0.98),
      y,
      str(row.get("coverage_label", "")),
      va="center",
      fontsize=8.5,
    )
  return work, labels


def build_coverage_figure(
  domain: pd.DataFrame,
  source: pd.DataFrame,
  cancer: pd.DataFrame,
  classes: pd.DataFrame,
  reasons: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
  """Build Figure S1 and return the exact combined source-data table."""
  required_frames = {
    "domain": domain,
    "source": source,
    "cancer": cancer,
    "mechanistic class": classes,
    "missingness reason": reasons,
  }
  missing = [
    name
    for name, frame in required_frames.items()
    if frame.empty
  ]
  if missing:
    raise RuntimeError(
      "Expanded Figure S1 requires evidence-audit tables; missing/empty: "
      + ", ".join(missing)
    )

  cancer_display = cancer.copy()
  cancer_display["cancer_label"] = (
    cancer_display["domain_family"].astype(str)
    + " | "
    + cancer_display["cancer"].astype(str).replace({
      "colon": "Colorectal",
      "stomach": "Gastric",
      "lung": "Lung",
    })
  )

  class_display = classes.copy()
  class_display["class_label"] = (
    class_display["domain_family"].astype(str)
    + " | "
    + class_display["mechanistic_class"]
    .fillna("Unclassified")
    .astype(str)
    .str.replace("_", " ", regex=False)
  )

  missing_only = reasons.loc[
    ~reasons["evidence_state"].astype(str).isin({
      "observed_evidence",
      "negative_evidence",
      "neutral_evidence",
    })
  ].copy()
  missing_only["reason_label"] = (
    missing_only["evidence_state"]
    .map(READABLE_STATUS)
    .fillna(missing_only["evidence_state"].astype(str))
    + " | "
    + missing_only["absence_reason"]
    .astype(str)
    .str.replace("_", " ", regex=False)
  )
  missing_only = (
    missing_only.groupby("reason_label", as_index=False)
    .agg(rows=("rows", "sum"))
    .sort_values("rows")
    .tail(20)
    .reset_index(drop=True)
  )

  domain_labels = _wrapped(
    _coverage_work(domain, "domain_label")["domain_label"],
    44,
  )
  source_labels = _wrapped(
    _coverage_work(source, "evidence_source")["evidence_source"],
    44,
  )
  cancer_labels = _wrapped(
    _coverage_work(cancer_display, "cancer_label")["cancer_label"],
    46,
  )
  class_labels = _wrapped(
    _coverage_work(class_display, "class_label", 20)["class_label"],
    48,
  )
  reason_labels = _wrapped(missing_only["reason_label"], 72)

  top_height = max(
    _label_burden_height(domain_labels, minimum=4.8),
    _label_burden_height(source_labels, minimum=4.8),
  )
  middle_height = max(
    _label_burden_height(cancer_labels, minimum=5.5),
    _label_burden_height(class_labels, minimum=7.0),
  )
  bottom_height = _label_burden_height(reason_labels, minimum=7.0)

  set_publication_style()
  fig = plt.figure(
    figsize=(22.0, top_height + middle_height + bottom_height),
    constrained_layout=True,
  )
  grid = fig.add_gridspec(
    3,
    2,
    height_ratios=[top_height, middle_height, bottom_height],
  )
  axes = [
    fig.add_subplot(grid[0, 0]),
    fig.add_subplot(grid[0, 1]),
    fig.add_subplot(grid[1, 0]),
    fig.add_subplot(grid[1, 1]),
    fig.add_subplot(grid[2, :]),
  ]

  _bar_coverage(
    axes[0],
    domain,
    "domain_label",
    wrap_width=44,
  )
  axes[0].set_title("Coverage by evidence domain")
  panel_label(axes[0], "A")

  _bar_coverage(
    axes[1],
    source,
    "evidence_source",
    wrap_width=44,
  )
  axes[1].set_title("Coverage by data source")
  panel_label(axes[1], "B")

  _bar_coverage(
    axes[2],
    cancer_display,
    "cancer_label",
    wrap_width=46,
  )
  axes[2].set_title("Coverage by cancer context")
  panel_label(axes[2], "C")

  _bar_coverage(
    axes[3],
    class_display,
    "class_label",
    limit=20,
    wrap_width=48,
  )
  axes[3].set_title("Coverage by mechanistic class")
  panel_label(axes[3], "D")

  positions = np.arange(len(missing_only), dtype=float)
  axes[4].barh(positions, missing_only["rows"])
  axes[4].set_yticks(positions, reason_labels)
  axes[4].tick_params(axis="y", labelsize=9.0, pad=5)
  for label in axes[4].get_yticklabels():
    label.set_linespacing(1.12)
  if len(missing_only):
    axes[4].set_ylim(-0.65, len(missing_only) - 0.35)
  axes[4].set_xlabel("Candidate–cancer–domain records")
  axes[4].set_title("Documented reasons for unavailable evidence")
  axes[4].grid(axis="x", alpha=0.25)
  for y, value in zip(positions, missing_only["rows"]):
    axes[4].text(
      float(value),
      y,
      f" {int(value)}",
      va="center",
      fontsize=8.5,
    )
  panel_label(axes[4], "E")

  source_data = pd.concat([
    domain.assign(panel="A_domain"),
    source.assign(panel="B_source"),
    cancer_display.assign(panel="C_cancer"),
    class_display.assign(panel="D_class"),
    missing_only.assign(panel="E_missingness_reason"),
  ], ignore_index=True, sort=False)
  return fig, source_data


def expanded_coverage_s1(
  module: Any,
  item: dict,
  ranking: pd.DataFrame,
  output_dir: Path,
  source_dir: Path,
  strict: bool,
  input_path: Path,
):
  """Replace Figure S1 with a data-adaptive strict-layout-safe figure."""
  del ranking
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

  fig, source_data = build_coverage_figure(
    domain,
    source,
    cancer,
    classes,
    reasons,
  )
  return module.save_record(
    fig,
    item,
    output_dir,
    source_dir,
    source_data,
    [
      input_path,
      domain_path,
      source_path,
      cancer_path,
      class_path,
      reason_path,
    ],
    strict,
  )
