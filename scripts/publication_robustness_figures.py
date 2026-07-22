#!/usr/bin/env python3
"""Comparability-aware robustness figures for RSES-Onco."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rses_onco.publication import FigureRecord, set_publication_style, wrap_label
from scripts.publication_audit_figures import read_tsv, save_record


def figure_s35(output_root: Path, strict: bool) -> FigureRecord:
  path = output_root / "tables" / "robustness" / "ranking_correlations.tsv"
  frame = read_tsv(path)
  if frame.empty:
    raise RuntimeError(f"Missing rank-correlation table: {path}")
  selected = frame.loc[
    frame["scenario"].astype(str).str.startswith("leave_out_")
    & frame["comparison_scope"].astype(str).eq("cancer_comparability_group")
  ].copy()
  if selected.empty:
    raise RuntimeError("No comparability-aware leave-one-domain-out correlations")
  preferred = "gene_pair"
  groups = selected["score_comparability_group"].dropna().astype(str).unique()
  comparability_group = preferred if preferred in groups else sorted(groups)[0]
  selected = selected.loc[
    selected["score_comparability_group"].astype(str).eq(comparability_group)
  ].copy()
  matrix = selected.pivot_table(
    index="scenario",
    columns="cancer",
    values="spearman_rho",
    aggfunc="first",
  )
  matrix = matrix[[
    column for column in ("colon", "stomach", "lung")
    if column in matrix
  ]]
  set_publication_style()
  fig, axis = plt.subplots(
    figsize=(11.5, max(6.5, 0.55 * len(matrix))),
    constrained_layout=True,
  )
  image = axis.imshow(
    matrix.to_numpy(dtype=float),
    aspect="auto",
    vmin=-1,
    vmax=1,
  )
  axis.set_xticks(
    np.arange(len(matrix.columns)),
    [str(value).title() for value in matrix.columns],
  )
  axis.set_yticks(
    np.arange(len(matrix)),
    [
      wrap_label(
        value.replace("leave_out_", "Without ").replace("_", " "),
        32,
      )
      for value in matrix.index
    ],
  )
  for y in range(len(matrix)):
    for x in range(len(matrix.columns)):
      value = matrix.iloc[y, x]
      axis.text(
        x,
        y,
        f"{value:.3f}" if pd.notna(value) else "NA",
        ha="center",
        va="center",
        fontsize=9,
      )
  colorbar = fig.colorbar(image, ax=axis, fraction=0.035, pad=0.03)
  colorbar.set_label("Spearman rank correlation")
  axis.set_xlabel(f"Cancer context | comparability group: {comparability_group}")
  return save_record(
    fig,
    "Figure_S35",
    "Figure_S35_leave_one_domain_out",
    "Leave-one-domain-out ranking robustness",
    "Baseline-versus-perturbed rank correlation within cancer and a common score-comparability group.",
    output_root,
    selected,
    [path],
    "scripts/make_audit_supplementary_figures.py",
    strict,
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
    evaluated_contexts=("score_comparability_group", "size"),
  ).sort_values("median_jaccard")
  set_publication_style()
  fig, axis = plt.subplots(
    figsize=(13.5, max(7.0, 0.35 * len(summary))),
    constrained_layout=True,
  )
  positions = np.arange(len(summary), dtype=float)
  lower = summary["median_jaccard"] - summary["minimum_jaccard"]
  upper = summary["maximum_jaccard"] - summary["median_jaccard"]
  axis.errorbar(
    summary["median_jaccard"],
    positions,
    xerr=np.vstack([lower, upper]),
    fmt="o",
    capsize=3,
  )
  axis.set_yticks(
    positions,
    [wrap_label(value.replace("_", " "), 42) for value in summary["scenario"]],
  )
  axis.set_xlim(0, 1.02)
  axis.set_xlabel(
    "Top-k Jaccard stability across cancer × score-comparability contexts"
  )
  axis.grid(axis="x", alpha=0.25)
  return save_record(
    fig,
    "Figure_S36",
    "Figure_S36_top_k_stability",
    "Stability of high-priority candidates",
    "Median and context-specific range of top-k overlap under leave-one-domain-out and controlled weight perturbations; rankings are never mixed across score-comparability groups.",
    output_root,
    frame,
    [path],
    "scripts/make_audit_supplementary_figures.py",
    strict,
  )
