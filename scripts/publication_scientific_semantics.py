#!/usr/bin/env python3
"""Scientific display helpers for main and supplementary publication figures.

These helpers preserve composite biomarkers as features rather than inventing gene
symbols, report the exact single FDR-supported all-target discovery when only one
passes, and distinguish compound-resolved pharmacology rows from target-only
evidence.
"""
from __future__ import annotations

from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rses_onco.publication import set_publication_style, wrap_label

NA_STRINGS = {"", "na", "nan", "none", "null", "<na>"}
CANCER_LABELS = {"colon": "Colorectal", "stomach": "Gastric", "lung": "Lung"}


def present_text(value: object) -> str:
  """Return a normalized non-missing scalar string or an empty string."""
  if value is None:
    return ""
  try:
    if pd.isna(value):
      return ""
  except (TypeError, ValueError):
    pass
  text = str(value).strip()
  return "" if text.casefold() in NA_STRINGS else text


def first_present(row: pd.Series | dict[str, Any], columns: Iterable[str]) -> str:
  """Resolve the first non-missing value in a documented precedence order."""
  for column in columns:
    try:
      value = row.get(column)
    except AttributeError:
      value = None
    text = present_text(value)
    if text:
      return text
  return ""


def resolved_pair_parts(row: pd.Series | dict[str, Any]) -> tuple[str, str, str]:
  """Resolve loss context, target and loss-entity type without inventing genes."""
  lost_gene = first_present(row, ("analysis_lost_gene", "lost_gene"))
  lost_feature = first_present(row, ("lost_feature",))
  target = first_present(
    row,
    ("analysis_target_gene", "target_gene", "target_feature"),
  ) or "Unresolved target"
  if lost_gene:
    return lost_gene, target, "gene"
  if lost_feature:
    return lost_feature, target, "composite_feature"
  return "Unresolved loss context", target, "unresolved"


def resolved_pair_label(
  row: pd.Series | dict[str, Any],
  width: int = 28,
  include_type_prefix: bool = False,
) -> str:
  """Build a traceable display label for gene loss or composite biomarkers."""
  lost, target, entity_type = resolved_pair_parts(row)
  if entity_type == "gene":
    label = f"{lost} → {target}"
  elif entity_type == "composite_feature":
    prefix = "Feature: " if include_type_prefix else ""
    label = f"{prefix}{lost} ⇒ {target}"
  else:
    label = f"{lost} ⇒ {target}"
  return wrap_label(label, width)


def pair_labels(frame: pd.DataFrame, width: int = 28) -> list[str]:
  return [resolved_pair_label(row, width) for _, row in frame.iterrows()]


def add_display_pair_columns(frame: pd.DataFrame) -> pd.DataFrame:
  output = frame.copy()
  if output.empty:
    output["display_lost_label"] = pd.Series(dtype=object)
    output["display_target_label"] = pd.Series(dtype=object)
    output["lost_entity_type"] = pd.Series(dtype=object)
    output["display_pair_label"] = pd.Series(dtype=object)
    return output
  parts = [resolved_pair_parts(row) for _, row in output.iterrows()]
  output["display_lost_label"] = [item[0] for item in parts]
  output["display_target_label"] = [item[1] for item in parts]
  output["lost_entity_type"] = [item[2] for item in parts]
  output["display_pair_label"] = [
    resolved_pair_label(row, width=80)
    for _, row in output.iterrows()
  ]
  return output


def compound_resolution(row: pd.Series | dict[str, Any]) -> tuple[str, str]:
  """Return display compound and a resolution class."""
  name = first_present(row, ("drug_name",))
  identifier = first_present(row, ("drug_id",))
  key = first_present(row, ("drug_key",))
  if name:
    display = f"{name} ({identifier})" if identifier and identifier != name else name
    return display, "named_compound"
  if identifier:
    return identifier, "compound_identifier"
  if key and key.upper() != "TARGET_ONLY":
    return key, "compound_identifier"
  return "Target-level evidence only", "target_only"


def figure_6(
  module: Any,
  item: dict,
  ranking: pd.DataFrame,
  discovery: pd.DataFrame,
  output_dir,
  source_dir,
  strict: bool,
  input_paths,
):
  """Render class distributions and exact FDR-supported discovery evidence."""
  set_publication_style()
  fig = plt.figure(figsize=(17.0, 8.8), constrained_layout=True)
  grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.30])
  axis_class = fig.add_subplot(grid[0, 0])
  axis_discovery = fig.add_subplot(grid[0, 1])

  class_summary = (
    ranking.assign(
      source_class=ranking.get("source_class", "unclassified")
        .fillna("unclassified")
        .astype(str)
    )
      .groupby("source_class", as_index=False)
      .agg(
        unique_directions=("pair_id", "nunique"),
        maximum_score=("coverage_adjusted_rses", "max"),
        median_score=("coverage_adjusted_rses", "median"),
      )
      .sort_values("median_score")
  )
  positions = np.arange(len(class_summary), dtype=float)
  axis_class.barh(
    positions,
    class_summary["median_score"],
    label="Median",
  )
  axis_class.scatter(
    class_summary["maximum_score"],
    positions,
    marker="D",
    s=54,
    label="Maximum",
    zorder=5,
  )
  axis_class.set_yticks(
    positions,
    [wrap_label(value, 24) for value in class_summary["source_class"]],
  )
  for y, row in zip(positions, class_summary.to_dict("records")):
    axis_class.annotate(
      f"n={int(row['unique_directions'])}",
      (float(row["maximum_score"]), y),
      xytext=(7, 0),
      textcoords="offset points",
      va="center",
      fontsize=8.5,
    )
  axis_class.set_xlabel("Coverage-adjusted RSES-Onco")
  axis_class.set_title("Class score distributions")
  axis_class.grid(axis="x", alpha=0.25)
  axis_class.legend(loc="lower right", frameon=False)
  module.panel_label(axis_class, "A")

  top_discovery = pd.DataFrame()
  q_column = (
    "q_value_bh_within_loss_cancer"
    if "q_value_bh_within_loss_cancer" in discovery
    else "q_value_bh"
  )
  if discovery.empty or "delta_effect" not in discovery or q_column not in discovery:
    module.placeholder(
      axis_discovery,
      "FDR-supported all-target discoveries",
      "No eligible all-target discovery rows were available.",
    )
  else:
    top_discovery = (
      discovery.sort_values([q_column, "delta_effect"], ascending=[True, True])
        .head(25)
        .copy()
    )
    top_discovery = add_display_pair_columns(top_discovery)
    supported_count = len(top_discovery)
    if supported_count <= 8:
      axis_discovery.set_axis_off()
      rows = []
      for row in top_discovery.to_dict("records"):
        p_value = pd.to_numeric(pd.Series([row.get("p_value")]), errors="coerce").iloc[0]
        q_value = pd.to_numeric(pd.Series([row.get(q_column)]), errors="coerce").iloc[0]
        delta = pd.to_numeric(pd.Series([row.get("delta_effect")]), errors="coerce").iloc[0]
        rows.append([
          present_text(row.get("display_pair_label")),
          CANCER_LABELS.get(str(row.get("cancer")), str(row.get("cancer"))),
          f"{delta:.3f}" if pd.notna(delta) else "NA",
          f"{p_value:.2e}" if pd.notna(p_value) else "NA",
          f"{q_value:.4f}" if pd.notna(q_value) else "NA",
          present_text(row.get("n_loss")) or "NA",
          present_text(row.get("n_intact")) or "NA",
        ])
      table = axis_discovery.table(
        cellText=rows,
        colLabels=["Pair/context", "Cancer", "Δ effect", "P", "FDR q", "n loss", "n intact"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.27, 0.12, 0.11, 0.13, 0.11, 0.10, 0.10],
      )
      table.auto_set_font_size(False)
      table.set_fontsize(8.2)
      table.scale(1.0, 1.65)
      axis_discovery.set_title(
        f"FDR-supported all-target discoveries (n={supported_count})",
        pad=12,
      )
    else:
      ordered = top_discovery.sort_values("delta_effect", ascending=False)
      labels = [
        wrap_label(
          f"{row['display_pair_label']} | "
          f"{CANCER_LABELS.get(str(row.get('cancer')), row.get('cancer'))} | "
          f"q={float(row[q_column]):.3g}",
          40,
        )
        for row in ordered.to_dict("records")
      ]
      positions_b = np.arange(len(ordered), dtype=float)
      axis_discovery.barh(positions_b, ordered["delta_effect"])
      axis_discovery.set_yticks(positions_b, labels)
      axis_discovery.axvline(0, linewidth=1, color="black")
      axis_discovery.set_xlabel(
        "Δ CRISPR effect (loss − intact); more negative is supportive"
      )
      axis_discovery.set_title(
        f"FDR-supported all-target discoveries (n={supported_count})"
      )
      axis_discovery.grid(axis="x", alpha=0.25)
  module.panel_label(axis_discovery, "B")
  fig.suptitle(str(item["title"]))

  source = pd.concat(
    [
      class_summary.assign(panel="A_class_summary"),
      top_discovery.assign(panel="B_all_target_discovery"),
    ],
    ignore_index=True,
    sort=False,
  )
  return module.save_record(
    fig=fig,
    item=item,
    output_dir=output_dir,
    source_dir=source_dir,
    source_data=source,
    inputs=input_paths,
    strict_layout=strict,
  )


def figure_7(
  module: Any,
  item: dict,
  pharmacology: pd.DataFrame,
  output_dir,
  source_dir,
  strict: bool,
  input_path,
  top_n: int,
):
  """Render density and compound-resolved hypotheses without NA display bugs."""
  set_publication_style()
  fig = plt.figure(figsize=(19.0, 9.4), constrained_layout=True)
  grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.55])
  axis_density = fig.add_subplot(grid[0, 0])
  axis_table = fig.add_subplot(grid[0, 1])
  axis_table.set_axis_off()

  required = {
    "coverage_adjusted_rses",
    "pharmacology_adjusted",
    "therapeutic_hypothesis_score",
  }
  eligible = (
    pharmacology.dropna(subset=list(required)).copy()
    if not pharmacology.empty and required.issubset(pharmacology.columns)
    else pd.DataFrame()
  )
  top = pd.DataFrame()
  if eligible.empty:
    module.placeholder(
      axis_density,
      "Pharmacological evidence density",
      "No pharmacology hypotheses were available.",
    )
    module.placeholder(
      axis_table,
      "Compound-resolved hypotheses",
      "Run pharmacology acquisition and prioritization.",
    )
  else:
    eligible = add_display_pair_columns(eligible)
    compound_parts = [compound_resolution(row) for _, row in eligible.iterrows()]
    eligible["compound_display"] = [part[0] for part in compound_parts]
    eligible["compound_resolution"] = [part[1] for part in compound_parts]

    density = axis_density.hexbin(
      pd.to_numeric(eligible["coverage_adjusted_rses"], errors="coerce"),
      pd.to_numeric(eligible["pharmacology_adjusted"], errors="coerce"),
      gridsize=42,
      mincnt=1,
    )
    colorbar = fig.colorbar(density, ax=axis_density, fraction=0.046, pad=0.04)
    colorbar.set_label("Hypothesis density")
    axis_density.set_xlabel("Coverage-adjusted vulnerability")
    axis_density.set_ylabel("Coverage-adjusted pharmacology")
    axis_density.set_title("All scored hypotheses")
    axis_density.grid(alpha=0.18)

    resolved = eligible.loc[
      eligible["compound_resolution"].ne("target_only")
    ].copy()
    dedup_columns = [
      column
      for column in ("cancer", "pair_id", "drug_key", "drug_id")
      if column in resolved.columns
    ]
    resolved = resolved.sort_values(
      "therapeutic_hypothesis_score",
      ascending=False,
    )
    if dedup_columns:
      resolved = resolved.drop_duplicates(dedup_columns)
    top = resolved.head(top_n).copy()

    offsets = [
      (7, 7), (7, -10), (-13, 7), (-13, -10),
      (14, 0), (-20, 0), (0, 14), (0, -16),
    ]
    for number, (_, row) in enumerate(top.iterrows(), start=1):
      x = float(row["coverage_adjusted_rses"])
      y = float(row["pharmacology_adjusted"])
      offset = offsets[(number - 1) % len(offsets)]
      axis_density.annotate(
        str(number),
        (x, y),
        xytext=offset,
        textcoords="offset points",
        ha="center",
        va="center",
        fontsize=8.5,
        fontweight="bold",
        bbox={"boxstyle": "circle,pad=0.20", "facecolor": "white", "edgecolor": "black"},
      )

    if top.empty:
      module.placeholder(
        axis_table,
        "Compound-resolved hypotheses",
        "No compound-resolved pharmacology rows were available; target-level evidence was retained only in supplementary source tables.",
      )
    else:
      table_rows = []
      for number, (_, row) in enumerate(top.iterrows(), start=1):
        cancer = CANCER_LABELS.get(
          str(row.get("cancer")),
          present_text(row.get("cancer")) or "NA",
        )
        sources = present_text(row.get("pharmacology_sources")) or "Unspecified"
        table_rows.append([
          number,
          wrap_label(cancer, 12),
          resolved_pair_label(row, width=22),
          wrap_label(row["compound_display"], 22),
          wrap_label(sources, 18),
          f"{float(row['therapeutic_hypothesis_score']):.3f}",
        ])
      table = axis_table.table(
        cellText=table_rows,
        colLabels=["#", "Cancer", "Vulnerability/context", "Compound / ID", "Sources", "Score"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.045, 0.10, 0.25, 0.23, 0.24, 0.075],
      )
      table.auto_set_font_size(False)
      table.set_fontsize(7.8)
      table.scale(1.0, 1.55)
      axis_table.set_title(
        "Top compound-resolved experimental hypotheses",
        pad=12,
      )

  module.panel_label(axis_density, "A")
  module.panel_label(axis_table, "B")
  fig.suptitle(
    "Pharmacological evidence and compound-resolved experimental hypotheses"
  )
  return module.save_record(
    fig=fig,
    item=item,
    output_dir=output_dir,
    source_dir=source_dir,
    source_data=top,
    inputs=[input_path],
    strict_layout=strict,
  )
