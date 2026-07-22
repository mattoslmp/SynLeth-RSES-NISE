#!/usr/bin/env python3
"""Publication-layout helpers for axis-off schematics and dense labels."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rses_onco.publication import (
  dynamic_height,
  save_figure_triplet as _save_figure_triplet,
  set_publication_style,
  wrap_label,
)


def save_figure_triplet(
  fig: plt.Figure,
  base_path: str | Path,
  figure_id: str,
  strict_layout: bool = True,
):
  """Clear latent ticks on axis-off panels before the strict layout audit.

  ``Axes.set_axis_off`` suppresses ticks in rendered output but Matplotlib keeps
  internal tick-label artists. Clearing only those latent artists prevents false
  clipping failures without weakening checks for visible axes.
  """
  for axis in fig.axes:
    if not axis.axison:
      axis.set_xticks([])
      axis.set_yticks([])
  return _save_figure_triplet(
    fig,
    base_path,
    figure_id,
    strict_layout=strict_layout,
  )


def dense_simple_bar(
  module: Any,
  item: dict,
  frame: pd.DataFrame,
  label: str,
  value: str,
  output_dir: Path,
  source_dir: Path,
  inputs: list[Path],
  strict: bool,
  xlabel: str,
):
  """Render horizontal bars with height based on wrapped label line count."""
  set_publication_style()
  ordered = frame.sort_values(value).copy() if not frame.empty else frame

  if ordered.empty:
    labels: list[str] = []
    visual_lines = 1
  else:
    labels = [
      wrap_label(item_value, 36)
      for item_value in ordered[label].astype(str)
    ]
    visual_lines = sum(
      max(1, current.count("\n") + 1)
      for current in labels
    )

  height = dynamic_height(
    visual_lines,
    minimum=6.5,
    per_row=0.42,
    maximum=30.0,
  )
  fig, axis = plt.subplots(
    figsize=(13.5, height),
    constrained_layout=True,
  )

  if ordered.empty:
    module.placeholder(
      axis,
      str(item["title"]),
      "No eligible observations were available in this release.",
    )
  else:
    positions = np.arange(len(ordered), dtype=float)
    axis.barh(
      positions,
      pd.to_numeric(ordered[value], errors="coerce"),
    )
    axis.set_yticks(positions, labels)
    axis.tick_params(axis="y", labelsize=9.0, pad=7)
    axis.margins(y=0.025)
    axis.set_xlabel(xlabel)
    axis.set_title(str(item["title"]), pad=12)
    axis.grid(axis="x", alpha=0.25)

  return module.save_record(
    fig,
    item,
    output_dir,
    source_dir,
    ordered,
    inputs,
    strict,
  )
