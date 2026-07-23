#!/usr/bin/env python3
"""Generate Figure S70 with a strict-layout-safe complete 25-ring legend."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import pandas as pd

import scripts.make_genomic_circos_figure as target


def compact_complete_legend(
  axis: plt.Axes,
  tracks: pd.DataFrame,
) -> None:
  """Show every ring and rendering rule without clipping or overlap."""
  axis.set_xticks([])
  axis.set_yticks([])
  axis.set_axis_off()
  y = 0.98
  axis.text(
    0.0,
    y,
    "Ring key",
    fontsize=12,
    fontweight="bold",
    va="top",
    transform=axis.transAxes,
  )
  y -= 0.045
  for panel, title in (
    ("A", "Panel A — top-level RSES-Onco"),
    ("B", "Panel B — microniche and internal layers"),
  ):
    axis.text(
      0.0,
      y,
      title,
      fontsize=9.0,
      fontweight="bold",
      va="top",
      transform=axis.transAxes,
    )
    y -= 0.030
    subset = tracks.loc[
      tracks["panel"].astype(str).eq(panel)
    ].sort_values("ring_order")
    for row in subset.to_dict("records"):
      axis.text(
        0.0,
        y,
        f"{row['track_id']}  {row['track_label']}",
        fontsize=6.8,
        va="top",
        transform=axis.transAxes,
      )
      y -= 0.020
    y -= 0.012
  y -= 0.005
  axis.text(
    0.0,
    y,
    "Rendering rules",
    fontsize=9.0,
    fontweight="bold",
    va="top",
    transform=axis.transAxes,
  )
  y -= 0.028
  rules = [
    "• every eligible NISE/paralog gene is a genomic tick",
    "• every coordinate-complete pair is a chord",
    "• ring values are maxima across pair × cancer rows",
    "• hollow markers are missing/non-eligible, never zero",
    "• chromosome positions use Ensembl/GRCh38 coordinates",
    "• exact inputs are Supplementary Tables S45–S52",
  ]
  for rule in rules:
    axis.text(
      0.0,
      y,
      rule,
      fontsize=6.8,
      va="top",
      transform=axis.transAxes,
    )
    y -= 0.023


target.legend_panel = compact_complete_legend
target.SCRIPT = "scripts/make_genomic_circos_figure_resilient.py"


if __name__ == "__main__":
  target.main()
