from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from scripts.make_supplementary_figures_resilient import simple_bar
from scripts.publication_layout_resilience import save_figure_triplet


def test_axis_off_panel_passes_strict_layout(tmp_path: Path) -> None:
  fig, axis = plt.subplots(figsize=(7, 5), constrained_layout=True)
  axis.set_axis_off()
  axis.text(
    0.5,
    0.5,
    "Axis-off publication schematic",
    transform=axis.transAxes,
    ha="center",
    va="center",
  )
  audit = save_figure_triplet(
    fig,
    tmp_path / "axis_off",
    "Figure_axis_off",
    strict_layout=True,
  )
  assert audit.status == "pass"
  assert not audit.warnings


def test_dense_multiline_bar_labels_pass_strict_layout(tmp_path: Path) -> None:
  frame = pd.DataFrame({
    "activity_label": [
      f"NISE group {index}: a deliberately long biological activity label "
      "requiring multiple wrapped lines"
      for index in range(1, 16)
    ],
    "directed_candidates": list(range(1, 16)),
  })
  item = {
    "id": "Figure_S2_test",
    "file": "Figure_S2_test",
    "title": "Complete curated NISE activity groups",
    "caption": "Dense-label regression test.",
  }
  record = simple_bar(
    item,
    frame,
    "activity_label",
    "directed_candidates",
    tmp_path / "figures",
    tmp_path / "source",
    [],
    True,
    "Directed NISE hypotheses",
  )
  assert record.layout_status == "pass"
  assert not record.layout_warnings
  for extension in ("png", "pdf", "svg"):
    assert (tmp_path / "figures" / f"Figure_S2_test.{extension}").exists()
