from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from rses_onco.publication import (
  FigureRecord,
  audit_figure_layout,
  save_figure_triplet,
  set_publication_style,
  write_figure_manifest,
  write_legends_markdown,
)


def test_scripted_figure_exports_png_pdf_svg(tmp_path: Path) -> None:
  set_publication_style()
  fig, axis = plt.subplots(figsize=(7, 5), constrained_layout=True)
  axis.plot([0, 1, 2], [0, 1, 0])
  axis.set_xlabel("Condition")
  axis.set_ylabel("Coverage-adjusted evidence")
  axis.set_title("Publication smoke figure")
  base = tmp_path / "Figure_smoke"
  audit = save_figure_triplet(fig, base, "Figure_smoke", strict_layout=True)
  assert audit.status == "pass"
  for extension in ("png", "pdf", "svg"):
    path = base.with_suffix(f".{extension}")
    assert path.exists()
    assert path.stat().st_size > 100
  payload = json.loads(base.with_suffix(".layout_audit.json").read_text())
  assert payload["status"] == "pass"
  assert payload["warnings"] == []


def test_axes_overlap_is_detected() -> None:
  set_publication_style()
  fig = plt.figure(figsize=(6, 4))
  first = fig.add_axes([0.1, 0.1, 0.7, 0.7])
  second = fig.add_axes([0.4, 0.4, 0.5, 0.5])
  first.plot([0, 1], [0, 1])
  second.plot([0, 1], [1, 0])
  audit = audit_figure_layout(fig, "Figure_overlap")
  plt.close(fig)
  assert audit.status == "warning"
  assert any(value.startswith("axes_overlap:") for value in audit.warnings)


def test_figure_manifest_and_legends(tmp_path: Path) -> None:
  record = FigureRecord(
    figure_id="Figure_1",
    category="main",
    title="Framework",
    caption="A fully scripted framework figure.",
    base_path="article_outputs/figures/main/Figure_1_framework",
    source_data_path="article_outputs/source_data/Figure_1.tsv",
    input_paths="input.tsv",
    formats="png;pdf;svg",
    layout_status="pass",
    layout_warnings="",
    script="scripts/make_main_figures.py",
  )
  manifest = write_figure_manifest([record], tmp_path / "manifest.tsv")
  legends = write_legends_markdown([record], tmp_path / "legends.md")
  frame = pd.read_csv(manifest, sep="\t")
  assert frame.loc[0, "figure_id"] == "Figure_1"
  assert "A fully scripted framework figure." in legends.read_text()


def test_strict_layout_rejects_overlap(tmp_path: Path) -> None:
  fig = plt.figure(figsize=(6, 4))
  fig.add_axes([0.1, 0.1, 0.7, 0.7])
  fig.add_axes([0.4, 0.4, 0.5, 0.5])
  with pytest.raises(RuntimeError, match="Layout audit failed"):
    save_figure_triplet(
      fig,
      tmp_path / "bad",
      "Figure_bad",
      strict_layout=True,
    )
  plt.close(fig)
