#!/usr/bin/env python3
"""Run main-figure generation with layout and scientific-display resilience."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.make_main_figures as target
from scripts.publication_layout_resilience import save_figure_triplet
from scripts.publication_scientific_semantics import (
  add_display_pair_columns,
  figure_6 as semantic_figure_6,
  figure_7 as semantic_figure_7,
  pair_labels,
)

_original_figure_3 = target.figure_3


def figure_3(
  item,
  ranking,
  output_dir,
  source_dir,
  strict,
  input_path,
  top_n,
):
  return _original_figure_3(
    item,
    add_display_pair_columns(ranking),
    output_dir,
    source_dir,
    strict,
    input_path,
    top_n,
  )


def figure_6(
  item,
  ranking,
  discovery,
  output_dir,
  source_dir,
  strict,
  input_paths,
):
  return semantic_figure_6(
    target,
    item,
    ranking,
    discovery,
    output_dir,
    source_dir,
    strict,
    input_paths,
  )


def figure_7(
  item,
  pharmacology,
  output_dir,
  source_dir,
  strict,
  input_path,
  top_n,
):
  return semantic_figure_7(
    target,
    item,
    pharmacology,
    output_dir,
    source_dir,
    strict,
    input_path,
    top_n,
  )


target.save_figure_triplet = save_figure_triplet
target.pair_labels = pair_labels
target.figure_3 = figure_3
target.figure_6 = figure_6
target.figure_7 = figure_7


if __name__ == "__main__":
  target.main()
