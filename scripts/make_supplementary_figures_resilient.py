#!/usr/bin/env python3
"""Run supplementary-figure generation with layout and semantic resilience."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.make_supplementary_figures as target
from scripts.publication_layout_resilience import (
  dense_simple_bar,
  save_figure_triplet,
)
from scripts.publication_scientific_semantics import pair_labels


def simple_bar(
  item,
  frame,
  label,
  value,
  output_dir,
  source_dir,
  inputs,
  strict,
  xlabel,
):
  return dense_simple_bar(
    target,
    item,
    frame,
    label,
    value,
    output_dir,
    source_dir,
    inputs,
    strict,
    xlabel,
  )


def pair_label(frame, width=32):
  return pair_labels(frame, width)


target.save_figure_triplet = save_figure_triplet
target.simple_bar = simple_bar
target.pair_label = pair_label


if __name__ == "__main__":
  target.main()
