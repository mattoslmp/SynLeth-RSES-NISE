#!/usr/bin/env python3
"""Run main-figure generation with layout, scientific-display and evidence-audit resilience."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.make_main_figures as target
from scripts.publication_layout_resilience import save_figure_triplet
from scripts.publication_main_figure_corrections import (
  figure_1 as corrected_figure_1,
  figure_2 as corrected_figure_2,
  figure_4 as corrected_figure_4,
  figure_5 as corrected_figure_5,
)
from scripts.publication_scientific_semantics import (
  add_display_pair_columns,
  figure_6 as semantic_figure_6,
  figure_7 as semantic_figure_7,
  pair_labels,
)

_original_figure_3 = target.figure_3


def figure_1(item, output_dir, source_dir, strict):
  return corrected_figure_1(target, item, output_dir, source_dir, strict)


def figure_2(item, candidates, output_dir, source_dir, strict, input_path):
  return corrected_figure_2(target, item, candidates, output_dir, source_dir, strict, input_path)


def figure_3(item, ranking, output_dir, source_dir, strict, input_path, top_n):
  enriched = add_display_pair_columns(ranking)
  component_columns = [
    column for column in enriched.columns
    if column.startswith("component_") or column.startswith("microniche_")
  ]
  enriched["eligible_domain_count"] = enriched[
    [column for column in component_columns if column.startswith("component_")]
  ].notna().sum(axis=1)
  enriched["evidence_domain_count"] = enriched[component_columns].notna().sum(axis=1)
  q_column = next((column for column in ("q_value_bh_within_loss_cancer", "q_value_bh") if column in enriched), None)
  p_column = "p_value" if "p_value" in enriched else None
  enriched["statistical_status"] = "prioritized_hypothesis"
  if p_column:
    p = enriched[p_column].apply(lambda value: float(value) if str(value) not in {"nan", "None"} else float("nan"))
    enriched.loc[p < 0.05, "statistical_status"] = "nominally_significant"
  if q_column:
    q = enriched[q_column].apply(lambda value: float(value) if str(value) not in {"nan", "None"} else float("nan"))
    enriched.loc[q < 0.05, "statistical_status"] = "fdr_supported"
  enriched["support_level"] = enriched.get("priority_class", "not_recorded")
  return _original_figure_3(
    item, enriched, output_dir, source_dir, strict, input_path, top_n,
  )


def figure_4(item, ranking, output_dir, source_dir, strict, input_path):
  return corrected_figure_4(target, item, ranking, output_dir, source_dir, strict, input_path)


def figure_5(item, ranking, output_dir, source_dir, strict, input_path, top_n):
  return corrected_figure_5(target, item, ranking, output_dir, source_dir, strict, input_path, top_n)


def figure_6(item, ranking, discovery, output_dir, source_dir, strict, input_paths):
  return semantic_figure_6(
    target, item, ranking, discovery, output_dir, source_dir, strict, input_paths,
  )


def figure_7(item, pharmacology, output_dir, source_dir, strict, input_path, top_n):
  return semantic_figure_7(
    target, item, pharmacology, output_dir, source_dir, strict, input_path, top_n,
  )


target.save_figure_triplet = save_figure_triplet
target.pair_labels = pair_labels
target.figure_1 = figure_1
target.figure_2 = figure_2
target.figure_3 = figure_3
target.figure_4 = figure_4
target.figure_5 = figure_5
target.figure_6 = figure_6
target.figure_7 = figure_7


if __name__ == "__main__":
  target.main()
