#!/usr/bin/env python3
"""Run main-figure generation with layout, scientific-display and evidence-audit resilience."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.make_main_figures as target
from scripts.publication_figure5_cancer_context import (
  figure_5 as corrected_figure_5,
)
from scripts.publication_layout_resilience import save_figure_triplet
from scripts.publication_main_figure_corrections import (
  figure_1 as corrected_figure_1,
  figure_2 as corrected_figure_2,
  figure_4 as corrected_figure_4,
)
from scripts.publication_scientific_semantics import (
  add_display_pair_columns,
  figure_6 as semantic_figure_6,
  figure_7 as semantic_figure_7,
  pair_labels,
)

_original_figure_3 = target.figure_3


def read_optional(path: Path) -> pd.DataFrame:
  """Read heterogeneous publication TSVs in one type-inference pass.

  The ranking and pharmacology tables intentionally mix identifiers, nullable
  values and text in some columns. ``low_memory=False`` prevents chunk-wise dtype
  inference and the associated DtypeWarning without coercing identifiers or
  missingness states to invented numeric values.
  """
  if not path.exists():
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t", low_memory=False)
  except EmptyDataError:
    return pd.DataFrame()


def figure_1(item, output_dir, source_dir, strict):
  return corrected_figure_1(
    target,
    item,
    output_dir,
    source_dir,
    strict,
  )


def figure_2(
  item,
  candidates,
  output_dir,
  source_dir,
  strict,
  input_path,
):
  return corrected_figure_2(
    target,
    item,
    candidates,
    output_dir,
    source_dir,
    strict,
    input_path,
  )


def figure_3(
  item,
  ranking,
  output_dir,
  source_dir,
  strict,
  input_path,
  top_n,
):
  enriched = add_display_pair_columns(ranking)
  audit_path = (
    ROOT
    / "article_outputs"
    / "tables"
    / "qc"
    / "candidate_domain_evidence_audit.tsv"
  )
  if not audit_path.exists() or audit_path.stat().st_size == 0:
    raise RuntimeError(
      f"Figure 3 requires the candidate-domain audit: {audit_path}"
    )
  audit = pd.read_csv(audit_path, sep="\t", low_memory=False)
  onco = audit.loc[
    audit["domain_family"].astype(str).eq("RSES-Onco")
  ].copy()
  counts = onco.groupby(
    ["cancer", "candidate_id"],
    as_index=False,
  ).agg(
    eligible_domain_count=("eligible", "sum"),
    evidence_domain_count=("evidence_present", "sum"),
    technical_failure_domain_count=(
      "evidence_state",
      lambda values: int(
        pd.Series(values).eq("technical_failure").sum()
      ),
    ),
    insufficient_sample_domain_count=(
      "evidence_state",
      lambda values: int(
        pd.Series(values).eq("insufficient_sample").sum()
      ),
    ),
  )
  enriched = enriched.merge(
    counts,
    left_on=["cancer", "pair_id"],
    right_on=["cancer", "candidate_id"],
    how="left",
  )
  enriched["eligible_domain_count"] = pd.to_numeric(
    enriched["eligible_domain_count"],
    errors="coerce",
  )
  enriched["evidence_domain_count"] = pd.to_numeric(
    enriched["evidence_domain_count"],
    errors="coerce",
  )
  enriched["statistical_status"] = "prioritized_hypothesis"
  p_column = "p_value" if "p_value" in enriched else None
  q_column = next(
    (
      column
      for column in (
        "q_value_bh_within_loss_cancer",
        "q_value_bh",
      )
      if column in enriched
    ),
    None,
  )
  if p_column:
    p_values = pd.to_numeric(
      enriched[p_column],
      errors="coerce",
    )
    enriched.loc[
      p_values < 0.05,
      "statistical_status",
    ] = "nominally_significant"
  if q_column:
    q_values = pd.to_numeric(
      enriched[q_column],
      errors="coerce",
    )
    enriched.loc[
      q_values < 0.05,
      "statistical_status",
    ] = "fdr_supported"
  status_text = enriched.get(
    "status",
    pd.Series("", index=enriched.index),
  ).fillna("").astype(str).str.casefold()
  known = status_text.str.contains(
    "known|validated|benchmark|control",
    regex=True,
  )
  enriched.loc[
    known
    & enriched["statistical_status"].eq(
      "prioritized_hypothesis"
    ),
    "statistical_status",
  ] = "known_or_benchmark_hypothesis"
  enriched["support_level"] = enriched.get(
    "priority_class",
    pd.Series("not_recorded", index=enriched.index),
  )
  return _original_figure_3(
    item,
    enriched,
    output_dir,
    source_dir,
    strict,
    input_path,
    top_n,
  )


def figure_4(
  item,
  ranking,
  output_dir,
  source_dir,
  strict,
  input_path,
):
  return corrected_figure_4(
    target,
    item,
    ranking,
    output_dir,
    source_dir,
    strict,
    input_path,
  )


def figure_5(
  item,
  ranking,
  output_dir,
  source_dir,
  strict,
  input_path,
  top_n,
):
  return corrected_figure_5(
    target,
    item,
    ranking,
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


target.read_optional = read_optional
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
