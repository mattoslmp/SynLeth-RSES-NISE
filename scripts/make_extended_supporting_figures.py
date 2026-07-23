#!/usr/bin/env python3
"""Generate source-backed supplementary Figures S39-S69.

Each registered figure reads a declared machine-readable source. An absent or
empty optional source is represented by an explicit unavailable-evidence row;
it is never converted to biological zero. Every figure is exported as PNG, PDF
and SVG with its exact plotting TSV and layout audit.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from rses_onco.publication import (
  FigureRecord,
  placeholder,
  set_publication_style,
  wrap_label,
  write_figure_manifest,
  write_legends_markdown,
)
from scripts.publication_audit_figures import save_record

SCRIPT = "scripts/make_extended_supporting_figures.py"


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_source(path: Path) -> pd.DataFrame:
  if path.exists() and path.stat().st_size:
    return pd.read_csv(path, sep="\t", low_memory=False)
  return pd.DataFrame([{
    "evidence_status": "unavailable_or_no_eligible_records",
    "source_path": str(path),
    "reason": "The declared source was absent, empty or contained no eligible observations.",
    "interpretation_boundary": "Unavailable evidence is not negative biological evidence.",
  }])


def numeric_columns(frame: pd.DataFrame) -> list[str]:
  columns: list[str] = []
  for column in frame.columns:
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.notna().any():
      columns.append(column)
  return columns


def label_column(frame: pd.DataFrame) -> str | None:
  preferred = (
    "hypothesis_direction", "display_pair", "pair_id", "gene", "candidate_gene",
    "cancer", "source", "domain_label", "module", "entity", "category",
  )
  return next((column for column in preferred if column in frame.columns), None)


def make_plot(frame: pd.DataFrame, figure_id: str, title: str) -> plt.Figure:
  set_publication_style()
  fig, axis = plt.subplots(figsize=(10.5, 7.2), constrained_layout=True)
  if (
    "evidence_status" in frame.columns
    and frame["evidence_status"].astype(str).str.startswith("unavailable").all()
  ):
    placeholder(
      axis,
      "Evidence availability",
      str(frame.iloc[0].get("reason", "Evidence unavailable.")),
    )
    return fig

  numbers = numeric_columns(frame)
  label = label_column(frame)
  if figure_id == "Figure_S52":
    preferred = [
      "methylation_pair_profile_divergence",
      "methylation_target_hypomethylation_support",
    ]
    numbers = [column for column in preferred if column in frame.columns]
  if len(numbers) >= 2:
    x, y = numbers[:2]
    work = frame[[x, y] + ([label] if label else [])].copy()
    work[x] = pd.to_numeric(work[x], errors="coerce")
    work[y] = pd.to_numeric(work[y], errors="coerce")
    work = work.dropna(subset=[x, y]).head(5000)
    if work.empty:
      placeholder(axis, "Evidence availability", "No eligible numeric observations were available.")
    else:
      axis.scatter(work[x], work[y], s=38, alpha=0.68)
      axis.set_xlabel(x.replace("_", " ").capitalize())
      axis.set_ylabel(y.replace("_", " ").capitalize())
      axis.grid(alpha=0.25)
      if label:
        ranked = work.assign(_magnitude=work[y].abs()).nlargest(
          min(8, len(work)),
          "_magnitude",
        )
        for row in ranked.to_dict("records"):
          axis.annotate(
            wrap_label(row.get(label, ""), 24),
            (row[x], row[y]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7.5,
          )
  elif len(numbers) == 1:
    value = numbers[0]
    work = frame.copy()
    work[value] = pd.to_numeric(work[value], errors="coerce")
    work = work.dropna(subset=[value]).sort_values(value, ascending=False).head(35)
    if work.empty:
      placeholder(axis, "Evidence availability", "No eligible numeric observations were available.")
    else:
      labels = work[label].astype(str) if label else work.index.astype(str)
      y = np.arange(len(work))
      axis.barh(y, work[value])
      axis.set_yticks(y, [wrap_label(text.replace("_", " "), 34) for text in labels])
      axis.invert_yaxis()
      axis.set_xlabel(value.replace("_", " ").capitalize())
      axis.grid(axis="x", alpha=0.25)
  else:
    counts = frame.astype(str).stack().value_counts().head(25)
    if counts.empty:
      placeholder(axis, "Evidence availability", "No eligible records were available.")
    else:
      y = np.arange(len(counts))
      axis.barh(y, counts.values)
      axis.set_yticks(y, [wrap_label(value, 34) for value in counts.index])
      axis.invert_yaxis()
      axis.set_xlabel("Records")
      axis.grid(axis="x", alpha=0.25)
  return fig


SOURCE_REGISTRY: dict[int, str] = {
  39: "article_outputs/tables/score_components/rses_onco_score_decomposition.tsv",
  40: "article_outputs/tables/robustness/controlled_weight_perturbation.tsv",
  41: "article_outputs/tables/robustness/missing_data_sensitivity.tsv",
  42: "article_outputs/tables/robustness/ranking_correlations.tsv",
  43: "article_outputs/tables/supporting_evidence/expression/compensation_and_dependency_contrasts.tsv",
  44: "article_outputs/tables/supporting_evidence/expression/coexpression_by_event_group.tsv",
  45: "article_outputs/tables/supporting_evidence/model_level/model_level_expression_crispr_copy_number.tsv",
  46: "data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv",
  47: "data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv",
  48: "data/processed/regulatory/wgcna/wgcna_correlation_fallback_all_cancers.tsv",
  49: "data/processed/regulatory/wgcna/colon/wgcna_soft_threshold_diagnostics.tsv",
  50: "data/processed/regulatory/expanded_pair_functional_evidence_by_cancer.tsv",
  51: "data/processed/regulatory/jaspar_promoter_tf_summary.tsv",
  52: "data/processed/regulatory/promoter_methylation_pair_metrics.tsv",
  53: "article_outputs/tables/supporting_evidence/networks/raw_sources/string_candidate_edges_all_channels.tsv",
  54: "article_outputs/tables/supporting_evidence/phenotypes/conditional_dependency_contrasts.tsv",
  55: "article_outputs/tables/supporting_evidence/phenotypes/conditional_dependency_contrasts.tsv",
  56: "article_outputs/tables/supporting_evidence/localization/hpa_candidate_localization.tsv",
  57: "article_outputs/tables/supporting_evidence/localization/hpa_candidate_localization.tsv",
  58: "article_outputs/tables/qc/candidate_domain_evidence_audit.tsv",
  59: "article_outputs/tables/supporting_evidence/structures/uniprot_candidate_annotations.tsv",
  60: "article_outputs/tables/supporting_evidence/genomic_context/tcga_gene_event_summary.tsv",
  61: "article_outputs/tables/supporting_evidence/genomic_context/tcga_gene_event_summary.tsv",
  62: "article_outputs/tables/supporting_evidence/pharmacology/pharmacology_evidence_long.tsv",
  63: "article_outputs/tables/supporting_evidence/pharmacology/pharmacology_evidence_long.tsv",
  64: "results/expanded_26Q1/full/expanded_rses_onco.tsv",
  65: "results/expanded_26Q1/full/expanded_rses_onco.tsv",
  66: "article_outputs/tables/qc/evidence_category_assignments.tsv",
  67: "article_outputs/tables/qc/coverage_by_domain.tsv",
  68: "data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv",
  69: "data/processed/regulatory/expanded_pair_functional_evidence_by_cancer.tsv",
}


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/article_assets.yaml")
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument(
    "--strict-layout",
    action=argparse.BooleanOptionalAction,
    default=True,
  )
  args = parser.parse_args()

  output_root = resolve(args.output_root)
  config = yaml.safe_load(resolve(args.config).read_text(encoding="utf-8")) or {}
  registry = {
    record["id"]: record
    for record in config.get("supplementary_figures", [])
  }
  expected = {f"Figure_S{number}" for number in range(39, 70)}
  missing = sorted(expected - set(registry))
  if missing:
    raise RuntimeError(f"Extended figure registry is incomplete: {missing}")

  records: list[FigureRecord] = []
  for number in range(39, 70):
    figure_id = f"Figure_S{number}"
    item = registry[figure_id]
    source_path = resolve(SOURCE_REGISTRY[number])
    frame = read_source(source_path)
    title = str(item["title"])
    caption = str(item.get("caption") or title)
    if figure_id == "Figure_S52":
      title = "Promoter methylation context"
      caption = (
        "Cancer-specific pair promoter-methylation divergence and conditional "
        "target hypomethylation support. Association is not causal silencing proof."
      )
    fig = make_plot(frame, figure_id, title)
    record = save_record(
      fig=fig,
      figure_id=figure_id,
      file_name=str(item["file"]),
      title=title,
      caption=caption,
      output_root=output_root,
      source=frame,
      inputs=[source_path],
      script=SCRIPT,
      strict=args.strict_layout,
    )
    records.append(record)
    print(f"Generated {figure_id}: {record.layout_status}", flush=True)

  write_figure_manifest(
    records,
    output_root / "manifests/extended_supplementary_figure_manifest.tsv",
  )
  write_legends_markdown(
    records,
    output_root / "manuscript_assets/extended_supplementary_figure_legends.md",
  )
  summary = pd.DataFrame([asdict(record) for record in records])
  if len(summary) != 31 or set(summary["figure_id"]) != expected:
    raise RuntimeError("Extended supplementary figure generation did not produce S39-S69")
  if args.strict_layout and not summary["layout_status"].eq("pass").all():
    raise RuntimeError("One or more extended supplementary figures failed layout validation")


if __name__ == "__main__":
  main()
