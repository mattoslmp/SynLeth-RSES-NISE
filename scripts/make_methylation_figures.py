#!/usr/bin/env python3
"""Generate supplementary Figures S70-S72 for methylation evidence."""
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

SCRIPT = "scripts/make_methylation_figures.py"


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if path.exists() and path.stat().st_size:
    return pd.read_csv(path, sep="\t", low_memory=False)
  return pd.DataFrame([{
    "evidence_status": "unavailable_technical_or_source_failure",
    "source_path": str(path),
    "reason": "No eligible TCGA/GDC methylation records were available.",
    "interpretation_boundary": (
      "Unavailable methylation is not negative biological evidence."
    ),
  }])


def unavailable(frame: pd.DataFrame) -> bool:
  return (
    "evidence_status" in frame.columns
    and frame["evidence_status"]
      .astype(str)
      .str.startswith("unavailable")
      .all()
  )


def coverage_figure(frame: pd.DataFrame) -> tuple[plt.Figure, pd.DataFrame]:
  set_publication_style()
  fig, axis = plt.subplots(figsize=(10.5, 6.8), constrained_layout=True)
  if unavailable(frame) or frame.empty:
    placeholder(axis, "Methylation coverage", "TCGA/GDC methylation was unavailable.")
    return fig, frame
  required = {"cancer", "gene", "probe_count", "primary_tumor_sample_count"}
  if not required.issubset(frame.columns):
    placeholder(axis, "Methylation coverage", "The gene summary lacks mandatory columns.")
    return fig, frame
  source = (
    frame.groupby("cancer", as_index=False)
    .agg(
      genes_with_methylation=("gene", "nunique"),
      median_probe_count=("probe_count", "median"),
      median_primary_tumor_samples=("primary_tumor_sample_count", "median"),
    )
  )
  x = np.arange(len(source))
  axis.bar(x, source["genes_with_methylation"])
  axis.set_xticks(x, [str(value).capitalize() for value in source["cancer"]])
  axis.set_ylabel("Candidate genes with methylation evidence")
  axis.set_xlabel("Cancer context")
  axis.grid(axis="y", alpha=0.25)
  for index, row in source.iterrows():
    axis.text(
      index,
      row["genes_with_methylation"],
      f"{int(row['genes_with_methylation'])} genes\n"
      f"median {row['median_primary_tumor_samples']:.0f} tumors",
      ha="center",
      va="bottom",
      fontsize=8,
    )
  return fig, source


def pair_figure(frame: pd.DataFrame) -> tuple[plt.Figure, pd.DataFrame]:
  set_publication_style()
  fig, axis = plt.subplots(figsize=(10.5, 7.0), constrained_layout=True)
  if unavailable(frame) or frame.empty:
    placeholder(axis, "Pair methylation context", "No eligible pair-level methylation was available.")
    return fig, frame
  source = frame.copy()
  for column in (
    "methylation_primary_tumor_overlap_n",
    "component_promoter_methylation_context",
  ):
    source[column] = pd.to_numeric(source.get(column), errors="coerce")
  plotted = source.dropna(subset=[
    "methylation_primary_tumor_overlap_n",
    "component_promoter_methylation_context",
  ])
  if plotted.empty:
    placeholder(axis, "Pair methylation context", "No pair satisfied probe and sample requirements.")
  else:
    axis.scatter(
      plotted["methylation_primary_tumor_overlap_n"],
      plotted["component_promoter_methylation_context"],
      s=38,
      alpha=0.68,
    )
    axis.set_xlabel("Overlapping primary-tumor cases")
    axis.set_ylabel("Coverage-adjusted methylation-context component")
    axis.set_ylim(-0.02, 1.02)
    axis.grid(alpha=0.25)
    label = "pair_id" if "pair_id" in plotted else None
    if label:
      for row in plotted.nlargest(
        min(10, len(plotted)),
        "component_promoter_methylation_context",
      ).to_dict("records"):
        axis.annotate(
          wrap_label(row[label], 26),
          (
            row["methylation_primary_tumor_overlap_n"],
            row["component_promoter_methylation_context"],
          ),
          xytext=(4, 4),
          textcoords="offset points",
          fontsize=7.2,
        )
  return fig, source


def integration_figure(frame: pd.DataFrame) -> tuple[plt.Figure, pd.DataFrame]:
  set_publication_style()
  fig, axis = plt.subplots(figsize=(10.5, 7.0), constrained_layout=True)
  if unavailable(frame) or frame.empty:
    placeholder(axis, "Methylation integration", "Integrated regulatory methylation was unavailable.")
    return fig, frame
  source = frame.copy()
  source["regulatory_promoter_methylation_context"] = pd.to_numeric(
    source.get("regulatory_promoter_methylation_context"),
    errors="coerce",
  )
  source["component_regulatory_network_composite"] = pd.to_numeric(
    source.get("component_regulatory_network_composite"),
    errors="coerce",
  )
  plotted = source.dropna(subset=[
    "regulatory_promoter_methylation_context",
    "component_regulatory_network_composite",
  ])
  if plotted.empty:
    placeholder(axis, "Methylation integration", "No pair had an observed methylation subcomponent.")
  else:
    axis.scatter(
      plotted["regulatory_promoter_methylation_context"],
      plotted["component_regulatory_network_composite"],
      s=38,
      alpha=0.68,
    )
    axis.set_xlabel("Methylation-context subcomponent")
    axis.set_ylabel("Integrated regulatory-network component")
    axis.set_xlim(-0.02, 1.02)
    axis.set_ylim(-0.02, 1.02)
    axis.grid(alpha=0.25)
  return fig, source


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
  config = yaml.safe_load(resolve(args.config).read_text(encoding="utf-8"))
  registry = {
    record["id"]: record
    for record in config.get("supplementary_figures", [])
  }
  required = {"Figure_S70", "Figure_S71", "Figure_S72"}
  if required - set(registry):
    raise RuntimeError(
      f"Methylation figure registry incomplete: {sorted(required - set(registry))}"
    )
  specifications = [
    (
      70,
      resolve("data/processed/epigenetics/methylation/tcga_nise_methylation_gene_summary.tsv"),
      coverage_figure,
    ),
    (
      71,
      resolve("data/processed/epigenetics/methylation/tcga_nise_methylation_pair_metrics.tsv"),
      pair_figure,
    ),
    (
      72,
      resolve("data/processed/regulatory/expanded_pair_functional_evidence_by_cancer.tsv"),
      integration_figure,
    ),
  ]
  records: list[FigureRecord] = []
  for number, source_path, builder in specifications:
    figure_id = f"Figure_S{number}"
    item = registry[figure_id]
    frame = read_optional(source_path)
    figure, exact_source = builder(frame)
    records.append(save_record(
      fig=figure,
      figure_id=figure_id,
      file_name=str(item["file"]),
      title=str(item["title"]),
      caption=str(item.get("caption") or item["title"]),
      output_root=output_root,
      source=exact_source,
      inputs=[source_path],
      script=SCRIPT,
      strict=args.strict_layout,
    ))
  write_figure_manifest(
    records,
    output_root / "manifests/methylation_supplementary_figure_manifest.tsv",
  )
  write_legends_markdown(
    records,
    output_root / "manuscript_assets/methylation_supplementary_figure_legends.md",
  )
  summary = pd.DataFrame([asdict(record) for record in records])
  if len(summary) != 3 or set(summary["figure_id"]) != required:
    raise RuntimeError("Methylation figures S70-S72 were not all generated")
  if args.strict_layout and not summary["layout_status"].eq("pass").all():
    raise RuntimeError("One or more methylation figures failed layout validation")


if __name__ == "__main__":
  main()
