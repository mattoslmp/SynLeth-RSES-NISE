#!/usr/bin/env python3
"""Create the mandatory manual 100%-zoom inspection checklist.

Automated layout checks cannot certify scientific correspondence or final-document
legibility. This script creates one pending review row for every registered figure
and records the exact rendered files and source-data table that must be inspected.
It never marks a figure as manually approved.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  manifest_path = article_root / "manifests/figure_manifest.tsv"
  if not manifest_path.exists() or manifest_path.stat().st_size == 0:
    raise FileNotFoundError(f"Missing figure manifest: {manifest_path}")
  figures = pd.read_csv(manifest_path, sep="\t", low_memory=False)
  rows = []
  for record in figures.to_dict("records"):
    base = Path(str(record["base_path"]))
    if not base.is_absolute():
      base = ROOT / base
    rows.append({
      "figure_id": record["figure_id"],
      "category": record["category"],
      "png": str(base.with_suffix(".png")),
      "pdf": str(base.with_suffix(".pdf")),
      "svg": str(base.with_suffix(".svg")),
      "source_data": record["source_data_path"],
      "inspect_at_zoom": "100%",
      "font_legible": "pending",
      "no_overlap": "pending",
      "no_clipping": "pending",
      "labels_complete": "pending",
      "legend_matches_figure": "pending",
      "figure_matches_source_data": "pending",
      "png_pdf_svg_consistent": "pending",
      "missing_values_rendered_correctly": "pending",
      "scientific_interpretation_appropriate": "pending",
      "document_page_checked": "pending_document_source",
      "reviewer": "",
      "review_date": "",
      "notes": "",
      "manual_review_status": "pending_manual_review",
    })
  checklist = pd.DataFrame(rows)
  review_dir = article_root / "review_records"
  review_dir.mkdir(parents=True, exist_ok=True)
  tsv_path = review_dir / "MANUAL_VISUAL_INSPECTION_CHECKLIST.tsv"
  checklist.to_csv(tsv_path, sep="\t", index=False)

  markdown = [
    "# Mandatory manual visual inspection", "",
    "Automated validation is not final approval. Open every PDF at 100% zoom and compare it with the exact source-data TSV.",
    "", "For every row confirm:", "",
    "- text and symbols are legible at final manuscript size;",
    "- no overlap, clipping, corrupted characters or off-canvas labels;",
    "- panel letters, legends, color bars and statistical annotations are correct;",
    "- missing, non-eligible and observed zero values are visually distinct;",
    "- the rendered values correspond to the exact source-data table;",
    "- PNG, PDF and SVG convey the same content;",
    "- no computational hypothesis is mislabeled as validated, clinical or efficacious;",
    "- the figure and caption match after insertion into the rendered manuscript/supplement.",
    "", f"Checklist TSV: `{tsv_path}`", "",
    "Packaging remains blocked until every row is manually completed and `manual_review_status` is changed to `passed_manual_review` with reviewer and date.",
  ]
  (review_dir / "MANUAL_VISUAL_INSPECTION_INSTRUCTIONS.md").write_text(
    "\n".join(markdown) + "\n",
    encoding="utf-8",
  )
  print(f"Created {len(checklist)} pending manual-inspection rows")
  print(f"Checklist: {tsv_path}")


if __name__ == "__main__":
  main()
