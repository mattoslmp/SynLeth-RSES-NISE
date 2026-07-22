#!/usr/bin/env python3
"""Block final packaging until every registered figure passes manual inspection."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REQUIRED_CHECKS = (
  "font_legible",
  "no_overlap",
  "no_clipping",
  "labels_complete",
  "legend_matches_figure",
  "figure_matches_source_data",
  "png_pdf_svg_consistent",
  "missing_values_rendered_correctly",
  "scientific_interpretation_appropriate",
  "document_page_checked",
)


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  root = resolve_path(args.article_root)
  checklist_path = (
    root / "review_records/MANUAL_VISUAL_INSPECTION_CHECKLIST.tsv"
  )
  manifest_path = root / "manifests/figure_manifest.tsv"
  if not checklist_path.exists() or checklist_path.stat().st_size == 0:
    raise SystemExit(
      "Manual-inspection checklist is absent. Generate and complete it before packaging."
    )
  checklist = pd.read_csv(checklist_path, sep="\t", low_memory=False)
  figures = pd.read_csv(manifest_path, sep="\t", low_memory=False)
  if len(checklist) != len(figures):
    raise SystemExit(
      f"Checklist has {len(checklist)} rows but figure manifest has {len(figures)}"
    )
  missing_columns = sorted(
    {"manual_review_status", "reviewer", "review_date", *REQUIRED_CHECKS}
    - set(checklist.columns)
  )
  if missing_columns:
    raise SystemExit(f"Checklist missing columns: {missing_columns}")
  failed = checklist.loc[
    ~checklist["manual_review_status"].astype(str).eq(
      "passed_manual_review"
    )
  ].copy()
  for column in REQUIRED_CHECKS:
    failed = pd.concat([
      failed,
      checklist.loc[
        ~checklist[column].astype(str).str.casefold().eq("pass")
      ],
    ]).drop_duplicates("figure_id")
  failed = pd.concat([
    failed,
    checklist.loc[
      checklist["reviewer"].fillna("").astype(str).str.strip().eq("")
      | checklist["review_date"].fillna("").astype(str).str.strip().eq("")
    ],
  ]).drop_duplicates("figure_id")
  if not failed.empty:
    raise SystemExit(
      "Final packaging is blocked. Incomplete manual inspection for:\n"
      + failed[["figure_id", "manual_review_status"]].to_string(index=False)
    )
  print(
    f"Manual inspection completion validated for {len(checklist)} figures."
  )


if __name__ == "__main__":
  main()
