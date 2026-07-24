#!/usr/bin/env python3
"""Validate the 86-figure/68-table v0.12.0 package around the stable core validator."""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

import pandas as pd


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def supplementary_number(value: object, prefix: str) -> int | None:
  match = re.fullmatch(rf"{re.escape(prefix)}(\d+)", str(value))
  return int(match.group(1)) if match else None


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument("--run-marker", default=None)
  parser.add_argument("--strict-sources", action="store_true")
  args = parser.parse_args()

  article_root = resolve(args.article_root)
  figure_manifest = article_root / "manifests/figure_manifest.tsv"
  table_manifest = article_root / "manifests/table_manifest.tsv"
  if not figure_manifest.exists() or not table_manifest.exists():
    raise FileNotFoundError("Publication manifests are absent")

  figures = pd.read_csv(figure_manifest, sep="\t", low_memory=False)
  tables = pd.read_csv(table_manifest, sep="\t", low_memory=False)

  expected_figures = {f"Figure_{index}" for index in range(1, 9)} | {
    f"Figure_S{index}" for index in range(1, 79)
  }
  observed_figures = set(figures["figure_id"].astype(str))
  if observed_figures != expected_figures:
    raise SystemExit(
      "Final figure set mismatch; missing="
      f"{sorted(expected_figures - observed_figures)}; extra="
      f"{sorted(observed_figures - expected_figures)}"
    )
  if len(figures) != 86:
    raise SystemExit(f"Expected 86 figures, observed {len(figures)}")
  if not figures["layout_status"].astype(str).eq("pass").all():
    failed = figures.loc[
      ~figures["layout_status"].astype(str).eq("pass")
    ]
    raise SystemExit("Layout failures:\n" + failed.to_string(index=False))

  main_tables = tables["category"].astype(str).eq("main")
  supplementary_tables = tables["category"].astype(str).eq("supplementary")
  if int(main_tables.sum()) != 4 or int(supplementary_tables.sum()) != 64:
    raise SystemExit(
      f"Expected 4 main and 64 supplementary tables; observed "
      f"{int(main_tables.sum())} and {int(supplementary_tables.sum())}"
    )
  if len(tables) != 68:
    raise SystemExit(f"Expected 68 tables, observed {len(tables)}")

  with tempfile.TemporaryDirectory(
    prefix="rses_v0120_core_validation_"
  ) as temporary:
    temporary_path = Path(temporary)
    figure_backup = temporary_path / "figure_manifest.tsv"
    table_backup = temporary_path / "table_manifest.tsv"
    shutil.copy2(figure_manifest, figure_backup)
    shutil.copy2(table_manifest, table_backup)
    try:
      core_figures = figures.loc[
        figures["figure_id"].astype(str).map(
          lambda value: not value.startswith("Figure_S")
          or (
            supplementary_number(value, "Figure_S") or 999
          ) <= 70
        )
      ]
      core_tables = tables.loc[
        tables["table_id"].astype(str).map(
          lambda value: (
            (match := re.match(r"^Table_S(\d+)", value)) is None
            or int(match.group(1)) <= 52
          )
        )
      ]
      core_figures.to_csv(figure_manifest, sep="\t", index=False)
      core_tables.to_csv(table_manifest, sep="\t", index=False)
      command = [
        sys.executable,
        "-u",
        str(ROOT / "scripts/validate_publication_outputs.py"),
        "--article-root",
        str(article_root),
      ]
      if args.run_marker:
        command.extend([
          "--run-marker",
          str(resolve(args.run_marker)),
        ])
      subprocess.run(command, cwd=ROOT, check=True)
    finally:
      shutil.copy2(figure_backup, figure_manifest)
      shutil.copy2(table_backup, table_manifest)

  extended_command = [
    sys.executable,
    "-u",
    str(ROOT / "scripts/validate_extended_multiomics_integrity.py"),
    "--article-root",
    str(article_root),
    "--require-publication-assets",
  ]
  if args.strict_sources:
    extended_command.append("--strict-sources")
  subprocess.run(extended_command, cwd=ROOT, check=True)

  image_count = sum(
    1
    for path in (article_root / "figures").rglob("*")
    if path.is_file()
    and path.suffix.lower() in {".png", ".pdf", ".svg"}
  )
  if image_count != 258:
    raise SystemExit(
      f"Expected 258 PNG/PDF/SVG files, observed {image_count}"
    )
  print("RSES-Onco v0.12.0 publication package validation passed.")
  print("Main figures: 8; supplementary figures: 78; image files: 258")
  print("Main tables: 4; supplementary tables: 64")


if __name__ == "__main__":
  main()
