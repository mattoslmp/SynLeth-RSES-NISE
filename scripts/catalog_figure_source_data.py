#!/usr/bin/env python3
"""Catalog and validate the exact machine-readable table used by every figure."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
  return pd.read_csv(path, sep="\t", low_memory=False)


def infer_category(figure_id: str) -> str:
  return (
    "main"
    if figure_id.startswith("Figure_")
    and not figure_id.startswith("Figure_S")
    else "supplementary"
  )


def command_for_script(script: str) -> str:
  if "make_genomic_circos_figure" in script:
    return (
      "MPLBACKEND=Agg python -u "
      "scripts/make_genomic_circos_figure.py "
      "--config config/genomic_circos_asset.yaml "
      "--output-root article_outputs --strict-layout"
    )
  if "make_main_figures" in script:
    return (
      "MPLBACKEND=Agg python -u "
      "scripts/make_main_figures_resilient.py "
      "--output-root article_outputs --strict-layout"
    )
  if "make_supplementary_figures" in script:
    return (
      "MPLBACKEND=Agg python -u "
      "scripts/make_supplementary_figures_resilient.py "
      "--output-root article_outputs --strict-layout"
    )
  if "make_audit_supplementary_figures" in script:
    return (
      "MPLBACKEND=Agg python -u "
      "scripts/make_audit_supplementary_figures.py "
      "--output-root article_outputs --strict-layout"
    )
  if "make_extended_supporting_figures" in script:
    return (
      "MPLBACKEND=Agg python -u "
      "scripts/make_extended_supporting_figures.py "
      "--output-root article_outputs --strict-layout"
    )
  if "structure" in script:
    return (
      "MPLBACKEND=Agg python -u "
      "scripts/make_nise_structure_figures.py "
      "--output-root article_outputs --strict-layout"
    )
  return (
    "MPLBACKEND=Agg bash scripts/run_publication_pipeline.sh figures"
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  manifest_path = article_root / "manifests/figure_manifest.tsv"
  if not manifest_path.exists() or manifest_path.stat().st_size == 0:
    raise FileNotFoundError(
      f"Missing or empty figure manifest: {manifest_path}"
    )
  manifest = read_tsv(manifest_path)
  required = {
    "figure_id",
    "source_data_path",
    "script",
    "input_paths",
    "base_path",
  }
  missing = sorted(required - set(manifest.columns))
  if missing:
    raise ValueError(f"Figure manifest missing columns: {missing}")

  rows = []
  column_rows = []
  missing_files = []
  for record in manifest.to_dict("records"):
    figure_id = str(record["figure_id"])
    category = infer_category(figure_id)
    source = Path(str(record["source_data_path"]))
    if not source.is_absolute():
      source = ROOT / source
    if not source.exists() or source.stat().st_size == 0:
      missing_files.append(str(source))
      continue
    frame = read_tsv(source)
    destination_dir = (
      article_root / "tables" / "figure_data" / category
    )
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / f"{figure_id}_source_data.tsv"
    shutil.copy2(source, destination)
    source_hash = sha256(source)
    destination_hash = sha256(destination)
    if source_hash != destination_hash:
      raise RuntimeError(
        f"Source-data copy checksum mismatch for {figure_id}"
      )

    script = str(record.get("script") or "")
    base = Path(str(record["base_path"]))
    if not base.is_absolute():
      base = ROOT / base
    formats = []
    for extension in ("png", "svg", "pdf"):
      path = base.with_suffix(f".{extension}")
      if not path.exists() or path.stat().st_size == 0:
        missing_files.append(str(path))
      else:
        formats.append(extension)

    rows.append({
      "figure_id": figure_id,
      "category": category,
      "source_table": str(destination),
      "source_table_sha256": destination_hash,
      "rows": len(frame),
      "columns": len(frame.columns),
      "generator_script": script,
      "input_paths": record.get("input_paths"),
      "reproduction_command": command_for_script(script),
      "rendered_formats": ";".join(formats),
      "figure_base_path": str(base),
      "catalogued_at_utc": datetime.now(timezone.utc).isoformat(),
    })
    for column in frame.columns:
      series = frame[column]
      column_rows.append({
        "figure_id": figure_id,
        "column": column,
        "dtype": str(series.dtype),
        "non_missing_values": int(series.notna().sum()),
        "missing_values": int(series.isna().sum()),
        "unit": "not_recorded_in_source_table",
        "meaning": (
          "Refer to the figure legend and generator script; "
          "no undocumented meaning is invented."
        ),
      })

  if missing_files:
    raise RuntimeError(
      "Missing mandatory figure assets/source data:\n"
      + "\n".join(missing_files)
    )
  inventory = pd.DataFrame(rows)
  columns = pd.DataFrame(column_rows)
  if len(inventory) != len(manifest):
    raise RuntimeError(
      f"Catalogued {len(inventory)} figures; "
      f"manifest contains {len(manifest)}"
    )

  catalog_dir = article_root / "tables/figure_data"
  inventory.to_csv(
    catalog_dir / "figure_source_data_inventory.tsv",
    sep="\t",
    index=False,
  )
  columns.to_csv(
    catalog_dir / "figure_source_data_column_dictionary.tsv",
    sep="\t",
    index=False,
  )
  readme_lines = [
    "# Figure source-data reproduction",
    "",
    (
      "Every registered figure has an exact TSV copy in `main/` "
      "or `supplementary/`."
    ),
    (
      "The TSV is copied byte-for-byte from the source table written "
      "by the figure generator before rendering."
    ),
    "",
    (
      "| Figure | Category | Generator | Exact source data | "
      "Reproduction command |"
    ),
    "|---|---|---|---|---|",
  ]
  for row in rows:
    readme_lines.append(
      f"| {row['figure_id']} | {row['category']} | "
      f"`{row['generator_script']}` | `{row['source_table']}` | "
      f"`{row['reproduction_command']}` |"
    )
  (catalog_dir / "README.md").write_text(
    "\n".join(readme_lines) + "\n",
    encoding="utf-8",
  )

  subprocess.run([
    sys.executable,
    "-u",
    "scripts/build_publication_methods_documentation.py",
    "--article-root",
    str(article_root),
  ], cwd=ROOT, check=True)
  print(f"Catalogued exact source data for {len(inventory)} figures.")
  print(
    "Inventory: "
    f"{catalog_dir / 'figure_source_data_inventory.tsv'}"
  )


if __name__ == "__main__":
  main()
