#!/usr/bin/env python3
"""Build the human-readable figure/table reproduction registry.

The registry is a supplementary table rather than an opaque technical manifest.
It records the responsible script, expected inputs, command, output and document
location for every registered publication figure and table.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def figure_script(figure_id: str) -> str:
  number_text = figure_id.removeprefix("Figure_S")
  if figure_id.startswith("Figure_") and not figure_id.startswith("Figure_S"):
    return "scripts/make_main_figures_resilient.py" if figure_id != "Figure_8" else "scripts/make_nise_structure_figures.py"
  number = int(number_text)
  if number <= 14:
    return "scripts/make_supplementary_figures_resilient.py"
  if number <= 32:
    return "scripts/make_nise_structure_figures.py"
  if number <= 38:
    return "scripts/make_audit_supplementary_figures.py"
  return "scripts/make_extended_supporting_figures.py"


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/article_assets.yaml")
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()

  config_path = resolve_path(args.config)
  config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
  article_root = resolve_path(args.article_root)
  rows: list[dict[str, object]] = []

  for category, key in (("main", "main_figures"), ("supplementary", "supplementary_figures")):
    for record in config.get(key, []):
      figure_id = str(record["id"])
      script = figure_script(figure_id)
      rows.append({
        "asset_id": figure_id,
        "asset_type": "figure",
        "category": category,
        "scientific_title": record.get("title"),
        "script": script,
        "primary_inputs": "Declared by the generator and exact figure-source inventory",
        "intermediate_data": f"{article_root}/source_data/figures/{category}/{record['file']}_source_data.tsv",
        "reproduction_command": "MPLBACKEND=Agg bash scripts/run_publication_pipeline.sh assets-only",
        "outputs": ";".join(
          f"{article_root}/figures/{category}/{record['file']}.{extension}"
          for extension in ("png", "pdf", "svg")
        ),
        "dependencies": "environment.yml; pyproject.toml",
        "document_location": "main manuscript" if category == "main" else "supplementary document",
      })

  for category, key in (("main", "main_tables"), ("supplementary", "supplementary_tables")):
    for name in config.get(key, []):
      table_id = Path(name).stem
      rows.append({
        "asset_id": table_id,
        "asset_type": "table",
        "category": category,
        "scientific_title": table_id.replace("_", " "),
        "script": "scripts/export_article_tables.py",
        "primary_inputs": "Declared by table manifest and source_paths field",
        "intermediate_data": f"{article_root}/tables/{category}/{name}",
        "reproduction_command": "bash scripts/run_publication_pipeline.sh assets-only",
        "outputs": f"{article_root}/tables/{category}/{name}",
        "dependencies": "environment.yml; pyproject.toml",
        "document_location": "main manuscript" if category == "main" else "supplementary data",
      })

  frame = pd.DataFrame(rows)
  output = article_root / "tables/supplementary/Table_S44_asset_reproduction_registry.tsv"
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)

  docs = article_root / "manuscript_assets/PUBLICATION_ASSET_REPRODUCTION.md"
  docs.parent.mkdir(parents=True, exist_ok=True)
  lines = [
    "# Publication asset reproduction", "",
    "| Asset | Type | Category | Script | Command | Output |",
    "|---|---|---|---|---|---|",
  ]
  for row in frame.to_dict("records"):
    lines.append(
      f"| {row['asset_id']} | {row['asset_type']} | {row['category']} | "
      f"`{row['script']}` | `{row['reproduction_command']}` | `{row['outputs']}` |"
    )
  docs.write_text("\n".join(lines) + "\n", encoding="utf-8")
  print(f"Wrote reproduction registry: {output} ({len(frame):,} assets)")


if __name__ == "__main__":
  main()
