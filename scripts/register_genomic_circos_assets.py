#!/usr/bin/env python3
"""Register Figure S70 and Supplementary Tables S45-S52 in publication manifests."""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

TABLE_SOURCES = {
  "Table_S45_genomic_circos_gene_coordinates.tsv": (
    "data/processed/circos/genomic_circos_gene_coordinates.tsv",
    "scripts/build_genomic_circos_inputs.py",
  ),
  "Table_S46_genomic_circos_pair_links.tsv": (
    "data/processed/circos/genomic_circos_pair_links.tsv",
    "scripts/build_genomic_circos_inputs.py",
  ),
  "Table_S47_genomic_circos_ring_values.tsv": (
    "data/processed/circos/genomic_circos_ring_values.tsv",
    "scripts/build_genomic_circos_inputs.py",
  ),
  "Table_S48_genomic_circos_track_definitions.tsv": (
    "data/processed/circos/genomic_circos_track_definitions.tsv",
    "scripts/build_genomic_circos_inputs.py",
  ),
  "Table_S49_genomic_circos_expression_summary.tsv": (
    "data/processed/circos/genomic_circos_expression_summary.tsv",
    "scripts/build_genomic_circos_inputs.py",
  ),
  "Table_S50_genomic_circos_expression_model_values.tsv": (
    "data/processed/circos/genomic_circos_expression_model_values.tsv",
    "scripts/build_genomic_circos_inputs.py",
  ),
  "Table_S51_pipeline_script_catalog.tsv": (
    "data/processed/documentation/pipeline_script_catalog.tsv",
    "scripts/build_script_documentation.py",
  ),
  "Table_S52_genomic_circos_source_provenance.tsv": (
    "data/processed/circos/genomic_circos_source_provenance.tsv",
    "scripts/build_genomic_circos_inputs.py",
  ),
}


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def require(path: Path) -> None:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(
      f"Mandatory genomic Circos asset is absent or empty: {path}"
    )


def atomic_copy(source: Path, destination: Path) -> None:
  destination.parent.mkdir(parents=True, exist_ok=True)
  temporary = destination.with_suffix(destination.suffix + ".tmp")
  shutil.copy2(source, temporary)
  temporary.replace(destination)


def append_figure_manifest(article_root: Path) -> None:
  combined_path = article_root / "manifests/figure_manifest.tsv"
  circos_path = (
    article_root / "manifests/genomic_circos_figure_manifest.tsv"
  )
  require(combined_path)
  require(circos_path)
  combined = pd.read_csv(combined_path, sep="\t", low_memory=False)
  circos = pd.read_csv(circos_path, sep="\t", low_memory=False)
  combined = pd.concat([
    combined.loc[
      ~combined["figure_id"].astype(str).eq("Figure_S70")
    ],
    circos,
  ], ignore_index=True, sort=False)
  if combined["figure_id"].duplicated().any():
    raise RuntimeError(
      "Duplicate figure IDs after registering Figure S70"
    )
  combined.to_csv(combined_path, sep="\t", index=False)


def register_tables(article_root: Path) -> None:
  manifest_path = article_root / "manifests/table_manifest.tsv"
  require(manifest_path)
  manifest = pd.read_csv(manifest_path, sep="\t", low_memory=False)
  supplementary = article_root / "tables/supplementary"
  source_copy = article_root / "source_data/tables"
  records = []
  for name, (source_value, script) in TABLE_SOURCES.items():
    source = resolve(source_value)
    require(source)
    destination = supplementary / name
    atomic_copy(source, destination)
    atomic_copy(source, source_copy / name)
    frame = pd.read_csv(source, sep="\t", low_memory=False)
    records.append({
      "table_id": Path(name).stem,
      "category": "supplementary",
      "path": str(destination),
      "rows": len(frame),
      "columns": len(frame.columns),
      "source_paths": str(source),
      "script": script,
      "status": (
        "ok" if not frame.empty else "empty_no_eligible_records"
      ),
    })
  new = pd.DataFrame(records)
  manifest = manifest.loc[
    ~manifest["table_id"].astype(str).isin(set(new["table_id"]))
  ]
  manifest = pd.concat([manifest, new], ignore_index=True, sort=False)
  if manifest["table_id"].duplicated().any():
    raise RuntimeError(
      "Duplicate table IDs after Circos table registration"
    )
  manifest.to_csv(manifest_path, sep="\t", index=False)


def update_reproduction_registry(article_root: Path) -> None:
  path = (
    article_root
    / "tables/supplementary/"
    "Table_S44_asset_reproduction_registry.tsv"
  )
  require(path)
  frame = pd.read_csv(path, sep="\t", low_memory=False)
  rows = [{
    "asset_id": "Figure_S70",
    "asset_type": "figure",
    "category": "supplementary",
    "scientific_title": (
      "Genomic Circos of NISE and homologous-paralog hypotheses"
    ),
    "script": "scripts/make_genomic_circos_figure.py",
    "primary_inputs": ";".join(
      value[0] for value in TABLE_SOURCES.values()
    ),
    "intermediate_data": (
      "article_outputs/source_data/figures/supplementary/"
      "Figure_S70_genomic_circos_rses_onco_source_data.tsv"
    ),
    "reproduction_command": (
      "MPLBACKEND=Agg python -u "
      "scripts/make_genomic_circos_figure.py "
      "--config config/genomic_circos_asset.yaml "
      "--output-root article_outputs --strict-layout"
    ),
    "outputs": ";".join(
      "article_outputs/figures/supplementary/"
      f"Figure_S70_genomic_circos_rses_onco.{extension}"
      for extension in ("png", "pdf", "svg")
    ),
    "dependencies": "environment.yml; pyproject.toml",
    "document_location": "supplementary document",
  }]
  for name, (source, script) in TABLE_SOURCES.items():
    table_id = Path(name).stem
    rows.append({
      "asset_id": table_id,
      "asset_type": "table",
      "category": "supplementary",
      "scientific_title": table_id.replace("_", " "),
      "script": script,
      "primary_inputs": source,
      "intermediate_data": source,
      "reproduction_command": (
        "bash scripts/run_publication_pipeline.sh assets-only"
      ),
      "outputs": f"article_outputs/tables/supplementary/{name}",
      "dependencies": "environment.yml; pyproject.toml",
      "document_location": "supplementary data",
    })
  additions = pd.DataFrame(rows)
  frame = frame.loc[
    ~frame["asset_id"].astype(str).isin(set(additions["asset_id"]))
  ]
  frame = pd.concat([frame, additions], ignore_index=True, sort=False)
  frame.to_csv(path, sep="\t", index=False)

  markdown = (
    article_root
    / "manuscript_assets/PUBLICATION_ASSET_REPRODUCTION.md"
  )
  lines = [
    "# Publication asset reproduction",
    "",
    "| Asset | Type | Category | Script | Command | Output |",
    "|---|---|---|---|---|---|",
  ]
  for row in frame.to_dict("records"):
    lines.append(
      f"| {row['asset_id']} | {row['asset_type']} | "
      f"{row['category']} | `{row['script']}` | "
      f"`{row['reproduction_command']}` | `{row['outputs']}` |"
    )
  markdown.parent.mkdir(parents=True, exist_ok=True)
  markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  article_root = resolve(args.article_root)
  append_figure_manifest(article_root)
  register_tables(article_root)
  update_reproduction_registry(article_root)
  print("Registered Figure S70 and Supplementary Tables S45-S52.")


if __name__ == "__main__":
  main()
