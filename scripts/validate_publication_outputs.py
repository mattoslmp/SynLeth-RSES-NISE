#!/usr/bin/env python3
"""Validate completeness and reproducibility of the publication asset package."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def require(path: Path, errors: list[str], minimum_size: int = 1) -> None:
  if not path.exists():
    errors.append(f"missing:{path}")
  elif not path.is_file():
    errors.append(f"not_file:{path}")
  elif path.stat().st_size < minimum_size:
    errors.append(f"too_small:{path}:{path.stat().st_size}")


def validate_structural_atlas(article_root: Path, errors: list[str]) -> None:
  proteins_path = ROOT / "data/curated/human_nise_bonafide_2017.tsv"
  structure_manifest_path = ROOT / "data/processed/structures/alphafold_structure_manifest.tsv"
  render_manifest_path = ROOT / "data/processed/structures/nise_structure_render_manifest.tsv"
  annotation_path = ROOT / "data/processed/structures/nise_structural_residue_annotations.tsv"
  for path in (proteins_path, structure_manifest_path, render_manifest_path, annotation_path):
    require(path, errors)
  if any(not path.exists() for path in (
    proteins_path, structure_manifest_path, render_manifest_path
  )):
    return

  proteins = pd.read_csv(proteins_path, sep="\t").drop_duplicates("uniprot_accession")
  structures = pd.read_csv(structure_manifest_path, sep="\t")
  renders = pd.read_csv(render_manifest_path, sep="\t")
  expected_accessions = set(proteins["uniprot_accession"].astype(str))
  observed_structures = set(
    structures.loc[structures["status"].astype(str).eq("ok"), "uniprot_accession"].astype(str)
  )
  if len(expected_accessions) != 70:
    errors.append(f"curated_nise_protein_count:{len(expected_accessions)}:expected_70")
  missing_structures = expected_accessions - observed_structures
  if missing_structures:
    errors.append(f"missing_alphafold_structures:{sorted(missing_structures)}")

  ok_renders = renders.loc[renders["status"].astype(str).eq("ok")].copy()
  for accession in sorted(expected_accessions):
    accession_rows = ok_renders.loc[
      ok_renders["uniprot_accession"].astype(str).eq(accession)
    ]
    views = set(accession_rows["view"].astype(str))
    missing_views = {"whole", "site"} - views
    if missing_views:
      errors.append(f"missing_structural_views:{accession}:{sorted(missing_views)}")
    for render_path in accession_rows.get("render_path", pd.Series(dtype=object)).dropna():
      path = resolve_path(str(render_path))
      require(path, errors, minimum_size=1000)
      if path.exists():
        try:
          with Image.open(path) as image:
            width, height = image.size
          if width < 2600 or height < 2000:
            errors.append(f"structural_render_too_small:{path}:{width}x{height}")
        except Exception as exc:
          errors.append(f"invalid_structural_image:{path}:{exc}")

  structure_dir = article_root / "structure_atlas" / "individual"
  if not structure_dir.is_dir():
    errors.append(f"missing_structure_directory:{structure_dir}")
  else:
    render_count = sum(1 for _ in structure_dir.rglob("*.png"))
    if render_count < 140:
      errors.append(f"individual_structure_render_count:{render_count}:expected_at_least_140")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()

  article_root = resolve_path(args.article_root)
  errors: list[str] = []
  manifest_path = article_root / "manifests" / "figure_manifest.tsv"
  table_manifest_path = article_root / "manifests" / "table_manifest.tsv"
  require(manifest_path, errors)
  require(table_manifest_path, errors)
  if errors:
    raise SystemExit("\n".join(errors))

  figures = pd.read_csv(manifest_path, sep="\t")
  tables = pd.read_csv(table_manifest_path, sep="\t")
  if len(figures) != 40:
    errors.append(f"figure_count:{len(figures)}:expected_40")
  if figures["figure_id"].duplicated().any():
    errors.append("duplicate_figure_ids")
  expected_main = {f"Figure_{index}" for index in range(1, 9)}
  expected_supplementary = {f"Figure_S{index}" for index in range(1, 33)}
  observed = set(figures["figure_id"].astype(str))
  missing_main = expected_main - observed
  missing_supplementary = expected_supplementary - observed
  if missing_main:
    errors.append(f"missing_main_figures:{sorted(missing_main)}")
  if missing_supplementary:
    errors.append(f"missing_supplementary_figures:{sorted(missing_supplementary)}")
  if "layout_status" not in figures or not figures["layout_status"].eq("pass").all():
    failed = figures.loc[
      ~figures.get("layout_status", pd.Series("missing", index=figures.index)).eq("pass"),
      [column for column in ("figure_id", "layout_status", "layout_warnings") if column in figures],
    ]
    errors.append("layout_audit_failed:" + failed.to_json(orient="records"))

  for record in figures.to_dict("records"):
    base = Path(str(record["base_path"]))
    if not base.is_absolute():
      base = ROOT / base
    for extension in ("png", "pdf", "svg"):
      require(base.with_suffix(f".{extension}"), errors, minimum_size=100)
    audit_path = base.with_suffix(".layout_audit.json")
    require(audit_path, errors)
    if audit_path.exists():
      payload = json.loads(audit_path.read_text(encoding="utf-8"))
      if payload.get("status") != "pass" or payload.get("warnings"):
        errors.append(f"audit_warning:{record['figure_id']}:{payload.get('warnings')}")
    source_path = Path(str(record["source_data_path"]))
    if not source_path.is_absolute():
      source_path = ROOT / source_path
    require(source_path, errors)

  if len(tables) != 22:
    errors.append(f"table_count:{len(tables)}:expected_22")
  main_count = int(tables["category"].astype(str).eq("main").sum())
  supplementary_count = int(tables["category"].astype(str).eq("supplementary").sum())
  if main_count != 4:
    errors.append(f"main_table_count:{main_count}:expected_4")
  if supplementary_count != 18:
    errors.append(f"supplementary_table_count:{supplementary_count}:expected_18")
  for path_value in tables["path"].astype(str):
    path = Path(path_value)
    if not path.is_absolute():
      path = ROOT / path
    require(path, errors)

  require(
    article_root / "workbooks" / "RSES_Onco_Article_Tables_and_Evidence.xlsx",
    errors,
    minimum_size=1000,
  )
  require(article_root / "manifests" / "publication_file_inventory.tsv", errors)
  require(article_root / "manifests" / "publication_provenance.json", errors)
  require(article_root / "manifests" / "SHA256SUMS.txt", errors)
  require(article_root / "manuscript_assets" / "all_figure_legends.md", errors)
  validate_structural_atlas(article_root, errors)

  if errors:
    raise SystemExit(
      "Publication package validation failed:\n" + "\n".join(f"- {error}" for error in errors)
    )
  print("Publication package validation passed.")
  print("Main figures: 8; supplementary figures: 32; exported image files: 120")
  print("Main tables: 4; supplementary tables: 18")
  print("All 70 curated NISE proteins have whole and site structural renders.")
  print("All registered figures passed automated layout audits.")


if __name__ == "__main__":
  main()
