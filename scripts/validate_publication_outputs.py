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


def require(
  path: Path,
  errors: list[str],
  minimum_size: int = 1,
) -> None:
  if not path.exists():
    errors.append(f"missing:{path}")
  elif not path.is_file():
    errors.append(f"not_file:{path}")
  elif path.stat().st_size < minimum_size:
    errors.append(f"too_small:{path}:{path.stat().st_size}")


def validate_structural_atlas(
  article_root: Path,
  errors: list[str],
) -> None:
  proteins_path = ROOT / "data/curated/human_nise_bonafide_2017.tsv"
  structure_manifest_path = (
    ROOT / "data/processed/structures/alphafold_structure_manifest.tsv"
  )
  render_manifest_path = (
    ROOT / "data/processed/structures/nise_structure_render_manifest.tsv"
  )
  annotation_path = (
    ROOT / "data/processed/structures/nise_structural_residue_annotations.tsv"
  )
  for path in (
    proteins_path,
    structure_manifest_path,
    render_manifest_path,
    annotation_path,
  ):
    require(path, errors)
  if any(
    not path.exists()
    for path in (
      proteins_path,
      structure_manifest_path,
      render_manifest_path,
    )
  ):
    return
  proteins = pd.read_csv(proteins_path, sep="\t").drop_duplicates(
    "uniprot_accession"
  )
  structures = pd.read_csv(structure_manifest_path, sep="\t")
  renders = pd.read_csv(render_manifest_path, sep="\t")
  expected_accessions = set(proteins["uniprot_accession"].astype(str))
  observed_structures = set(
    structures.loc[
      structures["status"].astype(str).eq("ok"),
      "uniprot_accession",
    ].astype(str)
  )
  if len(expected_accessions) != 70:
    errors.append(
      f"curated_nise_protein_count:{len(expected_accessions)}:expected_70"
    )
  missing_structures = expected_accessions - observed_structures
  if missing_structures:
    errors.append(
      f"missing_alphafold_structures:{sorted(missing_structures)}"
    )
  ok_renders = renders.loc[renders["status"].astype(str).eq("ok")].copy()
  for accession in sorted(expected_accessions):
    accession_rows = ok_renders.loc[
      ok_renders["uniprot_accession"].astype(str).eq(accession)
    ]
    missing_views = {"whole", "site"} - set(
      accession_rows["view"].astype(str)
    )
    if missing_views:
      errors.append(
        f"missing_structural_views:{accession}:{sorted(missing_views)}"
      )
    for render_path in accession_rows.get(
      "render_path",
      pd.Series(dtype=object),
    ).dropna():
      path = resolve_path(str(render_path))
      require(path, errors, minimum_size=1000)
      if path.exists():
        try:
          with Image.open(path) as image:
            width, height = image.size
          if width < 2600 or height < 2000:
            errors.append(
              f"structural_render_too_small:{path}:{width}x{height}"
            )
        except Exception as exc:  # noqa: BLE001
          errors.append(f"invalid_structural_image:{path}:{exc}")
  structure_dir = article_root / "structure_atlas/individual"
  if not structure_dir.is_dir():
    errors.append(f"missing_structure_directory:{structure_dir}")
  elif sum(1 for _ in structure_dir.rglob("*.png")) < 140:
    errors.append("individual_structure_render_count:expected_at_least_140")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument(
    "--run-marker",
    default=None,
    help=(
      "Optional file written immediately before the publication run; "
      "registered outputs must be at least as new."
    ),
  )
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  run_marker = resolve_path(args.run_marker) if args.run_marker else None
  marker_mtime = (
    run_marker.stat().st_mtime
    if run_marker and run_marker.exists()
    else None
  )
  errors: list[str] = []
  freshness_targets: list[Path] = []
  manifest_path = article_root / "manifests/figure_manifest.tsv"
  table_manifest_path = article_root / "manifests/table_manifest.tsv"
  require(manifest_path, errors)
  require(table_manifest_path, errors)
  if errors:
    raise SystemExit("\n".join(errors))

  figures = pd.read_csv(manifest_path, sep="\t")
  tables = pd.read_csv(table_manifest_path, sep="\t")
  if len(figures) != 77:
    errors.append(f"figure_count:{len(figures)}:expected_77")
  if figures["figure_id"].duplicated().any():
    errors.append("duplicate_figure_ids")
  expected_main = {f"Figure_{index}" for index in range(1, 9)}
  expected_supplementary = {
    f"Figure_S{index}" for index in range(1, 70)
  }
  observed = set(figures["figure_id"].astype(str))
  if expected_main - observed:
    errors.append(
      f"missing_main_figures:{sorted(expected_main - observed)}"
    )
  if expected_supplementary - observed:
    errors.append(
      "missing_supplementary_figures:"
      f"{sorted(expected_supplementary - observed)}"
    )
  if (
    "layout_status" not in figures
    or not figures["layout_status"].eq("pass").all()
  ):
    status = figures.get(
      "layout_status",
      pd.Series("missing", index=figures.index),
    )
    columns = [
      column
      for column in (
        "figure_id",
        "layout_status",
        "layout_warnings",
      )
      if column in figures
    ]
    failed = figures.loc[~status.eq("pass"), columns]
    errors.append(
      "layout_audit_failed:" + failed.to_json(orient="records")
    )

  for record in figures.to_dict("records"):
    base = Path(str(record["base_path"]))
    base = base if base.is_absolute() else ROOT / base
    for extension in ("png", "pdf", "svg"):
      rendered = base.with_suffix(f".{extension}")
      require(rendered, errors, minimum_size=100)
      freshness_targets.append(rendered)
    audit_path = base.with_suffix(".layout_audit.json")
    require(audit_path, errors)
    freshness_targets.append(audit_path)
    if audit_path.exists():
      payload = json.loads(audit_path.read_text(encoding="utf-8"))
      if payload.get("status") != "pass" or payload.get("warnings"):
        errors.append(
          f"audit_warning:{record['figure_id']}:"
          f"{payload.get('warnings')}"
        )
    source_path = Path(str(record["source_data_path"]))
    source_path = (
      source_path if source_path.is_absolute() else ROOT / source_path
    )
    require(source_path, errors)
    freshness_targets.append(source_path)

  if len(tables) != 51:
    errors.append(f"table_count:{len(tables)}:expected_51")
  main_count = int(tables["category"].astype(str).eq("main").sum())
  supplementary_count = int(
    tables["category"].astype(str).eq("supplementary").sum()
  )
  if main_count != 4:
    errors.append(f"main_table_count:{main_count}:expected_4")
  if supplementary_count != 47:
    errors.append(
      f"supplementary_table_count:{supplementary_count}:expected_47"
    )
  for path_value in tables["path"].astype(str):
    path = Path(path_value)
    path = path if path.is_absolute() else ROOT / path
    require(path, errors)
    freshness_targets.append(path)

  required_audit_outputs = [
    article_root / "tables/qc/candidate_domain_evidence_audit.tsv",
    article_root / "tables/qc/coverage_by_domain.tsv",
    article_root / "tables/qc/coverage_by_source.tsv",
    article_root / "tables/qc/coverage_by_cancer.tsv",
    article_root / "tables/qc/coverage_by_mechanistic_class.tsv",
    article_root / "tables/qc/missingness_reasons.tsv",
    article_root / "tables/qc/evidence_overlap_registry.tsv",
    article_root / "tables/score_components/rses_onco_score_decomposition.tsv",
    article_root / "tables/robustness/leave_one_domain_out.tsv",
    article_root / "tables/robustness/controlled_weight_perturbation.tsv",
    article_root / "tables/figure_data/figure_source_data_inventory.tsv",
    article_root / "tables/supporting_evidence/supporting_evidence_manifest.tsv",
    article_root / "tables/supplementary/Table_S44_asset_reproduction_registry.tsv",
    article_root / "tables/supplementary/Table_S45_pair_promoter_methylation_evidence.tsv",
    article_root / "tables/supplementary/Table_S46_gene_promoter_methylation_summary.tsv",
    article_root / "tables/supplementary/Table_S47_methylation_source_status.tsv",
    ROOT / "data/processed/regulatory/wgcna/wgcna_correlation_fallback_all_cancers.tsv",
    ROOT / "data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv",
    article_root / "manifests/scientific_integrity_validation.json",
  ]
  for path in required_audit_outputs:
    require(path, errors)
    if path.is_relative_to(article_root):
      freshness_targets.append(path)

  methylation_figure = figures.loc[
    figures["figure_id"].astype(str).eq("Figure_S69")
  ]
  if len(methylation_figure) != 1:
    errors.append("figure_s69_methylation_registration_missing")
  else:
    source_value = str(methylation_figure.iloc[0]["source_data_path"])
    if "methylation" not in source_value.casefold():
      errors.append(
        "figure_s69_does_not_use_methylation_source_data:"
        + source_value
      )

  workbook = (
    article_root / "workbooks/RSES_Onco_Article_Tables_and_Evidence.xlsx"
  )
  require(workbook, errors, minimum_size=1000)
  freshness_targets.append(workbook)
  require(
    article_root / "manifests/publication_file_inventory.tsv",
    errors,
  )
  require(
    article_root / "manifests/publication_provenance.json",
    errors,
  )
  require(article_root / "manifests/SHA256SUMS.txt", errors)
  require(article_root / "manuscript_assets/all_figure_legends.md", errors)
  validate_structural_atlas(article_root, errors)

  if marker_mtime is not None:
    for path in sorted(set(freshness_targets)):
      if (
        path.exists()
        and path.stat().st_mtime + 1.0 < marker_mtime
      ):
        errors.append(
          f"stale_output:{path}:mtime={path.stat().st_mtime}:"
          f"marker={marker_mtime}"
        )

  if errors:
    raise SystemExit(
      "Publication package validation failed:\n"
      + "\n".join(f"- {error}" for error in errors)
    )
  print("Publication package validation passed.")
  print(
    "Main figures: 8; supplementary figures: 69; "
    "exported image files: 231"
  )
  print("Main tables: 4; supplementary tables: 47")
  print(
    "All 70 curated NISE proteins have whole and site structural renders."
  )
  print("All registered figures passed automated layout audits.")
  print(
    "Coverage, missingness, overlap, methylation, robustness and exact "
    "figure-source tables are present."
  )


if __name__ == "__main__":
  main()
