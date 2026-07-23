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
    ROOT
    / "data/processed/structures/alphafold_structure_manifest.tsv"
  )
  render_manifest_path = (
    ROOT
    / "data/processed/structures/nise_structure_render_manifest.tsv"
  )
  annotation_path = (
    ROOT
    / "data/processed/structures/"
    "nise_structural_residue_annotations.tsv"
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
  proteins = pd.read_csv(
    proteins_path,
    sep="\t",
  ).drop_duplicates("uniprot_accession")
  structures = pd.read_csv(
    structure_manifest_path,
    sep="\t",
  )
  renders = pd.read_csv(render_manifest_path, sep="\t")
  expected_accessions = set(
    proteins["uniprot_accession"].astype(str)
  )
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
  ok_renders = renders.loc[
    renders["status"].astype(str).eq("ok")
  ].copy()
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
    errors.append(
      "individual_structure_render_count:expected_at_least_140"
    )


def validate_circos(
  article_root: Path,
  figures: pd.DataFrame,
  errors: list[str],
) -> None:
  coordinates_path = (
    article_root
    / "tables/supplementary/"
    "Table_S45_genomic_circos_gene_coordinates.tsv"
  )
  links_path = (
    article_root
    / "tables/supplementary/"
    "Table_S46_genomic_circos_pair_links.tsv"
  )
  rings_path = (
    article_root
    / "tables/supplementary/"
    "Table_S47_genomic_circos_ring_values.tsv"
  )
  tracks_path = (
    article_root
    / "tables/supplementary/"
    "Table_S48_genomic_circos_track_definitions.tsv"
  )
  expression_summary_path = (
    article_root
    / "tables/supplementary/"
    "Table_S49_genomic_circos_expression_summary.tsv"
  )
  expression_model_path = (
    article_root
    / "tables/supplementary/"
    "Table_S50_genomic_circos_expression_model_values.tsv"
  )
  script_catalog_path = (
    article_root
    / "tables/supplementary/"
    "Table_S51_pipeline_script_catalog.tsv"
  )
  provenance_path = (
    article_root
    / "tables/supplementary/"
    "Table_S52_genomic_circos_source_provenance.tsv"
  )
  required = [
    coordinates_path,
    links_path,
    rings_path,
    tracks_path,
    expression_summary_path,
    expression_model_path,
    script_catalog_path,
    provenance_path,
    ROOT / "docs/SCRIPT_CATALOG.md",
    ROOT / "docs/script_manifest.tsv",
    ROOT / "data/processed/circos/genomic_circos_status.json",
    article_root
    / "tables/figure_data/supplementary/Figure_S70_source_data.tsv",
  ]
  for path in required:
    require(path, errors)
  if any(not path.exists() for path in required[:8]):
    return

  coordinates = pd.read_csv(
    coordinates_path,
    sep="\t",
    low_memory=False,
  )
  links = pd.read_csv(links_path, sep="\t", low_memory=False)
  rings = pd.read_csv(rings_path, sep="\t", low_memory=False)
  tracks = pd.read_csv(tracks_path, sep="\t", low_memory=False)
  expression = pd.read_csv(
    expression_model_path,
    sep="\t",
    low_memory=False,
  )
  scripts = pd.read_csv(
    script_catalog_path,
    sep="\t",
    low_memory=False,
  )

  if coordinates.empty or coordinates["gene"].duplicated().any():
    errors.append("circos_gene_coordinates_empty_or_duplicated")
  if not coordinates.get(
    "coordinate_status",
    pd.Series("missing", index=coordinates.index),
  ).astype(str).eq("available").all():
    errors.append("circos_contains_gene_without_available_coordinate")
  expected_genes = set(coordinates["gene"].astype(str))
  linked_genes = set(links["lost_gene"].astype(str)) | set(
    links["target_gene"].astype(str)
  )
  if not linked_genes.issubset(expected_genes):
    errors.append("circos_links_reference_unmapped_genes")
  nise = links["pair_class"].astype(str).str.contains("NISE")
  paralog = links["pair_class"].astype(str).eq(
    "homologous_paralog"
  )
  if nise.any() and not links.loc[
    nise,
    "link_color",
  ].astype(str).eq("#C62828").all():
    errors.append("circos_nise_links_are_not_red")
  if paralog.any() and not links.loc[
    paralog,
    "link_color",
  ].astype(str).eq("#111111").all():
    errors.append("circos_paralog_links_are_not_black")

  expected_sources = {
    "coverage_adjusted_rses",
    "evidence_coverage",
    "component_tumor_event",
    "component_dependency",
    "component_selectivity",
    "component_expression_compensation",
    "component_functional_relation",
    "component_functional_microniche",
    "component_validation_tractability",
    "microniche_expression_context",
    "microniche_localization",
    "microniche_biochemical_structural",
    "microniche_genetic_phenotype",
    "microniche_interaction_network",
    "microniche_regulatory_network",
    "pairwise_expression_context",
    "wgcna_expression_network",
    "regulatory_tf_association_divergence",
    "regulatory_tf_expression_profile_divergence",
    "regulatory_promoter_motif_divergence",
    "component_promoter_methylation_context",
  }
  observed_sources = set(tracks["source_column"].astype(str))
  missing_sources = sorted(expected_sources - observed_sources)
  if missing_sources:
    errors.append(
      f"circos_missing_score_layers:{missing_sources}"
    )
  if set(rings["track_id"].astype(str)) != set(
    tracks["track_id"].astype(str)
  ):
    errors.append("circos_ring_track_definition_mismatch")
  expression_genes = set(expression["gene"].astype(str))
  if not expected_genes.issubset(expression_genes):
    errors.append("circos_model_expression_table_missing_genes")
  if not {
    "ModelID",
    "gene",
    "cancer",
    "expression_log2_tpm_plus_1",
    "source_file",
  }.issubset(expression.columns):
    errors.append("circos_expression_table_missing_required_columns")

  expected_script_paths = {
    path.relative_to(ROOT).as_posix()
    for directory in (ROOT / "scripts", ROOT / "src/rses_onco")
    for path in directory.rglob("*")
    if path.is_file()
    and path.suffix in {".py", ".sh", ".R", ".r"}
    and "__pycache__" not in path.parts
  }
  if set(scripts["script_path"].astype(str)) != expected_script_paths:
    errors.append("script_catalog_does_not_cover_every_pipeline_source")

  s70 = figures.loc[
    figures["figure_id"].astype(str).eq("Figure_S70")
  ]
  if len(s70) != 1:
    errors.append("figure_s70_registration_missing_or_duplicated")


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

  figures = pd.read_csv(manifest_path, sep="\t", low_memory=False)
  tables = pd.read_csv(
    table_manifest_path,
    sep="\t",
    low_memory=False,
  )
  if len(figures) != 78:
    errors.append(f"figure_count:{len(figures)}:expected_78")
  if figures["figure_id"].duplicated().any():
    errors.append("duplicate_figure_ids")
  expected_main = {f"Figure_{index}" for index in range(1, 9)}
  expected_supplementary = {
    f"Figure_S{index}" for index in range(1, 71)
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
    failed = figures.loc[
      ~figures.get(
        "layout_status",
        pd.Series("missing", index=figures.index),
      ).eq("pass"),
      [
        column
        for column in (
          "figure_id",
          "layout_status",
          "layout_warnings",
        )
        if column in figures
      ],
    ]
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

  if len(tables) != 56:
    errors.append(f"table_count:{len(tables)}:expected_56")
  main_count = int(
    tables["category"].astype(str).eq("main").sum()
  )
  supplementary_count = int(
    tables["category"].astype(str).eq("supplementary").sum()
  )
  if main_count != 4:
    errors.append(f"main_table_count:{main_count}:expected_4")
  if supplementary_count != 52:
    errors.append(
      "supplementary_table_count:"
      f"{supplementary_count}:expected_52"
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
    article_root
    / "tables/score_components/rses_onco_score_decomposition.tsv",
    article_root / "tables/robustness/leave_one_domain_out.tsv",
    article_root
    / "tables/robustness/controlled_weight_perturbation.tsv",
    article_root
    / "tables/figure_data/figure_source_data_inventory.tsv",
    article_root
    / "tables/supporting_evidence/supporting_evidence_manifest.tsv",
    article_root
    / "tables/supplementary/"
    "Table_S44_asset_reproduction_registry.tsv",
    ROOT
    / "data/processed/regulatory/wgcna/"
    "wgcna_correlation_fallback_all_cancers.tsv",
    ROOT
    / "data/processed/regulatory/wgcna/"
    "wgcna_run_diagnostics_all_cancers.tsv",
    article_root / "manifests/scientific_integrity_validation.json",
  ]
  for path in required_audit_outputs:
    require(path, errors)
    if path.is_relative_to(article_root):
      freshness_targets.append(path)

  validate_circos(article_root, figures, errors)

  workbook = (
    article_root
    / "workbooks/RSES_Onco_Article_Tables_and_Evidence.xlsx"
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
  require(
    article_root / "manuscript_assets/all_figure_legends.md",
    errors,
  )
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
    "Main figures: 8; supplementary figures: 70; "
    "exported image files: 234"
  )
  print("Main tables: 4; supplementary tables: 52")
  print(
    "All coordinate-complete NISE/paralog genes, links and RSES-Onco "
    "rings are registered in Figure S70."
  )
  print(
    "All model-level expression values used for Circos genes and every "
    "pipeline script/module are documented in supplementary tables."
  )
  print(
    "All registered figures passed automated layout audits."
  )


if __name__ == "__main__":
  main()
