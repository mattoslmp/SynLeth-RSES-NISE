from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_canonical_entrypoint_routes_to_complete_workflow() -> None:
  entrypoint = (
    ROOT / "scripts/run_publication_pipeline.sh"
  ).read_text(encoding="utf-8")
  assert "publication_pipeline_complete.sh" in entrypoint
  assert "assets-only" in entrypoint


def test_complete_wrapper_executes_extended_evidence_and_finalization() -> None:
  wrapper = (
    ROOT / "scripts/publication_pipeline_complete.sh"
  ).read_text(encoding="utf-8")
  required = {
    "build_model_level_supporting_evidence.py",
    "export_raw_functional_network_evidence.py",
    "export_wgcna_regulatory_supporting_evidence.py",
    "run_wgcna_regulatory_ablation.py",
    "validate_extended_supporting_evidence.py",
    "validate_wgcna_regulatory_evidence.py",
    "build_publication_methods_documentation.py",
    "build_wgcna_regulatory_methods.py",
    "create_manual_visual_inspection_checklist.py",
    'bash "$CORE" "$stage"',
    'bash "$CORE" workbook',
    'bash "$CORE" manifests',
    'bash "$CORE" validate',
  }
  missing = sorted(
    value
    for value in required
    if value not in wrapper
  )
  assert not missing, (
    f"complete publication wrapper missing: {missing}"
  )


def test_assets_only_core_includes_all_publication_stages() -> None:
  core = (
    ROOT / "scripts/publication_pipeline_steps.sh"
  ).read_text(encoding="utf-8")
  required = {
    "build_publication_evidence_audit_complete.py",
    "run_rses_robustness_analyses.py",
    "export_supporting_evidence_tables.py",
    "export_article_tables.py",
    "make_all_article_figures.py",
    "catalog_figure_source_data.py",
    "validate_publication_scientific_integrity.py",
    "build_article_workbook.py",
    "build_publication_manifest.py",
    "validate_publication_outputs.py",
    "materialize_extended_supplementary_tables.py",
    "build_asset_reproduction_registry.py",
    "build_publication_documents.py",
    "validate_publication_documents.py",
  }
  missing = sorted(
    value
    for value in required
    if value not in core
  )
  assert not missing, (
    f"publication core missing stages: {missing}"
  )


def test_publication_registry_matches_expanded_contract() -> None:
  config = yaml.safe_load(
    (ROOT / "config/article_assets.yaml").read_text(
      encoding="utf-8"
    )
  )
  main_ids = {
    str(record["id"])
    for record in config["main_figures"]
  }
  supplementary_ids = {
    str(record["id"])
    for record in config["supplementary_figures"]
  }
  assert main_ids == {
    f"Figure_{index}"
    for index in range(1, 9)
  }
  assert supplementary_ids == {
    f"Figure_S{index}"
    for index in range(1, 70)
  }
  assert len(config["main_tables"]) == 4
  assert len(config["supplementary_tables"]) == 44
