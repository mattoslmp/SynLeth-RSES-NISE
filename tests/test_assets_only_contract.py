from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_assets_only_entrypoint_includes_extended_evidence() -> None:
  entrypoint = (
    ROOT / "scripts/run_publication_pipeline.sh"
  ).read_text(encoding="utf-8")
  required = {
    "build_model_level_supporting_evidence.py",
    "export_raw_functional_network_evidence.py",
    "assets-only",
  }
  missing = sorted(value for value in required if value not in entrypoint)
  assert not missing, f"assets-only entrypoint missing: {missing}"


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
  }
  missing = sorted(value for value in required if value not in core)
  assert not missing, f"publication core missing stages: {missing}"


def test_publication_registry_matches_expanded_contract() -> None:
  config = yaml.safe_load(
    (ROOT / "config/article_assets.yaml").read_text(encoding="utf-8")
  )
  main_ids = {str(record["id"]) for record in config["main_figures"]}
  supplementary_ids = {
    str(record["id"]) for record in config["supplementary_figures"]
  }
  assert main_ids == {f"Figure_{index}" for index in range(1, 9)}
  assert supplementary_ids == {f"Figure_S{index}" for index in range(1, 39)}
  assert len(config["main_tables"]) == 4
  assert len(config["supplementary_tables"]) == 25
