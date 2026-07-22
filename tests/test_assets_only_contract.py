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
  assert len(config["main_figures"]) == 7
  assert len(config["supplementary_figures"]) == 20
  assert len(config["structural_main_figures"]) == 1
  assert len(config["structural_supplementary_figures"]) == 18
  assert len(config["main_tables"]) == 4
  assert len(config["supplementary_tables"]) == 25
  expected_supplementary = {f"Figure_S{index}" for index in range(1, 39)}
  observed_supplementary = {
    str(record["id"])
    for section in (
      "supplementary_figures",
      "structural_supplementary_figures",
    )
    for record in config[section]
  }
  assert observed_supplementary == expected_supplementary
