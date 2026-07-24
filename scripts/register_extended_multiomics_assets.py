#!/usr/bin/env python3
"""Register Figures S71-S78 and Supplementary Tables S53-S64."""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import yaml


TABLE_SOURCES = {
  "Table_S53_extended_multiomics_source_inventory.tsv": (
    "data/processed/extended_multiomics/extended_multiomics_source_inventory.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S54_extended_multiomics_source_status.tsv": (
    "data/processed/extended_multiomics/extended_multiomics_source_status.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S55_functional_loss_states.tsv": (
    "data/processed/extended_multiomics/functional_loss_states.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S56_extended_pair_evidence_by_cancer.tsv": (
    "data/processed/extended_multiomics/extended_pair_evidence_by_cancer.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S57_proteomics_pair_evidence_by_source.tsv": (
    "data/processed/extended_multiomics/proteomics_pair_evidence_by_source.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S58_extended_covariate_context.tsv": (
    "data/processed/extended_multiomics/extended_covariate_context.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S59_custom_drug_sensitivity_long.tsv": (
    "data/processed/extended_multiomics/custom_drug_sensitivity_long.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S60_gdsc_combination_evidence_long.tsv": (
    "data/processed/extended_multiomics/gdsc_combination_evidence_long.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S61_extended_multiomics_source_provenance.tsv": (
    "data/processed/extended_multiomics/extended_multiomics_source_provenance.tsv",
    "scripts/build_extended_multiomics_evidence.py",
  ),
  "Table_S62_extended_rses_onco_complete_ranking.tsv": (
    "results/expanded_26Q1/full/expanded_rses_onco.tsv",
    "scripts/recompute_rses_with_extended_multiomics.py",
  ),
}


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def require(path: Path) -> None:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(
      f"Mandatory extended multi-omics asset is absent or empty: {path}"
    )


def atomic_copy(source: Path, destination: Path) -> None:
  destination.parent.mkdir(parents=True, exist_ok=True)
  temporary = destination.with_suffix(destination.suffix + ".tmp")
  shutil.copy2(source, temporary)
  temporary.replace(destination)


def append_figure_manifest(article_root: Path) -> None:
  combined_path = article_root / "manifests/figure_manifest.tsv"
  extended_path = (
    article_root / "manifests/extended_multiomics_figure_manifest.tsv"
  )
  require(combined_path)
  require(extended_path)
  combined = pd.read_csv(combined_path, sep="\t", low_memory=False)
  extended = pd.read_csv(extended_path, sep="\t", low_memory=False)
  expected = {f"Figure_S{index}" for index in range(71, 79)}
  if set(extended["figure_id"].astype(str)) != expected:
    raise RuntimeError(
      "Extended figure manifest must contain exactly Figure S71-S78"
    )
  combined = combined.loc[
    ~combined["figure_id"].astype(str).isin(expected)
  ]
  combined = pd.concat([combined, extended], ignore_index=True, sort=False)
  if combined["figure_id"].duplicated().any():
    raise RuntimeError(
      "Duplicate figure IDs after extended multi-omics registration"
    )
  combined.to_csv(combined_path, sep="\t", index=False)


def derived_tables(ranking: pd.DataFrame) -> dict[str, pd.DataFrame]:
  comparison_columns = [
    column
    for column in (
      "pair_id",
      "cancer",
      "lost_gene",
      "target_gene",
      "source_class",
      "baseline_coverage_adjusted_rses",
      "coverage_adjusted_rses",
      "extended_score_delta",
      "baseline_rank_within_cancer",
      "extended_rank_within_cancer",
      "extended_rank_change",
      "extended_scored_layer_count",
      "extended_scored_layer_coverage",
      "score_version",
    )
    if column in ranking.columns
  ]
  comparison = ranking[comparison_columns].copy()
  definitions = pd.DataFrame([
    {
      "layer": "integrated_functional_loss_support",
      "domain": "tumor_event",
      "primary_score": True,
      "scientific_rule": (
        "Combines relative/absolute copy-number loss and LOH plus "
        "damaging mutation without double-counting overlapping events."
      ),
    },
    {
      "layer": "dependency_probability_support",
      "domain": "dependency",
      "primary_score": True,
      "scientific_rule": (
        "Uses CRISPR Gene Dependency as confidence within the existing "
        "Chronos dependency domain, not as an independent domain."
      ),
    },
    {
      "layer": "protein_compensation_support",
      "domain": "expression_compensation",
      "primary_score": True,
      "scientific_rule": (
        "Normalizes each protein platform separately and forms a "
        "coverage-aware cross-platform consensus."
      ),
    },
    {
      "layer": "rnai_orthogonal_support",
      "domain": "functional_microniche.genetic_phenotype",
      "primary_score": True,
      "scientific_rule": (
        "Uses DEMETER2 as orthogonal perturbation support while "
        "preventing a second full dependency-domain weight."
      ),
    },
    {
      "layer": (
        "metabolomics/miRNA/chromatin/ssGSEA/subtypes/signatures/MetMap"
      ),
      "domain": "context_and_robustness",
      "primary_score": False,
      "scientific_rule": (
        "Retained for interpretation, confounder control and validation; "
        "not directly scored without gene/reaction or causal mapping."
      ),
    },
    {
      "layer": "single and combination drug response",
      "domain": "translation_validation",
      "primary_score": False,
      "scientific_rule": (
        "Reported separately from discovery score to avoid "
        "pharmacology-driven circular prioritization."
      ),
    },
  ])
  return {
    "Table_S63_baseline_extended_score_comparison.tsv": comparison,
    "Table_S64_extended_layer_scientific_definitions.tsv": definitions,
  }


def register_tables(article_root: Path) -> None:
  manifest_path = article_root / "manifests/table_manifest.tsv"
  require(manifest_path)
  manifest = pd.read_csv(manifest_path, sep="\t", low_memory=False)
  supplementary = article_root / "tables/supplementary"
  source_copy = article_root / "source_data/tables"
  records = []
  table_frames: dict[str, pd.DataFrame] = {}
  for name, (source_value, script) in TABLE_SOURCES.items():
    source = resolve(source_value)
    require(source)
    frame = pd.read_csv(source, sep="\t", low_memory=False)
    table_frames[name] = frame
    destination = supplementary / name
    atomic_copy(source, destination)
    atomic_copy(source, source_copy / name)
    records.append({
      "table_id": Path(name).stem,
      "category": "supplementary",
      "path": str(destination),
      "rows": len(frame),
      "columns": len(frame.columns),
      "source_paths": str(source),
      "script": script,
      "status": "ok" if not frame.empty else "empty_no_evaluable_records",
    })
  ranking = table_frames["Table_S62_extended_rses_onco_complete_ranking.tsv"]
  for name, frame in derived_tables(ranking).items():
    destination = supplementary / name
    source = source_copy / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(destination, sep="\t", index=False)
    frame.to_csv(source, sep="\t", index=False)
    records.append({
      "table_id": Path(name).stem,
      "category": "supplementary",
      "path": str(destination),
      "rows": len(frame),
      "columns": len(frame.columns),
      "source_paths": str(resolve(
        "results/expanded_26Q1/full/expanded_rses_onco.tsv"
      )),
      "script": "scripts/register_extended_multiomics_assets.py",
      "status": "ok" if not frame.empty else "empty_no_evaluable_records",
    })
  new = pd.DataFrame(records)
  manifest = manifest.loc[
    ~manifest["table_id"].astype(str).isin(set(new["table_id"]))
  ]
  manifest = pd.concat([manifest, new], ignore_index=True, sort=False)
  if manifest["table_id"].duplicated().any():
    raise RuntimeError(
      "Duplicate table IDs after extended multi-omics registration"
    )
  manifest.to_csv(manifest_path, sep="\t", index=False)


def update_reproduction_registry(article_root: Path) -> None:
  path = (
    article_root
    / "tables/supplementary/Table_S44_asset_reproduction_registry.tsv"
  )
  require(path)
  frame = pd.read_csv(path, sep="\t", low_memory=False)
  rows = []
  asset_config = yaml.safe_load(
    resolve("config/extended_multiomics_asset.yaml")
      .read_text(encoding="utf-8")
  ) or {}
  figure_registry = {
    str(item["id"]): item
    for item in asset_config.get("supplementary_figures", [])
  }
  for index in range(71, 79):
    figure_id = f"Figure_S{index}"
    item = figure_registry[figure_id]
    file_name = str(item["file"])
    rows.append({
      "asset_id": figure_id,
      "asset_type": "figure",
      "category": "supplementary",
      "scientific_title": str(item["title"]),
      "script": "scripts/make_extended_multiomics_figures.py",
      "primary_inputs": (
        "data/processed/extended_multiomics;"
        "results/expanded_26Q1/full/expanded_rses_onco.tsv"
      ),
      "intermediate_data": (
        "article_outputs/source_data/figures/supplementary/"
        f"{file_name}_source_data.tsv"
      ),
      "reproduction_command": (
        "MPLBACKEND=Agg python -u "
        "scripts/make_extended_multiomics_figures.py --strict-layout"
      ),
      "outputs": ";".join(
        f"article_outputs/figures/supplementary/{file_name}." + extension
        for extension in ("png", "pdf", "svg")
      ),
      "dependencies": "environment.yml;pyproject.toml",
      "document_location": "supplementary document",
    })
  for name, (source, script) in TABLE_SOURCES.items():
    rows.append({
      "asset_id": Path(name).stem,
      "asset_type": "table",
      "category": "supplementary",
      "scientific_title": Path(name).stem.replace("_", " "),
      "script": script,
      "primary_inputs": source,
      "intermediate_data": source,
      "reproduction_command": (
        "bash scripts/run_extended_multiomics_pipeline.sh publication"
      ),
      "outputs": f"article_outputs/tables/supplementary/{name}",
      "dependencies": "environment.yml;pyproject.toml",
      "document_location": "supplementary data",
    })
  for name in (
    "Table_S63_baseline_extended_score_comparison.tsv",
    "Table_S64_extended_layer_scientific_definitions.tsv",
  ):
    rows.append({
      "asset_id": Path(name).stem,
      "asset_type": "table",
      "category": "supplementary",
      "scientific_title": Path(name).stem.replace("_", " "),
      "script": "scripts/register_extended_multiomics_assets.py",
      "primary_inputs": (
        "results/expanded_26Q1/full/expanded_rses_onco.tsv"
      ),
      "intermediate_data": (
        "results/expanded_26Q1/full/expanded_rses_onco.tsv"
      ),
      "reproduction_command": (
        "bash scripts/run_extended_multiomics_pipeline.sh publication"
      ),
      "outputs": f"article_outputs/tables/supplementary/{name}",
      "dependencies": "environment.yml;pyproject.toml",
      "document_location": "supplementary data",
    })
  additions = pd.DataFrame(rows)
  frame = frame.loc[
    ~frame["asset_id"].astype(str).isin(set(additions["asset_id"]))
  ]
  frame = pd.concat([frame, additions], ignore_index=True, sort=False)
  frame.to_csv(path, sep="\t", index=False)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  article_root = resolve(args.article_root)
  append_figure_manifest(article_root)
  register_tables(article_root)
  update_reproduction_registry(article_root)
  print("Registered Figures S71-S78 and Supplementary Tables S53-S64.")


if __name__ == "__main__":
  main()
