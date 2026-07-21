#!/usr/bin/env python3
"""Export the complete main and supplementary table package for the article."""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import shutil
from typing import Callable

import pandas as pd
from pandas.errors import EmptyDataError
import yaml

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class TableRecord:
  table_id: str
  category: str
  path: str
  rows: int
  columns: int
  source_paths: str
  script: str
  status: str


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists():
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t")
  except EmptyDataError:
    return pd.DataFrame()


def write_table(
  frame: pd.DataFrame,
  path: Path,
  category: str,
  source_paths: list[Path],
) -> TableRecord:
  path.parent.mkdir(parents=True, exist_ok=True)
  frame.to_csv(path, sep="\t", index=False)
  return TableRecord(
    table_id=path.stem,
    category=category,
    path=str(path),
    rows=len(frame),
    columns=len(frame.columns),
    source_paths=";".join(str(value) for value in source_paths),
    script="scripts/export_article_tables.py",
    status="ok" if not frame.empty else "empty_no_eligible_records",
  )


def significant_dependency(frame: pd.DataFrame, fdr: float) -> pd.DataFrame:
  if frame.empty:
    return frame
  q_column = next(
    (
      column for column in (
        "q_value_bh",
        "q_value_bh_within_loss_cancer",
      )
      if column in frame.columns
    ),
    None,
  )
  if q_column is None or "delta_effect" not in frame:
    return pd.DataFrame(columns=frame.columns)
  return frame.loc[
    (pd.to_numeric(frame[q_column], errors="coerce") < fdr)
    & (pd.to_numeric(frame["delta_effect"], errors="coerce") < 0)
  ].sort_values([q_column, "delta_effect"], ascending=[True, True])


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/article_assets.yaml")
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--candidates", default="data/processed/expanded_candidate_universe.tsv")
  parser.add_argument("--members", default="data/processed/expanded_class_member_inventory.tsv")
  parser.add_argument("--functional-evidence", default="data/processed/expanded_pair_functional_evidence.tsv")
  parser.add_argument("--dependency", default="results/expanded_26Q1/full/expanded_dependency_contrasts.tsv")
  parser.add_argument("--expression", default="results/expanded_26Q1/full/expanded_expression_compensation.tsv")
  parser.add_argument("--phenotype", default="results/expanded_26Q1/full/expanded_crispr_phenotype_profiles.tsv")
  parser.add_argument("--expression-context", default="results/expanded_26Q1/full/expanded_expression_context_profiles.tsv")
  parser.add_argument("--tcga-events", default="results/expanded_26Q1/full/tcga_gene_event_summary.tsv")
  parser.add_argument("--discovery", default="results/expanded_26Q1/discovery/all_target_dependency_screen.tsv")
  parser.add_argument("--pharmacology-evidence", default="data/processed/pharmacology/pharmacology_evidence_long.tsv")
  parser.add_argument("--pharmacology-ranking", default="results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv")
  parser.add_argument("--drug-sensitivity", default="data/processed/pharmacology/drug_response_selectivity.tsv")
  parser.add_argument("--pharmacology-source-status", default="data/processed/pharmacology/pharmacology_source_status.tsv")
  parser.add_argument("--pharmacology-source-coverage", default="results/expanded_26Q1/pharmacology/pharmacology_source_coverage.tsv")
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--top-n", type=int, default=20)
  parser.add_argument("--fdr", type=float, default=0.05)
  args = parser.parse_args()

  config = yaml.safe_load(resolve_path(args.config).read_text(encoding="utf-8")) or {}
  configured_main = list(config.get("main_tables") or [])
  configured_supplementary = list(config.get("supplementary_tables") or [])
  if len(configured_main) != 4 or len(configured_supplementary) != 15:
    raise ValueError(
      "config/article_assets.yaml must define 4 main and 15 supplementary tables"
    )

  paths = {
    "ranking": resolve_path(args.ranking),
    "candidates": resolve_path(args.candidates),
    "members": resolve_path(args.members),
    "functional": resolve_path(args.functional_evidence),
    "dependency": resolve_path(args.dependency),
    "expression": resolve_path(args.expression),
    "phenotype": resolve_path(args.phenotype),
    "expression_context": resolve_path(args.expression_context),
    "tcga": resolve_path(args.tcga_events),
    "discovery": resolve_path(args.discovery),
    "pharmacology_evidence": resolve_path(args.pharmacology_evidence),
    "pharmacology_ranking": resolve_path(args.pharmacology_ranking),
    "sensitivity": resolve_path(args.drug_sensitivity),
    "pharmacology_status": resolve_path(args.pharmacology_source_status),
    "pharmacology_coverage": resolve_path(args.pharmacology_source_coverage),
  }
  ranking = pd.read_csv(paths["ranking"], sep="\t")
  candidates = pd.read_csv(paths["candidates"], sep="\t")
  members = read_optional(paths["members"])
  functional = read_optional(paths["functional"])
  dependency = read_optional(paths["dependency"])
  expression = read_optional(paths["expression"])
  phenotype = read_optional(paths["phenotype"])
  expression_context = read_optional(paths["expression_context"])
  tcga = read_optional(paths["tcga"])
  discovery = read_optional(paths["discovery"])
  pharmacology_evidence = read_optional(paths["pharmacology_evidence"])
  pharmacology_ranking = read_optional(paths["pharmacology_ranking"])
  sensitivity = read_optional(paths["sensitivity"])
  pharmacology_status = read_optional(paths["pharmacology_status"])
  pharmacology_coverage = read_optional(paths["pharmacology_coverage"])

  output_root = resolve_path(args.output_root)
  main_dir = output_root / "tables" / "main"
  supplementary_dir = output_root / "tables" / "supplementary"
  source_dir = output_root / "source_data" / "tables"
  manifest_dir = output_root / "manifests"
  records: list[TableRecord] = []

  class_summary = (
    ranking.assign(
      source_class=ranking.get(
        "source_class", pd.Series("unclassified", index=ranking.index)
      ).fillna("unclassified").astype(str)
    )
    .groupby("source_class", as_index=False)
    .agg(
      scored_rows=("pair_id", "size"),
      unique_directions=("pair_id", "nunique"),
      unique_targets=("target_gene", "nunique"),
      maximum_adjusted_rses=("coverage_adjusted_rses", "max"),
      median_adjusted_rses=("coverage_adjusted_rses", "median"),
      median_coverage=("evidence_coverage", "median"),
    )
    .sort_values("maximum_adjusted_rses", ascending=False)
  )
  records.append(write_table(
    class_summary,
    main_dir / configured_main[0],
    "main",
    [paths["ranking"]],
  ))

  top_vulnerabilities = (
    ranking.sort_values(
      ["cancer", "coverage_adjusted_rses"],
      ascending=[True, False],
    )
    .groupby("cancer", group_keys=False)
    .head(args.top_n)
  )
  records.append(write_table(
    top_vulnerabilities,
    main_dir / configured_main[1],
    "main",
    [paths["ranking"]],
  ))

  significant = significant_dependency(dependency, args.fdr)
  records.append(write_table(
    significant,
    main_dir / configured_main[2],
    "main",
    [paths["dependency"]],
  ))

  top_pharmacology = (
    pharmacology_ranking.sort_values(
      "therapeutic_hypothesis_score", ascending=False
    ).head(args.top_n * 3)
    if not pharmacology_ranking.empty
    and "therapeutic_hypothesis_score" in pharmacology_ranking
    else pharmacology_ranking
  )
  records.append(write_table(
    top_pharmacology,
    main_dir / configured_main[3],
    "main",
    [paths["pharmacology_ranking"]],
  ))

  component_columns = [
    column for column in ranking.columns
    if column.startswith("component_")
    or column.startswith("microniche_")
    or column.startswith("functional_microniche_")
  ]
  component_keys = [
    column for column in (
      "cancer", "pair_id", "source_class", "lost_feature", "analysis_lost_gene",
      "lost_gene", "target_gene", "rses_onco", "evidence_coverage",
      "coverage_adjusted_rses", "priority_class",
    )
    if column in ranking
  ]
  component_matrix = ranking[component_keys + component_columns].copy()
  nise_directions = candidates.loc[
    candidates.get("source_class", pd.Series(index=candidates.index, dtype=object))
      .astype(str).eq("NISE")
  ].copy()
  network_columns = [
    column for column in functional.columns
    if column in {
      "pair_id", "lost_gene", "target_gene", "source_class",
      "component_localization", "component_biochemical_structural",
      "component_interaction_network", "string_direct_score",
      "string_neighbor_jaccard", "string_shared_neighbors",
      "component_regulatory_network", "regulator_jaccard",
      "shared_regulators",
    }
  ]
  network_evidence = functional[network_columns].copy() if network_columns else functional
  status_coverage = pd.concat([
    pharmacology_status.assign(record_type="source_status"),
    pharmacology_coverage.assign(record_type="source_coverage"),
  ], ignore_index=True, sort=False)

  supplementary_frames = [
    (candidates, [paths["candidates"]]),
    (nise_directions, [paths["candidates"]]),
    (members, [paths["members"]]),
    (component_matrix, [paths["ranking"]]),
    (dependency, [paths["dependency"]]),
    (expression, [paths["expression"]]),
    (phenotype, [paths["phenotype"]]),
    (expression_context, [paths["expression_context"]]),
    (network_evidence, [paths["functional"]]),
    (tcga, [paths["tcga"]]),
    (discovery, [paths["discovery"]]),
    (pharmacology_evidence, [paths["pharmacology_evidence"]]),
    (pharmacology_ranking, [paths["pharmacology_ranking"]]),
    (sensitivity, [paths["sensitivity"]]),
    (status_coverage, [paths["pharmacology_status"], paths["pharmacology_coverage"]]),
  ]
  for name, (frame, source_paths) in zip(configured_supplementary, supplementary_frames):
    records.append(write_table(frame, supplementary_dir / name, "supplementary", source_paths))

  source_dir.mkdir(parents=True, exist_ok=True)
  for record in records:
    table_path = Path(record.path)
    destination = source_dir / table_path.name
    shutil.copy2(table_path, destination)

  manifest_dir.mkdir(parents=True, exist_ok=True)
  manifest = pd.DataFrame([asdict(record) for record in records])
  manifest.to_csv(manifest_dir / "table_manifest.tsv", sep="\t", index=False)
  if len(records) != 19:
    raise RuntimeError(f"Expected 19 article tables; observed {len(records)}")
  print(manifest[["table_id", "category", "rows", "status"]].to_string(index=False))
  print(f"Wrote all article tables to {output_root / 'tables'}")


if __name__ == "__main__":
  main()
