#!/usr/bin/env python3
"""Export the complete main and supplementary table package for the article."""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import shutil

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
    return pd.read_csv(path, sep="\t", low_memory=False)
  except EmptyDataError:
    return pd.DataFrame()


def write_table(
  frame: pd.DataFrame,
  path: Path,
  category: str,
  source_paths: list[Path],
  script: str = "scripts/export_article_tables.py",
) -> TableRecord:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  return TableRecord(
    path.stem,
    category,
    str(path),
    len(frame),
    len(frame.columns),
    ";".join(str(value) for value in source_paths),
    script,
    "ok" if not frame.empty else "empty_no_eligible_records",
  )


def register_existing_table(
  path: Path,
  category: str,
  source_paths: list[Path],
  script: str,
) -> TableRecord:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(
      f"Mandatory supplementary table is missing or empty: {path}"
    )
  frame = pd.read_csv(path, sep="\t", low_memory=False)
  return TableRecord(
    path.stem,
    category,
    str(path),
    len(frame),
    len(frame.columns),
    ";".join(str(value) for value in source_paths),
    script,
    "ok" if not frame.empty else "empty_no_eligible_records",
  )


def significant_dependency(frame: pd.DataFrame, fdr: float) -> pd.DataFrame:
  if frame.empty:
    return frame
  q_column = next(
    (
      column
      for column in (
        "q_value_bh_within_loss_cancer",
        "q_value_bh",
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
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--members",
    default="data/processed/expanded_class_member_inventory.tsv",
  )
  parser.add_argument(
    "--functional-evidence",
    default="data/processed/expanded_pair_functional_evidence.tsv",
  )
  parser.add_argument(
    "--dependency",
    default="results/expanded_26Q1/full/expanded_dependency_contrasts.tsv",
  )
  parser.add_argument(
    "--expression",
    default="results/expanded_26Q1/full/expanded_expression_compensation.tsv",
  )
  parser.add_argument(
    "--phenotype",
    default="results/expanded_26Q1/full/expanded_crispr_phenotype_profiles.tsv",
  )
  parser.add_argument(
    "--expression-context",
    default="results/expanded_26Q1/full/expanded_expression_context_profiles.tsv",
  )
  parser.add_argument(
    "--tcga-events",
    default="results/expanded_26Q1/full/tcga_gene_event_summary.tsv",
  )
  parser.add_argument(
    "--discovery",
    default="results/expanded_26Q1/discovery/all_target_dependency_screen.tsv",
  )
  parser.add_argument(
    "--pharmacology-evidence",
    default="data/processed/pharmacology/pharmacology_evidence_long.tsv",
  )
  parser.add_argument(
    "--pharmacology-ranking",
    default="results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv",
  )
  parser.add_argument(
    "--drug-sensitivity",
    default="data/processed/pharmacology/drug_response_selectivity.tsv",
  )
  parser.add_argument(
    "--pharmacology-source-status",
    default="data/processed/pharmacology/pharmacology_source_status.tsv",
  )
  parser.add_argument(
    "--pharmacology-source-coverage",
    default="results/expanded_26Q1/pharmacology/pharmacology_source_coverage.tsv",
  )
  parser.add_argument(
    "--structure-manifest",
    default="data/processed/structures/alphafold_structure_manifest.tsv",
  )
  parser.add_argument(
    "--structural-annotations",
    default="data/processed/structures/nise_structural_residue_annotations.tsv",
  )
  parser.add_argument(
    "--structure-render-manifest",
    default="data/processed/structures/nise_structure_render_manifest.tsv",
  )
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--top-n", type=int, default=20)
  parser.add_argument("--fdr", type=float, default=0.05)
  args = parser.parse_args()

  config = yaml.safe_load(
    resolve_path(args.config).read_text(encoding="utf-8")
  ) or {}
  configured_main = list(config.get("main_tables") or [])
  configured_supplementary = list(config.get("supplementary_tables") or [])
  if len(configured_main) != 4 or len(configured_supplementary) != 44:
    raise ValueError(
      "article_assets.yaml must define 4 main and exactly 44 core "
      "supplementary tables"
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
    "structure_manifest": resolve_path(args.structure_manifest),
    "structural_annotations": resolve_path(args.structural_annotations),
    "structure_render_manifest": resolve_path(args.structure_render_manifest),
  }
  for required in ("ranking", "candidates"):
    if not paths[required].exists() or paths[required].stat().st_size == 0:
      raise FileNotFoundError(
        f"Missing or empty mandatory table input: {paths[required]}"
      )

  ranking = pd.read_csv(paths["ranking"], sep="\t", low_memory=False)
  candidates = pd.read_csv(paths["candidates"], sep="\t", low_memory=False)
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
  structure_manifest = read_optional(paths["structure_manifest"])
  structural_annotations = read_optional(paths["structural_annotations"])
  structure_render_manifest = read_optional(paths["structure_render_manifest"])

  output_root = resolve_path(args.output_root)
  main_dir = output_root / "tables/main"
  supplementary_dir = output_root / "tables/supplementary"
  source_dir = output_root / "source_data/tables"
  manifest_dir = output_root / "manifests"
  records: list[TableRecord] = []

  class_summary = (
    ranking.assign(
      source_class=ranking.get(
        "source_class",
        pd.Series("unclassified", index=ranking.index),
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
  records.append(
    write_table(
      class_summary,
      main_dir / configured_main[0],
      "main",
      [paths["ranking"]],
    )
  )
  top_vulnerabilities = (
    ranking.sort_values(
      ["cancer", "coverage_adjusted_rses"],
      ascending=[True, False],
    )
    .groupby("cancer", group_keys=False)
    .head(args.top_n)
  )
  records.append(
    write_table(
      top_vulnerabilities,
      main_dir / configured_main[1],
      "main",
      [paths["ranking"]],
    )
  )
  records.append(
    write_table(
      significant_dependency(dependency, args.fdr),
      main_dir / configured_main[2],
      "main",
      [paths["dependency"]],
    )
  )
  if (
    not pharmacology_ranking.empty
    and "therapeutic_hypothesis_score" in pharmacology_ranking
  ):
    top_pharmacology = pharmacology_ranking.sort_values(
      "therapeutic_hypothesis_score",
      ascending=False,
    ).head(args.top_n * 3)
  else:
    top_pharmacology = pharmacology_ranking
  records.append(
    write_table(
      top_pharmacology,
      main_dir / configured_main[3],
      "main",
      [paths["pharmacology_ranking"]],
    )
  )

  component_columns = [
    column
    for column in ranking.columns
    if column.startswith("component_")
    or column.startswith("microniche_")
    or column.startswith("functional_microniche_")
    or column.startswith("methylation_")
    or column.startswith("expression_methylation_")
  ]
  component_keys = [
    column
    for column in (
      "cancer",
      "pair_id",
      "source_class",
      "lost_feature",
      "analysis_lost_gene",
      "lost_gene",
      "target_gene",
      "rses_onco",
      "evidence_coverage",
      "coverage_adjusted_rses",
      "priority_class",
    )
    if column in ranking
  ]
  component_matrix = ranking[component_keys + component_columns].copy()
  nise_directions = candidates.loc[
    candidates.get(
      "source_class",
      pd.Series(index=candidates.index, dtype=object),
    ).astype(str).eq("NISE")
  ].copy()
  network_columns = [
    column
    for column in functional.columns
    if column in {
      "pair_id",
      "lost_gene",
      "target_gene",
      "source_class",
      "component_localization",
      "component_biochemical_structural",
      "component_interaction_network",
      "string_direct_score",
      "string_neighbor_jaccard",
      "string_shared_neighbors",
      "component_regulatory_network",
      "regulator_jaccard",
      "shared_regulators",
    }
  ]
  network_evidence = (
    functional[network_columns].copy() if network_columns else functional
  )
  status_coverage = pd.concat(
    [
      pharmacology_status.assign(record_type="source_status"),
      pharmacology_coverage.assign(record_type="source_coverage"),
    ],
    ignore_index=True,
    sort=False,
  )

  base_frames = [
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
    (
      status_coverage,
      [paths["pharmacology_status"], paths["pharmacology_coverage"]],
    ),
    (structure_manifest, [paths["structure_manifest"]]),
    (structural_annotations, [paths["structural_annotations"]]),
    (structure_render_manifest, [paths["structure_render_manifest"]]),
  ]
  for name, (frame, source_paths) in zip(
    configured_supplementary[:18],
    base_frames,
  ):
    records.append(
      write_table(
        frame,
        supplementary_dir / name,
        "supplementary",
        source_paths,
      )
    )

  external_sources = {
    "Table_S19_coverage_domain_source_cancer_class.tsv": (
      "scripts/build_publication_evidence_audit.py",
      [output_root / "tables/qc/coverage_by_domain_cancer_class.tsv"],
    ),
    "Table_S20_score_decomposition_and_coverage.tsv": (
      "scripts/build_publication_evidence_audit.py",
      [output_root / "tables/score_components/rses_onco_score_decomposition.tsv"],
    ),
    "Table_S21_evidence_overlap_and_roles.tsv": (
      "scripts/build_publication_evidence_audit.py",
      [output_root / "tables/qc/evidence_overlap_registry.tsv"],
    ),
    "Table_S22_candidate_domain_missingness_audit.tsv": (
      "scripts/build_publication_evidence_audit.py",
      [output_root / "tables/qc/candidate_domain_evidence_audit.tsv"],
    ),
    "Table_S23_rses_ranking_stability.tsv": (
      "scripts/run_rses_robustness_analyses.py",
      [output_root / "tables/robustness/top_k_stability.tsv"],
    ),
    "Table_S24_leave_one_domain_out.tsv": (
      "scripts/run_rses_robustness_analyses.py",
      [output_root / "tables/robustness/leave_one_domain_out.tsv"],
    ),
    "Table_S25_weight_sensitivity.tsv": (
      "scripts/run_rses_robustness_analyses.py",
      [output_root / "tables/robustness/controlled_weight_perturbation.tsv"],
    ),
    "Table_S26_missing_data_sensitivity.tsv": (
      "scripts/run_rses_robustness_analyses.py",
      [output_root / "tables/robustness/missing_data_sensitivity.tsv"],
    ),
    "Table_S27_evidence_category_assignments.tsv": (
      "scripts/build_publication_evidence_audit_complete.py",
      [output_root / "tables/qc/evidence_category_assignments.tsv"],
    ),
    "Table_S28_coexpression_by_event_group.tsv": (
      "scripts/build_model_level_supporting_evidence.py",
      [output_root / "tables/supporting_evidence/expression/coexpression_by_event_group.tsv"],
    ),
    "Table_S29_compensation_and_dependency_contrasts.tsv": (
      "scripts/build_model_level_supporting_evidence.py",
      [output_root / "tables/supporting_evidence/expression/compensation_and_dependency_contrasts.tsv"],
    ),
    "Table_S30_model_level_expression_crispr_copy_number.tsv": (
      "scripts/build_model_level_supporting_evidence.py",
      [output_root / "tables/supporting_evidence/model_level/model_level_expression_crispr_copy_number.tsv"],
    ),
    "Table_S31_string_candidate_edges_all_channels.tsv": (
      "scripts/export_raw_functional_network_evidence.py",
      [output_root / "tables/supporting_evidence/networks/raw_sources/string_candidate_edges_all_channels.tsv"],
    ),
    "Table_S32_dorothea_candidate_regulatory_edges.tsv": (
      "scripts/export_raw_functional_network_evidence.py",
      [output_root / "tables/supporting_evidence/networks/raw_sources/dorothea_candidate_regulatory_edges.tsv"],
    ),
    "Table_S33_hpa_candidate_localization.tsv": (
      "scripts/export_raw_functional_network_evidence.py",
      [output_root / "tables/supporting_evidence/localization/hpa_candidate_localization.tsv"],
    ),
    "Table_S34_uniprot_candidate_annotations.tsv": (
      "scripts/export_raw_functional_network_evidence.py",
      [output_root / "tables/supporting_evidence/structures/uniprot_candidate_annotations.tsv"],
    ),
    "Table_S35_wgcna_pair_metrics_all_cancers.tsv": (
      "scripts/build_wgcna_regulatory_layer.py",
      [ROOT / "data/processed/regulatory/wgcna/wgcna_pair_metrics_all_cancers.tsv"],
    ),
    "Table_S36_wgcna_correlation_fallback_audit.tsv": (
      "scripts/run_wgcna_expression_network.R",
      [ROOT / "data/processed/regulatory/wgcna/wgcna_correlation_fallback_all_cancers.tsv"],
    ),
    "Table_S37_wgcna_run_diagnostics.tsv": (
      "scripts/run_wgcna_expression_network.R",
      [ROOT / "data/processed/regulatory/wgcna/wgcna_run_diagnostics_all_cancers.tsv"],
    ),
    "Table_S38_jaspar_promoter_tf_summary.tsv": (
      "scripts/scan_promoter_motifs.py",
      [ROOT / "data/processed/regulatory/jaspar_promoter_tf_summary.tsv"],
    ),
    "Table_S39_conditional_dependency_contrasts.tsv": (
      "scripts/build_model_level_supporting_evidence.py",
      [output_root / "tables/supporting_evidence/phenotypes/conditional_dependency_contrasts.tsv"],
    ),
    "Table_S40_tcga_gene_event_summary.tsv": (
      "scripts/export_supporting_evidence_tables.py",
      [output_root / "tables/supporting_evidence/genomic_context/tcga_gene_event_summary.tsv"],
    ),
    "Table_S41_pharmacology_evidence_long_support.tsv": (
      "scripts/export_supporting_evidence_tables.py",
      [output_root / "tables/supporting_evidence/pharmacology/pharmacology_evidence_long.tsv"],
    ),
    "Table_S42_pharmacology_source_coverage.tsv": (
      "scripts/prioritize_pharmacology.py",
      [ROOT / "results/expanded_26Q1/pharmacology/pharmacology_source_coverage.tsv"],
    ),
    "Table_S43_evidence_category_definitions.tsv": (
      "scripts/materialize_extended_supplementary_tables.py",
      [output_root / "tables/supplementary/Table_S43_evidence_category_definitions.tsv"],
    ),
    "Table_S44_asset_reproduction_registry.tsv": (
      "scripts/build_asset_reproduction_registry.py",
      [output_root / "tables/supplementary/Table_S44_asset_reproduction_registry.tsv"],
    ),
  }
  for name in configured_supplementary[18:]:
    script, source_paths = external_sources.get(name, ("unknown", []))
    records.append(
      register_existing_table(
        supplementary_dir / name,
        "supplementary",
        source_paths,
        script,
      )
    )

  methylation_sources = {
    "Table_S45_pair_promoter_methylation_evidence.tsv": (
      "scripts/build_methylation_pair_evidence.py",
      [ROOT / "data/processed/methylation/pair_promoter_methylation_evidence.tsv"],
    ),
    "Table_S46_gene_promoter_methylation_summary.tsv": (
      "scripts/aggregate_gdc_methylation.py",
      [ROOT / "data/processed/methylation/gdc_promoter_methylation_gene_summary.tsv"],
    ),
    "Table_S47_methylation_source_status.tsv": (
      "scripts/aggregate_gdc_methylation.py",
      [ROOT / "data/processed/methylation/gdc_promoter_methylation_aggregation_status.tsv"],
    ),
  }
  for name, (script, source_paths) in methylation_sources.items():
    records.append(
      register_existing_table(
        supplementary_dir / name,
        "supplementary",
        source_paths,
        script,
      )
    )

  source_dir.mkdir(parents=True, exist_ok=True)
  for record in records:
    table_path = Path(record.path)
    shutil.copy2(table_path, source_dir / table_path.name)

  manifest_dir.mkdir(parents=True, exist_ok=True)
  manifest = pd.DataFrame([asdict(record) for record in records])
  manifest.to_csv(
    manifest_dir / "table_manifest.tsv",
    sep="\t",
    index=False,
  )
  expected = len(configured_main) + len(configured_supplementary) + 3
  if len(records) != expected:
    raise RuntimeError(
      f"Expected {expected} article tables; observed {len(records)}"
    )
  print(
    manifest[["table_id", "category", "rows", "status"]].to_string(
      index=False
    )
  )
  print(f"Wrote all article tables to {output_root / 'tables'}")


if __name__ == "__main__":
  main()
