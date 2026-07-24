#!/usr/bin/env python3
"""Build scored and exploratory extended multi-omics evidence for RSES-Onco."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
import yaml

from rses_onco.extended_multiomics import (
  DIRECT_SCORE_SOURCE_KEYS,
  SourceSpec,
  bh_adjust,
  build_functional_loss_table,
  build_gdsc_combination_table,
  build_global_context_evidence,
  build_pair_evidence,
  read_long_event_matrix,
  read_model_feature_matrix,
  read_table,
  sha256_file,
  source_inventory,
  standardize_matrix_long,
  write_status_json,
)


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  frame.to_csv(path, sep="\t", index=False)


def merge_binary(
  first: pd.DataFrame | None,
  second: pd.DataFrame | None,
) -> pd.DataFrame | None:
  if first is None or first.empty:
    return second
  if second is None or second.empty:
    return first
  rows = sorted(
    set(first.index.astype(str)) | set(second.index.astype(str))
  )
  columns = sorted(
    set(first.columns.astype(str)) | set(second.columns.astype(str))
  )
  left = first.reindex(index=rows, columns=columns)
  right = second.reindex(index=rows, columns=columns)
  return (
    pd.concat([left, right], axis=0)
      .groupby(level=0)
      .max()
      .reindex(columns=columns)
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--config", default="config/extended_multiomics_sources.yaml"
  )
  parser.add_argument("--input-dir", default="dmap_data")
  parser.add_argument("--models", default="data/raw/depmap/Model.csv")
  parser.add_argument(
    "--copy-number", default="data/raw/depmap/OmicsCNGeneWGS.csv"
  )
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--output-dir", default="data/processed/extended_multiomics"
  )
  parser.add_argument("--min-group-size", type=int, default=3)
  parser.add_argument("--loss-threshold", type=float, default=0.30)
  parser.add_argument("--strict", action="store_true")
  args = parser.parse_args()

  config_path = resolve(args.config)
  config = yaml.safe_load(
    config_path.read_text(encoding="utf-8")
  ) or {}
  input_dir = resolve(args.input_dir)
  output_dir = resolve(args.output_dir)
  output_dir.mkdir(parents=True, exist_ok=True)

  models_path = resolve(args.models)
  copy_number_path = resolve(args.copy_number)
  ranking_path = resolve(args.ranking)
  models = read_table(models_path)
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  candidates = ranking[[
    column
    for column in ("pair_id", "lost_gene", "target_gene", "cancer")
    if column in ranking.columns
  ]].drop_duplicates()

  source_config = config.get("sources") or {}
  specs = [
    SourceSpec(
      key=str(key),
      filename=str(value["filename"]),
      role=str(value.get("role") or "unspecified"),
      required=bool(value.get("required", False)),
    )
    for key, value in source_config.items()
  ]
  inventory = source_inventory(specs, input_dir)
  write_tsv(
    inventory,
    output_dir / "extended_multiomics_source_inventory.tsv",
  )

  statuses: list[dict[str, object]] = []
  matrices: dict[str, pd.DataFrame] = {}

  def load_matrix(
    key: str,
    *,
    gene_features: bool | None = None,
  ) -> pd.DataFrame | None:
    specification = source_config.get(key) or {}
    path = input_dir / str(specification.get("filename") or "")
    if not path.exists() or path.stat().st_size == 0:
      statuses.append({
        "source_key": key,
        "status": "not_found",
        "path": str(path),
      })
      if args.strict and specification.get("required"):
        raise FileNotFoundError(path)
      return None
    try:
      layout = str(specification.get("layout") or "matrix")
      if layout == "long_mutation":
        matrix = read_long_event_matrix(path, models, event="mutation")
      elif layout == "long_fusion":
        matrix = read_long_event_matrix(path, models, event="fusion")
      else:
        matrix = read_model_feature_matrix(
          path,
          models,
          gene_features=(
            bool(specification.get("gene_features", False))
            if gene_features is None
            else gene_features
          ),
        )
      if matrix.empty and args.strict and key in DIRECT_SCORE_SOURCE_KEYS:
        raise ValueError(
          f"Direct score source {key} is empty after standardization"
        )
      matrices[key] = matrix
      statuses.append({
        "source_key": key,
        "status": (
          "ok" if not matrix.empty else "empty_after_standardization"
        ),
        "path": str(path),
        "rows": len(matrix),
        "columns": len(matrix.columns),
        "sha256": sha256_file(path),
      })
      return matrix
    except Exception as exc:
      statuses.append({
        "source_key": key,
        "status": "failed",
        "path": str(path),
        "message": str(exc),
      })
      if args.strict and (
        key in DIRECT_SCORE_SOURCE_KEYS
        or bool(specification.get("required", False))
      ):
        raise
      return None

  relative_cn = read_model_feature_matrix(
    copy_number_path,
    models,
    gene_features=True,
  )
  dependency_probability = load_matrix("crispr_dependency")
  damaging = load_matrix("damaging_mutations")
  mutation_table = load_matrix("mutation_table")
  damaging = merge_binary(damaging, mutation_table)
  hotspot = load_matrix("hotspot_mutations")
  fusions = load_matrix("fusions")
  absolute_cn = load_matrix("absolute_cn")
  loh = load_matrix("loh")
  rnai = load_matrix("rnai_demeter2")

  proteomics = {}
  for key in (
    "proteomics_gygi",
    "proteomics_sanger",
    "proteomics_olink",
    "proteomics_rppa",
    "proteomics_rppa500",
  ):
    matrix = load_matrix(key)
    if matrix is not None and not matrix.empty:
      proteomics[key] = matrix

  loss_table = build_functional_loss_table(
    candidates,
    models,
    relative_cn=relative_cn,
    absolute_cn=absolute_cn,
    loh=loh,
    damaging=damaging,
    hotspot=hotspot,
    fusions=fusions,
    loss_threshold=args.loss_threshold,
  )
  write_tsv(loss_table, output_dir / "functional_loss_states.tsv")

  metabolomics = load_matrix("metabolomics")
  pair_evidence, protein_evidence, _ = build_pair_evidence(
    ranking,
    models,
    loss_table,
    dependency_probability=dependency_probability,
    proteomics=proteomics,
    rnai=rnai,
    metabolomics=metabolomics,
    min_group_size=args.min_group_size,
  )
  for p_column in (
    "dependency_probability_p_value",
    "rnai_p_value",
  ):
    if p_column in pair_evidence.columns:
      q_column = p_column.replace(
        "p_value", "q_value_bh_within_cancer"
      )
      pair_evidence[q_column] = (
        pair_evidence.groupby("cancer", group_keys=False)[p_column]
          .transform(lambda values: bh_adjust(values))
      )
  if not protein_evidence.empty and "p_value" in protein_evidence.columns:
    protein_evidence["q_value_bh_within_source_cancer"] = (
      protein_evidence.groupby(
        ["source", "cancer"], group_keys=False
      )["p_value"]
        .transform(lambda values: bh_adjust(values))
    )
  write_tsv(
    pair_evidence,
    output_dir / "extended_pair_evidence_by_cancer.tsv",
  )
  write_tsv(
    protein_evidence,
    output_dir / "proteomics_pair_evidence_by_source.tsv",
  )

  covariate_frames = []
  for key in (
    "ssgsea",
    "mirna",
    "global_chromatin",
    "omics_signatures",
    "molecular_subtypes",
    "subtype_matrix",
    "metmap_125",
    "metmap_500",
    "metmap_penetrance",
  ):
    matrix = load_matrix(key, gene_features=False)
    if matrix is None or matrix.empty:
      continue
    context = build_global_context_evidence(
      ranking,
      models,
      loss_table,
      matrix,
      source=key,
      min_group_size=args.min_group_size,
    )
    if not context.empty:
      covariate_frames.append(context)
    del matrix
    matrices.pop(key, None)
  covariates = (
    pd.concat(covariate_frames, ignore_index=True)
    if covariate_frames
    else pd.DataFrame()
  )
  write_tsv(
    covariates,
    output_dir / "extended_covariate_context.tsv",
  )

  drug_frames = []
  for key in (
    "prism_primary_viability",
    "prism_primary_compound_sample",
    "prism_secondary_auc",
    "prism_secondary_viability",
    "prism_secondary_log2_auc",
    "gdsc1_auc",
    "gdsc1_viability",
    "gdsc1_log2_auc",
    "gdsc2_auc",
    "gdsc2_viability",
    "gdsc2_log2_auc",
  ):
    matrix = load_matrix(key, gene_features=False)
    if matrix is None or matrix.empty:
      continue
    metric = key.replace("prism_", "").replace("gdsc", "GDSC")
    standardized_drug = standardize_matrix_long(
      matrix,
      source=key,
      metric=metric,
      source_file=input_dir / source_config[key]["filename"],
    ).rename(columns={"feature": "drug_name"})
    standardized_drug["drug_id"] = standardized_drug["drug_name"]
    standardized_drug["source_model_id"] = standardized_drug["model_id"]
    standardized_drug["lower_is_more_sensitive"] = True
    standardized_drug = standardized_drug[[
      "source",
      "model_id",
      "source_model_id",
      "drug_name",
      "drug_id",
      "response_value",
      "response_metric",
      "lower_is_more_sensitive",
      "source_file",
    ]]
    drug_frames.append(standardized_drug)
    del matrix
    matrices.pop(key, None)
  drug_long = (
    pd.concat(drug_frames, ignore_index=True)
    if drug_frames
    else pd.DataFrame()
  )
  write_tsv(
    drug_long,
    output_dir / "custom_drug_sensitivity_long.tsv",
  )

  combination_key_map = {
    "anchor_viability": "gdsc_combination_anchor_viability",
    "combination_auc": "gdsc_combination_combination_auc",
    "combination_viability": "gdsc_combination_combination_viability",
    "library_auc": "gdsc_combination_library_auc",
    "library_viability": "gdsc_combination_library_viability",
  }
  combination_matrices = {}
  combination_files = {}
  for short_key, source_key in combination_key_map.items():
    matrix = load_matrix(source_key, gene_features=False)
    if matrix is not None and not matrix.empty:
      combination_matrices[short_key] = matrix
      combination_files[short_key] = (
        input_dir / source_config[source_key]["filename"]
      )
  combination = build_gdsc_combination_table(
    combination_matrices,
    combination_files,
  )
  write_tsv(
    combination,
    output_dir / "gdsc_combination_evidence_long.tsv",
  )

  status = pd.DataFrame(statuses)
  write_tsv(
    status,
    output_dir / "extended_multiomics_source_status.tsv",
  )
  provenance = inventory.merge(
    status,
    on="source_key",
    how="left",
    suffixes=("_inventory", "_standardization"),
  )
  write_tsv(
    provenance,
    output_dir / "extended_multiomics_source_provenance.tsv",
  )

  payload = {
    "version": str(
      config.get("version")
      or "RSES-Onco-extended-multiomics-v0.12.0"
    ),
    "input_directory": str(input_dir),
    "ranking": str(ranking_path),
    "source_count": len(specs),
    "available_source_count": int(inventory["exists"].sum()),
    "standardized_source_count": (
      int((status.get("status") == "ok").sum())
      if not status.empty
      else 0
    ),
    "functional_loss_rows": len(loss_table),
    "pair_evidence_rows": len(pair_evidence),
    "protein_evidence_rows": len(protein_evidence),
    "covariate_context_rows": len(covariates),
    "drug_response_rows": len(drug_long),
    "combination_rows": len(combination),
    "primary_scored_layers": [
      "integrated functional loss",
      "CRISPR dependency probability",
      "multi-platform protein compensation",
      "RNAi orthogonal validation",
    ],
    "non_scored_context_layers": [
      "metabolomics without reaction mapping",
      "miRNA",
      "global chromatin",
      "ssGSEA",
      "molecular subtypes",
      "omics signatures",
      "MetMap",
      "hotspot mutations without loss-of-function annotation",
      "fusions without breakpoint-level disruption annotation",
    ],
  }
  write_status_json(
    output_dir / "extended_multiomics_status.json",
    payload,
  )
  print(json.dumps(payload, indent=2))


if __name__ == "__main__":
  main()
