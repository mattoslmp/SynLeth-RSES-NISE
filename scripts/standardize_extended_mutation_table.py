#!/usr/bin/env python3
"""Standardize the long DepMap mutation table to clear damaging/LoF events."""
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
  build_model_lookup,
  canonical_gene,
  find_first_column,
  map_model_id,
  read_table,
  sha256_file,
)


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def truthy_or_damaging(values: pd.Series) -> pd.Series:
  text = values.astype(str).str.casefold()
  return text.str.contains(
    r"(^|[^a-z0-9])(true|1|yes|y|damaging|deleterious|probably_damaging|"
    r"possibly_damaging|high|likely_lof|pathogenic)([^a-z0-9]|$)",
    regex=True,
    na=False,
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--config", default="config/extended_multiomics_sources.yaml"
  )
  parser.add_argument("--input-dir", default="dmap_data")
  parser.add_argument("--models", default="data/raw/depmap/Model.csv")
  parser.add_argument(
    "--output",
    default=(
      "data/processed/extended_multiomics/"
      "mutation_table_damaging_matrix.tsv"
    ),
  )
  parser.add_argument(
    "--event-output",
    default=(
      "data/processed/extended_multiomics/"
      "mutation_table_damaging_events.tsv"
    ),
  )
  parser.add_argument(
    "--status-output",
    default=(
      "data/processed/extended_multiomics/"
      "mutation_table_standardization_status.json"
    ),
  )
  parser.add_argument(
    "--runtime-config-output",
    default=(
      "data/processed/extended_multiomics/"
      "extended_multiomics_runtime_sources.yaml"
    ),
  )
  args = parser.parse_args()

  config_path = resolve(args.config)
  config = yaml.safe_load(
    config_path.read_text(encoding="utf-8")
  ) or {}
  source = config["sources"]["mutation_table"]
  input_path = resolve(args.input_dir) / str(source["filename"])
  models = read_table(resolve(args.models))
  frame = read_table(input_path)
  lookup = build_model_lookup(models)

  model_column = find_first_column(
    frame.columns,
    (
      "ModelID", "DepMap_ID", "depmap_id", "model_id", "CCLEName",
      "Tumor_Sample_Barcode", "cell_line",
    ),
  )
  gene_column = find_first_column(
    frame.columns,
    (
      "HugoSymbol", "Hugo_Symbol", "gene", "Gene", "GeneSymbol",
      "gene_symbol", "Symbol",
    ),
  )
  consequence_columns = [
    column
    for column in (
      "Variant_Classification", "Consequence", "consequence",
      "VEP_Consequence", "VariantInfo", "VariantType", "Variant_Type",
      "Protein_Change", "ProteinChange",
    )
    if column in frame.columns
  ]
  damaging_columns = [
    column
    for column in (
      "isDamaging", "damaging", "isDeleterious", "Deleterious",
      "LikelyLoF", "likely_lof", "PolyPhen", "SIFT", "VEP_Impact",
      "VepImpact", "IMPACT", "Oncogenic", "ClinicalSignificance",
    )
    if column in frame.columns
  ]
  if model_column is None or gene_column is None:
    raise ValueError(
      "Mutation table requires a model identifier and a gene-symbol column"
    )
  if not consequence_columns and not damaging_columns:
    raise ValueError(
      "Mutation table lacks functional consequence and damaging annotation fields"
    )

  combined_text = pd.Series("", index=frame.index, dtype="string")
  for column in consequence_columns:
    combined_text = combined_text.str.cat(
      frame[column].astype("string").fillna(""), sep=";"
    )
  normalized = combined_text.str.casefold()
  clear_lof = normalized.str.contains(
    r"frameshift|frame_shift|nonsense|stop_gained|stop_lost|"
    r"splice_acceptor|splice_donor|splice_site|start_lost|"
    r"translation_start_site|protein_truncating|\bfs\b|ter\d|\*",
    regex=True,
    na=False,
  )
  missense_or_inframe = normalized.str.contains(
    r"missense|in_frame|inframe", regex=True, na=False
  )
  annotated_damaging = pd.Series(False, index=frame.index)
  for column in damaging_columns:
    annotated_damaging |= truthy_or_damaging(frame[column])
  keep = clear_lof | (missense_or_inframe & annotated_damaging)

  events = pd.DataFrame({
    "ModelID": frame[model_column].map(
      lambda value: map_model_id(value, lookup)
    ),
    "gene": frame[gene_column].map(canonical_gene),
    "clear_loss_of_function": clear_lof,
    "annotated_damaging_missense_or_inframe": (
      missense_or_inframe & annotated_damaging
    ),
    "consequence_text": combined_text,
    "source_file": str(input_path),
  })
  events = events.loc[keep].dropna(subset=["ModelID"])
  events = events.loc[events["gene"].astype(bool)].drop_duplicates()
  if events.empty:
    raise ValueError(
      "No clear loss-of-function or annotated damaging variants remained"
    )
  events["value"] = 1.0
  matrix = events.pivot_table(
    index="ModelID",
    columns="gene",
    values="value",
    aggfunc="max",
    fill_value=0.0,
  ).astype(float)

  output = resolve(args.output)
  event_output = resolve(args.event_output)
  status_output = resolve(args.status_output)
  runtime_config_output = resolve(args.runtime_config_output)
  for path in (output, event_output, status_output, runtime_config_output):
    path.parent.mkdir(parents=True, exist_ok=True)
  matrix.reset_index().to_csv(output, sep="\t", index=False)
  events.to_csv(event_output, sep="\t", index=False)

  runtime = yaml.safe_load(
    config_path.read_text(encoding="utf-8")
  ) or {}
  runtime["sources"]["mutation_table"].update({
    "filename": str(output.resolve()),
    "layout": "matrix",
    "gene_features": True,
    "role": "scored_variant_level_functional_loss_standardized",
    "raw_source_file": str(input_path.resolve()),
  })
  runtime_config_output.write_text(
    yaml.safe_dump(runtime, sort_keys=False), encoding="utf-8"
  )
  status = {
    "status": "pass",
    "raw_source": str(input_path),
    "raw_sha256": sha256_file(input_path),
    "input_rows": len(frame),
    "retained_event_rows": len(events),
    "unique_models": events["ModelID"].nunique(),
    "unique_genes": events["gene"].nunique(),
    "clear_lof_rows": int(events["clear_loss_of_function"].sum()),
    "annotated_damaging_rows": int(
      events["annotated_damaging_missense_or_inframe"].sum()
    ),
    "matrix_output": str(output),
    "matrix_sha256": sha256_file(output),
    "runtime_config": str(runtime_config_output),
  }
  status_output.write_text(
    json.dumps(status, indent=2), encoding="utf-8"
  )
  print(json.dumps(status, indent=2))


if __name__ == "__main__":
  main()
