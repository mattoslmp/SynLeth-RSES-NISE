#!/usr/bin/env python3
"""Structural and overlap validation for real DepMap inputs."""
from __future__ import annotations

import argparse
import re
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rses_onco.depmap import cancer_model_ids, read_model_ids  # noqa: E402

SIMPLE_LOSS_PATTERNS = [
  re.compile(r"^([A-Za-z0-9-]+)$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+homozygous deletion$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss/low expression$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss or low expression$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss/low activity$", re.I),
]


def extract_single_lost_gene(feature: str) -> str | None:
  for pattern in SIMPLE_LOSS_PATTERNS:
    match = pattern.fullmatch(str(feature).strip())
    if match:
      return match.group(1).upper()
  return None


def resolve(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def genes(path: Path) -> dict[str, str]:
  columns = pd.read_csv(path, nrows=0).columns
  metadata = {
    "ModelID", "DepMap_ID", "DepMapID", "DepMap ID", "model_id",
    "ProfileID", "is_default_entry", "IsDefaultEntryForModel",
    "IsDefaultEntryForMC", "ModelConditionID", "SequencingID",
  }
  return {
    re.sub(r"\s*\([^)]*\)\s*$", "", str(c)).upper(): str(c)
    for c in columns
    if str(c) not in metadata and not str(c).startswith("Unnamed:")
  }


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument("--expression")
  parser.add_argument("--candidates", default="data/curated/synthetic_lethality_reference_pairs.tsv")
  args = parser.parse_args()

  effect_path = resolve(args.gene_effect)
  copy_path = resolve(args.copy_number)
  model_path = resolve(args.models)
  expr_path = resolve(args.expression) if args.expression else None

  effect_ids = read_model_ids(effect_path, "CRISPRGeneEffect")
  copy_ids = read_model_ids(copy_path, "CopyNumber")
  model_ids = read_model_ids(model_path, "Model.csv")
  model = pd.read_csv(model_path)
  if "ModelID" not in model.columns:
    raw = next((c for c in ("DepMap_ID", "DepMapID", "DepMap ID", "model_id") if c in model.columns), None)
    if raw:
      model = model.rename(columns={raw: "ModelID"})
  if "ModelID" not in model.columns:
    raise SystemExit("Model.csv lacks a ModelID-compatible column")
  if "OncotreeLineage" not in model.columns:
    raise SystemExit("Model.csv lacks OncotreeLineage")

  print(f"Gene effect models: {effect_ids.nunique():,}; duplicate rows: {effect_ids.duplicated().sum():,}")
  print(f"Copy-number models: {copy_ids.nunique():,}; duplicate rows: {copy_ids.duplicated().sum():,}")
  print(f"Metadata models: {model_ids.nunique():,}; duplicate rows: {model_ids.duplicated().sum():,}")
  common = set(effect_ids) & set(copy_ids) & set(model_ids)
  print(f"Common ModelID across effect/CN/metadata: {len(common):,}")
  for cancer, label in (("colon", "Colorectal"), ("stomach", "Gastric"), ("lung", "Lung")):
    selected = set(cancer_model_ids(model, cancer))
    count = len(selected & common)
    print(f"  {label}: {count:,} common models")

  print("Observed OncotreeLineage values among common models:")
  common_model = model.loc[model["ModelID"].astype(str).isin(common)].copy()
  observed = (
    common_model["OncotreeLineage"]
      .fillna("<missing>")
      .astype(str)
      .value_counts()
      .sort_index()
  )
  for value, count in observed.items():
    print(f"    {value}: {count:,}")

  effect_genes = genes(effect_path)
  copy_genes = genes(copy_path)
  expr_genes = genes(expr_path) if expr_path else {}
  candidates = pd.read_csv(resolve(args.candidates), sep="\t")
  simple = re.compile(r"^[A-Za-z0-9-]+$")
  lost = sorted({gene for x in candidates.lost_feature if (gene := extract_single_lost_gene(str(x)))})
  targets = sorted({str(x).upper() for x in candidates.target_gene if simple.fullmatch(str(x))})
  print("Missing simple lost genes in CN:", [g for g in lost if g not in copy_genes])
  print("Missing simple target genes in CRISPR:", [g for g in targets if g not in effect_genes])
  if expr_path:
    print("Missing simple target genes in expression:", [g for g in targets if g not in expr_genes])

  if not common:
    raise SystemExit("No shared ModelID values")


if __name__ == "__main__":
  main()
