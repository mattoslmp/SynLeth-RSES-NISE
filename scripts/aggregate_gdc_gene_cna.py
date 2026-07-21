#!/usr/bin/env python3
"""Aggregate downloaded GDC ASCAT3 gene-level integer CN into homdel matrices.

Output is deletion-only discrete coding: 0 total copies -> -2; >0 copies -> 0;
missing remains missing. This is compatible with the RSES-Onco homozygous-
deletion frequency function but must not be described as a GISTIC matrix.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANCER_PROJECTS = {
  "colon": ["TCGA-COAD", "TCGA-READ"],
  "stomach": ["TCGA-STAD"],
  "lung": ["TCGA-LUAD", "TCGA-LUSC"],
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sample_identifier(hit: dict) -> str:
  """Prefer the TCGA sample barcode; fall back to case ID and then file UUID."""
  for case in hit.get("cases") or []:
    for sample in case.get("samples") or []:
      value = sample.get("submitter_id")
      if value:
        return str(value)
  for case in hit.get("cases") or []:
    value = case.get("submitter_id")
    if value:
      return str(value)
  return str(hit["file_id"])


def choose_column(columns, candidates):
  lookup = {str(c).casefold(): c for c in columns}
  for candidate in candidates:
    if candidate.casefold() in lookup:
      return lookup[candidate.casefold()]
  return None


def read_gene_copy_number(path: Path) -> pd.Series:
  try:
    frame = pd.read_csv(path, sep="\t", comment="#")
  except Exception:
    frame = pd.read_csv(path, sep=None, engine="python", comment="#")
  gene_col = choose_column(frame.columns, [
    "gene_name", "gene_symbol", "hugo_symbol", "Gene Symbol", "gene",
  ])
  cn_col = choose_column(frame.columns, ["copy_number", "Copy_Number", "copy number"])
  if gene_col is None or cn_col is None:
    raise ValueError(f"Could not identify gene/copy_number columns in {path}; columns={list(frame.columns)}")
  genes = frame[gene_col].astype(str).str.replace(r"\..*$", "", regex=True).str.upper().str.strip()
  values = pd.to_numeric(frame[cn_col], errors="coerce")
  data = pd.DataFrame({"gene": genes, "copy_number": values}).dropna()
  data = data.loc[data["gene"].ne("")]
  # If a gene is represented more than once, retain the minimum total copy number.
  return data.groupby("gene", sort=False)["copy_number"].min()


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--raw-dir", default="data/raw/gdc")
  parser.add_argument("--manifest", default="data/raw/gdc/gdc_gene_level_copy_number_manifest.json")
  parser.add_argument("--output-dir", default="data/processed")
  args = parser.parse_args()
  raw_dir = resolve_path(args.raw_dir)
  manifest = json.loads(resolve_path(args.manifest).read_text(encoding="utf-8"))
  outdir = resolve_path(args.output_dir)
  outdir.mkdir(parents=True, exist_ok=True)

  project_matrices: dict[str, pd.DataFrame] = {}
  for project, hits in manifest.items():
    columns: dict[str, pd.Series] = {}
    for hit in hits:
      path = raw_dir / project / hit["file_name"]
      if not path.exists():
        raise FileNotFoundError(path)
      sample = sample_identifier(hit)
      if sample in columns:
        sample = f"{sample}__{str(hit['file_id'])[:8]}"
      columns[sample] = read_gene_copy_number(path)
    if not columns:
      raise ValueError(f"No files available for {project}")
    absolute = pd.concat(columns, axis=1).sort_index()
    discrete = pd.DataFrame(np.nan, index=absolute.index, columns=absolute.columns)
    discrete[absolute == 0] = -2
    discrete[absolute > 0] = 0
    discrete.index.name = "Hugo_Symbol"
    output = outdir / f"{project.replace('-', '_')}_homdel_discrete.tsv"
    discrete.to_csv(output, sep="\t", na_rep="NA")
    project_matrices[project] = discrete
    print(f"Wrote {output}: {discrete.shape[0]:,} genes x {discrete.shape[1]:,} cases", flush=True)

  for cancer, projects in CANCER_PROJECTS.items():
    missing = [project for project in projects if project not in project_matrices]
    if missing:
      raise ValueError(f"Missing projects for {cancer}: {missing}")
    combined = pd.concat([project_matrices[p] for p in projects], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()].sort_index()
    combined.index.name = "Hugo_Symbol"
    output = outdir / f"TCGA_{cancer.upper()}_homdel_discrete.tsv"
    combined.to_csv(output, sep="\t", na_rep="NA")
    print(f"Wrote {output}: {combined.shape[0]:,} genes x {combined.shape[1]:,} cases", flush=True)


if __name__ == "__main__":
  main()
