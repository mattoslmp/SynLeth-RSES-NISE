#!/usr/bin/env python3
"""Validate aggregated TCGA/GDC homozygous-deletion matrices and report coverage."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_DEFAULT = {
  "colon": 575,
  "stomach": 429,
  "lung": 993,
}
DEFAULT_FILES = {
  "colon": "data/processed/TCGA_COLON_homdel_discrete.tsv",
  "stomach": "data/processed/TCGA_STOMACH_homdel_discrete.tsv",
  "lung": "data/processed/TCGA_LUNG_homdel_discrete.tsv",
}
DEFAULT_GENES = [
  "MTAP", "ATM", "ARID1A", "SMARCA4", "ENO1", "ME2", "NTHL1",
  "OGG1", "APEX2", "SOD2", "PI4KA", "VPS4B", "STAG2",
]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def parse_mapping(items: list[str], defaults: dict[str, str]) -> dict[str, Path]:
  mapping = {key: resolve_path(value) for key, value in defaults.items()}
  for item in items:
    key, value = item.split("=", 1)
    mapping[key.strip().casefold()] = resolve_path(value)
  return mapping


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--matrix",
    action="append",
    default=[],
    help="Override cancer matrix as cancer=path/to/matrix.tsv",
  )
  parser.add_argument(
    "--expected-samples",
    action="append",
    default=[],
    help="Override expected sample count as cancer=integer; use 0 to skip exact-count enforcement",
  )
  parser.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
  parser.add_argument(
    "--output",
    default="results/empirical_26Q1/full/gdc_matrix_qc.tsv",
  )
  parser.add_argument(
    "--event-output",
    default="results/empirical_26Q1/full/tcga_gene_event_summary.tsv",
  )
  args = parser.parse_args()

  matrices = parse_mapping(args.matrix, DEFAULT_FILES)
  expected = dict(EXPECTED_DEFAULT)
  for item in args.expected_samples:
    key, value = item.split("=", 1)
    expected[key.strip().casefold()] = int(value)

  qc_rows: list[dict[str, object]] = []
  event_rows: list[dict[str, object]] = []
  failures: list[str] = []

  for cancer, path in matrices.items():
    if not path.exists():
      failures.append(f"Missing matrix for {cancer}: {path}")
      continue

    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str).str.upper().str.strip()
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    observed = set(pd.unique(numeric.to_numpy().ravel()))
    observed = {float(value) for value in observed if pd.notna(value)}
    invalid = sorted(observed - {-2.0, 0.0})
    duplicate_genes = int(matrix.index.duplicated().sum())
    duplicate_samples = int(matrix.columns.duplicated().sum())
    homdel_events = int(numeric.eq(-2).sum().sum())
    expected_samples = int(expected.get(cancer, 0))

    qc_rows.append({
      "cancer": cancer,
      "path": str(path),
      "genes": int(matrix.shape[0]),
      "samples": int(matrix.shape[1]),
      "expected_samples": expected_samples,
      "duplicate_genes": duplicate_genes,
      "duplicate_samples": duplicate_samples,
      "homdel_events": homdel_events,
      "invalid_values": ",".join(map(str, invalid)),
      "valid": (
        not invalid
        and duplicate_genes == 0
        and duplicate_samples == 0
        and (expected_samples == 0 or matrix.shape[1] == expected_samples)
      ),
    })

    if invalid:
      failures.append(f"{cancer}: invalid values {invalid}")
    if duplicate_genes:
      failures.append(f"{cancer}: {duplicate_genes} duplicated gene rows")
    if duplicate_samples:
      failures.append(f"{cancer}: {duplicate_samples} duplicated samples")
    if expected_samples and matrix.shape[1] != expected_samples:
      failures.append(
        f"{cancer}: expected {expected_samples} samples, observed {matrix.shape[1]}"
      )

    for gene in args.genes:
      gene_key = str(gene).upper().strip()
      if gene_key not in numeric.index:
        event_rows.append({
          "cancer": cancer,
          "gene": gene_key,
          "present": False,
          "homdel_n": pd.NA,
          "evaluable_n": pd.NA,
          "homdel_frequency": pd.NA,
        })
        continue
      values = pd.to_numeric(numeric.loc[gene_key], errors="coerce").dropna()
      deleted = int(values.eq(-2).sum())
      total = int(len(values))
      event_rows.append({
        "cancer": cancer,
        "gene": gene_key,
        "present": True,
        "homdel_n": deleted,
        "evaluable_n": total,
        "homdel_frequency": deleted / total if total else pd.NA,
      })

  output = resolve_path(args.output)
  event_output = resolve_path(args.event_output)
  output.parent.mkdir(parents=True, exist_ok=True)
  event_output.parent.mkdir(parents=True, exist_ok=True)
  pd.DataFrame(qc_rows).to_csv(output, sep="\t", index=False)
  pd.DataFrame(event_rows).to_csv(event_output, sep="\t", index=False)

  print(pd.DataFrame(qc_rows).to_string(index=False), flush=True)
  print(f"Wrote {output}", flush=True)
  print(f"Wrote {event_output}", flush=True)

  if failures:
    raise SystemExit("GDC matrix validation failed:\n- " + "\n- ".join(failures))
  print("GDC matrix validation passed.", flush=True)


if __name__ == "__main__":
  main()
