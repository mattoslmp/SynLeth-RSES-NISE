#!/usr/bin/env python3
"""Validate aggregated TCGA/GDC homozygous-deletion matrices and report coverage.

In addition to dimensions and deletion coding, this validator enforces that every
matrix column represents a TCGA Primary Tumor barcode (sample type code 01) and
matches the file-to-sample provenance emitted by ``aggregate_gdc_gene_cna.py``.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import re

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
CANCER_PROJECTS = {
  "colon": {"TCGA-COAD", "TCGA-READ"},
  "stomach": {"TCGA-STAD"},
  "lung": {"TCGA-LUAD", "TCGA-LUSC"},
}
DEFAULT_GENES = [
  "MTAP", "ATM", "ARID1A", "SMARCA4", "ENO1", "ME2", "NTHL1",
  "OGG1", "APEX2", "SOD2", "PI4KA", "VPS4B", "STAG2",
]
TCGA_SAMPLE = re.compile(
  r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-(?P<sample_code>\d{2})[A-Z]",
  flags=re.IGNORECASE,
)
PRIMARY_TUMOR_CODE = "01"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def parse_mapping(items: list[str], defaults: dict[str, str]) -> dict[str, Path]:
  mapping = {key: resolve_path(value) for key, value in defaults.items()}
  for item in items:
    key, value = item.split("=", 1)
    mapping[key.strip().casefold()] = resolve_path(value)
  return mapping


def tcga_sample_type_code(value: object) -> str | None:
  """Extract the two-digit TCGA sample type code from a matrix column."""
  match = TCGA_SAMPLE.match(str(value).strip())
  return match.group("sample_code") if match else None


def summarize_sample_types(columns: pd.Index) -> tuple[Counter, list[str], list[str]]:
  codes: Counter = Counter()
  non_primary: list[str] = []
  unparsed: list[str] = []
  for column in columns.astype(str):
    code = tcga_sample_type_code(column)
    if code is None:
      unparsed.append(column)
      continue
    codes[code] += 1
    if code != PRIMARY_TUMOR_CODE:
      non_primary.append(column)
  return codes, non_primary, unparsed


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
    help=(
      "Override expected sample count as cancer=integer; use 0 to skip "
      "exact-count enforcement"
    ),
  )
  parser.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
  parser.add_argument(
    "--provenance",
    default="data/processed/gdc_sample_provenance.tsv",
  )
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

  provenance_path = resolve_path(args.provenance)
  provenance: pd.DataFrame | None = None
  provenance_failures: list[str] = []
  if not provenance_path.exists() or provenance_path.stat().st_size == 0:
    provenance_failures.append(
      f"Missing or empty GDC sample provenance: {provenance_path}"
    )
  else:
    provenance = pd.read_csv(provenance_path, sep="\t", low_memory=False)
    required = {
      "project",
      "file_id",
      "file_name",
      "sample_id",
      "sample_type",
      "matrix_column",
      "status",
    }
    missing = sorted(required - set(provenance.columns))
    if missing:
      provenance_failures.append(
        f"GDC sample provenance lacks columns: {missing}"
      )
    elif not provenance["sample_type"].astype(str).str.casefold().eq(
      "primary tumor"
    ).all():
      invalid_types = sorted(
        set(
          provenance.loc[
            ~provenance["sample_type"]
              .astype(str)
              .str.casefold()
              .eq("primary tumor"),
            "sample_type",
          ].astype(str)
        )
      )
      provenance_failures.append(
        f"Non-primary sample types in provenance: {invalid_types}"
      )

  qc_rows: list[dict[str, object]] = []
  event_rows: list[dict[str, object]] = []
  failures: list[str] = list(provenance_failures)

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
    codes, non_primary, unparsed = summarize_sample_types(matrix.columns)

    provenance_column_match = False
    expected_columns: set[str] = set()
    if provenance is not None and {
      "project",
      "matrix_column",
    }.issubset(provenance.columns):
      projects = CANCER_PROJECTS.get(cancer, set())
      expected_columns = set(
        provenance.loc[
          provenance["project"].astype(str).isin(projects),
          "matrix_column",
        ].astype(str)
      )
      provenance_column_match = set(matrix.columns.astype(str)) == expected_columns

    valid = (
      not invalid
      and duplicate_genes == 0
      and duplicate_samples == 0
      and not non_primary
      and not unparsed
      and provenance_column_match
      and (expected_samples == 0 or matrix.shape[1] == expected_samples)
    )
    qc_rows.append({
      "cancer": cancer,
      "path": str(path),
      "genes": int(matrix.shape[0]),
      "samples": int(matrix.shape[1]),
      "expected_samples": expected_samples,
      "duplicate_genes": duplicate_genes,
      "duplicate_samples": duplicate_samples,
      "primary_tumor_samples": int(codes.get(PRIMARY_TUMOR_CODE, 0)),
      "sample_type_codes": ";".join(
        f"{code}={count}" for code, count in sorted(codes.items())
      ),
      "non_primary_tumor_samples": len(non_primary),
      "unparsed_sample_barcodes": len(unparsed),
      "non_primary_examples": ";".join(non_primary[:10]),
      "unparsed_examples": ";".join(unparsed[:10]),
      "provenance_expected_columns": len(expected_columns),
      "provenance_column_match": provenance_column_match,
      "homdel_events": homdel_events,
      "invalid_values": ",".join(map(str, invalid)),
      "valid": valid,
    })

    if invalid:
      failures.append(f"{cancer}: invalid values {invalid}")
    if duplicate_genes:
      failures.append(f"{cancer}: {duplicate_genes} duplicated gene rows")
    if duplicate_samples:
      failures.append(f"{cancer}: {duplicate_samples} duplicated samples")
    if non_primary:
      failures.append(
        f"{cancer}: {len(non_primary)} non-primary TCGA sample barcodes; "
        f"examples={non_primary[:10]}"
      )
    if unparsed:
      failures.append(
        f"{cancer}: {len(unparsed)} unparsed TCGA sample barcodes; "
        f"examples={unparsed[:10]}"
      )
    if provenance is not None and not provenance_column_match:
      failures.append(
        f"{cancer}: matrix columns do not exactly match GDC sample provenance"
      )
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
