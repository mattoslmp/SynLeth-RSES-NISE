#!/usr/bin/env python3
"""Complete the Circos expression summary with explicit unavailable rows.

The model-level expression table contains only observed values. This stage ensures
that the summary contains every Circos gene in every cancer context, distinguishing
missing gene columns or no observed values from biological zero.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANCERS = ("colon", "stomach", "lung")


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--coordinates",
    default=(
      "data/processed/circos/"
      "genomic_circos_gene_coordinates.tsv"
    ),
  )
  parser.add_argument(
    "--model-values",
    default=(
      "data/processed/circos/"
      "genomic_circos_expression_model_values.tsv"
    ),
  )
  parser.add_argument(
    "--summary",
    default=(
      "data/processed/circos/"
      "genomic_circos_expression_summary.tsv"
    ),
  )
  args = parser.parse_args()

  coordinates_path = resolve(args.coordinates)
  values_path = resolve(args.model_values)
  summary_path = resolve(args.summary)
  for path in (coordinates_path, values_path, summary_path):
    if not path.exists() or path.stat().st_size == 0:
      raise FileNotFoundError(path)

  coordinates = pd.read_csv(
    coordinates_path,
    sep="\t",
    low_memory=False,
  )
  values = pd.read_csv(values_path, sep="\t", low_memory=False)
  summary = pd.read_csv(summary_path, sep="\t", low_memory=False)
  genes = sorted(set(coordinates["gene"].astype(str)))
  source_file = (
    str(values["source_file"].dropna().iloc[0])
    if "source_file" in values
    and values["source_file"].notna().any()
    else "not_recorded"
  )

  complete = pd.MultiIndex.from_product(
    [CANCERS, genes],
    names=["cancer", "gene"],
  ).to_frame(index=False)
  result = complete.merge(
    summary,
    on=["cancer", "gene"],
    how="left",
  )
  observed = pd.to_numeric(
    result.get("observed_values"),
    errors="coerce",
  ).fillna(0)
  result["n_models"] = pd.to_numeric(
    result.get("n_models"),
    errors="coerce",
  ).fillna(0).astype(int)
  result["observed_values"] = observed.astype(int)
  result["unit"] = result.get(
    "unit",
    pd.Series(index=result.index, dtype=object),
  ).fillna("log2(TPM+1)")
  result["source_file"] = result.get(
    "source_file",
    pd.Series(index=result.index, dtype=object),
  ).fillna(source_file)
  result["evidence_status"] = "observed"
  result.loc[
    result["observed_values"].eq(0),
    "evidence_status",
  ] = "gene_or_context_unavailable_in_expression_matrix"
  result["absence_reason"] = ""
  result.loc[
    result["observed_values"].eq(0),
    "absence_reason",
  ] = "no_model_level_log2_tpm_plus_1_values_for_gene_and_cancer"
  result["missing_data_rule"] = (
    "missing_expression_remains_NA_and_is_not_numeric_zero"
  )
  result = result.sort_values(["cancer", "gene"])
  if len(result) != len(CANCERS) * len(genes):
    raise RuntimeError(
      "Circos expression summary does not cover every gene × cancer context"
    )
  atomic_tsv(result, summary_path)
  print(
    f"Completed Circos expression summary: {len(result):,} rows; "
    f"observed={int(result['observed_values'].gt(0).sum()):,}; "
    f"unavailable={int(result['observed_values'].eq(0).sum()):,}."
  )


if __name__ == "__main__":
  main()
