#!/usr/bin/env python3
"""Merge standardized custom-download drug responses into the canonical long table."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


COLUMNS = [
  "source",
  "model_id",
  "source_model_id",
  "drug_name",
  "drug_id",
  "response_value",
  "response_metric",
  "lower_is_more_sensitive",
  "source_file",
]


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    return pd.DataFrame(columns=COLUMNS)
  return pd.read_csv(path, sep="\t", low_memory=False)


def normalize(frame: pd.DataFrame) -> pd.DataFrame:
  result = frame.copy()
  for column in COLUMNS:
    if column not in result.columns:
      result[column] = pd.NA
  result = result[COLUMNS]
  result["response_value"] = pd.to_numeric(
    result["response_value"], errors="coerce"
  )
  result["model_id"] = result["model_id"].astype("string")
  result["drug_name"] = result["drug_name"].astype("string")
  result["source"] = result["source"].astype("string")
  result["response_metric"] = result["response_metric"].astype("string")
  result = result.dropna(
    subset=["model_id", "drug_name", "response_value"]
  )
  return result


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--canonical",
    default="data/processed/pharmacology/drug_sensitivity_long.tsv",
  )
  parser.add_argument(
    "--custom",
    default=(
      "data/processed/extended_multiomics/"
      "custom_drug_sensitivity_long.tsv"
    ),
  )
  parser.add_argument("--output", default=None)
  parser.add_argument(
    "--status-output",
    default=(
      "data/processed/pharmacology/"
      "extended_drug_sensitivity_merge_status.tsv"
    ),
  )
  args = parser.parse_args()

  canonical_path = Path(args.canonical)
  custom_path = Path(args.custom)
  output_path = Path(args.output) if args.output else canonical_path
  status_path = Path(args.status_output)

  canonical = normalize(read_optional(canonical_path))
  custom = normalize(read_optional(custom_path))
  combined = pd.concat([canonical, custom], ignore_index=True)
  key_columns = [
    "source",
    "model_id",
    "drug_name",
    "drug_id",
    "response_metric",
    "response_value",
  ]
  combined = combined.drop_duplicates(key_columns, keep="first")
  combined = combined.sort_values(
    ["source", "drug_name", "model_id", "response_metric"],
    na_position="last",
  )

  output_path.parent.mkdir(parents=True, exist_ok=True)
  temporary = output_path.with_suffix(output_path.suffix + ".tmp")
  combined.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output_path)

  status = pd.DataFrame([
    {
      "source_group": "canonical_preexisting",
      "input_path": str(canonical_path),
      "input_rows": len(canonical),
      "unique_sources": canonical["source"].nunique(),
      "unique_models": canonical["model_id"].nunique(),
      "unique_drugs": canonical["drug_name"].nunique(),
    },
    {
      "source_group": "custom_depmap_downloads",
      "input_path": str(custom_path),
      "input_rows": len(custom),
      "unique_sources": custom["source"].nunique(),
      "unique_models": custom["model_id"].nunique(),
      "unique_drugs": custom["drug_name"].nunique(),
    },
    {
      "source_group": "merged",
      "input_path": str(output_path),
      "input_rows": len(combined),
      "unique_sources": combined["source"].nunique(),
      "unique_models": combined["model_id"].nunique(),
      "unique_drugs": combined["drug_name"].nunique(),
    },
  ])
  status_path.parent.mkdir(parents=True, exist_ok=True)
  status.to_csv(status_path, sep="\t", index=False)
  print(f"Canonical rows before merge: {len(canonical):,}")
  print(f"Custom rows: {len(custom):,}")
  print(f"Merged rows: {len(combined):,}")
  print(f"Wrote {output_path}")
  print(f"Wrote {status_path}")


if __name__ == "__main__":
  main()
