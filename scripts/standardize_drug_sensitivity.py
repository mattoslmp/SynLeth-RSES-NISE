#!/usr/bin/env python3
"""Standardize PRISM, GDSC and CTRP response files into one long table.

The script is intentionally release-agnostic. It reads source definitions from a
YAML file, resolves common column names, maps cell-line names to DepMap ModelID
when possible, and records every skipped or ambiguous file in a status table.
Raw source files are never modified.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd
from pandas.errors import EmptyDataError
import yaml

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_table(path: Path) -> pd.DataFrame:
  separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
  return pd.read_csv(path, sep=separator, low_memory=False)


def first_existing(columns: list[str], candidates: list[str]) -> str | None:
  exact = {str(column): str(column) for column in columns}
  casefolded = {str(column).casefold(): str(column) for column in columns}
  for candidate in candidates:
    if candidate in exact:
      return exact[candidate]
    if candidate.casefold() in casefolded:
      return casefolded[candidate.casefold()]
  return None


def build_model_lookup(models: pd.DataFrame) -> dict[str, str]:
  if "ModelID" not in models:
    raise ValueError("Model.csv requires ModelID")
  candidate_columns = [
    column for column in (
      "ModelID",
      "CCLEName",
      "CellLineName",
      "StrippedCellLineName",
      "ModelConditionID",
      "COSMICID",
      "COSMIC_ID",
    )
    if column in models.columns
  ]
  lookup: dict[str, str] = {}
  for record in models[candidate_columns].to_dict("records"):
    model_id = str(record.get("ModelID") or "").strip()
    if not model_id:
      continue
    for column in candidate_columns:
      value = record.get(column)
      if pd.isna(value):
        continue
      key = str(value).strip().casefold()
      if key:
        lookup.setdefault(key, model_id)
  return lookup


def map_model_ids(values: pd.Series, lookup: dict[str, str]) -> pd.Series:
  def one(value: object) -> str | None:
    if pd.isna(value):
      return None
    text = str(value).strip()
    if text.upper().startswith("ACH-"):
      return text
    return lookup.get(text.casefold())
  return values.map(one)


def standardize_long(
  frame: pd.DataFrame,
  source_name: str,
  specification: dict[str, Any],
  source_file: Path,
  model_lookup: dict[str, str],
) -> pd.DataFrame:
  columns = [str(column) for column in frame.columns]
  model_column = first_existing(columns, list(specification.get("model_columns") or []))
  drug_column = first_existing(columns, list(specification.get("drug_columns") or []))
  drug_id_column = first_existing(columns, list(specification.get("drug_id_columns") or []))
  value_column = first_existing(columns, list(specification.get("value_columns") or []))
  missing = [
    name for name, value in (
      ("model", model_column),
      ("drug", drug_column),
      ("value", value_column),
    )
    if value is None
  ]
  if missing:
    raise ValueError(
      f"Could not resolve {missing}; first columns={columns[:25]}"
    )
  output = pd.DataFrame({
    "source": source_name,
    "model_id": map_model_ids(frame[model_column], model_lookup),
    "source_model_id": frame[model_column].astype(str),
    "drug_name": frame[drug_column].astype(str).str.strip(),
    "drug_id": (
      frame[drug_id_column].astype(str).str.strip()
      if drug_id_column is not None
      else frame[drug_column].astype(str).str.strip()
    ),
    "response_value": pd.to_numeric(frame[value_column], errors="coerce"),
    "response_metric": str(specification.get("metric") or value_column),
    "lower_is_more_sensitive": bool(specification.get("lower_is_more_sensitive", True)),
    "source_file": str(source_file),
  })
  return output.dropna(subset=["model_id", "drug_name", "response_value"])


def standardize_wide(
  frame: pd.DataFrame,
  source_name: str,
  specification: dict[str, Any],
  source_file: Path,
  model_lookup: dict[str, str],
) -> pd.DataFrame:
  row_id_column = specification.get("row_id_column")
  if not row_id_column or row_id_column not in frame.columns:
    raise ValueError("Wide source requires a valid row_id_column")
  value_columns = [column for column in frame.columns if column != row_id_column]
  melted = frame.melt(
    id_vars=[row_id_column],
    value_vars=value_columns,
    var_name="source_model_id",
    value_name="response_value",
  )
  output = pd.DataFrame({
    "source": source_name,
    "model_id": map_model_ids(melted["source_model_id"], model_lookup),
    "source_model_id": melted["source_model_id"].astype(str),
    "drug_name": melted[row_id_column].astype(str).str.strip(),
    "drug_id": melted[row_id_column].astype(str).str.strip(),
    "response_value": pd.to_numeric(melted["response_value"], errors="coerce"),
    "response_metric": str(specification.get("metric") or "response"),
    "lower_is_more_sensitive": bool(specification.get("lower_is_more_sensitive", True)),
    "source_file": str(source_file),
  })
  return output.dropna(subset=["model_id", "drug_name", "response_value"])


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--config",
    default="config/drug_sensitivity_sources.yaml",
  )
  parser.add_argument(
    "--models",
    default="data/raw/depmap/Model.csv",
  )
  parser.add_argument(
    "--output",
    default="data/processed/pharmacology/drug_sensitivity_long.tsv",
  )
  parser.add_argument(
    "--status-output",
    default="data/processed/pharmacology/drug_sensitivity_source_status.tsv",
  )
  parser.add_argument("--strict", action="store_true")
  args = parser.parse_args()

  config = yaml.safe_load(resolve_path(args.config).read_text(encoding="utf-8")) or {}
  sources = config.get("sources") or {}
  models = pd.read_csv(resolve_path(args.models), low_memory=False)
  model_lookup = build_model_lookup(models)
  standardized: list[pd.DataFrame] = []
  statuses: list[dict[str, object]] = []

  for source_name, specification in sources.items():
    if not specification.get("enabled", True):
      statuses.append({"source": source_name, "status": "disabled"})
      continue
    paths: list[Path] = []
    for pattern in specification.get("path_globs") or []:
      paths.extend(sorted(ROOT.glob(pattern)))
    paths = sorted(set(path for path in paths if path.is_file()))
    if not paths:
      statuses.append({
        "source": source_name,
        "status": "not_found",
        "message": ";".join(specification.get("path_globs") or []),
      })
      if args.strict:
        raise FileNotFoundError(f"No files found for {source_name}")
      continue
    for path in paths:
      try:
        frame = read_table(path)
        if frame.empty:
          raise EmptyDataError("empty table")
        layout = str(specification.get("layout") or "long").casefold()
        if layout == "wide":
          result = standardize_wide(
            frame,
            source_name,
            specification,
            path,
            model_lookup,
          )
        else:
          result = standardize_long(
            frame,
            source_name,
            specification,
            path,
            model_lookup,
          )
        standardized.append(result)
        statuses.append({
          "source": source_name,
          "source_file": str(path),
          "status": "ok",
          "input_rows": len(frame),
          "output_rows": len(result),
          "unique_models": result["model_id"].nunique(),
          "unique_drugs": result["drug_name"].nunique(),
        })
      except Exception as exc:
        statuses.append({
          "source": source_name,
          "source_file": str(path),
          "status": "failed",
          "message": str(exc),
        })
        if args.strict:
          raise

  columns = [
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
  if standardized:
    result = pd.concat(standardized, ignore_index=True)
    result = result[columns].drop_duplicates()
  else:
    result = pd.DataFrame(columns=columns)
  output = resolve_path(args.output)
  status_output = resolve_path(args.status_output)
  output.parent.mkdir(parents=True, exist_ok=True)
  status_output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)
  pd.DataFrame(statuses).to_csv(status_output, sep="\t", index=False)
  print(f"Standardized sensitivity rows: {len(result):,}")
  print(f"Unique models: {result['model_id'].nunique():,}")
  print(f"Unique drugs: {result['drug_name'].nunique():,}")
  print(f"Wrote {output}")
  print(f"Wrote {status_output}")


if __name__ == "__main__":
  main()
