#!/usr/bin/env python3
"""Validate local DepMap extension schemas before any scientific calculation."""
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

from rses_onco.extended_multiomics import find_first_column, read_table, sha256_file


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/extended_multiomics_sources.yaml")
  parser.add_argument("--input-dir", default="dmap_data")
  parser.add_argument(
    "--output",
    default="data/processed/extended_multiomics/extended_input_schema_validation.tsv",
  )
  parser.add_argument("--strict", action="store_true")
  args = parser.parse_args()

  config = yaml.safe_load(resolve(args.config).read_text(encoding="utf-8")) or {}
  data_dir = resolve(args.input_dir)
  rows: list[dict[str, object]] = []
  errors: list[str] = []
  hashes: dict[str, list[str]] = {}

  for key, specification in (config.get("sources") or {}).items():
    path = data_dir / str(specification["filename"])
    record: dict[str, object] = {
      "source_key": key,
      "path": str(path),
      "role": specification.get("role", ""),
      "exists": path.exists(),
      "size_bytes": path.stat().st_size if path.exists() else 0,
      "schema_status": "not_found",
      "n_columns_preview": 0,
      "message": "",
    }
    if not path.exists() or path.stat().st_size == 0:
      if args.strict:
        errors.append(f"missing_or_empty:{key}:{path}")
      rows.append(record)
      continue
    digest = sha256_file(path)
    record["sha256"] = digest
    hashes.setdefault(digest, []).append(key)
    try:
      preview = read_table(path, nrows=20)
      record["n_columns_preview"] = len(preview.columns)
      if len(preview.columns) < 2:
        raise ValueError("fewer_than_two_columns")
      if key == "mutation_table":
        model_column = find_first_column(
          preview.columns,
          (
            "ModelID",
            "DepMap_ID",
            "depmap_id",
            "model_id",
            "CCLEName",
            "Tumor_Sample_Barcode",
            "cell_line",
          ),
        )
        gene_column = find_first_column(
          preview.columns,
          (
            "HugoSymbol",
            "Hugo_Symbol",
            "gene",
            "Gene",
            "GeneSymbol",
            "gene_symbol",
            "Symbol",
          ),
        )
        consequence_column = find_first_column(
          preview.columns,
          (
            "Variant_Classification",
            "Consequence",
            "consequence",
            "Protein_Change",
            "VariantType",
            "Variant_Type",
          ),
        )
        missing = [
          name
          for name, value in (
            ("model", model_column),
            ("gene", gene_column),
            ("consequence", consequence_column),
          )
          if value is None
        ]
        if missing:
          raise ValueError(
            "mutation_table_missing_required_fields:" + ",".join(missing)
          )
      record["schema_status"] = "pass"
    except Exception as exc:
      record["schema_status"] = "fail"
      record["message"] = str(exc)
      errors.append(f"schema_failure:{key}:{exc}")
    rows.append(record)

  for keys in hashes.values():
    if len(keys) > 1:
      for record in rows:
        if record.get("source_key") in keys:
          record["identical_content_group"] = ";".join(sorted(keys))

  output = resolve(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
  payload = {
    "status": "pass" if not errors else "fail",
    "errors": errors,
    "sources": len(rows),
    "schema_pass": sum(row["schema_status"] == "pass" for row in rows),
    "identical_content_groups": [
      keys for keys in hashes.values() if len(keys) > 1
    ],
  }
  output.with_suffix(".json").write_text(
    json.dumps(payload, indent=2), encoding="utf-8"
  )
  if errors:
    raise SystemExit(
      "Extended input schema validation failed:\n" + "\n".join(errors)
    )
  print(json.dumps(payload, indent=2))


if __name__ == "__main__":
  main()
