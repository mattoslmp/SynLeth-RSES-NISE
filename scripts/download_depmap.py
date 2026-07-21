#!/usr/bin/env python3
"""Validate locally downloaded DepMap release files.

The validator accepts current and legacy filenames and normalizes historical
row-identifier encodings (ModelID, DepMap_ID, or an unnamed ACH-* index).
It does not rewrite the large source matrices.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rses_onco.depmap import detect_model_id_column  # noqa: E402

REQUIRED = {
  "CRISPRGeneEffect": ["CRISPRGeneEffect.csv", "CRISPR_gene_effect.csv", "CRISPRGeneEffect.csv.gz"],
  "CopyNumber": [
    "OmicsCNGeneWGS.csv",
    "OmicsCNGeneWGS.csv.gz",
    "OmicsCNGene.csv",
    "OmicsCNGene.csv.gz",
  ],
  "Model": ["Model.csv", "Model.csv.gz"],
  "Expression": [
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv.gz",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv.gz",
    "OmicsExpression.csv",
  ],
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(block)
  return digest.hexdigest()


def looks_like_html(path: Path) -> bool:
  opener = __import__("gzip").open if path.suffix == ".gz" else open
  try:
    with opener(path, "rb") as handle:
      prefix = handle.read(512).lstrip().lower()
  except OSError:
    return False
  return prefix.startswith(b"<!doctype html") or prefix.startswith(b"<html")


def main() -> None:
  parser = argparse.ArgumentParser(description="Validate locally downloaded DepMap release files.")
  parser.add_argument("--input-dir", default="data/raw/depmap")
  parser.add_argument("--allow-missing-expression", action="store_true")
  parser.add_argument("--write-checksums", action="store_true")
  args = parser.parse_args()

  directory = resolve_path(args.input_dir)
  directory.mkdir(parents=True, exist_ok=True)
  resolved: dict[str, Path] = {}
  missing: list[tuple[str, list[str]]] = []

  for logical, alternatives in REQUIRED.items():
    found = next((directory / name for name in alternatives if (directory / name).exists()), None)
    if found:
      if found.stat().st_size == 0:
        raise SystemExit(f"Empty file: {found}")
      if looks_like_html(found):
        raise SystemExit(f"Downloaded HTML instead of a dataset: {found}")
      resolved[logical] = found
      print(f"OK  {logical}: {found} ({found.stat().st_size:,} bytes)")
    elif args.allow_missing_expression and logical == "Expression":
      print("OPTIONAL  Expression absent; expression-compensation coverage will be reduced")
    else:
      missing.append((logical, alternatives))

  if missing:
    print("\nMissing DepMap inputs. Download them from the official current-release page:")
    print("  https://depmap.org/portal/data_page/?tab=currentRelease")
    for logical, alternatives in missing:
      print(f"  - {logical} (accepted names: {', '.join(alternatives)})")
    raise SystemExit(2)

  for logical, path in resolved.items():
    header = pd.read_csv(path, nrows=0)
    try:
      raw_identifier, inferred = detect_model_id_column(path, logical)
    except ValueError as exc:
      raise SystemExit(str(exc)) from exc
    qualifier = "normalized to ModelID" if inferred else "ModelID present"
    print(
      f"SCHEMA  {logical}: {len(header.columns):,} columns; "
      f"identifier={raw_identifier!r} ({qualifier})"
    )

  if args.write_checksums:
    output = directory / "SHA256SUMS.txt"
    output.write_text(
      "".join(f"{sha256(path)}  {path.name}\n" for path in resolved.values()),
      encoding="utf-8",
    )
    print(f"Wrote {output}")

  print("All required DepMap inputs are present and have readable CSV schemas.")


if __name__ == "__main__":
  main()
