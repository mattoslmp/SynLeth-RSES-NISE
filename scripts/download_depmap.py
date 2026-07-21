#!/usr/bin/env python3
"""Prepare and validate the DepMap input directory.

DepMap directs bulk users to the official downloads section and may protect the
portal with browser verification. This script deliberately does not scrape the
portal. Download the current public release (schema tested against 26Q1), then
point --input-dir to those files.
"""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

REQUIRED = {
  "CRISPRGeneEffect.csv": ["CRISPRGeneEffect.csv", "CRISPR_gene_effect.csv", "CRISPRGeneEffect.csv.gz"],
  "OmicsCNGene.csv": ["OmicsCNGene.csv", "OmicsCNGene.csv.gz"],
  "Model.csv": ["Model.csv", "Model.csv.gz"],
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv": [
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv.gz",
    "OmicsExpression.csv",
  ],
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser(description="Validate locally downloaded DepMap release files.")
  parser.add_argument("--input-dir", default="data/raw/depmap")
  parser.add_argument("--allow-missing-expression", action="store_true")
  args = parser.parse_args()
  directory = resolve_path(args.input_dir)
  directory.mkdir(parents=True, exist_ok=True)
  missing: list[tuple[str, list[str]]] = []
  for logical, alternatives in REQUIRED.items():
    found = next((directory / name for name in alternatives if (directory / name).exists()), None)
    if found:
      print(f"OK  {logical}: {found}")
    elif args.allow_missing_expression and logical.startswith("OmicsExpression"):
      print(f"OPTIONAL  {logical}: absent; expression-compensation coverage will be reduced")
    else:
      missing.append((logical, alternatives))
  if missing:
    print("\nMissing DepMap inputs. Download them from the official current-release page:")
    print("  https://depmap.org/portal/data_page/?tab=currentRelease")
    for logical, alternatives in missing:
      print(f"  - {logical} (accepted names: {', '.join(alternatives)})")
    raise SystemExit(2)
  print("All required DepMap inputs are present.")


if __name__ == "__main__":
  main()
