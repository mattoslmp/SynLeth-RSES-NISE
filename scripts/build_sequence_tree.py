#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from rses_onco.tree import build_sequence_tree, catalog_accessions, download_uniprot_fasta

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser(description="Build activity-aware sequence trees for the human NISE catalogue.")
  parser.add_argument("--catalog", default="data/curated/human_nise_bonafide_2017.tsv")
  parser.add_argument("--output-prefix", default="results/human_nise_all")
  parser.add_argument("--threads", type=int, default=4)
  args = parser.parse_args()
  catalog = resolve_path(args.catalog)
  prefix = resolve_path(args.output_prefix)
  fasta = prefix.with_suffix(".fasta")
  accessions = catalog_accessions(catalog)
  download_uniprot_fasta(accessions, fasta)
  alignment, tree = build_sequence_tree(fasta, prefix, args.threads)
  print(f"Alignment: {alignment}\nTree: {tree}")
  print("Caution: whole-catalog alignments are exploratory; interpret activity-specific trees whenever possible.")


if __name__ == "__main__":
  main()
