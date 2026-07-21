#!/usr/bin/env python3
"""Download the current Swiss-Prot NISE dataset and extract Homo sapiens rows.

The 2026 dataset uses UniProtKB release 2025_01. Column names are detected
case-insensitively because future versions may adjust the exact header labels.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd
import requests

URL = "https://zenodo.org/records/18008936/files/SwissProt_NISE.tsv.gz"
ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--raw", default="data/raw/SwissProt_NISE.tsv.gz")
  parser.add_argument("--output", default="data/processed/human_nise_2026.tsv")
  args = parser.parse_args()
  raw = resolve_path(args.raw)
  raw.parent.mkdir(parents=True, exist_ok=True)
  if not raw.exists():
    with requests.get(URL, stream=True, timeout=300) as response:
      response.raise_for_status()
      with raw.open("wb") as handle:
        for chunk in response.iter_content(1024 * 1024):
          handle.write(chunk)
  frame = pd.read_csv(raw, sep="\t", compression="gzip", low_memory=False)
  organism_col = next(c for c in frame.columns if "organism" in c.casefold())
  human = frame[frame[organism_col].astype(str).str.contains("Homo sapiens", case=False, na=False)].copy()
  out = resolve_path(args.output)
  out.parent.mkdir(parents=True, exist_ok=True)
  human.to_csv(out, sep="\t", index=False)
  print(f"Extracted {len(human):,} Homo sapiens records to {out}")


if __name__ == "__main__":
  main()
