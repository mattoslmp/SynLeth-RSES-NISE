#!/usr/bin/env python3
"""Preserve the eligibility-aware v1 contract while recording the v0.10.8 sublayer."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERNAL_VERSION = "eligibility-aware-wgcna-regulatory-v2"
PUBLIC_ELIGIBILITY_VERSION = "eligibility-aware-v1"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--ranking", required=True)
  args = parser.parse_args()
  path = resolve_path(args.ranking)
  frame = pd.read_csv(path, sep="\t", low_memory=False)
  versions = set(frame["scoring_semantics_version"].dropna().astype(str))
  if versions != {INTERNAL_VERSION}:
    raise ValueError(
      f"Expected internal WGCNA semantics {INTERNAL_VERSION}; observed {sorted(versions)}"
    )
  frame["expression_regulatory_semantics_version"] = INTERNAL_VERSION
  frame["scoring_semantics_version"] = PUBLIC_ELIGIBILITY_VERSION
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  print(
    f"Preserved {PUBLIC_ELIGIBILITY_VERSION} eligibility contract and recorded "
    f"{INTERNAL_VERSION}: {path}"
  )


if __name__ == "__main__":
  main()
