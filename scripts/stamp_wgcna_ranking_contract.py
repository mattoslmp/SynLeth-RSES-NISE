#!/usr/bin/env python3
"""Record WGCNA/regulatory extension while retaining eligibility-aware-v1 semantics.

The missingness/eligibility algebra is unchanged from eligibility-aware-v1. WGCNA and
promoter-aware regulation refine existing functional-microniche domains rather than
adding independent top-level weights, so the extension receives a separate version
field instead of masquerading as a new missing-data contract.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--ranking", required=True)
  args = parser.parse_args()
  path = resolve_path(args.ranking)
  frame = pd.read_csv(path, sep="\t", low_memory=False)
  required = {
    "component_wgcna_expression_network",
    "regulatory_promoter_motif_divergence",
    "coverage_adjusted_rses",
  }
  missing = sorted(required - set(frame.columns))
  if missing:
    raise ValueError(
      "Cannot stamp WGCNA/regulatory ranking; missing fields: "
      + ", ".join(missing)
    )
  frame["scoring_semantics_version"] = "eligibility-aware-v1"
  frame["scoring_extension_version"] = "wgcna-promoter-regulatory-v1"
  frame["score_version"] = "RSES-Onco-expanded-v0.10.8"
  frame["direct_promoter_binding_claim"] = False
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  print(
    "Stamped ranking contract: eligibility-aware-v1 + "
    "wgcna-promoter-regulatory-v1"
  )


if __name__ == "__main__":
  main()
