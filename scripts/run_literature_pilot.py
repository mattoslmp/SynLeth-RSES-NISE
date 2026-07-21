#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from rses_onco.integrate import load_reference_candidates, score_candidate_table

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser(description="Score the transparent literature-anchored RSES-Onco benchmark.")
  parser.add_argument("--input", default="data/curated/synthetic_lethality_reference_pairs.tsv")
  parser.add_argument("--output", default="results/literature_anchored_candidates.tsv")
  args = parser.parse_args()
  candidates = load_reference_candidates(resolve_path(args.input))
  result = score_candidate_table(candidates)
  out = resolve_path(args.output)
  out.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(out, sep="\t", index=False)
  print(result[["pair_id", "lost_feature", "target_gene", "rses_onco", "evidence_coverage", "priority_class"]].to_string(index=False))
  print(f"\nWrote {len(result)} candidates to {out}")


if __name__ == "__main__":
  main()
