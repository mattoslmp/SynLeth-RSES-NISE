from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .integrate import load_reference_candidates, score_candidate_table


def build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(prog="rses-onco")
  sub = parser.add_subparsers(dest="command", required=True)
  score = sub.add_parser("score-literature", help="Score the bundled literature-anchored candidate table")
  score.add_argument("--input", default="data/curated/synthetic_lethality_reference_pairs.tsv")
  score.add_argument("--output", default="results/literature_anchored_candidates.tsv")
  return parser


def main() -> None:
  args = build_parser().parse_args()
  if args.command == "score-literature":
    frame = load_reference_candidates(args.input)
    result = score_candidate_table(frame)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)
    print(f"Wrote {len(result)} candidates to {output}")


if __name__ == "__main__":
  main()
