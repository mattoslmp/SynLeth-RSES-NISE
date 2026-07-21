#!/usr/bin/env python3
"""Build the expanded RSES-Onco candidate universe.

The mandatory universe contains every directed cross-cluster human NISE pair plus
all curated benchmark vulnerabilities. Optional standardized tables can add
Ensembl paralogs, pathway backups, collateral-deletion pairs, downstream
vulnerabilities or other explicitly sourced classes.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from rses_onco.expanded import (
  build_directed_nise_candidates,
  class_member_inventory,
  load_optional_table,
  merge_candidate_sources,
)

ROOT = Path(__file__).resolve().parents[1]
SIMPLE_LOSS_PATTERNS = [
  re.compile(r"^([A-Za-z0-9-]+)$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+homozygous deletion$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss/low expression$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss or low expression$", re.I),
  re.compile(r"^([A-Za-z0-9-]+)\s+loss/low activity$", re.I),
]


def resolve_path(value: str | None) -> Path | None:
  if value is None:
    return None
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def simple_lost_gene(value: object) -> str | None:
  text = str(value).strip()
  for pattern in SIMPLE_LOSS_PATTERNS:
    match = pattern.fullmatch(text)
    if match:
      return match.group(1).upper()
  return None


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--nise-pairs",
    default="data/processed/human_nise_all_within_activity_pairs_2017.tsv",
  )
  parser.add_argument(
    "--benchmarks",
    default="data/curated/synthetic_lethality_reference_pairs.tsv",
  )
  parser.add_argument(
    "--additional",
    action="append",
    default=[],
    help=(
      "Optional standardized TSV/CSV with pair_id, lost_feature or lost_gene, "
      "target_gene, source_class/relation_type and evidence fields. Repeatable."
    ),
  )
  parser.add_argument(
    "--output",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--members-output",
    default="data/processed/expanded_class_member_inventory.tsv",
  )
  args = parser.parse_args()

  nise_pairs = pd.read_csv(resolve_path(args.nise_pairs), sep="\t")
  benchmarks = pd.read_csv(resolve_path(args.benchmarks), sep="\t")
  if "source_class" not in benchmarks:
    benchmarks["source_class"] = benchmarks["relation_type"].astype(str)
  if "lost_gene" not in benchmarks:
    benchmarks["lost_gene"] = benchmarks["lost_feature"].map(simple_lost_gene)

  nise_candidates = build_directed_nise_candidates(nise_pairs)
  additional = []
  for value in args.additional:
    table = load_optional_table(resolve_path(value))
    if table is not None:
      additional.append(table)

  universe = merge_candidate_sources(nise_candidates, benchmarks, additional)
  members = class_member_inventory(universe)

  output = resolve_path(args.output)
  members_output = resolve_path(args.members_output)
  assert output is not None and members_output is not None
  output.parent.mkdir(parents=True, exist_ok=True)
  members_output.parent.mkdir(parents=True, exist_ok=True)
  universe.to_csv(output, sep="\t", index=False)
  members.to_csv(members_output, sep="\t", index=False)

  directed_nise = int((universe["source_class"] == "NISE").sum())
  nise_genes = set(
    members.loc[members["source_class"] == "NISE", "gene"].astype(str)
  )
  complex_benchmarks = int(benchmarks["lost_gene"].isna().sum())
  print(f"Wrote {len(universe):,} candidate directions to {output}")
  print(f"Directed NISE candidates: {directed_nise:,}")
  print(f"Unique NISE genes represented: {len(nise_genes):,}")
  print(f"Complex benchmark biomarkers retained without fake lost genes: {complex_benchmarks:,}")
  print(f"Wrote class-member inventory to {members_output}")


if __name__ == "__main__":
  main()
