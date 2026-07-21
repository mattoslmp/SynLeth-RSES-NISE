#!/usr/bin/env python3
"""Build the expanded RSES-Onco candidate universe.

The mandatory universe contains every directed cross-cluster human NISE pair plus
all curated benchmark vulnerabilities. Optional standardized tables can add
Ensembl paralogs, pathway backups, collateral-deletion pairs, downstream
vulnerabilities or other explicitly sourced classes.

Composite biological loss states are preserved as composite features. In contrast,
an explicitly gene-valued target field such as ``CDK4/CDK6`` is expanded into one
atomic target hypothesis per gene while retaining the original target feature and
pair-group provenance.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

from rses_onco.expanded import (
  build_directed_nise_candidates,
  class_member_inventory,
  load_optional_table,
  merge_candidate_sources,
)
from rses_onco.utils import atomic_gene_symbols, canonical_gene_name

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


def expand_explicit_target_gene_lists(frame: pd.DataFrame) -> pd.DataFrame:
  """Expand delimiter-separated target genes without altering composite loss states.

  A target field is expanded only when every delimiter-separated component is a
  valid atomic HGNC-style symbol. The original target text and original pair ID are
  retained in ``target_feature`` and ``target_group_id``. Prose or ambiguous fields
  are left unchanged so the pipeline never invents a partial target set.
  """
  if frame.empty or "target_gene" not in frame:
    return frame.copy()

  rows: list[dict[str, object]] = []
  for record in frame.to_dict("records"):
    raw_target = record.get("target_gene")
    components = atomic_gene_symbols(raw_target)
    original_pair_id = str(record.get("pair_id") or "")

    if len(components) <= 1:
      row = dict(record)
      if components:
        row["target_gene"] = components[0]
      else:
        row["target_gene"] = canonical_gene_name(raw_target)
      row.setdefault("target_feature", str(raw_target or ""))
      row.setdefault("target_group_id", original_pair_id)
      row.setdefault("target_component_index", 1)
      row.setdefault("target_component_count", 1)
      rows.append(row)
      continue

    for component_index, component in enumerate(components, start=1):
      row = dict(record)
      row["target_feature"] = str(raw_target)
      row["target_group_id"] = original_pair_id
      row["target_component_index"] = component_index
      row["target_component_count"] = len(components)
      row["target_gene"] = component
      row["pair_id"] = (
        f"{original_pair_id}__TARGET_{component}"
        if original_pair_id
        else f"TARGET_{component_index}_{component}"
      )
      rows.append(row)

  return pd.DataFrame(rows)


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
  benchmarks = expand_explicit_target_gene_lists(benchmarks)

  nise_candidates = build_directed_nise_candidates(nise_pairs)
  additional = []
  for value in args.additional:
    try:
      table = load_optional_table(resolve_path(value))
    except EmptyDataError:
      table = None
    if table is not None and not table.empty:
      additional.append(expand_explicit_target_gene_lists(table))

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
  expanded_target_groups = int(
    benchmarks.get("target_component_count", pd.Series(dtype=int))
      .fillna(1).astype(int).gt(1).sum()
  )
  print(f"Wrote {len(universe):,} candidate directions to {output}")
  print(f"Directed NISE candidates: {directed_nise:,}")
  print(f"Unique NISE genes represented: {len(nise_genes):,}")
  print(f"Complex benchmark biomarkers retained without fake lost genes: {complex_benchmarks:,}")
  print(f"Atomic rows derived from composite target lists: {expanded_target_groups:,}")
  print(f"Additional non-empty source catalogs: {len(additional):,}")
  print(f"Wrote class-member inventory to {members_output}")


if __name__ == "__main__":
  main()
