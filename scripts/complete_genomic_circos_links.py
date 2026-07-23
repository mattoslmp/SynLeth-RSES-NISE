#!/usr/bin/env python3
"""Guarantee one genomic Circos chord for every simple NISE/paralog pair.

The primary Circos builder derives link score summaries from the final ranking.
This completion stage preserves pairs that are present in the candidate universe
but have no score row, recording their score as missing instead of silently
omitting the biological relationship.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from rses_onco.utils import canonical_gene_name
from scripts.build_genomic_circos_inputs import relation_flags


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def selected_pairs(path: Path) -> pd.DataFrame:
  frame = relation_flags(
    pd.read_csv(path, sep="\t", low_memory=False)
  )
  frame["lost_gene"] = frame["lost_gene"].map(canonical_gene_name)
  frame["target_gene"] = frame["target_gene"].map(
    canonical_gene_name
  )
  return frame.loc[
    (frame["is_nise"] | frame["is_paralog"])
    & frame["pair_id"].notna()
    & frame["lost_gene"].ne("")
    & frame["target_gene"].ne("")
  ].drop_duplicates("pair_id")


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--coordinates",
    default=(
      "data/processed/circos/"
      "genomic_circos_gene_coordinates.tsv"
    ),
  )
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--links",
    default="data/processed/circos/genomic_circos_pair_links.tsv",
  )
  args = parser.parse_args()

  candidates_path = resolve(args.candidates)
  coordinates_path = resolve(args.coordinates)
  ranking_path = resolve(args.ranking)
  links_path = resolve(args.links)
  for path in (
    candidates_path,
    coordinates_path,
    ranking_path,
    links_path,
  ):
    if not path.exists() or path.stat().st_size == 0:
      raise FileNotFoundError(path)

  pairs = selected_pairs(candidates_path)
  coordinates = pd.read_csv(
    coordinates_path,
    sep="\t",
    low_memory=False,
  ).set_index("gene")
  ranking = pd.read_csv(
    ranking_path,
    sep="\t",
    low_memory=False,
  )
  links = pd.read_csv(links_path, sep="\t", low_memory=False)
  existing = set(links["pair_id"].astype(str))
  additions = []
  for pair in pairs.to_dict("records"):
    pair_id = str(pair["pair_id"])
    if pair_id in existing:
      continue
    lost = str(pair["lost_gene"])
    target = str(pair["target_gene"])
    if lost not in coordinates.index or target not in coordinates.index:
      raise RuntimeError(
        f"Pair {pair_id} has no complete Circos coordinates: "
        f"{lost} -> {target}"
      )
    subset = ranking.loc[
      ranking["pair_id"].astype(str).eq(pair_id)
    ]
    score = pd.to_numeric(
      subset.get(
        "coverage_adjusted_rses",
        pd.Series(dtype=float),
      ),
      errors="coerce",
    ).dropna()
    coverage = pd.to_numeric(
      subset.get("evidence_coverage", pd.Series(dtype=float)),
      errors="coerce",
    ).dropna()
    maximum = float(score.max()) if not score.empty else np.nan
    median = float(score.median()) if not score.empty else np.nan
    lost_coordinate = coordinates.loc[lost]
    target_coordinate = coordinates.loc[target]
    pair_class = str(pair["pair_class"])
    additions.append({
      "pair_id": pair_id,
      "lost_gene": lost,
      "target_gene": target,
      "pair_class": pair_class,
      "lost_chromosome": lost_coordinate["chromosome"],
      "lost_position": lost_coordinate["genomic_position"],
      "target_chromosome": target_coordinate["chromosome"],
      "target_position": target_coordinate["genomic_position"],
      "cancers": (
        ";".join(sorted(set(subset["cancer"].dropna().astype(str))))
        if "cancer" in subset
        else ""
      ),
      "maximum_coverage_adjusted_rses": maximum,
      "median_coverage_adjusted_rses": median,
      "maximum_evidence_coverage": (
        float(coverage.max()) if not coverage.empty else np.nan
      ),
      "link_width": (
        0.25 + 2.75 * maximum
        if np.isfinite(maximum)
        else 0.25
      ),
      "link_alpha": (
        0.06 + 0.44 * maximum
        if np.isfinite(maximum)
        else 0.06
      ),
      "link_color": (
        "#C62828" if "NISE" in pair_class else "#111111"
      ),
      "link_status": (
        "available" if np.isfinite(maximum) else "score_missing"
      ),
    })

  result = pd.concat(
    [links, pd.DataFrame(additions)],
    ignore_index=True,
    sort=False,
  )
  result = result.drop_duplicates("pair_id").sort_values(
    ["maximum_coverage_adjusted_rses", "pair_id"],
    ascending=[False, True],
    na_position="last",
  )
  expected = set(pairs["pair_id"].astype(str))
  observed = set(result["pair_id"].astype(str))
  if expected != observed:
    raise RuntimeError(
      "Circos link completion failed; missing="
      f"{sorted(expected - observed)[:50]}, extra="
      f"{sorted(observed - expected)[:50]}"
    )
  atomic_tsv(result, links_path)
  print(
    f"Circos links complete: {len(result):,}/{len(expected):,} "
    f"NISE/paralog pairs; added {len(additions):,}."
  )


if __name__ == "__main__":
  main()
