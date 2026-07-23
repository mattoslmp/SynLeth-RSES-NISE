#!/usr/bin/env python3
"""Add score-internal WGCNA, methylation and validation layers to Circos rings.

The base Circos table contains every top-level RSES-Onco and functional-microniche
component. This stage adds the internal values that are combined to create WGCNA,
promoter methylation and validation/tractability components so the visualization
and Supplementary Table S47 expose every calculational layer available to the final
score.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]

ADDITIONAL_TRACKS = [
  ("A00", "Observed-domain RSES", "rses_onco", "A", "RSES-Onco", "global", 0),
  ("A10", "Genetic-screen validation", "genetic_screen", "A", "Validation/tractability", "validation_tractability", 10),
  ("A11", "Isogenic validation", "isogenic_validation", "A", "Validation/tractability", "validation_tractability", 11),
  ("A12", "In vivo validation", "in_vivo", "A", "Validation/tractability", "validation_tractability", 12),
  ("A13", "Clinical tractability", "clinical_tractability", "A", "Validation/tractability", "validation_tractability", 13),
  ("B17", "WGCNA TOM divergence", "wgcna_tom_divergence", "B", "Expression context", "wgcna_expression_network", 17),
  ("B18", "WGCNA module divergence", "wgcna_module_divergence", "B", "Expression context", "wgcna_expression_network", 18),
  ("B19", "WGCNA kME divergence", "wgcna_kME_divergence", "B", "Expression context", "wgcna_expression_network", 19),
  ("B20", "Methylation-profile divergence", "methylation_pair_profile_divergence", "B", "Regulatory network", "promoter_methylation", 20),
  ("B21", "Conditional target hypomethylation", "methylation_target_hypomethylation_support", "B", "Regulatory network", "promoter_methylation", 21),
]


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if np.isfinite(result) else None


def pair_mapping(candidates: pd.DataFrame) -> pd.DataFrame:
  result = candidates[["pair_id", "lost_gene", "target_gene"]].copy()
  result["lost_gene"] = result["lost_gene"].map(canonical_gene_name)
  result["target_gene"] = result["target_gene"].map(canonical_gene_name)
  return result.drop_duplicates("pair_id")


def aggregate_gene_track(
  data: pd.DataFrame,
  coordinates: pd.DataFrame,
  track: dict[str, object],
) -> pd.DataFrame:
  rows = []
  column = str(track["source_column"])
  for gene in coordinates["gene"].astype(str):
    subset = data.loc[
      data["lost_gene"].astype(str).eq(gene)
      | data["target_gene"].astype(str).eq(gene)
    ]
    values = (
      pd.to_numeric(subset[column], errors="coerce")
      if column in subset
      else pd.Series(dtype=float)
    ).dropna()
    rows.append({
      "gene": gene,
      **track,
      "value": float(values.max()) if not values.empty else np.nan,
      "median_value": (
        float(values.median()) if not values.empty else np.nan
      ),
      "minimum_value": (
        float(values.min()) if not values.empty else np.nan
      ),
      "observed_pair_cancer_rows": int(len(values)),
      "eligible_pair_cancer_rows": int(len(subset)),
      "evidence_status": (
        "observed" if not values.empty else "missing_or_not_eligible"
      ),
      "aggregation": (
        "maximum_across_all_associated_pair_cancer_rows"
      ),
    })
  return pd.DataFrame(rows).merge(
    coordinates[[
      "gene",
      "gene_class",
      "chromosome",
      "genomic_position",
    ]],
    on="gene",
    how="left",
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--wgcna",
    default=(
      "data/processed/regulatory/wgcna/"
      "wgcna_pair_metrics_all_cancers.tsv"
    ),
  )
  parser.add_argument(
    "--coordinates",
    default=(
      "data/processed/circos/"
      "genomic_circos_gene_coordinates.tsv"
    ),
  )
  parser.add_argument(
    "--rings",
    default=(
      "data/processed/circos/"
      "genomic_circos_ring_values.tsv"
    ),
  )
  parser.add_argument(
    "--tracks",
    default=(
      "data/processed/circos/"
      "genomic_circos_track_definitions.tsv"
    ),
  )
  args = parser.parse_args()

  paths = {
    name: resolve(value)
    for name, value in {
      "ranking": args.ranking,
      "candidates": args.candidates,
      "wgcna": args.wgcna,
      "coordinates": args.coordinates,
      "rings": args.rings,
      "tracks": args.tracks,
    }.items()
  }
  for path in paths.values():
    if not path.exists() or path.stat().st_size == 0:
      raise FileNotFoundError(path)

  ranking = pd.read_csv(paths["ranking"], sep="\t", low_memory=False)
  candidates = pd.read_csv(
    paths["candidates"],
    sep="\t",
    low_memory=False,
  )
  mapping = pair_mapping(candidates)
  validation_columns = [
    column
    for column in (
      "pair_id",
      "genetic_screen",
      "isogenic_validation",
      "in_vivo",
      "clinical_tractability",
    )
    if column in candidates
  ]
  candidate_values = candidates[validation_columns].drop_duplicates(
    "pair_id"
  ).merge(mapping, on="pair_id", how="left")

  ranking = ranking.drop(
    columns=[
      column
      for column in ("lost_gene", "target_gene")
      if column in ranking
    ],
    errors="ignore",
  ).merge(mapping, on="pair_id", how="left")
  wgcna = pd.read_csv(paths["wgcna"], sep="\t", low_memory=False)
  if not {"lost_gene", "target_gene"}.issubset(wgcna.columns):
    wgcna = wgcna.merge(mapping, on="pair_id", how="left")
  wgcna["lost_gene"] = wgcna["lost_gene"].map(canonical_gene_name)
  wgcna["target_gene"] = wgcna["target_gene"].map(
    canonical_gene_name
  )

  coordinates = pd.read_csv(
    paths["coordinates"],
    sep="\t",
    low_memory=False,
  )
  rings = pd.read_csv(paths["rings"], sep="\t", low_memory=False)
  tracks = pd.read_csv(paths["tracks"], sep="\t", low_memory=False)
  definitions = pd.DataFrame(
    ADDITIONAL_TRACKS,
    columns=[
      "track_id",
      "track_label",
      "source_column",
      "panel",
      "domain_family",
      "parent_domain",
      "ring_order",
    ],
  )
  definitions["aggregation"] = (
    "maximum_across_all_associated_pair_cancer_rows"
  )
  definitions["value_range"] = "0_to_1"
  definitions["missing_data_rule"] = (
    "missing_remains_NA_and_is_rendered_as_hollow_marker"
  )

  additions = []
  for track in definitions.to_dict("records"):
    source_column = str(track["source_column"])
    if source_column.startswith("wgcna_"):
      source = wgcna
    elif source_column in {
      "genetic_screen",
      "isogenic_validation",
      "in_vivo",
      "clinical_tractability",
    }:
      source = candidate_values
    else:
      source = ranking
    additions.append(
      aggregate_gene_track(source, coordinates, track)
    )
  new_rings = pd.concat(additions, ignore_index=True, sort=False)
  rings = pd.concat([
    rings.loc[
      ~rings["track_id"].astype(str).isin(
        set(definitions["track_id"].astype(str))
      )
    ],
    new_rings,
  ], ignore_index=True, sort=False)
  tracks = pd.concat([
    tracks.loc[
      ~tracks["track_id"].astype(str).isin(
        set(definitions["track_id"].astype(str))
      )
    ],
    definitions,
  ], ignore_index=True, sort=False)
  tracks = tracks.sort_values(["panel", "ring_order"])
  rings = rings.sort_values([
    "panel",
    "ring_order",
    "chromosome",
    "genomic_position",
  ])

  expected = set(tracks["track_id"].astype(str))
  observed = set(rings["track_id"].astype(str))
  if expected != observed:
    raise RuntimeError(
      "Circos internal-layer enrichment failed; missing tracks="
      f"{sorted(expected - observed)}"
    )
  atomic_tsv(tracks, paths["tracks"])
  atomic_tsv(rings, paths["rings"])
  print(
    f"Circos enriched to {len(tracks):,} score/domain/internal tracks "
    f"and {len(rings):,} gene-track rows."
  )


if __name__ == "__main__":
  main()
