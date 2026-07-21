from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

from .utils import canonical_gene_name


@dataclass(frozen=True)
class PairSetMetrics:
  gene_a: str
  gene_b: str
  shared_count: int
  exclusive_a_count: int
  exclusive_b_count: int
  jaccard: float | None
  divergence: float | None


def _split_locations(value: object) -> set[str]:
  if value is None or (isinstance(value, float) and np.isnan(value)):
    return set()
  text = str(value).strip()
  if not text:
    return set()
  normalized = text.replace(";", ",")
  return {
    item.strip().casefold()
    for item in normalized.split(",")
    if item.strip()
  }


def set_metrics(gene_a: str, gene_b: str, values_a: Iterable[str], values_b: Iterable[str]) -> PairSetMetrics:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  set_a = {str(value) for value in values_a if str(value)}
  set_b = {str(value) for value in values_b if str(value)}
  union = set_a | set_b
  if not union:
    return PairSetMetrics(gene_a, gene_b, 0, 0, 0, None, None)
  intersection = set_a & set_b
  jaccard = len(intersection) / len(union)
  return PairSetMetrics(
    gene_a=gene_a,
    gene_b=gene_b,
    shared_count=len(intersection),
    exclusive_a_count=len(set_a - set_b),
    exclusive_b_count=len(set_b - set_a),
    jaccard=float(jaccard),
    divergence=float(1.0 - jaccard),
  )


def string_neighbor_map(edges: pd.DataFrame, minimum_score: float = 0.7) -> dict[str, set[str]]:
  aliases = {
    "gene_a": ["preferredName_A", "preferred_name_a", "gene_a", "protein1"],
    "gene_b": ["preferredName_B", "preferred_name_b", "gene_b", "protein2"],
    "score": ["score", "combined_score", "combinedScore"],
  }

  def resolve(options: list[str]) -> str:
    for option in options:
      if option in edges.columns:
        return option
    raise ValueError(f"Could not resolve any of {options} in STRING columns")

  a_column = resolve(aliases["gene_a"])
  b_column = resolve(aliases["gene_b"])
  score_column = resolve(aliases["score"])
  neighbors: dict[str, set[str]] = defaultdict(set)
  for row in edges[[a_column, b_column, score_column]].itertuples(index=False, name=None):
    gene_a, gene_b, score = row
    try:
      score_value = float(score)
    except (TypeError, ValueError):
      continue
    if score_value > 1.0:
      score_value /= 1000.0
    if score_value < minimum_score:
      continue
    gene_a = canonical_gene_name(gene_a)
    gene_b = canonical_gene_name(gene_b)
    if not gene_a or not gene_b or gene_a == gene_b:
      continue
    neighbors[gene_a].add(gene_b)
    neighbors[gene_b].add(gene_a)
  return dict(neighbors)


def string_pair_metrics(edges: pd.DataFrame, gene_a: str, gene_b: str) -> PairSetMetrics:
  neighbors = string_neighbor_map(edges)
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  return set_metrics(gene_a, gene_b, neighbors.get(gene_a, set()), neighbors.get(gene_b, set()))


def regulatory_target_map(interactions: pd.DataFrame) -> dict[str, set[str]]:
  source_column = next(
    (column for column in ("source_genesymbol", "source", "tf") if column in interactions.columns),
    None,
  )
  target_column = next(
    (column for column in ("target_genesymbol", "target", "gene") if column in interactions.columns),
    None,
  )
  if source_column is None or target_column is None:
    raise ValueError("Regulatory table requires source and target gene-symbol columns")
  regulators: dict[str, set[str]] = defaultdict(set)
  for source, target in interactions[[source_column, target_column]].itertuples(index=False, name=None):
    source = canonical_gene_name(source)
    target = canonical_gene_name(target)
    if source and target:
      regulators[target].add(source)
  return dict(regulators)


def regulatory_pair_metrics(interactions: pd.DataFrame, gene_a: str, gene_b: str) -> PairSetMetrics:
  regulators = regulatory_target_map(interactions)
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  return set_metrics(gene_a, gene_b, regulators.get(gene_a, set()), regulators.get(gene_b, set()))


def hpa_location_map(localization: pd.DataFrame) -> dict[str, set[str]]:
  gene_column = next(
    (column for column in ("Gene name", "Gene_name", "gene", "gene_symbol") if column in localization.columns),
    None,
  )
  if gene_column is None:
    raise ValueError("HPA localization table lacks a gene-symbol column")
  location_columns = [
    column for column in (
      "Enhanced", "Supported", "Approved", "Validated", "Uncertain",
      "Main location", "Additional location",
    )
    if column in localization.columns
  ]
  if not location_columns:
    raise ValueError("HPA localization table lacks recognized location columns")
  result: dict[str, set[str]] = defaultdict(set)
  for record in localization[[gene_column, *location_columns]].to_dict("records"):
    gene = canonical_gene_name(record[gene_column])
    if not gene:
      continue
    for column in location_columns:
      result[gene].update(_split_locations(record.get(column)))
  return dict(result)


def localization_pair_metrics(localization: pd.DataFrame, gene_a: str, gene_b: str) -> PairSetMetrics:
  locations = hpa_location_map(localization)
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  return set_metrics(gene_a, gene_b, locations.get(gene_a, set()), locations.get(gene_b, set()))


def direct_string_score(edges: pd.DataFrame, gene_a: str, gene_b: str) -> float | None:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  if not gene_a or not gene_b:
    return None
  aliases = {
    "gene_a": next((x for x in ("preferredName_A", "preferred_name_a", "gene_a", "protein1") if x in edges.columns), None),
    "gene_b": next((x for x in ("preferredName_B", "preferred_name_b", "gene_b", "protein2") if x in edges.columns), None),
    "score": next((x for x in ("score", "combined_score", "combinedScore") if x in edges.columns), None),
  }
  if any(value is None for value in aliases.values()):
    return None
  matches = []
  for a_value, b_value, score in edges[
    [aliases["gene_a"], aliases["gene_b"], aliases["score"]]
  ].itertuples(index=False, name=None):
    a_name = canonical_gene_name(a_value)
    b_name = canonical_gene_name(b_value)
    if {a_name, b_name} != {gene_a, gene_b}:
      continue
    try:
      score_value = float(score)
    except (TypeError, ValueError):
      continue
    if score_value > 1.0:
      score_value /= 1000.0
    matches.append(score_value)
  return float(max(matches)) if matches else None
