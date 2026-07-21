from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .utils import clamp01


DEFAULT_WEIGHTS = {
  "tumor_event": 0.18,
  "dependency": 0.25,
  "selectivity": 0.17,
  "expression_compensation": 0.10,
  "functional_relation": 0.12,
  "validation_tractability": 0.18,
}


@dataclass(frozen=True)
class ScoreResult:
  observed_score: float
  coverage: float
  adjusted_score: float
  n_domains: int
  interpretation: str


def rses_onco_score(
  components: Mapping[str, float | None],
  weights: Mapping[str, float] | None = None,
  min_domains_for_high: int = 4,
  min_coverage_for_high: float = 0.67,
) -> ScoreResult:
  """Calculate a coverage-aware RSES-Onco score.

  Missing values are omitted from the observed-domain denominator and reduce
  coverage; they are never converted to zero. This mirrors the conservative
  missing-data rule in the original RSES framework.
  """
  weights = dict(weights or DEFAULT_WEIGHTS)
  total_weight = float(sum(weights.values()))
  numerator = 0.0
  observed_weight = 0.0
  n_domains = 0
  for domain, weight in weights.items():
    value = clamp01(components.get(domain))
    if value is None:
      continue
    numerator += weight * value
    observed_weight += weight
    n_domains += 1
  observed = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / total_weight if total_weight else float("nan")
  adjusted = observed * coverage if np.isfinite(observed) else float("nan")
  if not np.isfinite(observed):
    label = "insufficient evidence"
  elif observed >= 0.72 and coverage >= min_coverage_for_high and n_domains >= min_domains_for_high:
    label = "high priority"
  elif observed >= 0.48:
    label = "moderate priority"
  else:
    label = "exploratory"
  return ScoreResult(observed, coverage, adjusted, n_domains, label)


def literature_prior_components(row: Mapping[str, object]) -> dict[str, float]:
  """Map a literature-evidence row to transparent RSES-Onco pilot domains."""
  validation = np.mean([
    float(row["genetic_screen"]),
    float(row["isogenic_validation"]),
    float(row["in_vivo"]),
    float(row["clinical_tractability"]),
  ])
  lineage_flags = [float(row.get(k, 0) or 0) for k in ("colon", "stomach", "lung")]
  return {
    "tumor_event": float(row["lineage_relevance"]),
    "dependency": float(row["genetic_screen"]),
    "selectivity": float(row["isogenic_validation"]),
    "expression_compensation": None,
    "functional_relation": float(row["relation_confidence"]),
    "validation_tractability": float(validation),
  }
