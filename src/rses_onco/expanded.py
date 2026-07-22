from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .depmap import cancer_model_ids
from .utils import canonical_gene_name, clamp01


FUNCTIONAL_MICRONICHE_WEIGHTS = {
  "expression_context": 0.20,
  "localization": 0.15,
  "biochemical_structural": 0.15,
  "genetic_phenotype": 0.20,
  "interaction_network": 0.15,
  "regulatory_network": 0.15,
}

EXPANDED_ONCO_WEIGHTS = {
  "tumor_event": 0.16,
  "dependency": 0.22,
  "selectivity": 0.14,
  "expression_compensation": 0.08,
  "functional_relation": 0.06,
  "functional_microniche": 0.16,
  "validation_tractability": 0.18,
}


@dataclass(frozen=True)
class CoverageAwareResult:
  observed_score: float
  coverage: float
  adjusted_score: float
  n_domains: int
  eligible_domains: int = 0
  observed_weight: float = float("nan")
  eligible_weight: float = float("nan")


@dataclass(frozen=True)
class PairProfileMetrics:
  gene_a: str
  gene_b: str
  cancer: str
  n_models: int
  spearman_rho: float | None
  median_absolute_difference: float | None
  set_jaccard: float | None
  divergence: float | None


def coverage_aware_score(
  components: Mapping[str, float | None],
  weights: Mapping[str, float],
  eligible_domains: Iterable[str] | None = None,
) -> CoverageAwareResult:
  """Calculate an eligibility-aware coverage-adjusted score.

  Missing eligible components are omitted from the observed-score denominator and
  penalized through coverage. Ineligible domains are excluded from both the observed
  and eligible denominators. An observed value of zero remains a real zero.
  """
  eligible = (
    set(weights)
    if eligible_domains is None
    else {str(domain) for domain in eligible_domains if str(domain) in weights}
  )
  eligible_weight = float(sum(weights[domain] for domain in eligible))
  numerator = 0.0
  observed_weight = 0.0
  n_domains = 0
  for domain, weight in weights.items():
    if domain not in eligible:
      continue
    value = clamp01(components.get(domain))
    if value is None:
      continue
    numerator += float(weight) * float(value)
    observed_weight += float(weight)
    n_domains += 1
  observed = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / eligible_weight if eligible_weight else float("nan")
  adjusted = observed * coverage if np.isfinite(observed) and np.isfinite(coverage) else float("nan")
  return CoverageAwareResult(
    observed_score=observed,
    coverage=coverage,
    adjusted_score=adjusted,
    n_domains=n_domains,
    eligible_domains=len(eligible),
    observed_weight=observed_weight,
    eligible_weight=eligible_weight,
  )


def functional_microniche_score(
  components: Mapping[str, float | None],
  eligible_domains: Iterable[str] | None = None,
) -> CoverageAwareResult:
  return coverage_aware_score(
    components,
    FUNCTIONAL_MICRONICHE_WEIGHTS,
    eligible_domains=eligible_domains,
  )


def expanded_onco_score(
  components: Mapping[str, float | None],
  eligible_domains: Iterable[str] | None = None,
) -> CoverageAwareResult:
  return coverage_aware_score(
    components,
    EXPANDED_ONCO_WEIGHTS,
    eligible_domains=eligible_domains,
  )


def _normalized_pair_table(
  frame: pd.DataFrame,
  models: pd.DataFrame,
  gene_a: str,
  gene_b: str,
  cancer: str,
) -> pd.DataFrame:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  if gene_a not in frame or gene_b not in frame:
    return pd.DataFrame(columns=["ModelID", gene_a, gene_b])
  selected = set(cancer_model_ids(models, cancer))
  if not selected:
    return pd.DataFrame(columns=["ModelID", gene_a, gene_b])
  table = frame[["ModelID", gene_a, gene_b]].copy()
  table = table.loc[table["ModelID"].astype(str).isin(selected)]
  table[gene_a] = pd.to_numeric(table[gene_a], errors="coerce")
  table[gene_b] = pd.to_numeric(table[gene_b], errors="coerce")
  return table.dropna(subset=[gene_a, gene_b])


def expression_profile_metrics(
  expression: pd.DataFrame,
  models: pd.DataFrame,
  gene_a: str,
  gene_b: str,
  cancer: str,
  separation_saturation: float = 2.0,
) -> PairProfileMetrics:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  table = _normalized_pair_table(expression, models, gene_a, gene_b, cancer)
  if len(table) < 3:
    return PairProfileMetrics(gene_a, gene_b, cancer, len(table), None, None, None, None)
  rho = float(spearmanr(table[gene_a], table[gene_b], nan_policy="omit").statistic)
  if not np.isfinite(rho):
    rho = None
  median_difference = float(np.median(np.abs(table[gene_a] - table[gene_b])))
  correlation_divergence = None if rho is None else (1.0 - rho) / 2.0
  separation = float(np.clip(median_difference / separation_saturation, 0, 1))
  available = [value for value in (correlation_divergence, separation) if value is not None]
  divergence = float(np.mean(available)) if available else None
  return PairProfileMetrics(
    gene_a, gene_b, cancer, len(table), rho, median_difference, None, divergence,
  )


def phenotype_profile_metrics(
  gene_effect: pd.DataFrame,
  models: pd.DataFrame,
  gene_a: str,
  gene_b: str,
  cancer: str,
  dependency_threshold: float = -0.5,
  separation_saturation: float = 1.0,
) -> PairProfileMetrics:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  table = _normalized_pair_table(gene_effect, models, gene_a, gene_b, cancer)
  if len(table) < 3:
    return PairProfileMetrics(gene_a, gene_b, cancer, len(table), None, None, None, None)
  rho = float(spearmanr(table[gene_a], table[gene_b], nan_policy="omit").statistic)
  if not np.isfinite(rho):
    rho = None
  median_difference = float(np.median(np.abs(table[gene_a] - table[gene_b])))
  set_a = set(table.loc[table[gene_a] < dependency_threshold, "ModelID"].astype(str))
  set_b = set(table.loc[table[gene_b] < dependency_threshold, "ModelID"].astype(str))
  union = set_a | set_b
  jaccard = float(len(set_a & set_b) / len(union)) if union else None
  correlation_divergence = None if rho is None else (1.0 - rho) / 2.0
  separation = float(np.clip(median_difference / separation_saturation, 0, 1))
  set_divergence = None if jaccard is None else 1.0 - jaccard
  available = [
    value for value in (correlation_divergence, separation, set_divergence)
    if value is not None
  ]
  divergence = float(np.mean(available)) if available else None
  return PairProfileMetrics(
    gene_a, gene_b, cancer, len(table), rho, median_difference, jaccard, divergence,
  )
