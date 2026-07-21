from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np
import pandas as pd

from .utils import clamp01


PHARMACOLOGY_WEIGHTS = {
  "target_tractability": 0.18,
  "direct_interaction": 0.18,
  "compound_potency": 0.18,
  "clinical_maturity": 0.14,
  "cancer_relevance": 0.12,
  "biomarker_selective_response": 0.20,
}


@dataclass(frozen=True)
class PharmacologyScore:
  observed_score: float
  coverage: float
  adjusted_score: float
  n_domains: int
  interpretation: str


def coverage_aware_pharmacology_score(
  components: Mapping[str, float | None],
  weights: Mapping[str, float] | None = None,
) -> PharmacologyScore:
  weights = dict(weights or PHARMACOLOGY_WEIGHTS)
  total_weight = float(sum(weights.values()))
  observed_weight = 0.0
  numerator = 0.0
  n_domains = 0
  for domain, weight in weights.items():
    value = clamp01(components.get(domain))
    if value is None:
      continue
    observed_weight += float(weight)
    numerator += float(weight) * float(value)
    n_domains += 1
  observed = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / total_weight if total_weight else float("nan")
  adjusted = observed * coverage if np.isfinite(observed) else float("nan")
  if not np.isfinite(observed):
    interpretation = "insufficient pharmacology evidence"
  elif observed >= 0.75 and coverage >= 0.65 and n_domains >= 4:
    interpretation = "high experimental actionability"
  elif observed >= 0.50:
    interpretation = "moderate experimental actionability"
  else:
    interpretation = "early pharmacology hypothesis"
  return PharmacologyScore(observed, coverage, adjusted, n_domains, interpretation)


def normalize_pharos_tdl(value: object) -> float | None:
  mapping = {
    "TCLIN": 1.0,
    "TCHEM": 0.82,
    "TBIO": 0.42,
    "TDARK": 0.12,
  }
  key = str(value).strip().upper()
  return mapping.get(key)


def normalize_open_targets_tractability(values: object) -> float | None:
  if values is None:
    return None
  if isinstance(values, Mapping):
    values = [values]
  scores = []
  if isinstance(values, (list, tuple)):
    for item in values:
      if not isinstance(item, Mapping):
        continue
      raw_value = item.get("value")
      if isinstance(raw_value, bool):
        scores.append(1.0 if raw_value else 0.0)
      else:
        try:
          scores.append(float(raw_value))
        except (TypeError, ValueError):
          continue
  return max(scores) if scores else None


def normalize_dgidb_interaction_score(value: object) -> float | None:
  try:
    score = max(0.0, float(value))
  except (TypeError, ValueError):
    return None
  return float(1.0 - math.exp(-score / 5.0))


def normalize_pchembl(value: object) -> float | None:
  try:
    pchembl = float(value)
  except (TypeError, ValueError):
    return None
  return float(np.clip((pchembl - 5.0) / 4.0, 0, 1))


def normalize_phase(value: object) -> float | None:
  try:
    phase = float(value)
  except (TypeError, ValueError):
    return None
  return float(np.clip(phase / 4.0, 0, 1))


def normalize_source_count(value: object, saturation: float = 4.0) -> float | None:
  try:
    count = max(0.0, float(value))
  except (TypeError, ValueError):
    return None
  return float(np.clip(count / saturation, 0, 1))


def normalize_sensitivity_selectivity(
  delta_response: object,
  lower_is_more_sensitive: bool = True,
  saturation: float = 1.0,
) -> float | None:
  try:
    delta = float(delta_response)
  except (TypeError, ValueError):
    return None
  supportive = -delta if lower_is_more_sensitive else delta
  return float(np.clip(supportive / saturation, 0, 1))


def therapeutic_hypothesis_score(
  vulnerability_adjusted: object,
  pharmacology_adjusted: object,
) -> float | None:
  """Combine vulnerability and pharmacology by geometric concordance.

  A geometric combination prevents a highly druggable but biologically weak
  target, or a strong vulnerability without any actionable compound evidence,
  from being presented as a mature therapeutic hypothesis.
  """
  try:
    vulnerability = float(vulnerability_adjusted)
    pharmacology = float(pharmacology_adjusted)
  except (TypeError, ValueError):
    return None
  if not np.isfinite(vulnerability) or not np.isfinite(pharmacology):
    return None
  return float(math.sqrt(max(0.0, vulnerability) * max(0.0, pharmacology)))


def evidence_rows_to_components(rows: pd.DataFrame) -> dict[str, float | None]:
  if rows.empty:
    return {domain: None for domain in PHARMACOLOGY_WEIGHTS}

  tractability_values = []
  direct_values = []
  potency_values = []
  phase_values = []
  cancer_values = []
  sensitivity_values = []

  for record in rows.to_dict("records"):
    source = str(record.get("source", "")).casefold()
    if source == "pharos":
      value = normalize_pharos_tdl(record.get("target_development_level"))
      if value is not None:
        tractability_values.append(value)
    if source == "open_targets":
      try:
        value = float(record.get("tractability_score"))
      except (TypeError, ValueError):
        value = None
      if value is not None and np.isfinite(value):
        tractability_values.append(value)
      phase = normalize_phase(record.get("max_phase"))
      if phase is not None:
        phase_values.append(phase)
      if record.get("disease_name"):
        cancer_values.append(1.0)
    if source == "dgidb":
      value = normalize_dgidb_interaction_score(record.get("interaction_score"))
      if value is not None:
        direct_values.append(value)
      if str(record.get("approved", "")).casefold() in {"1", "true", "yes"}:
        phase_values.append(1.0)
    if source == "chembl":
      value = normalize_pchembl(record.get("pchembl_value"))
      if value is not None:
        potency_values.append(value)
      phase = normalize_phase(record.get("max_phase"))
      if phase is not None:
        phase_values.append(phase)
      if record.get("action_type") or record.get("mechanism_of_action"):
        direct_values.append(1.0)
    if source == "civic":
      if str(record.get("civic_gene_record", "")).casefold() in {"1", "true", "yes"}:
        cancer_values.append(0.60)
    if source in {"prism", "gdsc", "ctrp"}:
      value = normalize_sensitivity_selectivity(
        record.get("delta_response"),
        lower_is_more_sensitive=str(record.get("lower_is_more_sensitive", "true")).casefold()
        in {"1", "true", "yes"},
      )
      if value is not None:
        sensitivity_values.append(value)

  def maximum(values: list[float]) -> float | None:
    return max(values) if values else None

  return {
    "target_tractability": maximum(tractability_values),
    "direct_interaction": maximum(direct_values),
    "compound_potency": maximum(potency_values),
    "clinical_maturity": maximum(phase_values),
    "cancer_relevance": maximum(cancer_values),
    "biomarker_selective_response": maximum(sensitivity_values),
  }
