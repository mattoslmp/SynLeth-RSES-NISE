from __future__ import annotations

from datetime import datetime, timezone
from typing import Mapping

import numpy as np
import pandas as pd

from .audit import (
  DomainSpec,
  candidate_identity,
  domain_specs,
  evidence_state,
  metadata_for_source,
  numeric,
  present,
  score_family_summary,
  source_metadata,
)


def explicit_boolean(value: object) -> bool | None:
  if not present(value):
    return None
  if isinstance(value, bool):
    return value
  text = str(value).strip().casefold()
  if text in {"1", "true", "yes", "eligible"}:
    return True
  if text in {"0", "false", "no", "not_eligible", "ineligible"}:
    return False
  try:
    return bool(int(float(value)))
  except (TypeError, ValueError):
    return None


def domain_eligibility(
  record: Mapping[str, object],
  spec: DomainSpec,
) -> tuple[bool, str]:
  explicit_columns = (
    f"eligible_component_{spec.domain}",
    f"component_{spec.domain}_eligible",
    f"eligible_microniche_{spec.domain}",
    f"microniche_{spec.domain}_eligible",
    f"eligible_{spec.domain}",
    f"{spec.domain}_eligible",
  )
  for column in explicit_columns:
    if column in record:
      parsed = explicit_boolean(record.get(column))
      if parsed is not None:
        return (
          parsed,
          "explicit_pipeline_eligibility"
          if parsed
          else "explicitly_not_eligible_for_hypothesis_type",
        )

  identity = candidate_identity(record)
  gene_pair = identity["hypothesis_type"] == "gene"
  if spec.family == "Functional microniche" and not gene_pair:
    return False, "composite_event_not_executable_as_two_protein_microniche"
  if spec.domain in {
    "tumor_event",
    "dependency",
    "selectivity",
    "expression_compensation",
  } and not gene_pair:
    return False, "composite_event_not_executable_as_single_gene_analysis"
  return True, "eligible_by_candidate_and_domain_definition"


def build_candidate_domain_audit(
  ranking: pd.DataFrame,
  metadata_root,
) -> pd.DataFrame:
  metadata = source_metadata(metadata_root)
  rows: list[dict[str, object]] = []
  specs = domain_specs()
  for record in ranking.to_dict("records"):
    identity = candidate_identity(record)
    eligibility = {
      (spec.family, spec.domain): domain_eligibility(record, spec)
      for spec in specs
    }
    eligible_weights = {
      family: float(sum(
        spec.weight
        for spec in specs
        if spec.family == family
        and eligibility[(spec.family, spec.domain)][0]
      ))
      for family in {spec.family for spec in specs}
    }
    eligible_counts = {
      family: int(sum(
        eligibility[(spec.family, spec.domain)][0]
        for spec in specs
        if spec.family == family
      ))
      for family in {spec.family for spec in specs}
    }

    for spec in specs:
      eligible, eligibility_reason = eligibility[(spec.family, spec.domain)]
      original = numeric(record.get(spec.column))
      normalized = None if original is None else float(np.clip(original, 0, 1))
      state, reason = evidence_state(record, spec, eligible, normalized)
      if not eligible:
        state = "not_eligible"
        reason = eligibility_reason
        normalized = None
        original = None
      observed, coverage, adjusted, n_domains = score_family_summary(
        record,
        spec.family,
      )
      eligible_weight = eligible_weights[spec.family]
      contribution = (
        spec.weight * normalized / eligible_weight
        if eligible and normalized is not None and eligible_weight > 0
        else None
      )
      source_info = metadata_for_source(spec.source, metadata)
      evidence_count = record.get(f"{spec.domain}_evidence_count")
      independent_count = record.get(
        f"{spec.domain}_independent_evidence_count"
      )
      rows.append({
        **identity,
        "cancer": record.get("cancer"),
        "mechanistic_class": record.get("source_class") or record.get("relation_type"),
        "score_comparability_group": record.get("score_comparability_group"),
        "scoring_semantics_version": record.get("scoring_semantics_version"),
        "domain_family": spec.family,
        "domain": spec.domain,
        "domain_label": spec.label,
        "internal_component_column": spec.column,
        "eligible": bool(eligible),
        "eligibility_reason": eligibility_reason,
        "family_eligible_domains": eligible_counts[spec.family],
        "family_eligible_weight": eligible_weight,
        "evidence_state": state,
        "evidence_present": normalized is not None,
        "absence_reason": "" if normalized is not None else reason,
        "evidence_source": spec.source,
        "evidence_role": spec.role,
        "source_version": source_info.get("source_version", "not_recorded"),
        "access_date": source_info.get("access_date", "not_recorded"),
        "source_status": source_info.get("source_status", "not_recorded"),
        "source_metadata_file": source_info.get("source_metadata_file", ""),
        "source_url": source_info.get("source_url", "not_recorded"),
        "evidence_count": evidence_count if present(evidence_count) else np.nan,
        "independent_evidence_count": (
          independent_count if present(independent_count) else np.nan
        ),
        "component_original": original,
        "component_normalized": normalized,
        "component_coverage_adjusted": contribution,
        "domain_weight": spec.weight,
        "weighted_numerator_contribution": (
          spec.weight * normalized
          if normalized is not None
          else None
        ),
        "final_score_contribution": contribution,
        "family_observed_score": observed,
        "family_coverage": coverage,
        "family_adjusted_score": adjusted,
        "family_observed_domains": n_domains,
        "audit_generated_at_utc": datetime.now(timezone.utc).isoformat(),
      })
  return pd.DataFrame(rows)


def score_decomposition(audit: pd.DataFrame) -> pd.DataFrame:
  keys = [
    "candidate_id",
    "hypothesis_direction",
    "origin_entity",
    "target_entity",
    "hypothesis_type",
    "cancer",
    "mechanistic_class",
    "score_comparability_group",
    "domain_family",
  ]
  summaries: list[dict[str, object]] = []
  for values, group in audit.groupby(keys, dropna=False):
    row = dict(zip(keys, values))
    eligible = group.loc[group["eligible"].fillna(False)].copy()
    observed = eligible.loc[eligible["component_normalized"].notna()].copy()
    numerator = float(
      (observed["domain_weight"] * observed["component_normalized"]).sum()
    )
    observed_weight = float(observed["domain_weight"].sum())
    eligible_weight = float(eligible["domain_weight"].sum())
    raw = numerator / observed_weight if observed_weight else float("nan")
    coverage = observed_weight / eligible_weight if eligible_weight else float("nan")
    adjusted = numerator / eligible_weight if eligible_weight else float("nan")
    first = group.iloc[0]
    row.update({
      "recomputed_weighted_numerator": numerator,
      "recomputed_observed_weight": observed_weight,
      "recomputed_eligible_weight": eligible_weight,
      "recomputed_observed_score": raw,
      "recomputed_coverage": coverage,
      "recomputed_adjusted_score": adjusted,
      "recomputed_n_domains": int(len(observed)),
      "recomputed_eligible_domains": int(len(eligible)),
      "pipeline_observed_score": first["family_observed_score"],
      "pipeline_coverage": first["family_coverage"],
      "pipeline_adjusted_score": first["family_adjusted_score"],
      "pipeline_n_domains": first["family_observed_domains"],
      "pipeline_eligible_domains": first["family_eligible_domains"],
      "formula": (
        "adjusted_score = sum(w_i*x_i for observed eligible i) / "
        "sum(w_i for eligible i)"
      ),
      "observed_score_formula": (
        "observed_score = sum(w_i*x_i for observed eligible i) / "
        "sum(w_i for observed eligible i)"
      ),
      "coverage_formula": (
        "coverage = sum(w_i for observed eligible i) / "
        "sum(w_i for eligible i)"
      ),
      "missing_data_rule": (
        "Missing eligible domains are omitted from the observed denominator and "
        "penalized through coverage; non-eligible domains enter neither denominator."
      ),
    })
    summaries.append(row)
  return pd.DataFrame(summaries)
