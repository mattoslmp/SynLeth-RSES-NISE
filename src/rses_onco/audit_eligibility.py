from __future__ import annotations

from datetime import datetime, timezone
from typing import Mapping

import numpy as np
import pandas as pd

from .audit import (
  DomainSpec,
  candidate_identity,
  domain_specs,
  metadata_for_source,
  numeric,
  present,
  score_family_summary,
  source_metadata,
)


TECHNICAL_STATUS_TOKENS = {
  "failed",
  "failure",
  "error",
  "timeout",
  "unavailable",
  "request_failed",
  "mapping_request_failed",
  "invalid_cache",
  "technical_failure",
}
INSUFFICIENT_REASON_TOKENS = {
  "insufficient",
  "fewer_than",
  "too_few",
  "no_models_in_loss_group",
  "no_models_in_intact_group",
  "no_compatible_cancer_models",
  "no_evaluable_models",
}


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


def explicit_absence_reason(
  record: Mapping[str, object],
  spec: DomainSpec,
) -> str:
  columns: tuple[str, ...]
  if spec.domain == "tumor_event":
    columns = ("tumor_event_absence_reason",)
  elif spec.domain in {"dependency", "selectivity"}:
    columns = ("dependency_absence_reason",)
  elif spec.domain == "expression_compensation":
    columns = ("expression_compensation_absence_reason",)
  elif spec.domain == "expression_context":
    columns = ("expression_context_absence_reason",)
  elif spec.domain == "genetic_phenotype":
    columns = ("phenotype_profile_absence_reason",)
  else:
    columns = (
      f"{spec.domain}_absence_reason",
      f"component_{spec.domain}_absence_reason",
      f"microniche_{spec.domain}_absence_reason",
    )
  for column in columns:
    value = record.get(column)
    if present(value):
      return str(value).strip()

  if spec.domain == "expression_context":
    n_models = numeric(record.get("expression_context_n_models"))
    if n_models is not None and n_models < 3:
      return "fewer_than_three_compatible_expression_models"
  if spec.domain == "genetic_phenotype":
    n_models = numeric(record.get("phenotype_profile_n_models"))
    if n_models is not None and n_models < 3:
      return "fewer_than_three_compatible_crispr_models"
  if spec.domain == "functional_microniche":
    observed = numeric(record.get("functional_microniche_n_domains"))
    if observed is not None and observed == 0:
      return "no_observed_eligible_functional_microniche_domains"
  return ""


def fallback_absence_reason(spec: DomainSpec) -> str:
  reasons = {
    "tumor_event": "no_evaluable_tcga_event_for_gene_and_cancer",
    "dependency": "conditional_dependency_analysis_not_executable_or_no_observed_contrast",
    "selectivity": "loss_selectivity_analysis_not_executable_or_no_observed_contrast",
    "expression_compensation": "expression_compensation_analysis_not_executable_or_no_observed_contrast",
    "functional_relation": "no_curated_functional_relation_value",
    "functional_microniche": "no_observed_eligible_functional_microniche_domains",
    "validation_tractability": "no_curated_validation_or_tractability_value",
    "expression_context": "no_evaluable_expression_context_profile",
    "localization": "no_valid_localization_annotation_for_one_or_both_proteins",
    "biochemical_structural": "no_traceable_biochemical_or_structural_evidence",
    "genetic_phenotype": "no_evaluable_crispr_phenotype_profile",
    "interaction_network": "no_string_mapping_or_no_interaction_evidence",
    "regulatory_network": "no_regulatory_edges_for_one_or_both_genes",
  }
  return reasons.get(spec.domain, "evidence_absent_reason_not_recorded_by_upstream_stage")


def source_failed(source_info: Mapping[str, object]) -> bool:
  status = str(source_info.get("source_status") or "").casefold()
  available = source_info.get("source_available")
  if available is False:
    return True
  return any(token in status for token in TECHNICAL_STATUS_TOKENS)


def evidence_state_and_reason(
  *,
  eligible: bool,
  eligibility_reason: str,
  normalized: float | None,
  explicit_reason: str,
  source_info: Mapping[str, object],
  spec: DomainSpec,
) -> tuple[str, str]:
  if not eligible:
    return "not_eligible", eligibility_reason
  if normalized is not None:
    if np.isclose(normalized, 0.0, atol=1e-12):
      return "negative_evidence", "observed_component_equals_zero"
    if np.isclose(normalized, 0.5, atol=1e-12):
      return "neutral_evidence", "observed_component_equals_neutral_midpoint"
    return "observed_evidence", "observed_component_value"

  reason = explicit_reason or fallback_absence_reason(spec)
  reason_cf = reason.casefold()
  if source_failed(source_info):
    return "technical_failure", (
      explicit_reason or "technical_failure_or_source_unavailable"
    )
  if any(token in reason_cf for token in INSUFFICIENT_REASON_TOKENS):
    return "insufficient_sample", reason
  if any(token in reason_cf for token in TECHNICAL_STATUS_TOKENS):
    return "technical_failure", reason
  return "missing", reason


def evidence_count_details(
  record: Mapping[str, object],
  spec: DomainSpec,
  normalized: float | None,
) -> tuple[float | None, float | None, str]:
  explicit_count = numeric(record.get(f"{spec.domain}_evidence_count"))
  explicit_independent = numeric(
    record.get(f"{spec.domain}_independent_evidence_count")
  )
  if explicit_count is not None or explicit_independent is not None:
    return (
      explicit_count,
      explicit_independent,
      "explicit_upstream_evidence_count",
    )

  if spec.domain in {"dependency", "selectivity"}:
    n_loss = numeric(
      record.get("dependency_n_loss_observed")
      or record.get("dependency_n_loss")
    )
    n_intact = numeric(
      record.get("dependency_n_intact_observed")
      or record.get("dependency_n_intact")
    )
    if n_loss is not None and n_intact is not None:
      return n_loss + n_intact, 1.0, "models_in_one_independent_loss_vs_intact_contrast"
  if spec.domain == "expression_compensation":
    n_loss = numeric(
      record.get("expression_n_loss_observed")
      or record.get("expression_n_loss")
    )
    n_intact = numeric(
      record.get("expression_n_intact_observed")
      or record.get("expression_n_intact")
    )
    if n_loss is not None and n_intact is not None:
      return n_loss + n_intact, 1.0, "models_in_one_independent_expression_contrast"
  if spec.domain == "tumor_event":
    evaluable = numeric(record.get("tcga_evaluable_n"))
    if evaluable is not None:
      return evaluable, 1.0, "evaluable_tcga_samples_in_one_cancer_event_summary"
  if spec.domain == "expression_context":
    n_models = numeric(record.get("expression_context_n_models"))
    if n_models is not None:
      return n_models, 1.0 if normalized is not None else None, "models_in_one_expression_profile_comparison"
  if spec.domain == "genetic_phenotype":
    n_models = numeric(record.get("phenotype_profile_n_models"))
    if n_models is not None:
      return n_models, 1.0 if normalized is not None else None, "models_in_one_crispr_profile_comparison"
  if spec.domain == "functional_microniche":
    count = numeric(record.get("functional_microniche_n_domains"))
    if count is not None:
      return count, count, "observed_independent_microniche_domains"
  if spec.domain == "validation_tractability":
    count = sum(
      present(record.get(column))
      for column in (
        "genetic_screen",
        "isogenic_validation",
        "in_vivo",
        "clinical_tractability",
      )
    )
    if count:
      return float(count), float(count), "observed_curated_validation_fields"
  if normalized is not None:
    return 1.0, 1.0, "one_derived_pair_level_component_raw_counts_in_support_table"
  return None, None, "no_observed_evidence_count"


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
      if not eligible:
        original = None
        normalized = None
      source_info = metadata_for_source(spec.source, metadata)
      explicit_reason = explicit_absence_reason(record, spec)
      state, reason = evidence_state_and_reason(
        eligible=eligible,
        eligibility_reason=eligibility_reason,
        normalized=normalized,
        explicit_reason=explicit_reason,
        source_info=source_info,
        spec=spec,
      )
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
      evidence_count, independent_count, count_basis = evidence_count_details(
        record,
        spec,
        normalized,
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
        "evidence_count": evidence_count,
        "independent_evidence_count": independent_count,
        "evidence_count_basis": count_basis,
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
