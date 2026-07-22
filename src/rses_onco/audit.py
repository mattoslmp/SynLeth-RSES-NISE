from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .expanded import EXPANDED_ONCO_WEIGHTS, FUNCTIONAL_MICRONICHE_WEIGHTS
from .utils import clamp01


MISSING_STRINGS = {"", "na", "nan", "none", "null", "<na>"}

DOMAIN_LABELS = {
  "tumor_event": "Tumor event frequency",
  "dependency": "Conditional dependency",
  "selectivity": "Loss-associated selectivity",
  "expression_compensation": "Transcriptional compensation",
  "functional_relation": "Functional relationship",
  "functional_microniche": "Functional microniche",
  "validation_tractability": "Validation and experimental tractability",
  "expression_context": "Expression-context divergence",
  "localization": "Subcellular localization relationship",
  "biochemical_structural": "Biochemical and structural evidence",
  "genetic_phenotype": "CRISPR phenotype divergence",
  "interaction_network": "Protein interaction network",
  "regulatory_network": "Transcriptional regulatory network",
}

DOMAIN_SOURCES = {
  "tumor_event": "TCGA/GDC ASCAT3",
  "dependency": "DepMap CRISPR + WGS copy number",
  "selectivity": "DepMap CRISPR + WGS copy number",
  "expression_compensation": "DepMap expression + WGS copy number",
  "functional_relation": "Curated candidate evidence",
  "functional_microniche": "Integrated functional-microniche domains",
  "validation_tractability": "Curated validation and tractability evidence",
  "expression_context": "DepMap expression",
  "localization": "Human Protein Atlas",
  "biochemical_structural": "UniProtKB/PDB/curated structural evidence",
  "genetic_phenotype": "DepMap CRISPR phenotype profiles",
  "interaction_network": "STRING",
  "regulatory_network": "DoRothEA/OmniPath",
}

DOMAIN_ROLES = {
  "tumor_event": "score",
  "dependency": "score_and_statistical_support",
  "selectivity": "score_and_statistical_support",
  "expression_compensation": "score",
  "functional_relation": "score",
  "functional_microniche": "score",
  "validation_tractability": "score_and_prioritization",
  "expression_context": "microniche_score",
  "localization": "microniche_score",
  "biochemical_structural": "microniche_score",
  "genetic_phenotype": "microniche_score",
  "interaction_network": "microniche_score",
  "regulatory_network": "microniche_score",
}

TECHNICAL_FAILURE_TOKENS = {
  "failed", "failure", "error", "timeout", "unavailable", "request_failed",
  "mapping_request_failed", "invalid_cache", "technical_failure",
}


@dataclass(frozen=True)
class DomainSpec:
  domain: str
  column: str
  family: str
  weight: float
  label: str
  source: str
  role: str


def domain_specs() -> list[DomainSpec]:
  specs: list[DomainSpec] = []
  for domain, weight in EXPANDED_ONCO_WEIGHTS.items():
    specs.append(DomainSpec(
      domain=domain,
      column=f"component_{domain}",
      family="RSES-Onco",
      weight=float(weight),
      label=DOMAIN_LABELS[domain],
      source=DOMAIN_SOURCES[domain],
      role=DOMAIN_ROLES[domain],
    ))
  for domain, weight in FUNCTIONAL_MICRONICHE_WEIGHTS.items():
    specs.append(DomainSpec(
      domain=domain,
      column=f"microniche_{domain}",
      family="Functional microniche",
      weight=float(weight),
      label=DOMAIN_LABELS[domain],
      source=DOMAIN_SOURCES[domain],
      role=DOMAIN_ROLES[domain],
    ))
  return specs


def present(value: object) -> bool:
  if value is None:
    return False
  try:
    if pd.isna(value):
      return False
  except (TypeError, ValueError):
    pass
  return str(value).strip().casefold() not in MISSING_STRINGS


def first_present(record: Mapping[str, object], columns: Sequence[str]) -> tuple[str, str]:
  for column in columns:
    value = record.get(column)
    if present(value):
      return str(value).strip(), column
  return "", ""


def candidate_identity(record: Mapping[str, object]) -> dict[str, str]:
  lost_gene, lost_gene_source = first_present(
    record, ("analysis_lost_gene", "lost_gene")
  )
  lost_feature, lost_feature_source = first_present(record, ("lost_feature",))
  target, target_source = first_present(
    record, ("analysis_target_gene", "target_gene", "target_feature")
  )
  if lost_gene:
    lost = lost_gene
    entity_type = "gene"
    arrow = "→"
    lost_source = lost_gene_source
  elif lost_feature:
    lost = lost_feature
    entity_type = "composite_feature"
    arrow = "⇒"
    lost_source = lost_feature_source
  else:
    lost = "Unresolved loss context"
    entity_type = "unresolved"
    arrow = "⇒"
    lost_source = ""
  target = target or "Unresolved target"
  return {
    "candidate_id": str(record.get("pair_id") or ""),
    "hypothesis_direction": f"{lost} {arrow} {target}",
    "origin_entity": lost,
    "target_entity": target,
    "hypothesis_type": entity_type,
    "origin_display_source_column": lost_source,
    "target_display_source_column": target_source,
  }


def _status_text(record: Mapping[str, object], domain: str) -> str:
  candidates = (
    f"{domain}_status",
    f"status_{domain}",
    f"component_{domain}_status",
    f"microniche_{domain}_status",
    "source_status",
    "status",
  )
  values = []
  for column in candidates:
    value = record.get(column)
    if present(value):
      values.append(str(value).strip())
  return ";".join(values)


def technical_failure(record: Mapping[str, object], domain: str) -> bool:
  status = _status_text(record, domain).casefold()
  return any(token in status for token in TECHNICAL_FAILURE_TOKENS)


def domain_eligibility(record: Mapping[str, object], spec: DomainSpec) -> tuple[bool, str]:
  identity = candidate_identity(record)
  gene_pair = identity["hypothesis_type"] == "gene"
  cancer = str(record.get("cancer") or "")

  if spec.domain in {
    "dependency", "selectivity", "expression_compensation",
    "expression_context", "genetic_phenotype", "tumor_event",
  } and not gene_pair:
    return False, "composite_event_not_executable_as_single_gene_analysis"

  explicit_columns = (
    f"eligible_{spec.domain}",
    f"{spec.domain}_eligible",
    f"component_{spec.domain}_eligible",
    f"microniche_{spec.domain}_eligible",
  )
  for column in explicit_columns:
    value = record.get(column)
    if present(value):
      try:
        eligible = bool(int(float(value)))
      except (TypeError, ValueError):
        eligible = str(value).strip().casefold() in {"true", "yes", "eligible"}
      return eligible, "explicit_pipeline_eligibility" if eligible else "explicitly_not_eligible"

  lineage_flag = record.get(cancer) if cancer else None
  if present(lineage_flag):
    try:
      if int(float(lineage_flag)) == 0:
        return False, "candidate_not_applicable_to_cancer_context"
    except (TypeError, ValueError):
      pass

  return True, "eligible_by_candidate_and_domain_definition"


def missing_reason(record: Mapping[str, object], spec: DomainSpec, eligible: bool) -> str:
  if not eligible:
    return domain_eligibility(record, spec)[1]
  if technical_failure(record, spec.domain):
    return "technical_failure_or_source_unavailable"

  if spec.domain in {"dependency", "selectivity"}:
    if record.get("has_empirical_dependency") is False:
      return "analysis_not_executable_or_insufficient_group_size"
  if spec.domain == "expression_compensation":
    if record.get("has_empirical_expression_compensation") is False:
      return "analysis_not_executable_or_insufficient_group_size"
  if spec.domain == "expression_context":
    if record.get("has_empirical_expression_context") is False:
      return "fewer_than_three_compatible_models_or_unmapped_gene"
  if spec.domain == "genetic_phenotype":
    if record.get("has_empirical_phenotype_profile") is False:
      return "fewer_than_three_compatible_models_or_unmapped_gene"
  if spec.domain == "tumor_event":
    if record.get("has_empirical_tcga") is False:
      return "no_evaluable_tcga_event_for_gene_and_cancer"
  if spec.domain == "regulatory_network":
    return "no_regulatory_edges_or_source_unavailable"
  if spec.domain == "interaction_network":
    return "no_string_mapping_or_no_interaction_evidence"
  if spec.domain == "localization":
    return "no_valid_localization_annotation_for_one_or_both_proteins"
  if spec.domain == "biochemical_structural":
    return "no_traceable_biochemical_or_structural_evidence"
  if spec.domain == "functional_relation":
    return "no_curated_functional_relation_value"
  if spec.domain == "validation_tractability":
    return "no_curated_validation_or_tractability_value"
  if spec.domain == "functional_microniche":
    return "no_observed_microniche_domains"
  return "evidence_absent_reason_not_recorded_by_upstream_stage"


def evidence_state(
  record: Mapping[str, object],
  spec: DomainSpec,
  eligible: bool,
  value: float | None,
) -> tuple[str, str]:
  if not eligible:
    return "not_eligible", domain_eligibility(record, spec)[1]
  if value is None:
    reason = missing_reason(record, spec, eligible)
    if reason.startswith("technical_failure") or "source_unavailable" in reason:
      return "technical_failure", reason
    if "insufficient" in reason or "fewer_than" in reason:
      return "insufficient_sample", reason
    return "missing", reason
  if value == 0.0:
    return "negative_evidence", "observed_component_equals_zero"
  if np.isclose(value, 0.5, atol=1e-12):
    return "neutral_evidence", "observed_component_equals_neutral_midpoint"
  return "observed_evidence", "observed_component_value"


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if np.isfinite(result) else None


def source_metadata(root: Path) -> dict[str, dict[str, object]]:
  """Collect source version/date/status without inventing missing metadata."""
  metadata: dict[str, dict[str, object]] = {}
  candidates = list(root.rglob("*status*.json")) + list(root.rglob("*metadata*.json"))
  for path in sorted(set(candidates)):
    try:
      payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
      continue
    key = path.stem
    metadata[key] = {
      "source_metadata_file": str(path),
      "source_version": payload.get("version") or payload.get("release") or "not_recorded",
      "access_date": payload.get("retrieved_at_utc") or payload.get("accessed_at") or payload.get("date") or "not_recorded",
      "source_status": payload.get("status") or "not_recorded",
      "source_available": payload.get("available"),
      "source_url": payload.get("source_url") or payload.get("source") or "not_recorded",
    }
  return metadata


def metadata_for_source(
  source: str,
  metadata: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
  tokens = {
    "STRING": ("string",),
    "DoRothEA/OmniPath": ("dorothea", "omnipath"),
    "Human Protein Atlas": ("hpa", "subcellular"),
    "UniProtKB/PDB/curated structural evidence": ("uniprot", "structure"),
    "TCGA/GDC ASCAT3": ("gdc", "tcga", "ascat"),
    "DepMap CRISPR + WGS copy number": ("depmap",),
    "DepMap expression + WGS copy number": ("depmap",),
    "DepMap expression": ("depmap",),
    "DepMap CRISPR phenotype profiles": ("depmap",),
  }.get(source, ())
  matches = [
    value for key, value in metadata.items()
    if any(token in key.casefold() for token in tokens)
  ]
  if not matches:
    return {
      "source_version": "not_recorded",
      "access_date": "not_recorded",
      "source_status": "not_recorded",
      "source_metadata_file": "",
      "source_url": "not_recorded",
    }
  merged = dict(matches[0])
  return merged


def score_family_summary(
  record: Mapping[str, object],
  family: str,
) -> tuple[float | None, float | None, float | None, int | None]:
  if family == "RSES-Onco":
    return (
      numeric(record.get("rses_onco")),
      numeric(record.get("evidence_coverage")),
      numeric(record.get("coverage_adjusted_rses")),
      int(record["n_domains"]) if present(record.get("n_domains")) else None,
    )
  return (
    numeric(record.get("functional_microniche_rses")),
    numeric(record.get("functional_microniche_coverage")),
    numeric(record.get("functional_microniche_adjusted")),
    int(record["functional_microniche_n_domains"])
    if present(record.get("functional_microniche_n_domains")) else None,
  )


def build_candidate_domain_audit(
  ranking: pd.DataFrame,
  metadata_root: Path,
) -> pd.DataFrame:
  metadata = source_metadata(metadata_root)
  rows: list[dict[str, object]] = []
  total_weights = {
    "RSES-Onco": float(sum(EXPANDED_ONCO_WEIGHTS.values())),
    "Functional microniche": float(sum(FUNCTIONAL_MICRONICHE_WEIGHTS.values())),
  }
  for record in ranking.to_dict("records"):
    identity = candidate_identity(record)
    for spec in domain_specs():
      eligible, eligibility_reason = domain_eligibility(record, spec)
      original = numeric(record.get(spec.column))
      normalized = clamp01(original)
      normalized = float(normalized) if normalized is not None else None
      state, reason = evidence_state(record, spec, eligible, normalized)
      observed, coverage, adjusted, n_domains = score_family_summary(record, spec.family)
      contribution = (
        spec.weight * normalized / total_weights[spec.family]
        if normalized is not None else None
      )
      source_info = metadata_for_source(spec.source, metadata)
      evidence_count = record.get(f"{spec.domain}_evidence_count")
      independent_count = record.get(f"{spec.domain}_independent_evidence_count")
      rows.append({
        **identity,
        "cancer": record.get("cancer"),
        "mechanistic_class": record.get("source_class") or record.get("relation_type"),
        "domain_family": spec.family,
        "domain": spec.domain,
        "domain_label": spec.label,
        "internal_component_column": spec.column,
        "eligible": bool(eligible),
        "eligibility_reason": eligibility_reason,
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
        "independent_evidence_count": independent_count if present(independent_count) else np.nan,
        "component_original": original,
        "component_normalized": normalized,
        "component_coverage_adjusted": contribution,
        "domain_weight": spec.weight,
        "weighted_numerator_contribution": spec.weight * normalized if normalized is not None else None,
        "final_score_contribution": contribution,
        "family_observed_score": observed,
        "family_coverage": coverage,
        "family_adjusted_score": adjusted,
        "family_observed_domains": n_domains,
        "audit_generated_at_utc": datetime.now(timezone.utc).isoformat(),
      })
  return pd.DataFrame(rows)


def coverage_summary(
  audit: pd.DataFrame,
  group_columns: Sequence[str],
) -> pd.DataFrame:
  eligible = audit.loc[audit["eligible"].fillna(False)].copy()
  grouped = (
    eligible.groupby(list(group_columns), dropna=False)
      .agg(
        eligible_hypotheses=("candidate_id", "size"),
        hypotheses_with_evidence=("evidence_present", "sum"),
        observed_evidence_count=("evidence_count", "sum"),
        independent_evidence_count=("independent_evidence_count", "sum"),
        median_component=("component_normalized", "median"),
        mean_component=("component_normalized", "mean"),
      )
      .reset_index()
  )
  grouped["coverage_fraction"] = np.where(
    grouped["eligible_hypotheses"] > 0,
    grouped["hypotheses_with_evidence"] / grouped["eligible_hypotheses"],
    np.nan,
  )
  grouped["coverage_label"] = (
    grouped["hypotheses_with_evidence"].astype(int).astype(str)
    + "/"
    + grouped["eligible_hypotheses"].astype(int).astype(str)
  )
  return grouped


def missingness_summary(audit: pd.DataFrame) -> pd.DataFrame:
  return (
    audit.groupby(
      ["domain_family", "domain", "domain_label", "evidence_state", "absence_reason"],
      dropna=False,
    )
      .agg(rows=("candidate_id", "size"), unique_candidates=("candidate_id", "nunique"))
      .reset_index()
      .sort_values(["domain_family", "domain", "rows"], ascending=[True, True, False])
  )


def recompute_family(
  group: pd.DataFrame,
  total_weight: float,
) -> dict[str, float | int]:
  observed = group.loc[group["component_normalized"].notna()].copy()
  numerator = float((observed["domain_weight"] * observed["component_normalized"]).sum())
  observed_weight = float(observed["domain_weight"].sum())
  raw = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / total_weight if total_weight else float("nan")
  adjusted = raw * coverage if np.isfinite(raw) else float("nan")
  return {
    "recomputed_weighted_numerator": numerator,
    "recomputed_observed_weight": observed_weight,
    "recomputed_observed_score": raw,
    "recomputed_coverage": coverage,
    "recomputed_adjusted_score": adjusted,
    "recomputed_n_domains": int(len(observed)),
  }


def score_decomposition(audit: pd.DataFrame) -> pd.DataFrame:
  totals = {
    "RSES-Onco": float(sum(EXPANDED_ONCO_WEIGHTS.values())),
    "Functional microniche": float(sum(FUNCTIONAL_MICRONICHE_WEIGHTS.values())),
  }
  keys = [
    "candidate_id", "hypothesis_direction", "origin_entity", "target_entity",
    "hypothesis_type", "cancer", "mechanistic_class", "domain_family",
  ]
  summaries = []
  for values, group in audit.groupby(keys, dropna=False):
    row = dict(zip(keys, values))
    row.update(recompute_family(group, totals[str(row["domain_family"]) ]))
    first = group.iloc[0]
    row.update({
      "pipeline_observed_score": first["family_observed_score"],
      "pipeline_coverage": first["family_coverage"],
      "pipeline_adjusted_score": first["family_adjusted_score"],
      "pipeline_n_domains": first["family_observed_domains"],
      "formula": "adjusted_score = sum(w_i*x_i_observed) / sum(w_i_observed) * sum(w_i_observed) / sum(w_i_all)",
      "missing_data_rule": "missing and non-eligible domains are not converted to zero; only observed eligible domains enter the observed-score denominator",
    })
    summaries.append(row)
  return pd.DataFrame(summaries)


def assign_evidence_category(record: Mapping[str, object]) -> dict[str, object]:
  q_candidates = (
    record.get("q_value_bh_within_loss_cancer"),
    record.get("q_value_bh"),
  )
  q_value = next((numeric(value) for value in q_candidates if numeric(value) is not None), None)
  p_value = numeric(record.get("p_value"))
  adjusted = numeric(record.get("coverage_adjusted_rses"))
  microniche = numeric(record.get("functional_microniche_adjusted"))
  pharmacology = numeric(record.get("pharmacology_adjusted"))
  clinical_phase = numeric(record.get("maximum_clinical_phase_normalized"))
  dependency = bool(record.get("has_empirical_dependency"))

  categories = ["candidate_universe"]
  categories.append("computational_hypothesis")
  if adjusted is not None and adjusted >= 0.48:
    categories.append("prioritized_hypothesis")
  if microniche is not None:
    categories.append("microniche_supported_hypothesis")
  if dependency:
    categories.append("conditional_dependency_supported_hypothesis")
  if p_value is not None and p_value < 0.05:
    categories.append("nominally_significant_result")
  if q_value is not None and q_value < 0.05:
    categories.append("fdr_supported_result")
  if pharmacology is not None:
    categories.append("experimentally_tractable_candidate")
  if clinical_phase is not None and clinical_phase > 0:
    categories.append("candidate_with_clinical_development_evidence")
  return {
    "evidence_categories": ";".join(categories),
    "highest_evidence_category": categories[-1],
    "terminology_boundary": (
      "Priority scores are computational hypotheses. FDR support refers only to the "
      "specified statistical family. Tractability is not clinical efficacy."
    ),
  }


def evidence_category_table(ranking: pd.DataFrame) -> pd.DataFrame:
  rows = []
  for record in ranking.to_dict("records"):
    rows.append({**candidate_identity(record), "cancer": record.get("cancer"), **assign_evidence_category(record)})
  return pd.DataFrame(rows)
