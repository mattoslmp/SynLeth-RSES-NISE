#!/usr/bin/env python3
"""Run eligibility-aware expanded RSES-Onco across all candidate directions.

The workflow separates human functional-microniche evidence from cancer-specific
therapeutic prioritization. Missing eligible evidence remains missing and lowers
explicit coverage. Domains that cannot be evaluated for a hypothesis type are marked
non-eligible and excluded from the eligible denominator. Technical or source failure
is never converted into biological negative evidence.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.depmap import (
  cancer_model_ids,
  dependency_components,
  dependency_contrast,
  expression_component,
  expression_contrast,
  read_depmap_inputs,
)
from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  FUNCTIONAL_MICRONICHE_WEIGHTS,
  expression_profile_metrics,
  expanded_onco_score,
  functional_microniche_score,
  phenotype_profile_metrics,
)
from rses_onco.tcga import (
  event_component,
  homozygous_deletion_frequency,
  read_gistic_matrix,
)
from rses_onco.utils import bh_adjust, canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
LINEAGES = {"colon": "colon", "stomach": "stomach", "lung": "lung"}
SIMPLE_GENE = re.compile(r"^[A-Za-z0-9-]+$")
SCORING_SEMANTICS_VERSION = "eligibility-aware-v1"


def resolve_path(value: str | None) -> Path | None:
  if value is None:
    return None
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if np.isfinite(result) else None


def validation_tractability(record: dict[str, object]) -> float | None:
  values = [
    numeric(record.get(column))
    for column in (
      "genetic_screen",
      "isogenic_validation",
      "in_vivo",
      "clinical_tractability",
    )
  ]
  observed = [value for value in values if value is not None]
  return float(np.mean(observed)) if observed else None


def priority_label(
  observed: float,
  coverage: float,
  n_domains: int,
  eligible_domains: int,
) -> str:
  if not np.isfinite(observed) or eligible_domains == 0:
    return "insufficient evidence"
  if observed >= 0.72 and coverage >= 0.70 and n_domains >= 5:
    return "high priority"
  if observed >= 0.48 and coverage >= 0.50 and n_domains >= 3:
    return "moderate priority"
  return "exploratory"


def group_diagnostics(
  matrix: pd.DataFrame,
  copy_number: pd.DataFrame,
  models: pd.DataFrame,
  lost_gene: str,
  target_gene: str,
  cancer: str,
  loss_threshold: float,
) -> dict[str, object]:
  """Return mapping, sample-size and group-balance diagnostics without scoring."""
  result: dict[str, object] = {
    "lost_gene_in_copy_number": lost_gene in copy_number,
    "target_gene_in_matrix": target_gene in matrix,
    "cancer_model_count": 0,
    "evaluable_model_count": 0,
    "n_loss": 0,
    "n_intact": 0,
    "absence_reason": "",
  }
  selected = set(cancer_model_ids(models, cancer).astype(str))
  result["cancer_model_count"] = len(selected)
  if not selected:
    result["absence_reason"] = "no_compatible_cancer_models"
    return result
  if lost_gene not in copy_number:
    result["absence_reason"] = "lost_gene_unmapped_in_copy_number"
    return result
  if target_gene not in matrix:
    result["absence_reason"] = "target_gene_unmapped_in_analysis_matrix"
    return result
  table = copy_number[["ModelID", lost_gene]].merge(
    matrix[["ModelID", target_gene]],
    on="ModelID",
    how="inner",
  )
  table = table.loc[table["ModelID"].astype(str).isin(selected)].copy()
  table[lost_gene] = pd.to_numeric(table[lost_gene], errors="coerce")
  table[target_gene] = pd.to_numeric(table[target_gene], errors="coerce")
  table = table.dropna(subset=[lost_gene, target_gene])
  result["evaluable_model_count"] = len(table)
  result["n_loss"] = int((table[lost_gene] < loss_threshold).sum())
  result["n_intact"] = int((table[lost_gene] >= loss_threshold).sum())
  if table.empty:
    result["absence_reason"] = "no_evaluable_models_after_mapping_and_join"
  elif result["n_loss"] == 0:
    result["absence_reason"] = "no_models_in_loss_group"
  elif result["n_intact"] == 0:
    result["absence_reason"] = "no_models_in_intact_group"
  return result


def eligibility_for_candidate(simple_pair: bool) -> tuple[set[str], set[str]]:
  if simple_pair:
    return set(FUNCTIONAL_MICRONICHE_WEIGHTS), set(EXPANDED_ONCO_WEIGHTS)
  return set(), {"functional_relation", "validation_tractability"}


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument("--expression", required=True)
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--functional-evidence",
    default="data/processed/expanded_pair_functional_evidence.tsv",
  )
  parser.add_argument("--tcga", action="append", default=[])
  parser.add_argument("--loss-threshold", type=float, default=0.30)
  parser.add_argument("--min-group-size", type=int, default=3)
  parser.add_argument("--dependency-threshold", type=float, default=-0.5)
  parser.add_argument(
    "--output",
    default="results/expanded_26Q1/expanded_rses_onco.tsv",
  )
  args = parser.parse_args()

  effect, copy_number, models, expression = read_depmap_inputs(
    resolve_path(args.gene_effect),
    resolve_path(args.copy_number),
    resolve_path(args.models),
    resolve_path(args.expression),
  )
  if expression is None:
    raise ValueError("Expanded RSES-Onco requires an expression matrix")

  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t")
  evidence_path = resolve_path(args.functional_evidence)
  if evidence_path is not None and evidence_path.exists():
    functional_evidence = pd.read_csv(evidence_path, sep="\t")
  else:
    functional_evidence = pd.DataFrame()
  evidence_by_id = {
    str(record["pair_id"]): record
    for record in functional_evidence.to_dict("records")
    if pd.notna(record.get("pair_id"))
  }

  tcga: dict[str, pd.DataFrame] = {}
  for item in args.tcga:
    cancer, path = item.split("=", 1)
    tcga[cancer] = read_gistic_matrix(resolve_path(path))

  scored_rows: list[dict[str, object]] = []
  dependency_rows: list[dict[str, object]] = []
  compensation_rows: list[dict[str, object]] = []
  expression_profile_rows: list[dict[str, object]] = []
  phenotype_profile_rows: list[dict[str, object]] = []

  for record in candidates.to_dict("records"):
    pair_id = str(record["pair_id"])
    lost = canonical_gene_name(record.get("lost_gene"))
    target = canonical_gene_name(record.get("target_gene"))
    pair_evidence = evidence_by_id.get(pair_id, {})

    for cancer, lineage in LINEAGES.items():
      flag = record.get(cancer, 1)
      try:
        included = int(float(flag)) == 1
      except (TypeError, ValueError):
        included = True
      if not included:
        continue

      dependency = None
      compensation = None
      empirical_components: dict[str, float | None] = {
        "tumor_event": None,
        "dependency": None,
        "selectivity": None,
        "expression_compensation": None,
      }
      tcga_values: dict[str, object] = {}
      simple_pair = bool(
        lost
        and target
        and SIMPLE_GENE.fullmatch(lost)
        and SIMPLE_GENE.fullmatch(target)
      )
      microniche_eligible, onco_eligible = eligibility_for_candidate(simple_pair)

      dependency_diagnostics = {
        "cancer_model_count": 0,
        "evaluable_model_count": 0,
        "n_loss": 0,
        "n_intact": 0,
        "absence_reason": "composite_event_not_executable_as_single_gene_analysis",
      }
      expression_diagnostics = dict(dependency_diagnostics)

      if simple_pair:
        dependency_diagnostics = group_diagnostics(
          effect,
          copy_number,
          models,
          lost,
          target,
          lineage,
          args.loss_threshold,
        )
        dependency = dependency_contrast(
          effect,
          copy_number,
          models,
          lost,
          target,
          lineage,
          loss_threshold=args.loss_threshold,
          min_group_size=args.min_group_size,
        )
        if dependency is not None:
          dependency_rows.append({
            "pair_id": pair_id,
            "cancer": cancer,
            **dependency.__dict__,
          })
          empirical_components.update(dependency_components(dependency))
        elif not dependency_diagnostics["absence_reason"]:
          dependency_diagnostics["absence_reason"] = (
            "insufficient_loss_or_intact_group_size"
          )

        expression_diagnostics = group_diagnostics(
          expression,
          copy_number,
          models,
          lost,
          target,
          lineage,
          args.loss_threshold,
        )
        compensation = expression_contrast(
          expression,
          copy_number,
          models,
          lost,
          target,
          lineage,
          loss_threshold=args.loss_threshold,
          min_group_size=args.min_group_size,
        )
        if compensation is not None:
          compensation_rows.append({
            "pair_id": pair_id,
            "cancer": cancer,
            **compensation.__dict__,
          })
          empirical_components["expression_compensation"] = (
            expression_component(compensation)
          )
        elif not expression_diagnostics["absence_reason"]:
          expression_diagnostics["absence_reason"] = (
            "insufficient_loss_or_intact_group_size"
          )

        if cancer in tcga:
          frequency = homozygous_deletion_frequency(tcga[cancer], lost)
          if frequency is not None:
            empirical_components["tumor_event"] = event_component(frequency[2])
            tcga_values = {
              "tcga_homdel_n": frequency[0],
              "tcga_evaluable_n": frequency[1],
              "tcga_homdel_frequency": frequency[2],
              "tumor_event_absence_reason": "",
            }
          else:
            tcga_values["tumor_event_absence_reason"] = (
              "lost_gene_unavailable_or_no_evaluable_tcga_samples"
            )
        else:
          tcga_values["tumor_event_absence_reason"] = (
            "tcga_source_not_provided_for_cancer"
          )

      expression_profile = (
        expression_profile_metrics(expression, models, lost, target, lineage)
        if simple_pair
        else None
      )
      phenotype_profile = (
        phenotype_profile_metrics(
          effect,
          models,
          lost,
          target,
          lineage,
          dependency_threshold=args.dependency_threshold,
        )
        if simple_pair
        else None
      )

      if expression_profile is not None:
        expression_profile_rows.append({
          "pair_id": pair_id,
          **expression_profile.__dict__,
        })
      if phenotype_profile is not None:
        phenotype_profile_rows.append({
          "pair_id": pair_id,
          **phenotype_profile.__dict__,
        })

      microniche_components = {
        "expression_context": (
          expression_profile.divergence
          if expression_profile is not None
          else None
        ),
        "localization": numeric(pair_evidence.get("component_localization")),
        "biochemical_structural": numeric(
          pair_evidence.get("component_biochemical_structural")
        ),
        "genetic_phenotype": (
          phenotype_profile.divergence
          if phenotype_profile is not None
          else None
        ),
        "interaction_network": numeric(
          pair_evidence.get("component_interaction_network")
        ),
        "regulatory_network": numeric(
          pair_evidence.get("component_regulatory_network")
        ),
      }
      microniche = functional_microniche_score(
        microniche_components,
        eligible_domains=microniche_eligible,
      )

      onco_components = {
        **empirical_components,
        "functional_relation": numeric(record.get("relation_confidence")),
        "functional_microniche": (
          microniche.adjusted_score
          if np.isfinite(microniche.adjusted_score)
          else None
        ),
        "validation_tractability": validation_tractability(record),
      }
      onco = expanded_onco_score(
        onco_components,
        eligible_domains=onco_eligible,
      )

      output_record = dict(record)
      output_record.update({
        "cancer": cancer,
        "analysis_lost_gene": lost,
        "analysis_target_gene": target,
        "hypothesis_type": "gene" if simple_pair else "composite_feature",
        "score_comparability_group": (
          "gene_pair" if simple_pair else "composite_event_to_gene"
        ),
        "scoring_semantics_version": SCORING_SEMANTICS_VERSION,
        **{
          f"component_{key}": value
          for key, value in onco_components.items()
        },
        **{
          f"eligible_component_{key}": key in onco_eligible
          for key in EXPANDED_ONCO_WEIGHTS
        },
        **{
          f"microniche_{key}": value
          for key, value in microniche_components.items()
        },
        **{
          f"eligible_microniche_{key}": key in microniche_eligible
          for key in FUNCTIONAL_MICRONICHE_WEIGHTS
        },
        "functional_microniche_rses": microniche.observed_score,
        "functional_microniche_coverage": microniche.coverage,
        "functional_microniche_adjusted": microniche.adjusted_score,
        "functional_microniche_n_domains": microniche.n_domains,
        "functional_microniche_eligible_domains": microniche.eligible_domains,
        "functional_microniche_observed_weight": microniche.observed_weight,
        "functional_microniche_eligible_weight": microniche.eligible_weight,
        "rses_onco": onco.observed_score,
        "evidence_coverage": onco.coverage,
        "coverage_adjusted_rses": onco.adjusted_score,
        "n_domains": onco.n_domains,
        "eligible_domains": onco.eligible_domains,
        "observed_domain_weight": onco.observed_weight,
        "eligible_domain_weight": onco.eligible_weight,
        "priority_class": priority_label(
          onco.observed_score,
          onco.coverage,
          onco.n_domains,
          onco.eligible_domains,
        ),
        "has_empirical_dependency": dependency is not None,
        "has_empirical_expression_compensation": compensation is not None,
        "has_empirical_expression_context": (
          expression_profile is not None
          and expression_profile.divergence is not None
        ),
        "has_empirical_phenotype_profile": (
          phenotype_profile is not None
          and phenotype_profile.divergence is not None
        ),
        "has_empirical_tcga": empirical_components["tumor_event"] is not None,
        "dependency_cancer_model_count": dependency_diagnostics["cancer_model_count"],
        "dependency_evaluable_model_count": dependency_diagnostics["evaluable_model_count"],
        "dependency_n_loss": dependency_diagnostics["n_loss"],
        "dependency_n_intact": dependency_diagnostics["n_intact"],
        "dependency_absence_reason": (
          "" if dependency is not None else dependency_diagnostics["absence_reason"]
        ),
        "expression_cancer_model_count": expression_diagnostics["cancer_model_count"],
        "expression_evaluable_model_count": expression_diagnostics["evaluable_model_count"],
        "expression_n_loss": expression_diagnostics["n_loss"],
        "expression_n_intact": expression_diagnostics["n_intact"],
        "expression_compensation_absence_reason": (
          "" if compensation is not None else expression_diagnostics["absence_reason"]
        ),
        "expression_context_n_models": (
          expression_profile.n_models if expression_profile is not None else 0
        ),
        "phenotype_profile_n_models": (
          phenotype_profile.n_models if phenotype_profile is not None else 0
        ),
        "score_version": "RSES-Onco-expanded-v0.10.7",
        "score_domain_weights": ";".join(
          f"{key}={value}"
          for key, value in EXPANDED_ONCO_WEIGHTS.items()
        ),
        "eligible_score_domains": ";".join(sorted(onco_eligible)),
        "eligible_microniche_domains": ";".join(
          sorted(microniche_eligible)
        ),
        **tcga_values,
        "string_direct_score": pair_evidence.get("string_direct_score"),
        "string_neighbor_jaccard": pair_evidence.get(
          "string_neighbor_jaccard"
        ),
        "regulator_jaccard": pair_evidence.get("regulator_jaccard"),
      })
      if dependency is not None:
        output_record.update({
          "dependency_n_loss_observed": dependency.n_loss,
          "dependency_n_intact_observed": dependency.n_intact,
          "dependency_delta_effect": dependency.delta_effect,
          "dependency_p_value": dependency.p_value,
        })
      if compensation is not None:
        output_record.update({
          "expression_n_loss_observed": compensation.n_loss,
          "expression_n_intact_observed": compensation.n_intact,
          "expression_delta": compensation.delta_expression,
          "expression_p_value": compensation.p_value,
        })
      scored_rows.append(output_record)

  dependency_table = pd.DataFrame(dependency_rows)
  if not dependency_table.empty:
    dependency_table["q_value_bh"] = bh_adjust(dependency_table["p_value"])
    dependency_table["q_value_bh_within_cancer"] = (
      dependency_table.groupby("cancer", group_keys=False)["p_value"]
        .transform(lambda values: bh_adjust(values))
    )

  compensation_table = pd.DataFrame(compensation_rows)
  if not compensation_table.empty:
    compensation_table["q_value_bh"] = bh_adjust(compensation_table["p_value"])
    compensation_table["q_value_bh_within_cancer"] = (
      compensation_table.groupby("cancer", group_keys=False)["p_value"]
        .transform(lambda values: bh_adjust(values))
    )

  result = pd.DataFrame(scored_rows)
  if not dependency_table.empty:
    dependency_merge = dependency_table[[
      "pair_id",
      "cancer",
      "delta_effect",
      "p_value",
      "q_value_bh",
      "q_value_bh_within_cancer",
      "n_loss",
      "n_intact",
    ]].rename(columns={
      "delta_effect": "dependency_delta_effect",
      "p_value": "dependency_p_value",
      "q_value_bh": "dependency_q_value_bh",
      "q_value_bh_within_cancer": "dependency_q_value_bh_within_cancer",
      "n_loss": "dependency_n_loss_observed",
      "n_intact": "dependency_n_intact_observed",
    })
    duplicate_columns = [
      column
      for column in dependency_merge.columns
      if column in result.columns and column not in {"pair_id", "cancer"}
    ]
    result = result.drop(columns=duplicate_columns).merge(
      dependency_merge,
      on=["pair_id", "cancer"],
      how="left",
    )
  if not compensation_table.empty:
    compensation_merge = compensation_table[[
      "pair_id",
      "cancer",
      "delta_expression",
      "p_value",
      "q_value_bh",
      "q_value_bh_within_cancer",
      "n_loss",
      "n_intact",
    ]].rename(columns={
      "delta_expression": "expression_delta",
      "p_value": "expression_p_value",
      "q_value_bh": "expression_q_value_bh",
      "q_value_bh_within_cancer": "expression_q_value_bh_within_cancer",
      "n_loss": "expression_n_loss_observed",
      "n_intact": "expression_n_intact_observed",
    })
    duplicate_columns = [
      column
      for column in compensation_merge.columns
      if column in result.columns and column not in {"pair_id", "cancer"}
    ]
    result = result.drop(columns=duplicate_columns).merge(
      compensation_merge,
      on=["pair_id", "cancer"],
      how="left",
    )

  result = result.sort_values(
    [
      "cancer",
      "score_comparability_group",
      "coverage_adjusted_rses",
      "functional_microniche_adjusted",
    ],
    ascending=[True, True, False, False],
  )

  output = resolve_path(args.output)
  assert output is not None
  output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)
  dependency_table.to_csv(
    output.with_name("expanded_dependency_contrasts.tsv"),
    sep="\t",
    index=False,
  )
  compensation_table.to_csv(
    output.with_name("expanded_expression_compensation.tsv"),
    sep="\t",
    index=False,
  )
  pd.DataFrame(expression_profile_rows).drop_duplicates().to_csv(
    output.with_name("expanded_expression_context_profiles.tsv"),
    sep="\t",
    index=False,
  )
  pd.DataFrame(phenotype_profile_rows).drop_duplicates().to_csv(
    output.with_name("expanded_crispr_phenotype_profiles.tsv"),
    sep="\t",
    index=False,
  )

  nise_rows = int(
    (result.get("source_class", pd.Series(dtype=str)) == "NISE").sum()
  )
  print(f"Wrote {len(result):,} cancer-specific scored rows to {output}")
  print(f"Cancer-specific NISE rows: {nise_rows:,}")
  print(f"Unique candidate directions: {result['pair_id'].nunique():,}")
  print(f"Scoring semantics: {SCORING_SEMANTICS_VERSION}")


if __name__ == "__main__":
  main()
