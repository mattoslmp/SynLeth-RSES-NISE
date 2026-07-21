#!/usr/bin/env python3
"""Run expanded RSES-Onco across all curated human NISE directions and benchmarks.

The expanded workflow separates two levels:
1. Human functional-microniche RSES: expression context, localization,
   biochemical/structural implementation, CRISPR phenotype profile, STRING
   neighborhood and transcriptional-regulatory neighborhood.
2. Cancer therapeutic prioritization: tumor event, conditional dependency,
   selectivity, expression compensation, functional relation, the adjusted
   functional-microniche score and validation/tractability.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.depmap import (
  dependency_components,
  dependency_contrast,
  expression_component,
  expression_contrast,
  read_depmap_inputs,
)
from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  expression_profile_metrics,
  expanded_onco_score,
  functional_microniche_score,
  phenotype_profile_metrics,
)
from rses_onco.tcga import event_component, homozygous_deletion_frequency, read_gistic_matrix
from rses_onco.utils import bh_adjust, canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
LINEAGES = {"colon": "colon", "stomach": "stomach", "lung": "lung"}
SIMPLE_GENE = re.compile(r"^[A-Za-z0-9-]+$")


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
    for column in ("genetic_screen", "isogenic_validation", "in_vivo", "clinical_tractability")
  ]
  observed = [value for value in values if value is not None]
  return float(np.mean(observed)) if observed else None


def priority_label(observed: float, coverage: float, n_domains: int) -> str:
  if not np.isfinite(observed):
    return "insufficient evidence"
  if observed >= 0.72 and coverage >= 0.70 and n_domains >= 5:
    return "high priority"
  if observed >= 0.48:
    return "moderate priority"
  return "exploratory"


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
  parser.add_argument("--output", default="results/expanded_26Q1/expanded_rses_onco.tsv")
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

      simple_pair = bool(lost and target and SIMPLE_GENE.fullmatch(lost) and SIMPLE_GENE.fullmatch(target))
      if simple_pair:
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
          dependency_rows.append({"pair_id": pair_id, "cancer": cancer, **dependency.__dict__})
          empirical_components.update(dependency_components(dependency))

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
          compensation_rows.append({"pair_id": pair_id, "cancer": cancer, **compensation.__dict__})
          empirical_components["expression_compensation"] = expression_component(compensation)

        if cancer in tcga:
          frequency = homozygous_deletion_frequency(tcga[cancer], lost)
          if frequency is not None:
            empirical_components["tumor_event"] = event_component(frequency[2])
            tcga_values = {
              "tcga_homdel_n": frequency[0],
              "tcga_evaluable_n": frequency[1],
              "tcga_homdel_frequency": frequency[2],
            }

      expression_profile = expression_profile_metrics(
        expression, models, lost, target, lineage,
      ) if simple_pair else None
      phenotype_profile = phenotype_profile_metrics(
        effect,
        models,
        lost,
        target,
        lineage,
        dependency_threshold=args.dependency_threshold,
      ) if simple_pair else None

      if expression_profile is not None:
        expression_profile_rows.append({"pair_id": pair_id, **expression_profile.__dict__})
      if phenotype_profile is not None:
        phenotype_profile_rows.append({"pair_id": pair_id, **phenotype_profile.__dict__})

      microniche_components = {
        "expression_context": (
          expression_profile.divergence if expression_profile is not None else None
        ),
        "localization": numeric(pair_evidence.get("component_localization")),
        "biochemical_structural": numeric(
          pair_evidence.get("component_biochemical_structural")
        ),
        "genetic_phenotype": (
          phenotype_profile.divergence if phenotype_profile is not None else None
        ),
        "interaction_network": numeric(
          pair_evidence.get("component_interaction_network")
        ),
        "regulatory_network": numeric(
          pair_evidence.get("component_regulatory_network")
        ),
      }
      microniche = functional_microniche_score(microniche_components)

      onco_components = {
        **empirical_components,
        "functional_relation": numeric(record.get("relation_confidence")),
        "functional_microniche": (
          microniche.adjusted_score if np.isfinite(microniche.adjusted_score) else None
        ),
        "validation_tractability": validation_tractability(record),
      }
      onco = expanded_onco_score(onco_components)

      output_record = dict(record)
      output_record.update({
        "cancer": cancer,
        "analysis_lost_gene": lost,
        "analysis_target_gene": target,
        **{f"component_{key}": value for key, value in onco_components.items()},
        **{f"microniche_{key}": value for key, value in microniche_components.items()},
        "functional_microniche_rses": microniche.observed_score,
        "functional_microniche_coverage": microniche.coverage,
        "functional_microniche_adjusted": microniche.adjusted_score,
        "functional_microniche_n_domains": microniche.n_domains,
        "rses_onco": onco.observed_score,
        "evidence_coverage": onco.coverage,
        "coverage_adjusted_rses": onco.adjusted_score,
        "n_domains": onco.n_domains,
        "priority_class": priority_label(onco.observed_score, onco.coverage, onco.n_domains),
        "has_empirical_dependency": dependency is not None,
        "has_empirical_expression_compensation": compensation is not None,
        "has_empirical_expression_context": (
          expression_profile is not None and expression_profile.divergence is not None
        ),
        "has_empirical_phenotype_profile": (
          phenotype_profile is not None and phenotype_profile.divergence is not None
        ),
        "has_empirical_tcga": empirical_components["tumor_event"] is not None,
        "score_version": "RSES-Onco-expanded-v0.8",
        "score_domain_weights": ";".join(
          f"{key}={value}" for key, value in EXPANDED_ONCO_WEIGHTS.items()
        ),
        **tcga_values,
        "string_direct_score": pair_evidence.get("string_direct_score"),
        "string_neighbor_jaccard": pair_evidence.get("string_neighbor_jaccard"),
        "regulator_jaccard": pair_evidence.get("regulator_jaccard"),
      })
      scored_rows.append(output_record)

  output = resolve_path(args.output)
  assert output is not None
  output.parent.mkdir(parents=True, exist_ok=True)
  result = pd.DataFrame(scored_rows)
  result = result.sort_values(
    ["cancer", "coverage_adjusted_rses", "functional_microniche_adjusted"],
    ascending=[True, False, False],
  )
  result.to_csv(output, sep="\t", index=False)

  dependency_table = pd.DataFrame(dependency_rows)
  if not dependency_table.empty:
    dependency_table["q_value_bh"] = bh_adjust(dependency_table["p_value"])
  dependency_table.to_csv(output.with_name("expanded_dependency_contrasts.tsv"), sep="\t", index=False)

  compensation_table = pd.DataFrame(compensation_rows)
  if not compensation_table.empty:
    compensation_table["q_value_bh"] = bh_adjust(compensation_table["p_value"])
  compensation_table.to_csv(
    output.with_name("expanded_expression_compensation.tsv"), sep="\t", index=False
  )

  pd.DataFrame(expression_profile_rows).drop_duplicates().to_csv(
    output.with_name("expanded_expression_context_profiles.tsv"), sep="\t", index=False
  )
  pd.DataFrame(phenotype_profile_rows).drop_duplicates().to_csv(
    output.with_name("expanded_crispr_phenotype_profiles.tsv"), sep="\t", index=False
  )

  nise_rows = int((result.get("source_class", pd.Series(dtype=str)) == "NISE").sum())
  print(f"Wrote {len(result):,} cancer-specific scored rows to {output}")
  print(f"Cancer-specific NISE rows: {nise_rows:,}")
  print(f"Unique candidate directions: {result['pair_id'].nunique():,}")


if __name__ == "__main__":
  main()
