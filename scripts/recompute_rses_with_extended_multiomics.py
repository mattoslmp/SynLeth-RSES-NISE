#!/usr/bin/env python3
"""Recompute RSES-Onco v0.12.0 with causal extended multi-omics layers."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  FUNCTIONAL_MICRONICHE_WEIGHTS,
  coverage_aware_score,
)
from rses_onco.extended_multiomics import clamp01, coverage_consensus


ABLATIONS = {
  "integrated_functional_loss": "integrated_functional_loss_support",
  "dependency_probability": "dependency_probability_support",
  "protein_compensation": "protein_compensation_support",
  "rnai_orthogonal_support": "rnai_orthogonal_support",
}


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def numeric(record: pd.Series, column: str) -> float | None:
  if column not in record:
    return None
  return clamp01(record[column])


def combine(
  baseline: float | None,
  extension: float | None,
  *,
  baseline_weight: float,
  extension_eligible: bool,
) -> tuple[float | None, float]:
  if not extension_eligible:
    return baseline, 1.0 if baseline is not None else 0.0
  return coverage_consensus(
    {"baseline": baseline, "extension": extension},
    {"baseline": baseline_weight, "extension": 1.0 - baseline_weight},
  )


def score_row(
  row: pd.Series,
  *,
  exclude: str | None = None,
) -> dict[str, object]:
  extension_values = {
    key: None if exclude == key else numeric(row, column)
    for key, column in ABLATIONS.items()
  }
  tumor_event, tumor_coverage = combine(
    numeric(row, "baseline_component_tumor_event"),
    extension_values["integrated_functional_loss"],
    baseline_weight=0.65,
    extension_eligible=bool(
      row.get("eligible_integrated_functional_loss", False)
    ),
  )
  dependency, dependency_coverage = combine(
    numeric(row, "baseline_component_dependency"),
    extension_values["dependency_probability"],
    baseline_weight=0.75,
    extension_eligible=bool(
      row.get("eligible_dependency_probability", False)
    ),
  )
  expression_compensation, expression_coverage = combine(
    numeric(row, "baseline_component_expression_compensation"),
    extension_values["protein_compensation"],
    baseline_weight=0.65,
    extension_eligible=bool(
      row.get("eligible_protein_compensation", False)
    ),
  )
  genetic_phenotype, genetic_coverage = combine(
    numeric(row, "baseline_microniche_genetic_phenotype"),
    extension_values["rnai_orthogonal_support"],
    baseline_weight=0.70,
    extension_eligible=bool(
      row.get("eligible_rnai_orthogonal_support", False)
    ),
  )

  microniche_components = {
    "expression_context": numeric(
      row, "baseline_microniche_expression_context"
    ),
    "localization": numeric(row, "baseline_microniche_localization"),
    "biochemical_structural": numeric(
      row, "baseline_microniche_biochemical_structural"
    ),
    "genetic_phenotype": genetic_phenotype,
    "interaction_network": numeric(
      row, "baseline_microniche_interaction_network"
    ),
    "regulatory_network": numeric(
      row, "baseline_microniche_regulatory_network"
    ),
  }
  eligible_microniche = {
    key
    for key in FUNCTIONAL_MICRONICHE_WEIGHTS
    if bool(row.get(f"eligible_microniche_{key}", True))
  }
  microniche = coverage_aware_score(
    microniche_components,
    FUNCTIONAL_MICRONICHE_WEIGHTS,
    eligible_domains=eligible_microniche,
  )

  onco_components = {
    "tumor_event": tumor_event,
    "dependency": dependency,
    "selectivity": numeric(row, "baseline_component_selectivity"),
    "expression_compensation": expression_compensation,
    "functional_relation": numeric(
      row, "baseline_component_functional_relation"
    ),
    "functional_microniche": (
      microniche.adjusted_score
      if np.isfinite(microniche.adjusted_score)
      else None
    ),
    "validation_tractability": numeric(
      row, "baseline_component_validation_tractability"
    ),
  }
  eligible_onco = {
    key
    for key in EXPANDED_ONCO_WEIGHTS
    if bool(row.get(f"eligible_component_{key}", True))
  }
  onco = coverage_aware_score(
    onco_components,
    EXPANDED_ONCO_WEIGHTS,
    eligible_domains=eligible_onco,
  )
  return {
    "microniche_components": microniche_components,
    "onco_components": onco_components,
    "microniche": microniche,
    "onco": onco,
    "tumor_coverage": tumor_coverage,
    "dependency_coverage": dependency_coverage,
    "expression_coverage": expression_coverage,
    "genetic_coverage": genetic_coverage,
  }


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--evidence",
    default=(
      "data/processed/extended_multiomics/"
      "extended_pair_evidence_by_cancer.tsv"
    ),
  )
  parser.add_argument("--output", default=None)
  args = parser.parse_args()

  ranking_path = resolve(args.ranking)
  evidence_path = resolve(args.evidence)
  output_path = resolve(args.output) if args.output else ranking_path

  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  evidence = pd.read_csv(evidence_path, sep="\t", low_memory=False)
  keys = ["pair_id", "cancer"]
  if evidence.duplicated(keys).any():
    raise ValueError("Extended evidence contains duplicate pair_id/cancer rows")
  evidence_columns = [
    column
    for column in evidence.columns
    if column not in ranking.columns or column in keys
  ]
  merged = ranking.merge(
    evidence[evidence_columns],
    on=keys,
    how="left",
    validate="many_to_one",
  )

  baseline_columns = [
    "rses_onco",
    "evidence_coverage",
    "coverage_adjusted_rses",
    "functional_microniche_rses",
    "functional_microniche_coverage",
    "functional_microniche_adjusted",
    *[f"component_{key}" for key in EXPANDED_ONCO_WEIGHTS],
    *[f"microniche_{key}" for key in FUNCTIONAL_MICRONICHE_WEIGHTS],
  ]
  for column in baseline_columns:
    if column in merged.columns and f"baseline_{column}" not in merged.columns:
      merged[f"baseline_{column}"] = merged[column]

  output_rows = []
  for _, row in merged.iterrows():
    full = score_row(row)
    microniche = full["microniche"]
    onco = full["onco"]
    updated = row.to_dict()
    for key, value in full["microniche_components"].items():
      updated[f"microniche_{key}"] = value
    for key, value in full["onco_components"].items():
      updated[f"component_{key}"] = value
    updated.update({
      "functional_microniche_rses": microniche.observed_score,
      "functional_microniche_coverage": microniche.coverage,
      "functional_microniche_adjusted": microniche.adjusted_score,
      "functional_microniche_n_domains": microniche.n_domains,
      "functional_microniche_observed_weight": microniche.observed_weight,
      "functional_microniche_eligible_weight": microniche.eligible_weight,
      "rses_onco": onco.observed_score,
      "evidence_coverage": onco.coverage,
      "coverage_adjusted_rses": onco.adjusted_score,
      "n_domains": onco.n_domains,
      "observed_domain_weight": onco.observed_weight,
      "eligible_domain_weight": onco.eligible_weight,
      "extended_tumor_event_internal_coverage": full["tumor_coverage"],
      "extended_dependency_internal_coverage": full["dependency_coverage"],
      "extended_expression_internal_coverage": full["expression_coverage"],
      "extended_genetic_phenotype_internal_coverage": full["genetic_coverage"],
      "score_version": "RSES-Onco-expanded-v0.12.0",
      "extended_multiomics_semantics_version": (
        "causal-layers-scored-context-layers-separated-v1"
      ),
      "extended_internal_weights": (
        "tumor_event:baseline=0.65,functional_loss=0.35;"
        "dependency:Chronos=0.75,gene_dependency=0.25;"
        "expression_compensation:RNA=0.65,protein=0.35;"
        "genetic_phenotype:CRISPR=0.70,RNAi=0.30"
      ),
      "extended_primary_score_includes": (
        "integrated_functional_loss;crispr_dependency_probability;"
        "multi_platform_protein_compensation;rnai_orthogonal_support"
      ),
      "extended_context_not_directly_scored": (
        "metabolomics_without_reaction_mapping;miRNA;global_chromatin;"
        "ssGSEA;subtypes;omics_signatures;MetMap;unresolved_fusions;"
        "hotspot_mutations_without_LOF_annotation"
      ),
    })
    for ablation in ABLATIONS:
      ablated = score_row(row, exclude=ablation)
      ablated_onco = ablated["onco"]
      updated[
        f"ablation_without_{ablation}_coverage_adjusted_rses"
      ] = ablated_onco.adjusted_score
    output_rows.append(updated)

  result = pd.DataFrame(output_rows)
  result["extended_rank_within_cancer"] = result.groupby("cancer")[
    "coverage_adjusted_rses"
  ].rank(method="min", ascending=False, na_option="bottom")
  if "baseline_coverage_adjusted_rses" in result.columns:
    result["extended_score_delta"] = (
      pd.to_numeric(result["coverage_adjusted_rses"], errors="coerce")
      - pd.to_numeric(
        result["baseline_coverage_adjusted_rses"], errors="coerce"
      )
    )
    result["baseline_rank_within_cancer"] = result.groupby("cancer")[
      "baseline_coverage_adjusted_rses"
    ].rank(method="min", ascending=False, na_option="bottom")
    result["extended_rank_change"] = (
      result["baseline_rank_within_cancer"]
      - result["extended_rank_within_cancer"]
    )
  for ablation in ABLATIONS:
    score_column = f"ablation_without_{ablation}_coverage_adjusted_rses"
    rank_column = f"ablation_without_{ablation}_rank_within_cancer"
    result[rank_column] = result.groupby("cancer")[score_column].rank(
      method="min",
      ascending=False,
      na_option="bottom",
    )

  output_path.parent.mkdir(parents=True, exist_ok=True)
  temporary = output_path.with_suffix(output_path.suffix + ".tmp")
  result.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output_path)
  print(f"Extended ranking rows: {len(result):,}")
  print(f"Extended layer ablations: {len(ABLATIONS)}")
  print(f"Wrote {output_path}")


if __name__ == "__main__":
  main()
