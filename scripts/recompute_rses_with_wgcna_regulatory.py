#!/usr/bin/env python3
"""Recompute RSES-Onco with WGCNA and promoter-aware regulatory sublayers.

The total functional-microniche weights are unchanged. Pairwise expression and WGCNA
share the existing expression-context weight; DoRothEA, TF-expression consistency and
promoter motif support share the existing regulatory-network weight. This prevents
expression-derived evidence from receiving multiple independent full-domain weights.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  FUNCTIONAL_MICRONICHE_WEIGHTS,
  coverage_aware_score,
  expanded_onco_score,
  functional_microniche_score,
)

ROOT = Path(__file__).resolve().parents[1]
SCORING_SEMANTICS_VERSION = "eligibility-aware-wgcna-regulatory-v2"
EXPRESSION_SUBWEIGHTS = {
  "pairwise_expression_context": 0.50,
  "wgcna_expression_network": 0.50,
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if np.isfinite(result) else None


def explicit_bool(value: object) -> bool:
  if isinstance(value, bool):
    return value
  text = str(value).strip().casefold()
  return text in {"1", "true", "yes", "eligible"}


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


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--ranking", required=True)
  parser.add_argument("--functional-evidence", required=True)
  parser.add_argument("--output", required=True)
  args = parser.parse_args()

  ranking = pd.read_csv(resolve_path(args.ranking), sep="\t", low_memory=False)
  evidence = pd.read_csv(
    resolve_path(args.functional_evidence),
    sep="\t",
    low_memory=False,
  )
  required_evidence = {
    "pair_id",
    "component_wgcna_expression_network",
    "component_regulatory_network_composite",
  }
  missing = sorted(required_evidence - set(evidence.columns))
  if missing:
    raise ValueError(
      "Functional evidence lacks WGCNA/regulatory fields: "
      + ", ".join(missing)
    )
  evidence_by_pair = {
    str(record["pair_id"]): record
    for record in evidence.to_dict("records")
  }

  rows: list[dict[str, object]] = []
  for record in ranking.to_dict("records"):
    pair_id = str(record.get("pair_id"))
    pair_evidence = evidence_by_pair.get(pair_id, {})
    pairwise_expression = numeric(
      record.get("microniche_expression_context")
    )
    wgcna_expression = numeric(
      pair_evidence.get("component_wgcna_expression_network")
    )
    expression_score = coverage_aware_score(
      {
        "pairwise_expression_context": pairwise_expression,
        "wgcna_expression_network": wgcna_expression,
      },
      EXPRESSION_SUBWEIGHTS,
    )
    regulatory_composite = numeric(
      pair_evidence.get("component_regulatory_network_composite")
    )

    microniche_components = {
      domain: numeric(record.get(f"microniche_{domain}"))
      for domain in FUNCTIONAL_MICRONICHE_WEIGHTS
    }
    microniche_components["expression_context"] = (
      expression_score.adjusted_score
      if np.isfinite(expression_score.adjusted_score)
      else None
    )
    microniche_components["regulatory_network"] = regulatory_composite
    eligible_microniche = {
      domain
      for domain in FUNCTIONAL_MICRONICHE_WEIGHTS
      if explicit_bool(record.get(f"eligible_microniche_{domain}", True))
    }
    microniche = functional_microniche_score(
      microniche_components,
      eligible_domains=eligible_microniche,
    )

    onco_components = {
      domain: numeric(record.get(f"component_{domain}"))
      for domain in EXPANDED_ONCO_WEIGHTS
    }
    onco_components["functional_microniche"] = (
      microniche.adjusted_score
      if np.isfinite(microniche.adjusted_score)
      else None
    )
    eligible_onco = {
      domain
      for domain in EXPANDED_ONCO_WEIGHTS
      if explicit_bool(record.get(f"eligible_component_{domain}", True))
    }
    onco = expanded_onco_score(
      onco_components,
      eligible_domains=eligible_onco,
    )

    updated = dict(record)
    for domain, value in microniche_components.items():
      updated[f"microniche_{domain}"] = value
    for domain, value in onco_components.items():
      updated[f"component_{domain}"] = value
    updated.update({
      "pairwise_expression_context": pairwise_expression,
      "wgcna_expression_network": wgcna_expression,
      "expression_context_raw": expression_score.observed_score,
      "expression_context_subcoverage": expression_score.coverage,
      "expression_context_adjusted": expression_score.adjusted_score,
      "expression_context_observed_subcomponents": expression_score.n_domains,
      "regulatory_tf_association_divergence": pair_evidence.get(
        "regulatory_tf_association_divergence"
      ),
      "regulatory_tf_expression_profile_divergence": pair_evidence.get(
        "regulatory_tf_expression_profile_divergence"
      ),
      "regulatory_promoter_motif_divergence": pair_evidence.get(
        "regulatory_promoter_motif_divergence"
      ),
      "regulatory_network_raw": pair_evidence.get("regulatory_network_raw"),
      "regulatory_network_subcoverage": pair_evidence.get(
        "regulatory_network_coverage"
      ),
      "promoter_evidence_type": pair_evidence.get("promoter_evidence_type"),
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
      "scoring_semantics_version": SCORING_SEMANTICS_VERSION,
      "score_version": "RSES-Onco-expanded-v0.10.8",
      "expression_context_formula": (
        "0.5*pairwise_expression_context + 0.5*WGCNA_context, "
        "coverage-adjusted within the existing expression-context domain"
      ),
      "regulatory_network_formula": (
        "0.40*DoRothEA_regulator_divergence + "
        "0.35*TF_expression_profile_divergence + "
        "0.25*JASPAR_promoter_motif_divergence, coverage-adjusted within "
        "the existing regulatory-network domain"
      ),
      "direct_promoter_binding_claim": False,
    })
    rows.append(updated)

  result = pd.DataFrame(rows).sort_values(
    [
      "cancer",
      "score_comparability_group",
      "coverage_adjusted_rses",
      "functional_microniche_adjusted",
    ],
    ascending=[True, True, False, False],
  )
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  result.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)
  print(
    f"Recomputed WGCNA/promoter-aware RSES-Onco: {output} "
    f"({len(result):,} rows)"
  )
  print(f"Scoring semantics: {SCORING_SEMANTICS_VERSION}")


if __name__ == "__main__":
  main()
