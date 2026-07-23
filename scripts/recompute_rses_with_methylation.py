#!/usr/bin/env python3
"""Integrate promoter methylation into the expression-compensation score domain.

The global RSES-Onco domain weights remain unchanged. Promoter methylation is an
internal subcomponent of the existing `expression_compensation` domain, avoiding
double counting with tumor-event and functional-microniche domains. The integration
is applied only to a TCGA-integrated ranking; a DepMap-only ranking remains free of
TCGA methylation evidence.
"""
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
  coverage_aware_score,
  expanded_onco_score,
)

METHYLATION_COMPENSATION_SUBWEIGHTS = {
  "expression_compensation": 0.70,
  "promoter_methylation_context": 0.30,
}
SCORING_SEMANTICS_VERSION = "eligibility-aware-v1"
EXPRESSION_REGULATORY_SEMANTICS_VERSION = (
  "eligibility-aware-wgcna-regulatory-methylation-v4"
)
SCORE_VERSION = "RSES-Onco-expanded-v0.11.1"


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
  return str(value).strip().casefold() in {
    "1",
    "true",
    "yes",
    "eligible",
  }


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


def has_tcga_event_evidence(ranking: pd.DataFrame) -> bool:
  if "component_tumor_event" not in ranking.columns:
    return False
  values = pd.to_numeric(
    ranking["component_tumor_event"],
    errors="coerce",
  )
  return bool(np.isfinite(values).any())


def write_atomic(frame: pd.DataFrame, output: Path) -> None:
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--ranking", required=True)
  parser.add_argument("--methylation-evidence", required=True)
  parser.add_argument("--output", required=True)
  args = parser.parse_args()

  ranking_path = resolve_path(args.ranking)
  output = resolve_path(args.output)
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)

  if not has_tcga_event_evidence(ranking):
    if output.resolve() != ranking_path.resolve():
      write_atomic(ranking, output)
    print(
      "Skipped GDC methylation integration because the ranking contains no "
      "observed TCGA tumor-event component; DepMap-only semantics were preserved."
    )
    return

  evidence = pd.read_csv(
    resolve_path(args.methylation_evidence),
    sep="\t",
    low_memory=False,
  )
  required = {
    "pair_id",
    "cancer",
    "evidence_status",
    "promoter_methylation_context_score",
  }
  missing = sorted(required - set(evidence.columns))
  if missing:
    raise ValueError(
      "Methylation evidence lacks required columns: " + ", ".join(missing)
    )
  evidence_by_context = {
    (str(row["pair_id"]), str(row["cancer"])): row
    for row in evidence.to_dict("records")
  }

  rows: list[dict[str, object]] = []
  for record in ranking.to_dict("records"):
    pair_id = str(record.get("pair_id"))
    cancer = str(record.get("cancer"))
    methylation = evidence_by_context.get((pair_id, cancer), {})
    expression_only = numeric(
      record.get("component_expression_compensation")
    )
    methylation_score = numeric(
      methylation.get("promoter_methylation_context_score")
    )
    eligible_expression = explicit_bool(
      record.get("eligible_component_expression_compensation", True)
    )

    if eligible_expression:
      composite = coverage_aware_score(
        {
          "expression_compensation": expression_only,
          "promoter_methylation_context": methylation_score,
        },
        METHYLATION_COMPENSATION_SUBWEIGHTS,
      )
      integrated_expression = (
        composite.adjusted_score
        if np.isfinite(composite.adjusted_score)
        else None
      )
    else:
      composite = coverage_aware_score(
        {},
        METHYLATION_COMPENSATION_SUBWEIGHTS,
        eligible_domains=(),
      )
      integrated_expression = None

    onco_components = {
      domain: numeric(record.get(f"component_{domain}"))
      for domain in EXPANDED_ONCO_WEIGHTS
    }
    onco_components["expression_compensation"] = integrated_expression
    eligible_onco = {
      domain
      for domain in EXPANDED_ONCO_WEIGHTS
      if explicit_bool(
        record.get(f"eligible_component_{domain}", True)
      )
    }
    onco = expanded_onco_score(
      onco_components,
      eligible_domains=eligible_onco,
    )

    updated = dict(record)
    for domain, value in onco_components.items():
      updated[f"component_{domain}"] = value
    updated.update({
      "expression_compensation_expression_only": expression_only,
      "promoter_methylation_context_score": methylation_score,
      "expression_methylation_raw": composite.observed_score,
      "expression_methylation_subcoverage": composite.coverage,
      "expression_methylation_adjusted": composite.adjusted_score,
      "expression_methylation_observed_subcomponents": composite.n_domains,
      "expression_methylation_eligible_subcomponents": (
        composite.eligible_domains
      ),
      "methylation_evidence_status": methylation.get(
        "evidence_status",
        "missing",
      ),
      "methylation_absence_reason": methylation.get(
        "absence_reason",
        "pair_or_cancer_not_present_in_methylation_evidence",
      ),
      "methylation_n_samples": methylation.get("n_samples"),
      "lost_gene_median_promoter_beta": methylation.get(
        "lost_median_beta"
      ),
      "target_gene_median_promoter_beta": methylation.get(
        "target_median_beta"
      ),
      "methylation_median_delta_beta": methylation.get(
        "median_delta_beta"
      ),
      "methylation_spearman_rho": methylation.get(
        "methylation_spearman_rho"
      ),
      "methylation_p_value": methylation.get(
        "paired_wilcoxon_p_value"
      ),
      "methylation_q_value_bh": methylation.get(
        "paired_wilcoxon_q_value_bh"
      ),
      "methylation_q_value_bh_within_cancer": methylation.get(
        "paired_wilcoxon_q_value_bh_within_cancer"
      ),
      "methylation_source": methylation.get("methylation_source"),
      "methylation_workflow": methylation.get("methylation_workflow"),
      "methylation_score_role": "expression_compensation_sublayer",
      "methylation_direct_silencing_claim": False,
      "expression_methylation_formula": (
        "0.70*expression_compensation + "
        "0.30*(lost_promoter_beta*(1-target_promoter_beta)), "
        "coverage-adjusted within the existing expression-compensation domain"
      ),
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
      "expression_regulatory_semantics_version": (
        EXPRESSION_REGULATORY_SEMANTICS_VERSION
      ),
      "score_version": SCORE_VERSION,
      "score_domain_weights": ";".join(
        f"{key}={value}"
        for key, value in EXPANDED_ONCO_WEIGHTS.items()
      ),
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
  write_atomic(result, output)
  print(
    f"Recomputed methylation-aware RSES-Onco: {output} "
    f"({len(result):,} rows)"
  )
  print(f"Score version: {SCORE_VERSION}")
  print(
    "Expression/regulatory semantics: "
    f"{EXPRESSION_REGULATORY_SEMANTICS_VERSION}"
  )


if __name__ == "__main__":
  main()
