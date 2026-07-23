from __future__ import annotations

import pandas as pd


EVIDENCE_CATEGORY_DEFINITIONS: tuple[dict[str, str], ...] = (
  {
    "category": "candidate_universe",
    "criterion": "Included by a documented curated or systematic candidate-source rule.",
    "allowed_claim": "Candidate in the evaluated universe.",
  },
  {
    "category": "computational_hypothesis",
    "criterion": "Candidate received a source-bounded computational score.",
    "allowed_claim": "Computational hypothesis.",
  },
  {
    "category": "prioritized_hypothesis",
    "criterion": "Coverage-adjusted score met the declared prioritization rule.",
    "allowed_claim": "Prioritized hypothesis, not a discovery.",
  },
  {
    "category": "microniche_supported_hypothesis",
    "criterion": "At least one traceable functional-microniche domain was observed.",
    "allowed_claim": "Hypothesis with functional-microniche support.",
  },
  {
    "category": "conditional_dependency_supported_hypothesis",
    "criterion": "A loss-versus-intact DepMap contrast was executable and observed.",
    "allowed_claim": "Observed conditional-dependency contrast.",
  },
  {
    "category": "nominally_significant_result",
    "criterion": "P < 0.05 within the declared test family before multiple-testing correction.",
    "allowed_claim": "Nominally significant result.",
  },
  {
    "category": "fdr_supported_result",
    "criterion": "Benjamini-Hochberg q < 0.05 within the declared family and supportive effect direction.",
    "allowed_claim": "FDR-supported result within the declared family.",
  },
  {
    "category": "externally_validated_result",
    "criterion": "Independent data not used to construct the score reproduce the result.",
    "allowed_claim": "Externally validated only when the independent source is traceable.",
  },
  {
    "category": "experimentally_tractable_candidate",
    "criterion": "A compound, assay or experimental tool is traceably available.",
    "allowed_claim": "Experimentally tractable; not clinical efficacy.",
  },
  {
    "category": "clinical_evidence",
    "criterion": "Traceable clinical evidence exists for the candidate and context.",
    "allowed_claim": "Clinical evidence only when directly supported by the cited clinical source.",
  },
)


def evidence_category_definitions() -> pd.DataFrame:
  return pd.DataFrame(EVIDENCE_CATEGORY_DEFINITIONS)
