from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.audit import (
  build_candidate_domain_audit,
  candidate_identity,
  coverage_summary,
  score_decomposition,
)
from rses_onco.expanded import EXPANDED_ONCO_WEIGHTS, expanded_onco_score
from scripts.run_rses_robustness_analyses import recompute


def synthetic_ranking() -> pd.DataFrame:
  components = {
    "tumor_event": 0.50,
    "dependency": 0.80,
    "selectivity": 0.70,
    "expression_compensation": None,
    "functional_relation": 1.00,
    "functional_microniche": 0.60,
    "validation_tractability": 0.25,
  }
  scored = expanded_onco_score(components)
  return pd.DataFrame([{
    "pair_id": "PAIR1",
    "cancer": "lung",
    "source_class": "NISE",
    "lost_gene": "GENEA",
    "target_gene": "GENEB",
    **{f"component_{key}": value for key, value in components.items()},
    "microniche_expression_context": 0.4,
    "microniche_localization": None,
    "microniche_biochemical_structural": 0.8,
    "microniche_genetic_phenotype": 0.7,
    "microniche_interaction_network": 0.6,
    "microniche_regulatory_network": None,
    "functional_microniche_rses": 0.625,
    "functional_microniche_coverage": 0.70,
    "functional_microniche_adjusted": 0.4375,
    "functional_microniche_n_domains": 4,
    "rses_onco": scored.observed_score,
    "evidence_coverage": scored.coverage,
    "coverage_adjusted_rses": scored.adjusted_score,
    "n_domains": scored.n_domains,
    "has_empirical_dependency": True,
    "has_empirical_expression_compensation": False,
    "has_empirical_expression_context": True,
    "has_empirical_phenotype_profile": True,
    "has_empirical_tcga": True,
  }])


def test_missing_component_is_not_zero_or_contribution(tmp_path: Path) -> None:
  audit = build_candidate_domain_audit(synthetic_ranking(), tmp_path)
  row = audit.loc[
    (audit["domain_family"] == "RSES-Onco")
    & (audit["domain"] == "expression_compensation")
  ].iloc[0]
  assert row["evidence_state"] in {"missing", "insufficient_sample"}
  assert pd.isna(row["component_normalized"])
  assert pd.isna(row["final_score_contribution"])


def test_score_decomposition_reproduces_pipeline(tmp_path: Path) -> None:
  audit = build_candidate_domain_audit(synthetic_ranking(), tmp_path)
  decomposition = score_decomposition(audit)
  onco = decomposition.loc[decomposition["domain_family"] == "RSES-Onco"].iloc[0]
  assert np.isclose(onco["pipeline_observed_score"], onco["recomputed_observed_score"])
  assert np.isclose(onco["pipeline_coverage"], onco["recomputed_coverage"])
  assert np.isclose(onco["pipeline_adjusted_score"], onco["recomputed_adjusted_score"])


def test_coverage_summary_has_numerator_and_denominator(tmp_path: Path) -> None:
  audit = build_candidate_domain_audit(synthetic_ranking(), tmp_path)
  summary = coverage_summary(audit, ["domain_family", "domain", "domain_label"])
  dependency = summary.loc[summary["domain"] == "dependency"].iloc[0]
  assert dependency["eligible_hypotheses"] == 1
  assert dependency["hypotheses_with_evidence"] == 1
  assert dependency["coverage_label"] == "1/1"


def test_composite_feature_remains_feature_and_is_not_gene_eligible(tmp_path: Path) -> None:
  frame = synthetic_ranking().copy()
  frame["pair_id"] = "COMPOSITE"
  frame["lost_gene"] = np.nan
  frame["analysis_lost_gene"] = np.nan
  frame["lost_feature"] = "BRCA1/BRCA2 or HRD"
  frame["component_dependency"] = np.nan
  identity = candidate_identity(frame.iloc[0].to_dict())
  assert identity["hypothesis_type"] == "composite_feature"
  assert identity["hypothesis_direction"] == "BRCA1/BRCA2 or HRD ⇒ GENEB"
  audit = build_candidate_domain_audit(frame, tmp_path)
  dependency = audit.loc[audit["domain"] == "dependency"].iloc[0]
  assert dependency["evidence_state"] == "not_eligible"
  assert dependency["eligibility_reason"] == "composite_event_not_executable_as_single_gene_analysis"


def test_robustness_recompute_preserves_missing_values() -> None:
  frame = synthetic_ranking()
  recomputed = recompute(frame, dict(EXPANDED_ONCO_WEIGHTS))
  assert np.isclose(recomputed.loc[0, "recomputed_adjusted_score"], frame.loc[0, "coverage_adjusted_rses"])
  assert recomputed.loc[0, "recomputed_observed_domains"] == 6
