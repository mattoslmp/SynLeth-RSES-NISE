from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.audit_eligibility import (
  build_candidate_domain_audit,
  score_decomposition,
)
from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  expanded_onco_score,
)
from scripts.run_rses_robustness_analyses import recompute


def test_noneligible_domains_do_not_enter_eligible_denominator() -> None:
  components = {
    "tumor_event": None,
    "dependency": None,
    "selectivity": None,
    "expression_compensation": None,
    "functional_relation": 0.80,
    "functional_microniche": None,
    "validation_tractability": 0.40,
  }
  eligible = {"functional_relation", "validation_tractability"}
  score = expanded_onco_score(components, eligible_domains=eligible)
  expected = (
    EXPANDED_ONCO_WEIGHTS["functional_relation"] * 0.80
    + EXPANDED_ONCO_WEIGHTS["validation_tractability"] * 0.40
  ) / (
    EXPANDED_ONCO_WEIGHTS["functional_relation"]
    + EXPANDED_ONCO_WEIGHTS["validation_tractability"]
  )
  assert np.isclose(score.observed_score, expected)
  assert np.isclose(score.coverage, 1.0)
  assert np.isclose(score.adjusted_score, expected)
  assert score.eligible_domains == 2


def test_missing_eligible_domain_lowers_coverage_without_zero_imputation() -> None:
  components = {
    "functional_relation": 0.80,
    "validation_tractability": None,
  }
  eligible = {"functional_relation", "validation_tractability"}
  score = expanded_onco_score(components, eligible_domains=eligible)
  expected_coverage = (
    EXPANDED_ONCO_WEIGHTS["functional_relation"]
    / (
      EXPANDED_ONCO_WEIGHTS["functional_relation"]
      + EXPANDED_ONCO_WEIGHTS["validation_tractability"]
    )
  )
  assert np.isclose(score.observed_score, 0.80)
  assert np.isclose(score.coverage, expected_coverage)
  assert np.isclose(score.adjusted_score, 0.80 * expected_coverage)


def composite_ranking() -> pd.DataFrame:
  eligible = {"functional_relation", "validation_tractability"}
  components = {
    "tumor_event": None,
    "dependency": None,
    "selectivity": None,
    "expression_compensation": None,
    "functional_relation": 0.75,
    "functional_microniche": None,
    "validation_tractability": 0.50,
  }
  score = expanded_onco_score(components, eligible_domains=eligible)
  row = {
    "pair_id": "COMPOSITE1",
    "cancer": "lung",
    "source_class": "curated_composite_biomarker",
    "lost_feature": "BRCA1/BRCA2 or HRD",
    "lost_gene": np.nan,
    "analysis_lost_gene": np.nan,
    "target_gene": "PARP1",
    "analysis_target_gene": "PARP1",
    "hypothesis_type": "composite_feature",
    "score_comparability_group": "composite_event_to_gene",
    "scoring_semantics_version": "eligibility-aware-v1",
    "rses_onco": score.observed_score,
    "evidence_coverage": score.coverage,
    "coverage_adjusted_rses": score.adjusted_score,
    "n_domains": score.n_domains,
    "eligible_domains": score.eligible_domains,
    **{f"component_{key}": value for key, value in components.items()},
  }
  for domain in EXPANDED_ONCO_WEIGHTS:
    row[f"eligible_component_{domain}"] = domain in eligible
  for domain in (
    "expression_context",
    "localization",
    "biochemical_structural",
    "genetic_phenotype",
    "interaction_network",
    "regulatory_network",
  ):
    row[f"microniche_{domain}"] = None
    row[f"eligible_microniche_{domain}"] = False
  row.update({
    "functional_microniche_rses": np.nan,
    "functional_microniche_coverage": np.nan,
    "functional_microniche_adjusted": np.nan,
    "functional_microniche_n_domains": 0,
  })
  return pd.DataFrame([row])


def test_composite_audit_marks_empirical_domains_noneligible(tmp_path: Path) -> None:
  audit = build_candidate_domain_audit(composite_ranking(), tmp_path)
  empirical = audit.loc[
    (audit["domain_family"] == "RSES-Onco")
    & audit["domain"].isin(
      {
        "tumor_event",
        "dependency",
        "selectivity",
        "expression_compensation",
        "functional_microniche",
      }
    )
  ]
  assert not empirical.empty
  assert empirical["evidence_state"].eq("not_eligible").all()
  assert empirical["component_normalized"].isna().all()
  assert empirical["final_score_contribution"].isna().all()


def test_composite_score_decomposition_reproduces_pipeline(tmp_path: Path) -> None:
  audit = build_candidate_domain_audit(composite_ranking(), tmp_path)
  decomposition = score_decomposition(audit)
  row = decomposition.loc[
    decomposition["domain_family"].eq("RSES-Onco")
  ].iloc[0]
  assert row["recomputed_eligible_domains"] == 2
  assert np.isclose(
    row["pipeline_adjusted_score"],
    row["recomputed_adjusted_score"],
  )


def test_robustness_recompute_uses_row_specific_eligibility() -> None:
  frame = composite_ranking()
  recomputed = recompute(frame, dict(EXPANDED_ONCO_WEIGHTS))
  assert recomputed.loc[0, "recomputed_eligible_domains"] == 2
  assert np.isclose(
    recomputed.loc[0, "recomputed_adjusted_score"],
    frame.loc[0, "coverage_adjusted_rses"],
  )
