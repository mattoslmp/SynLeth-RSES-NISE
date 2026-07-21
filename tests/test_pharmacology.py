from __future__ import annotations

import math

import pandas as pd

from rses_onco.pharmacology import (
  coverage_aware_pharmacology_score,
  evidence_rows_to_components,
  normalize_pchembl,
  normalize_pharos_tdl,
  therapeutic_hypothesis_score,
)
from scripts.prioritize_pharmacology import candidate_drug_rows


def test_pharmacology_missingness_reduces_coverage() -> None:
  score = coverage_aware_pharmacology_score({
    "target_tractability": 1.0,
    "direct_interaction": 0.8,
    "compound_potency": None,
    "clinical_maturity": None,
    "cancer_relevance": 0.6,
    "biomarker_selective_response": None,
  })
  assert score.n_domains == 3
  assert 0 < score.coverage < 1
  assert math.isclose(score.adjusted_score, score.observed_score * score.coverage)


def test_pharmacology_normalizers() -> None:
  assert normalize_pharos_tdl("Tclin") == 1.0
  assert normalize_pharos_tdl("Tdark") == 0.12
  assert normalize_pchembl(9) == 1.0
  assert normalize_pchembl(5) == 0.0


def test_evidence_rows_to_components() -> None:
  evidence = pd.DataFrame([
    {
      "source": "chembl",
      "pchembl_value": 8.0,
      "action_type": "INHIBITOR",
      "max_phase": 2,
    },
    {
      "source": "pharos",
      "target_development_level": "Tchem",
    },
    {
      "source": "open_targets",
      "disease_name": "colorectal carcinoma",
      "tractability_score": 1.0,
      "max_phase": 3,
    },
  ])
  components = evidence_rows_to_components(evidence)
  assert components["target_tractability"] == 1.0
  assert components["direct_interaction"] == 1.0
  assert components["compound_potency"] > 0
  assert components["clinical_maturity"] == 0.75
  assert components["cancer_relevance"] == 1.0
  assert components["biomarker_selective_response"] is None


def test_candidate_drug_rows_preserves_research_boundary() -> None:
  ranking = pd.DataFrame([
    {
      "cancer": "colon",
      "pair_id": "TEST_A_TO_B",
      "analysis_lost_gene": "A",
      "target_gene": "B",
      "coverage_adjusted_rses": 0.8,
    }
  ])
  evidence = pd.DataFrame([
    {
      "source": "chembl",
      "target_gene": "B",
      "drug_name": "Test inhibitor",
      "drug_id": "CHEMBLTEST",
      "action_type": "INHIBITOR",
      "pchembl_value": 8.5,
      "max_phase": 2,
    }
  ])
  result = candidate_drug_rows(ranking, evidence, pd.DataFrame())
  assert len(result) == 1
  assert bool(result.iloc[0]["research_only"])
  assert result.iloc[0]["therapeutic_hypothesis_score"] > 0
  assert "not evidence of clinical efficacy" in result.iloc[0]["interpretation_boundary"]


def test_geometric_concordance_requires_both_layers() -> None:
  assert therapeutic_hypothesis_score(0.81, 0.49) == math.sqrt(0.81 * 0.49)
  assert therapeutic_hypothesis_score(0.81, float("nan")) is None
