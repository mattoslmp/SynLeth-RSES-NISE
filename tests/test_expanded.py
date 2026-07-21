from __future__ import annotations

import math

import pandas as pd

from rses_onco.expanded import (
  build_directed_nise_candidates,
  expression_profile_metrics,
  functional_microniche_score,
  phenotype_profile_metrics,
)
from rses_onco.networks import set_metrics, string_pair_metrics
from rses_onco.utils import canonical_gene_name


def test_missing_gene_symbols_remain_empty() -> None:
  assert canonical_gene_name(None) == ""
  assert canonical_gene_name(float("nan")) == ""
  assert canonical_gene_name("SOD2 (6648)") == "SOD2"


def test_directed_nise_universe_contains_both_directions() -> None:
  source = pd.DataFrame([
    {
      "group_id": "HNISE02",
      "ec_number": "1.15.1.1",
      "activity": "superoxide dismutase activity",
      "gene_a": "SOD2",
      "uniprot_a": "P04179",
      "cluster_a": 1,
      "gene_b": "SOD1",
      "uniprot_b": "P00441",
      "cluster_b": 2,
    }
  ])
  result = build_directed_nise_candidates(source)
  assert len(result) == 2
  assert set(zip(result["lost_gene"], result["target_gene"])) == {
    ("SOD2", "SOD1"),
    ("SOD1", "SOD2"),
  }
  assert set(result["source_class"]) == {"NISE"}


def test_functional_microniche_preserves_missingness() -> None:
  result = functional_microniche_score({
    "expression_context": 0.8,
    "localization": None,
    "biochemical_structural": 1.0,
    "genetic_phenotype": None,
    "interaction_network": 0.6,
    "regulatory_network": None,
  })
  assert result.n_domains == 3
  assert 0 < result.coverage < 1
  assert math.isclose(result.adjusted_score, result.observed_score * result.coverage)


def test_expression_and_phenotype_profiles_are_cancer_specific() -> None:
  models = pd.DataFrame({
    "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
    "OncotreeLineage": ["Lung", "Lung", "Lung", "Bowel"],
    "OncotreePrimaryDisease": ["Lung Cancer", "Lung Cancer", "Lung Cancer", "Colorectal Adenocarcinoma"],
    "OncotreeSubtype": ["LUAD", "LUAD", "LUSC", "COAD"],
  })
  expression = pd.DataFrame({
    "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
    "A": [1.0, 2.0, 3.0, 5.0],
    "B": [3.0, 2.0, 1.0, 5.0],
  })
  effect = pd.DataFrame({
    "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
    "A": [-1.0, -0.8, 0.1, 0.0],
    "B": [0.0, -0.1, -1.0, 0.0],
  })
  expression_metrics = expression_profile_metrics(expression, models, "A", "B", "lung")
  phenotype_metrics = phenotype_profile_metrics(effect, models, "A", "B", "lung")
  assert expression_metrics.n_models == 3
  assert expression_metrics.divergence is not None
  assert phenotype_metrics.n_models == 3
  assert phenotype_metrics.divergence is not None


def test_string_neighborhood_divergence() -> None:
  edges = pd.DataFrame({
    "preferredName_A": ["A", "A", "B"],
    "preferredName_B": ["X", "Y", "X"],
    "score": [0.9, 0.8, 0.95],
  })
  metrics = string_pair_metrics(edges, "A", "B")
  assert metrics.shared_count == 1
  assert metrics.exclusive_a_count == 1
  assert metrics.exclusive_b_count == 0
  assert math.isclose(metrics.jaccard or 0, 0.5)
  direct = set_metrics("A", "B", {"X", "Y"}, {"X"})
  assert math.isclose(direct.divergence or 0, 0.5)
