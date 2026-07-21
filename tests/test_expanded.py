from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.expanded import (
  build_directed_nise_candidates,
  expression_profile_metrics,
  functional_microniche_score,
  phenotype_profile_metrics,
)
from rses_onco.networks import set_metrics, string_pair_metrics
from rses_onco.utils import atomic_gene_symbols, canonical_gene_name
from scripts.build_expanded_candidate_universe import expand_explicit_target_gene_lists
from scripts.discover_conditional_dependencies import vectorized_target_tests

ROOT = Path(__file__).resolve().parents[1]


def test_missing_gene_symbols_remain_empty() -> None:
  assert canonical_gene_name(None) == ""
  assert canonical_gene_name(float("nan")) == ""
  assert canonical_gene_name("SOD2 (6648)") == "SOD2"


def test_atomic_gene_list_parser_is_strict() -> None:
  assert atomic_gene_symbols("CDK4/CDK6") == ["CDK4", "CDK6"]
  assert atomic_gene_symbols("BRCA1 + BRCA2") == ["BRCA1", "BRCA2"]
  assert atomic_gene_symbols("SOD2 (6648)") == ["SOD2"]
  assert atomic_gene_symbols("BRCA1/BRCA2 or HRD") == []


def test_composite_target_is_expanded_without_splitting_loss_state() -> None:
  source = pd.DataFrame([
    {
      "pair_id": "SL011",
      "lost_feature": "SMARCA4 loss",
      "lost_gene": "SMARCA4",
      "target_gene": "CDK4/CDK6",
      "source_class": "downstream_synthetic_lethality",
    },
    {
      "pair_id": "SL025",
      "lost_feature": "SMARCA4/SMARCA2 dual loss",
      "lost_gene": None,
      "target_gene": "MCL1",
      "source_class": "downstream_apoptotic_dependency",
    },
  ])
  expanded = expand_explicit_target_gene_lists(source)
  sl011 = expanded.loc[expanded["target_group_id"].eq("SL011")]
  assert set(sl011["target_gene"]) == {"CDK4", "CDK6"}
  assert set(sl011["pair_id"]) == {"SL011__TARGET_CDK4", "SL011__TARGET_CDK6"}
  assert set(sl011["target_feature"]) == {"CDK4/CDK6"}
  assert set(sl011["lost_gene"]) == {"SMARCA4"}
  sl025 = expanded.loc[expanded["pair_id"].eq("SL025")].iloc[0]
  assert pd.isna(sl025["lost_gene"])
  assert sl025["lost_feature"] == "SMARCA4/SMARCA2 dual loss"
  assert sl025["target_gene"] == "MCL1"


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


def test_complete_curated_nise_catalogue_is_included() -> None:
  source = pd.read_csv(
    ROOT / "data/processed/human_nise_all_within_activity_pairs_2017.tsv",
    sep="\t",
  )
  result = build_directed_nise_candidates(source)
  genes = set(result["lost_gene"]) | set(result["target_gene"])
  assert len(source) == 101
  assert len(result) == 202
  assert len(genes) == 70
  assert result["group_id"].nunique() == 15


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
    "OncotreePrimaryDisease": [
      "Lung Cancer", "Lung Cancer", "Lung Cancer", "Colorectal Adenocarcinoma"
    ],
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


def test_all_target_tests_adjust_complete_target_family() -> None:
  loss = pd.DataFrame({
    "TARGET_A": [-1.2, -1.1, -1.0, -0.9],
    "TARGET_B": [-0.2, -0.1, -0.1, -0.2],
    "TARGET_C": [np.nan, -0.8, -0.7, -0.9],
  })
  intact = pd.DataFrame({
    "TARGET_A": [-0.1, 0.0, -0.2, -0.1],
    "TARGET_B": [-0.2, -0.2, -0.1, -0.1],
    "TARGET_C": [-0.1, -0.2, -0.1, -0.2],
  })
  tested = vectorized_target_tests(loss, intact, minimum_group_size=3)
  assert len(tested) == 3
  assert set(tested["target_gene"]) == {"TARGET_A", "TARGET_B", "TARGET_C"}
  assert tested["target_family_size"].nunique() == 1
  assert int(tested["target_family_size"].iloc[0]) == 3
  target_a = tested.set_index("target_gene").loc["TARGET_A"]
  assert target_a["delta_effect"] < 0
  assert 0 <= target_a["q_value_bh_within_loss_cancer"] <= 1
