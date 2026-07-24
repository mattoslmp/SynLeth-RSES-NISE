from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from rses_onco.extended_multiomics import (
  build_functional_loss_table,
  build_gdsc_combination_table,
  build_pair_evidence,
  coverage_consensus,
  read_model_feature_matrix,
)


ROOT = Path(__file__).resolve().parents[1]


def load_recompute_module():
  path = ROOT / "scripts/recompute_rses_with_extended_multiomics.py"
  specification = importlib.util.spec_from_file_location(
    "recompute_rses_with_extended_multiomics",
    path,
  )
  assert specification is not None and specification.loader is not None
  module = importlib.util.module_from_spec(specification)
  specification.loader.exec_module(module)
  return module


def models() -> pd.DataFrame:
  return pd.DataFrame({
    "ModelID": [f"ACH-{index:06d}" for index in range(1, 9)],
    "OncotreeLineage": ["Colon Adenocarcinoma"] * 8,
    "CCLEName": [f"CELL{index}" for index in range(1, 9)],
  })


def test_read_model_feature_matrix_model_by_gene(
  tmp_path: Path,
) -> None:
  frame = pd.DataFrame({
    "ModelID": models()["ModelID"],
    "GENEA (1)": np.arange(8),
    "GENEB (2)": np.arange(8) + 1,
  })
  path = tmp_path / "matrix.csv"
  frame.to_csv(path, index=False)
  matrix = read_model_feature_matrix(path, models(), gene_features=True)
  assert matrix.shape == (8, 2)
  assert list(matrix.columns) == ["GENEA", "GENEB"]
  assert matrix.index[0] == "ACH-000001"


def test_read_model_feature_matrix_feature_by_model(
  tmp_path: Path,
) -> None:
  frame = pd.DataFrame({
    "feature": ["GENEA (1)", "GENEB (2)"],
    **{
      model_id: [index, index + 1]
      for index, model_id in enumerate(models()["ModelID"])
    },
  })
  path = tmp_path / "transposed.csv"
  frame.to_csv(path, index=False)
  matrix = read_model_feature_matrix(path, models(), gene_features=True)
  assert matrix.shape == (8, 2)
  assert {"GENEA", "GENEB"} == set(matrix.columns)


def test_functional_loss_and_pair_evidence() -> None:
  model_table = models()
  ids = model_table["ModelID"].tolist()
  candidates = pd.DataFrame({
    "pair_id": ["P1"],
    "lost_gene": ["GENEA"],
    "target_gene": ["GENEB"],
    "cancer": ["colon"],
  })
  relative_cn = pd.DataFrame({
    "GENEA": [0.0, 0.1, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0],
  }, index=ids)
  absolute_cn = pd.DataFrame({
    "GENEA": [0.0, 0.2, 0.4, 2.0, 2.0, 2.0, 2.0, 2.0],
  }, index=ids)
  loh = pd.DataFrame({
    "GENEA": [1, 1, 1, 0, 0, 0, 0, 0]
  }, index=ids)
  damaging = pd.DataFrame({
    "GENEA": [1, 0, 0, 0, 0, 0, 0, 0]
  }, index=ids)
  loss = build_functional_loss_table(
    candidates,
    model_table,
    relative_cn=relative_cn,
    absolute_cn=absolute_cn,
    loh=loh,
    damaging=damaging,
  )
  assert len(loss) == 16
  gene_a = loss.loc[loss["gene"].eq("GENEA")]
  assert int(
    gene_a["functional_loss_state"]
      .eq("biallelic_or_homdel")
      .sum()
  ) >= 3

  dependency = pd.DataFrame({
    "GENEB": [0.95, 0.90, 0.85, 0.10, 0.15, 0.20, 0.10, 0.05],
  }, index=ids)
  protein = pd.DataFrame({
    "GENEB": [3.0, 3.2, 3.1, 1.0, 1.1, 0.9, 1.0, 1.2],
  }, index=ids)
  rnai = pd.DataFrame({
    "GENEB": [-1.2, -1.1, -1.0, -0.1, -0.2, -0.1, -0.1, -0.2],
  }, index=ids)
  evidence, protein_evidence, _ = build_pair_evidence(
    candidates,
    model_table,
    loss,
    dependency_probability=dependency,
    proteomics={"proteomics_test": protein},
    rnai=rnai,
    min_group_size=3,
  )
  row = evidence.iloc[0]
  assert row["integrated_functional_loss_support"] > 0
  assert row["dependency_probability_support"] > 0.5
  assert row["protein_compensation_support"] > 0.5
  assert row["rnai_orthogonal_support"] > 0.5
  assert protein_evidence.iloc[0]["status"] == "ok"


def test_gdsc_combination_alignment_and_bliss() -> None:
  index = ["ACH-000001", "ACH-000002"]
  matrices = {
    "anchor_viability": pd.DataFrame({
      "DrugA + DrugB": [0.8, 0.9]
    }, index=index),
    "library_viability": pd.DataFrame({
      "DrugA + DrugB": [0.7, 0.8]
    }, index=index),
    "combination_viability": pd.DataFrame({
      "DrugA + DrugB": [0.3, 0.5]
    }, index=index),
    "library_auc": pd.DataFrame({
      "DrugA + DrugB": [0.7, 0.8]
    }, index=index),
    "combination_auc": pd.DataFrame({
      "DrugA + DrugB": [0.4, 0.6]
    }, index=index),
  }
  result = build_gdsc_combination_table(
    matrices,
    {key: f"{key}.csv" for key in matrices},
  )
  assert len(result) == 2
  assert np.allclose(result["auc_sensitization"], [0.3, 0.2])
  assert (result["bliss_excess_effect"] > 0).all()
  assert result.iloc[0]["anchor_drug"] == "DrugA"
  assert result.iloc[0]["library_drug"] == "DrugB"


def test_coverage_consensus_penalizes_missing_eligible_extension() -> None:
  complete, complete_coverage = coverage_consensus(
    {"baseline": 0.8, "extension": 0.6},
    {"baseline": 0.7, "extension": 0.3},
  )
  missing, missing_coverage = coverage_consensus(
    {"baseline": 0.8, "extension": None},
    {"baseline": 0.7, "extension": 0.3},
  )
  assert complete is not None and missing is not None
  assert complete_coverage == 1.0
  assert missing_coverage == 0.7
  assert missing < 0.8


def test_noneligible_extension_does_not_penalize_baseline() -> None:
  module = load_recompute_module()
  score, coverage = module.combine(
    0.8,
    None,
    baseline_weight=0.7,
    extension_eligible=False,
  )
  assert score == 0.8
  assert coverage == 1.0


def test_missing_eligible_extension_reduces_internal_coverage() -> None:
  module = load_recompute_module()
  score, coverage = module.combine(
    0.8,
    None,
    baseline_weight=0.7,
    extension_eligible=True,
  )
  assert score is not None and score < 0.8
  assert coverage == 0.7
