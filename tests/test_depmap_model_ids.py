from pathlib import Path

import pandas as pd
import pytest

from rses_onco.depmap import (
  detect_model_id_column,
  normalize_model_id_column,
  read_model_ids,
)


def test_normalize_unnamed_ach_index() -> None:
  frame = pd.DataFrame({"Unnamed: 0": ["ACH-000001", "ACH-000002"], "MTAP (4507)": [0.1, 0.2]})
  result = normalize_model_id_column(frame, "test")
  assert result.columns[0] == "ModelID"
  assert result["ModelID"].tolist() == ["ACH-000001", "ACH-000002"]


def test_normalize_legacy_depmap_id() -> None:
  frame = pd.DataFrame({"DepMap_ID": ["ACH-000001"], "PRMT5 (10419)": [-1.2]})
  result = normalize_model_id_column(frame, "test")
  assert "ModelID" in result
  assert "DepMap_ID" not in result


def test_detect_and_read_unnamed_index(tmp_path: Path) -> None:
  path = tmp_path / "matrix.csv"
  pd.DataFrame(
    {"MTAP (4507)": [0.1, 0.2]},
    index=pd.Index(["ACH-000001", "ACH-000002"]),
  ).to_csv(path)
  raw, inferred = detect_model_id_column(path, "test")
  assert raw.startswith("Unnamed:")
  assert inferred is True
  assert read_model_ids(path, "test").tolist() == ["ACH-000001", "ACH-000002"]


def test_reject_non_model_index() -> None:
  frame = pd.DataFrame({"Unnamed: 0": ["PR-001", "PR-002"], "GENE": [1, 2]})
  with pytest.raises(ValueError):
    normalize_model_id_column(frame, "test")
