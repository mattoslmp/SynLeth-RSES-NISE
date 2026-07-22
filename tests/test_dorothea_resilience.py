from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.download_human_functional_evidence_resilient import (
  acquire_dorothea,
  load_complete_string,
)


def valid_dorothea_frame() -> pd.DataFrame:
  return pd.DataFrame({
    "source_genesymbol": ["TP53"],
    "target_genesymbol": ["CDKN1A"],
    "dorothea_level": ["A"],
    "sources": ["DoRothEA"],
    "references": ["12345678"],
  })


def test_dorothea_reuses_valid_cache(tmp_path: Path) -> None:
  output = tmp_path / "dorothea.tsv"
  status_output = tmp_path / "dorothea_status.json"
  valid_dorothea_frame().to_csv(output, sep="\t", index=False)

  frame, summary = acquire_dorothea(
    output,
    status_output,
    local_file=None,
    retries=1,
    refresh=False,
    strict=True,
  )

  assert len(frame) == 1
  assert summary["status"] == "cache"
  assert summary["available"] is True
  assert json.loads(status_output.read_text())["status"] == "cache"


def test_local_dorothea_file_is_copied_and_validated(tmp_path: Path) -> None:
  local_file = tmp_path / "local.tsv"
  output = tmp_path / "dorothea.tsv"
  status_output = tmp_path / "dorothea_status.json"
  valid_dorothea_frame().to_csv(local_file, sep="\t", index=False)

  frame, summary = acquire_dorothea(
    output,
    status_output,
    local_file=local_file,
    retries=1,
    refresh=True,
    strict=True,
  )

  assert len(frame) == 1
  assert output.exists()
  assert summary["status"] == "local_file"


def test_dorothea_outage_becomes_explicit_missingness(
  tmp_path: Path,
  monkeypatch,
) -> None:
  output = tmp_path / "dorothea.tsv"
  status_output = tmp_path / "dorothea_status.json"

  def fail(*args, **kwargs):
    raise RuntimeError("HTTP 502")

  monkeypatch.setattr(
    "scripts.download_human_functional_evidence_resilient.request_with_retries",
    fail,
  )

  frame, summary = acquire_dorothea(
    output,
    status_output,
    local_file=None,
    retries=1,
    refresh=True,
    strict=False,
  )

  assert frame.empty
  assert {"source_genesymbol", "target_genesymbol"}.issubset(frame.columns)
  assert summary["status"] == "unavailable"
  assert summary["available"] is False
  assert len(summary["attempts"]) == 8
  persisted = json.loads(status_output.read_text())
  assert persisted["scientific_interpretation"].startswith(
    "Regulatory-network evidence is missing"
  )


def test_strict_dorothea_outage_fails_after_status_write(
  tmp_path: Path,
  monkeypatch,
) -> None:
  output = tmp_path / "dorothea.tsv"
  status_output = tmp_path / "dorothea_status.json"

  monkeypatch.setattr(
    "scripts.download_human_functional_evidence_resilient.request_with_retries",
    lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("HTTP 502")),
  )

  with pytest.raises(RuntimeError, match="DoRothEA remained unavailable"):
    acquire_dorothea(
      output,
      status_output,
      local_file=None,
      retries=1,
      refresh=True,
      strict=True,
    )

  assert output.exists()
  assert status_output.exists()


def test_complete_string_output_is_reused(tmp_path: Path) -> None:
  string_path = tmp_path / "string.tsv"
  status_path = tmp_path / "status.tsv"
  pd.DataFrame({
    "preferredName_A": ["TP53"],
    "preferredName_B": ["CDKN1A"],
    "score": [0.9],
    "query_gene": ["TP53"],
    "query_string_id": ["9606.ENSP00000269305"],
  }).to_csv(string_path, sep="\t", index=False)
  pd.DataFrame({
    "query_gene": ["TP53", "CDKN1A"],
    "status": ["ok", "no_interactions"],
  }).to_csv(status_path, sep="\t", index=False)

  edges, summary = load_complete_string(
    ["CDKN1A", "TP53"],
    string_path,
    status_path,
  )

  assert edges is not None
  assert summary is not None
  assert summary["discovery_status"] == "complete_existing_acquisition"
  assert len(edges) == 1
