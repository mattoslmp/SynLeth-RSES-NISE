from __future__ import annotations

import pickle

import pandas as pd
import requests

from rses_onco.dorothea import acquire_dorothea, normalize_dorothea_frame


class DummySession:
  def __init__(self, fallback_content: bytes):
    self.fallback_content = fallback_content
    self.headers: dict[str, str] = {}
    self.calls: list[tuple[str, str]] = []

  def request(self, method, url, timeout=None, **kwargs):
    self.calls.append((method, url))
    response = requests.Response()
    response.url = url
    if "omnipathdb.org" in url:
      response.status_code = 502
      response._content = b"bad gateway"
      response.headers = {}
      return response
    response.status_code = 200
    response._content = self.fallback_content
    response.headers = {}
    return response


class NoNetworkSession:
  def __init__(self):
    self.headers: dict[str, str] = {}

  def request(self, *args, **kwargs):
    raise AssertionError("network should not be used when normalized cache exists")


def fallback_frame() -> pd.DataFrame:
  return pd.DataFrame({
    "tf": ["TP53", "MYC", "SP1", "BAD"],
    "target": ["CDKN1A", "TERT", "ALDOA", "DROP"],
    "confidence": ["A", "B", "C", "D"],
    "mor": [1.0, 1.0, 1.0, -1.0],
  })


def test_normalize_official_pickle_columns_and_levels() -> None:
  result = normalize_dorothea_frame(
    fallback_frame(), acquisition_source="unit_test"
  )
  assert set(result["source_genesymbol"]) == {"TP53", "MYC", "SP1"}
  assert set(result["target_genesymbol"]) == {"CDKN1A", "TERT", "ALDOA"}
  assert set(result["dorothea_level"]) == {"A", "B", "C"}
  assert set(result["acquisition_source"]) == {"unit_test"}


def test_omnipath_502_uses_pinned_official_fallback(tmp_path) -> None:
  content = pickle.dumps(fallback_frame())
  session = DummySession(content)
  output = tmp_path / "omnipath_dorothea.tsv"
  metadata = tmp_path / "dorothea_metadata.json"
  fallback_pickle = tmp_path / "dorothea_hs_official.pkl"

  result, summary = acquire_dorothea(
    output,
    metadata_output=metadata,
    fallback_pickle=fallback_pickle,
    retries=1,
    session=session,
    sleep=lambda _: None,
  )

  assert len(result) == 3
  assert summary["status"] == "fallback"
  assert summary["fallback_used"] is True
  assert output.exists()
  assert metadata.exists()
  assert fallback_pickle.read_bytes() == content
  assert any("omnipathdb.org" in url for _, url in session.calls)
  assert any("raw.githubusercontent.com" in url for _, url in session.calls)


def test_normalized_tsv_cache_is_reused_without_network(tmp_path) -> None:
  output = tmp_path / "omnipath_dorothea.tsv"
  normalized = normalize_dorothea_frame(
    fallback_frame(), acquisition_source="first_run"
  )
  normalized.to_csv(output, sep="\t", index=False)

  result, summary = acquire_dorothea(
    output,
    metadata_output=tmp_path / "metadata.json",
    session=NoNetworkSession(),
  )

  assert len(result) == 3
  assert summary["status"] == "cache"
  assert summary["source"] == "cached_normalized_tsv"
