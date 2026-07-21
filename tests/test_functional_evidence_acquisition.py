from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import requests

from scripts.download_human_functional_evidence import (
  RequestFailure,
  cached_string_result,
  discover_string_api_root,
  map_string_ids,
  parse_string_mapping,
  request_with_retries,
)


class DummyResponse:
  def __init__(
    self,
    *,
    status_code: int = 200,
    text: str = "",
    json_payload=None,
    headers: dict[str, str] | None = None,
  ):
    self.status_code = status_code
    self.text = text
    self._json_payload = json_payload
    self.headers = headers or {}

  def raise_for_status(self) -> None:
    if self.status_code >= 400:
      raise requests.HTTPError(
        f"HTTP {self.status_code}",
        response=self,
      )

  def json(self):
    return self._json_payload


class DummySession:
  def __init__(self, responses: list[DummyResponse]):
    self.responses = list(responses)
    self.calls: list[dict[str, object]] = []
    self.headers: dict[str, str] = {}

  def request(
    self,
    method,
    url,
    timeout=None,
    **kwargs,
  ):
    self.calls.append({
      "method": method,
      "url": url,
      "timeout": timeout,
      **kwargs,
    })
    if not self.responses:
      raise AssertionError("Unexpected request")
    return self.responses.pop(0)


def test_parse_string_mapping_uses_echoed_query_items() -> None:
  text = (
    "queryItem\tqueryIndex\tstringId\tncbiTaxonId\ttaxonName\tpreferredName\tannotation\n"
    "ACP1\t0\t9606.ENSP00000357741\t9606\tHomo sapiens\tACP1\tannotation\n"
    "TP53\t1\t9606.ENSP00000269305\t9606\tHomo sapiens\tTP53\tannotation\n"
  )
  mapping = parse_string_mapping(text, ["ACP1", "TP53"])
  assert mapping["ACP1"]["string_id"] == "9606.ENSP00000357741"
  assert mapping["TP53"]["preferred_name"] == "TP53"


def test_non_retryable_404_is_not_retried(monkeypatch) -> None:
  session = DummySession([DummyResponse(status_code=404)])
  monkeypatch.setattr("time.sleep", lambda _: None)
  with pytest.raises(RequestFailure) as error:
    request_with_retries(
      session,
      "POST",
      "https://example.invalid/api",
      retries=7,
    )
  assert error.value.status_code == 404
  assert len(session.calls) == 1


def test_retryable_500_is_retried(monkeypatch) -> None:
  session = DummySession([
    DummyResponse(status_code=500),
    DummyResponse(status_code=200, text="ok"),
  ])
  monkeypatch.setattr("time.sleep", lambda _: None)
  response = request_with_retries(
    session,
    "GET",
    "https://example.invalid/api",
    retries=3,
  )
  assert response.text == "ok"
  assert len(session.calls) == 2


def test_string_mapping_records_unmapped_404_without_crashing(monkeypatch) -> None:
  session = DummySession([DummyResponse(status_code=404)])
  monkeypatch.setattr("time.sleep", lambda _: None)
  result = map_string_ids(
    session,
    ["UNKNOWN1", "UNKNOWN2"],
    api_root="https://version-12-0.string-db.org/api",
    caller_identity="rses-onco-tests",
    retries=7,
    chunk_size=200,
  )
  assert set(result["mapping_status"]) == {"unmapped"}
  assert result["string_id"].eq("").all()
  assert len(session.calls) == 1


def test_string_version_discovery_uses_stable_address() -> None:
  session = DummySession([
    DummyResponse(
      json_payload=[{
        "string_version": "12.0",
        "string_stable_address": "https://version-12-0.string-db.org/",
      }]
    )
  ])
  root, metadata = discover_string_api_root(session, retries=1)
  assert root == "https://version-12-0.string-db.org/api"
  assert metadata["discovery_status"] == "ok"


def test_cached_string_result_reuses_complete_gene_cache(tmp_path: Path) -> None:
  cache_tsv = tmp_path / "ACP1.tsv"
  cache_json = tmp_path / "ACP1.json"
  frame = pd.DataFrame({
    "preferredName_A": ["ACP1"],
    "preferredName_B": ["SRC"],
    "score": [0.9],
    "query_gene": ["ACP1"],
    "query_string_id": ["9606.ENSP00000357741"],
  })
  frame.to_csv(cache_tsv, sep="\t", index=False)
  cache_json.write_text(
    """{
      "api_root": "https://version-12-0.string-db.org/api",
      "string_id": "9606.ENSP00000357741",
      "required_score": 700,
      "limit": 100,
      "network_type": "functional",
      "status": "ok",
      "edge_rows": 1
    }""",
    encoding="utf-8",
  )
  expected = {
    "api_root": "https://version-12-0.string-db.org/api",
    "string_id": "9606.ENSP00000357741",
    "required_score": 700,
    "limit": 100,
    "network_type": "functional",
  }
  cached, metadata = cached_string_result(cache_tsv, cache_json, expected)
  assert cached is not None
  assert len(cached) == 1
  assert metadata is not None
  assert metadata["status"] == "ok"
