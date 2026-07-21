from __future__ import annotations

from scripts.download_ensembl_paralogs import (
  extract_symbol_from_lookup,
  extract_symbol_from_xrefs,
  resolve_symbols_batch,
)


class DummyResponse:
  def __init__(self, payload, status_code: int = 200):
    self._payload = payload
    self.status_code = status_code
    self.headers: dict[str, str] = {}

  def raise_for_status(self) -> None:
    if self.status_code >= 400:
      raise RuntimeError(f"HTTP {self.status_code}")

  def json(self):
    return self._payload


class DummySession:
  def __init__(self, post_payload):
    self.post_payload = post_payload
    self.calls = []

  def request(
    self,
    method,
    url,
    headers=None,
    params=None,
    json=None,
    timeout=None,
  ):
    self.calls.append({
      "method": method,
      "url": url,
      "headers": headers,
      "params": params,
      "json": json,
      "timeout": timeout,
    })
    return DummyResponse(self.post_payload)


def test_extract_symbol_from_lookup() -> None:
  assert extract_symbol_from_lookup({"display_name": "TP53"}) == "TP53"
  assert extract_symbol_from_lookup({"display_name": None}) is None


def test_extract_symbol_from_xrefs_prefers_hgnc_display_id() -> None:
  payload = [
    {"dbname": "EntrezGene", "display_id": "7157"},
    {"dbname": "HGNC", "display_id": "TP53"},
  ]
  assert extract_symbol_from_xrefs(payload) == "TP53"


def test_resolve_symbols_batch_uses_post_lookup() -> None:
  session = DummySession({
    "ENSG00000141510": {"display_name": "TP53"},
    "ENSG00000146648": {"display_name": "EGFR"},
  })
  cache: dict[str, str] = {}
  resolved, failures = resolve_symbols_batch(
    session,
    ["ENSG00000141510", "ENSG00000146648"],
    cache,
    batch_size=1000,
    retries=1,
  )
  assert resolved == {
    "ENSG00000141510": "TP53",
    "ENSG00000146648": "EGFR",
  }
  assert failures == []
  assert session.calls[0]["method"] == "POST"
  assert session.calls[0]["json"] == {
    "ids": ["ENSG00000141510", "ENSG00000146648"]
  }


def test_resolve_symbols_batch_reuses_cache_without_request() -> None:
  session = DummySession({})
  cache = {"ENSG00000141510": "TP53"}
  resolved, failures = resolve_symbols_batch(
    session,
    ["ENSG00000141510"],
    cache,
    retries=1,
  )
  assert resolved == {"ENSG00000141510": "TP53"}
  assert failures == []
  assert session.calls == []
