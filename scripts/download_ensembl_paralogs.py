#!/usr/bin/env python3
"""Expand the analyzed human gene universe with Ensembl Compara paralogs.

This is an explicit source-bounded definition of the homologous-paralog class:
all human paralogs returned by the selected Ensembl release for the supplied seed
genes. It does not claim that all human functional backups are paralogs.

The downloader is deliberately resumable and tolerant of transient Ensembl REST
failures. Homology responses are cached per seed gene, target identifiers are
resolved in batches through POST /lookup/id, and unresolved identifiers are
reported instead of aborting the complete workflow on one temporary HTTP 5xx.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import re
import time
from typing import Any, Iterable

import pandas as pd
import requests

from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
ENSEMBL = "https://rest.ensembl.org"
SPECIES = "homo_sapiens"
RETRYABLE_STATUS = {408, 425, 429, 500, 502, 503, 504}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def safe_name(value: str) -> str:
  return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def chunked(values: Iterable[str], size: int) -> Iterable[list[str]]:
  chunk: list[str] = []
  for value in values:
    chunk.append(value)
    if len(chunk) >= size:
      yield chunk
      chunk = []
  if chunk:
    yield chunk


def retry_delay(response: requests.Response | None, attempt: int) -> float:
  if response is not None:
    header = response.headers.get("Retry-After")
    if header:
      try:
        return max(0.0, min(120.0, float(header)))
      except ValueError:
        pass
  return min(60.0, (2.0 ** attempt) + random.uniform(0.0, 0.75))


def request_json(
  session: requests.Session,
  method: str,
  url: str,
  *,
  retries: int = 7,
  params: dict[str, object] | None = None,
  json_body: dict[str, object] | None = None,
  timeout: int = 180,
) -> Any:
  """Request JSON with Retry-After handling and exponential backoff."""
  error: Exception | None = None
  for attempt in range(1, retries + 1):
    response: requests.Response | None = None
    try:
      response = session.request(
        method,
        url,
        headers={
          "Content-Type": "application/json",
          "Accept": "application/json",
        },
        params=params,
        json=json_body,
        timeout=timeout,
      )
      if response.status_code in RETRYABLE_STATUS:
        raise requests.HTTPError(
          f"retryable HTTP {response.status_code}", response=response
        )
      response.raise_for_status()
      return response.json()
    except Exception as exc:
      error = exc
      if attempt == retries:
        break
      delay = retry_delay(response, attempt)
      print(
        f"Ensembl request attempt {attempt}/{retries} failed for {url}: "
        f"{exc}; retrying in {delay:.1f}s",
        flush=True,
      )
      time.sleep(delay)
  raise RuntimeError(
    f"Ensembl request failed after {retries} attempts: {method} {url}: {error}"
  )


def get_json(
  session: requests.Session,
  url: str,
  retries: int = 7,
  params: dict[str, object] | None = None,
) -> Any:
  return request_json(
    session,
    "GET",
    url,
    retries=retries,
    params=params,
  )


def post_json(
  session: requests.Session,
  url: str,
  payload: dict[str, object],
  retries: int = 7,
  params: dict[str, object] | None = None,
) -> Any:
  return request_json(
    session,
    "POST",
    url,
    retries=retries,
    params=params,
    json_body=payload,
  )


def extract_symbol_from_lookup(payload: Any) -> str | None:
  if not isinstance(payload, dict):
    return None
  symbol = canonical_gene_name(payload.get("display_name"))
  return symbol or None


def extract_symbol_from_xrefs(payload: Any) -> str | None:
  if not isinstance(payload, list):
    return None
  preferred = []
  fallback = []
  for record in payload:
    if not isinstance(record, dict):
      continue
    symbol = canonical_gene_name(
      record.get("display_id") or record.get("primary_id")
    )
    if not symbol:
      continue
    database = str(record.get("dbname") or record.get("db_display_name") or "")
    if database.upper().startswith("HGNC"):
      preferred.append(symbol)
    else:
      fallback.append(symbol)
  return preferred[0] if preferred else (fallback[0] if fallback else None)


def lookup_symbol_fallback(
  session: requests.Session,
  ensembl_id: str,
  *,
  retries: int,
) -> tuple[str | None, str | None]:
  """Resolve one ID using GET lookup then HGNC xrefs without raising upstream."""
  errors: list[str] = []
  try:
    payload = get_json(
      session,
      f"{ENSEMBL}/lookup/id/{ensembl_id}",
      retries=retries,
      params={"expand": 0, "species": SPECIES},
    )
    symbol = extract_symbol_from_lookup(payload)
    if symbol:
      return symbol, "lookup_get"
  except Exception as exc:
    errors.append(f"lookup_get={exc}")

  try:
    payload = get_json(
      session,
      f"{ENSEMBL}/xrefs/id/{ensembl_id}",
      retries=retries,
      params={
        "external_db": "HGNC",
        "object_type": "gene",
        "species": SPECIES,
      },
    )
    symbol = extract_symbol_from_xrefs(payload)
    if symbol:
      return symbol, "hgnc_xref"
  except Exception as exc:
    errors.append(f"hgnc_xref={exc}")

  return None, "; ".join(errors) or "no display_name or HGNC xref returned"


def resolve_symbols_batch(
  session: requests.Session,
  identifiers: Iterable[str],
  cache: dict[str, str],
  *,
  batch_size: int = 500,
  retries: int = 7,
) -> tuple[dict[str, str], list[dict[str, str]]]:
  """Resolve Ensembl IDs with POST /lookup/id and per-ID fallback.

  The official POST endpoint accepts up to 1000 identifiers. A smaller default
  batch limits the amount of work repeated after a transient server failure.
  """
  unique_ids = sorted({str(value) for value in identifiers if str(value)})
  resolved = {
    identifier: canonical_gene_name(cache.get(identifier))
    for identifier in unique_ids
    if canonical_gene_name(cache.get(identifier))
  }
  pending = [identifier for identifier in unique_ids if identifier not in resolved]
  failures: list[dict[str, str]] = []

  for batch_number, batch in enumerate(chunked(pending, batch_size), start=1):
    batch_payload: Any = {}
    batch_error: str | None = None
    try:
      batch_payload = post_json(
        session,
        f"{ENSEMBL}/lookup/id",
        {"ids": batch},
        retries=retries,
        params={
          "expand": 0,
          "format": "condensed",
          "object_type": "gene",
          "species": SPECIES,
        },
      )
    except Exception as exc:
      batch_error = str(exc)
      batch_payload = {}

    for identifier in batch:
      symbol = extract_symbol_from_lookup(
        batch_payload.get(identifier) if isinstance(batch_payload, dict) else None
      )
      method = "lookup_post"
      detail = batch_error
      if not symbol:
        symbol, fallback_detail = lookup_symbol_fallback(
          session,
          identifier,
          retries=max(2, min(retries, 4)),
        )
        method = fallback_detail or "unresolved"
        detail = fallback_detail or batch_error
      if symbol:
        resolved[identifier] = symbol
        cache[identifier] = symbol
      else:
        failures.append({
          "ensembl_id": identifier,
          "status": "unresolved",
          "detail": detail or "No symbol returned by batch or fallback endpoints",
        })

    print(
      f"[Ensembl symbol batch {batch_number}] resolved "
      f"{sum(identifier in resolved for identifier in batch)}/{len(batch)}",
      flush=True,
    )

  return resolved, failures


def load_or_fetch_homology(
  session: requests.Session,
  seed: str,
  cache_dir: Path,
  *,
  retries: int,
  refresh: bool,
) -> Any:
  cache_path = cache_dir / f"{safe_name(seed)}.json"
  if cache_path.exists() and not refresh:
    return json.loads(cache_path.read_text(encoding="utf-8"))
  payload = get_json(
    session,
    f"{ENSEMBL}/homology/symbol/{SPECIES}/{seed}",
    retries=retries,
    params={
      "type": "paralogues",
      "target_species": SPECIES,
      "sequence": "none",
    },
  )
  cache_path.parent.mkdir(parents=True, exist_ok=True)
  cache_path.write_text(
    json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
  )
  return payload


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--output",
    default="data/raw/ensembl/human_seed_paralogs.tsv",
  )
  parser.add_argument(
    "--cache",
    default="data/raw/ensembl/ensembl_symbol_cache.json",
  )
  parser.add_argument(
    "--homology-cache-dir",
    default="data/raw/ensembl/homology_cache",
  )
  parser.add_argument(
    "--unresolved-output",
    default="data/raw/ensembl/ensembl_unresolved_identifiers.tsv",
  )
  parser.add_argument(
    "--metadata-output",
    default="data/raw/ensembl/ensembl_acquisition_metadata.json",
  )
  parser.add_argument("--sleep", type=float, default=0.15)
  parser.add_argument("--retries", type=int, default=7)
  parser.add_argument("--batch-size", type=int, default=500)
  parser.add_argument("--refresh", action="store_true")
  parser.add_argument(
    "--strict-completeness",
    action="store_true",
    help="Fail after writing outputs when any seed or target identifier remains unresolved.",
  )
  args = parser.parse_args()

  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t")
  genes = sorted({
    canonical_gene_name(value)
    for column in ("lost_gene", "target_gene")
    if column in candidates
    for value in candidates[column].dropna()
    if canonical_gene_name(value)
  })
  output = resolve_path(args.output)
  cache_path = resolve_path(args.cache)
  homology_cache_dir = resolve_path(args.homology_cache_dir)
  unresolved_path = resolve_path(args.unresolved_output)
  metadata_path = resolve_path(args.metadata_output)
  for path in (
    output.parent,
    cache_path.parent,
    homology_cache_dir,
    unresolved_path.parent,
    metadata_path.parent,
  ):
    path.mkdir(parents=True, exist_ok=True)

  if cache_path.exists():
    cache = json.loads(cache_path.read_text(encoding="utf-8"))
  else:
    cache = {}

  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.10 Ensembl Compara acquisition",
    "Accept": "application/json",
  })
  info = get_json(session, f"{ENSEMBL}/info/data", retries=args.retries)

  pending_homologies: list[dict[str, Any]] = []
  seed_failures: list[dict[str, str]] = []
  target_ids: set[str] = set()
  for index, seed in enumerate(genes, start=1):
    try:
      payload = load_or_fetch_homology(
        session,
        seed,
        homology_cache_dir,
        retries=args.retries,
        refresh=args.refresh,
      )
    except Exception as exc:
      seed_failures.append({
        "seed_gene": seed,
        "status": "homology_request_failed",
        "detail": str(exc),
      })
      print(
        f"[Ensembl {index}/{len(genes)}] {seed}: request failed after retries; "
        "recorded for resumption",
        flush=True,
      )
      continue

    homologies: list[dict[str, Any]] = []
    if isinstance(payload, dict):
      for item in payload.get("data", []):
        if isinstance(item, dict):
          homologies.extend(
            value for value in item.get("homologies", []) if isinstance(value, dict)
          )
    for homology in homologies:
      target = homology.get("target", {})
      if not isinstance(target, dict):
        continue
      target_id = str(target.get("id") or "")
      if not target_id:
        continue
      target_ids.add(target_id)
      pending_homologies.append({
        "seed": seed,
        "target_id": target_id,
        "target_percent_identity": target.get("perc_id"),
        "homology_type": homology.get("type"),
      })

    print(
      f"[Ensembl {index}/{len(genes)}] {seed}: {len(homologies)} paralog records",
      flush=True,
    )
    if args.sleep:
      time.sleep(args.sleep)

  symbols, symbol_failures = resolve_symbols_batch(
    session,
    target_ids,
    cache,
    batch_size=args.batch_size,
    retries=args.retries,
  )
  cache_path.write_text(
    json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8"
  )

  rows: list[dict[str, object]] = []
  for record in pending_homologies:
    seed = str(record["seed"])
    target_id = str(record["target_id"])
    target_symbol = symbols.get(target_id)
    if not target_symbol or target_symbol == seed:
      continue
    rows.append({
      "pair_id": f"ENSEMBL_PARALOG_{seed}_TO_{target_symbol}",
      "lost_feature": f"{seed} loss",
      "lost_gene": seed,
      "target_gene": target_symbol,
      "source_class": "homologous_paralog",
      "relation_type": "ensembl_compara_paralog",
      "mechanism": (
        f"Ensembl Compara identifies {seed} and {target_symbol} as human paralogs; "
        "the pair is evaluated as a candidate homologous backup, not assumed synthetic lethal."
      ),
      "colon": 1,
      "stomach": 1,
      "lung": 1,
      "relation_confidence": 0.85,
      "genetic_screen": 0.0,
      "isogenic_validation": 0.0,
      "in_vivo": 0.0,
      "clinical_tractability": 0.20,
      "lineage_relevance": 0.0,
      "evidence_stage": "systematic discovery",
      "primary_doi": "",
      "supporting_doi": "",
      "status": "unvalidated Ensembl paralog direction",
      "ensembl_target_id": target_id,
      "ensembl_homology_type": record.get("homology_type"),
      "ensembl_target_percent_identity": record.get("target_percent_identity"),
    })

  result = pd.DataFrame(rows)
  if not result.empty:
    result = result.drop_duplicates(["lost_gene", "target_gene", "source_class"])
    result = result.sort_values(["lost_gene", "target_gene"])
  result.to_csv(output, sep="\t", index=False)

  unresolved = pd.DataFrame([
    *seed_failures,
    *symbol_failures,
  ])
  unresolved.to_csv(unresolved_path, sep="\t", index=False)
  complete = not seed_failures and not symbol_failures
  metadata_path.write_text(
    json.dumps(
      {
        "endpoint": ENSEMBL,
        "species": SPECIES,
        "release_information": info,
        "seed_gene_count": len(genes),
        "successful_seed_homology_queries": len(genes) - len(seed_failures),
        "failed_seed_homology_queries": len(seed_failures),
        "target_identifier_count": len(target_ids),
        "resolved_target_identifier_count": len(symbols),
        "unresolved_target_identifier_count": len(symbol_failures),
        "directed_paralog_count": len(result),
        "complete": complete,
        "lookup_strategy": "POST /lookup/id batches with GET lookup and HGNC xref fallback",
        "homology_cache_dir": str(homology_cache_dir),
      },
      indent=2,
      sort_keys=True,
    ),
    encoding="utf-8",
  )
  print(f"Seed genes: {len(genes):,}")
  print(f"Successful homology queries: {len(genes) - len(seed_failures):,}")
  print(f"Resolved target identifiers: {len(symbols):,}/{len(target_ids):,}")
  print(f"Directed Ensembl paralog candidates: {len(result):,}")
  print(f"Complete acquisition: {complete}")
  print(f"Wrote {output}")
  print(f"Wrote {unresolved_path}")
  print(f"Wrote {metadata_path}")

  if args.strict_completeness and not complete:
    raise SystemExit(
      "Ensembl acquisition was incomplete after retries. Outputs and resume caches "
      f"were written; inspect {unresolved_path} and rerun the same command."
    )


if __name__ == "__main__":
  main()
