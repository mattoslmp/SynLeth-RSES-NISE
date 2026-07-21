#!/usr/bin/env python3
"""Acquire human functional-microniche evidence for the expanded universe.

Sources:
- STRING functional interaction partners (taxon 9606);
- OmniPath/DoRothEA transcription-factor target interactions;
- Human Protein Atlas subcellular localization;
- UniProtKB reviewed annotations for biochemical/structural traceability.

Raw responses are cached and pair-level divergence metrics are written separately.
Missing resources remain missing and never become score zero.
"""
from __future__ import annotations

import argparse
import io
import json
import random
import re
import time
import zipfile
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd
import requests

from rses_onco.networks import (
  direct_string_score,
  localization_pair_metrics,
  regulatory_pair_metrics,
  string_pair_metrics,
)
from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
STRING_DISCOVERY_ROOT = "https://string-db.org/api"
STRING_VERSION_ENDPOINT = f"{STRING_DISCOVERY_ROOT}/json/version"
STRING_FALLBACK_STABLE_ROOT = "https://version-12-0.string-db.org/api"
OMNIPATH_ENDPOINT = "https://omnipathdb.org/interactions"
HPA_URL = "https://www.proteinatlas.org/download/subcellular_location.tsv.zip"
UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"
RETRYABLE_STATUS = {408, 425, 429, 500, 502, 503, 504}
STRING_REQUIRED_COLUMNS = [
  "stringId_A",
  "stringId_B",
  "preferredName_A",
  "preferredName_B",
  "ncbiTaxonId",
  "score",
  "nscore",
  "fscore",
  "pscore",
  "ascore",
  "escore",
  "dscore",
  "tscore",
  "query_gene",
  "query_string_id",
]


class RequestFailure(RuntimeError):
  """HTTP/network failure with an optional status code."""

  def __init__(self, message: str, status_code: int | None = None):
    super().__init__(message)
    self.status_code = status_code


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
    retry_after = response.headers.get("Retry-After")
    if retry_after:
      try:
        return max(0.0, min(120.0, float(retry_after)))
      except ValueError:
        pass
  return min(60.0, (2.0 ** attempt) + random.uniform(0.0, 0.75))


def request_with_retries(
  session: requests.Session,
  method: str,
  url: str,
  *,
  retries: int = 7,
  timeout: int = 180,
  **kwargs: object,
) -> requests.Response:
  """Request with retry only for network failures, rate limits and HTTP 5xx."""
  error: Exception | None = None
  status_code: int | None = None
  attempts_used = 0
  for attempt in range(1, retries + 1):
    attempts_used = attempt
    response: requests.Response | None = None
    status_code = None
    try:
      response = session.request(method, url, timeout=timeout, **kwargs)
      status_code = response.status_code
      if response.status_code in RETRYABLE_STATUS:
        raise requests.HTTPError(
          f"retryable HTTP {response.status_code}",
          response=response,
        )
      response.raise_for_status()
      return response
    except Exception as exc:
      error = exc
      if isinstance(exc, requests.HTTPError) and exc.response is not None:
        status_code = exc.response.status_code
      retryable = (
        status_code in RETRYABLE_STATUS
        or not isinstance(exc, requests.HTTPError)
      )
      if attempt == retries or not retryable:
        break
      delay = retry_delay(response, attempt)
      print(
        f"Request attempt {attempt}/{retries} failed for {url}: "
        f"{exc}; retrying in {delay:.1f}s",
        flush=True,
      )
      time.sleep(delay)
  raise RequestFailure(
    f"Request failed after {attempts_used}/{retries} attempts: "
    f"{method} {url}: {error}",
    status_code=status_code,
  )


def discover_string_api_root(
  session: requests.Session,
  retries: int,
) -> tuple[str, dict[str, object]]:
  """Resolve the current stable STRING version, with a pinned fallback."""
  try:
    response = request_with_retries(
      session,
      "GET",
      STRING_VERSION_ENDPOINT,
      retries=retries,
    )
    payload: Any = response.json()
    records = payload if isinstance(payload, list) else [payload]
    record = next(
      (item for item in records if isinstance(item, dict)),
      {},
    )
    stable_address = str(record.get("string_stable_address") or "").strip()
    if not stable_address:
      raise ValueError("STRING version response lacks string_stable_address")
    api_root = stable_address.rstrip("/")
    if not api_root.endswith("/api"):
      api_root += "/api"
    return api_root, {
      "discovery_status": "ok",
      "string_version": record.get("string_version"),
      "string_stable_address": stable_address,
      "api_root": api_root,
    }
  except Exception as exc:
    return STRING_FALLBACK_STABLE_ROOT, {
      "discovery_status": "fallback",
      "discovery_error": str(exc),
      "string_version": "12.0 fallback",
      "string_stable_address": STRING_FALLBACK_STABLE_ROOT.removesuffix("/api"),
      "api_root": STRING_FALLBACK_STABLE_ROOT,
    }


def parse_string_mapping(
  text: str,
  submitted_genes: list[str],
) -> dict[str, dict[str, str]]:
  """Parse STRING get_string_ids output with echo_query/queryIndex fallback."""
  if not text.strip():
    return {}
  frame = pd.read_csv(io.StringIO(text), sep="\t")
  if "stringId" not in frame.columns:
    raise ValueError(
      f"STRING mapping response lacks stringId: {frame.columns.tolist()}"
    )
  mappings: dict[str, dict[str, str]] = {}
  for record in frame.to_dict("records"):
    query_item = canonical_gene_name(record.get("queryItem"))
    if not query_item and pd.notna(record.get("queryIndex")):
      try:
        query_index = int(record["queryIndex"])
      except (TypeError, ValueError):
        query_index = -1
      if 0 <= query_index < len(submitted_genes):
        query_item = canonical_gene_name(submitted_genes[query_index])
    string_id = str(record.get("stringId") or "").strip()
    if not query_item or not string_id:
      continue
    mappings[query_item] = {
      "string_id": string_id,
      "preferred_name": canonical_gene_name(record.get("preferredName")),
    }
  return mappings


def map_string_ids(
  session: requests.Session,
  genes: list[str],
  *,
  api_root: str,
  caller_identity: str,
  retries: int,
  chunk_size: int,
) -> pd.DataFrame:
  """Map gene symbols to exact STRING IDs in bounded POST batches."""
  rows: list[dict[str, object]] = []
  for batch_number, batch in enumerate(chunked(genes, chunk_size), start=1):
    batch_mappings: dict[str, dict[str, str]] = {}
    batch_status = "mapped"
    detail = ""
    try:
      response = request_with_retries(
        session,
        "POST",
        f"{api_root}/tsv/get_string_ids",
        retries=retries,
        data={
          "identifiers": "\r".join(batch),
          "species": 9606,
          "echo_query": 1,
          "caller_identity": caller_identity,
        },
      )
      batch_mappings = parse_string_mapping(response.text, batch)
    except RequestFailure as exc:
      if exc.status_code == 404:
        batch_status = "unmapped"
      else:
        batch_status = "mapping_request_failed"
      detail = str(exc)
    except Exception as exc:
      batch_status = "mapping_request_failed"
      detail = str(exc)

    mapped_count = 0
    for gene in batch:
      mapped = batch_mappings.get(gene)
      if mapped:
        mapped_count += 1
        rows.append({
          "query_gene": gene,
          "string_id": mapped["string_id"],
          "preferred_name": mapped["preferred_name"],
          "mapping_status": "mapped",
          "detail": "",
        })
      else:
        status = batch_status if batch_status != "mapped" else "unmapped"
        rows.append({
          "query_gene": gene,
          "string_id": "",
          "preferred_name": "",
          "mapping_status": status,
          "detail": detail or "No STRING identifier returned",
        })
    print(
      f"[STRING mapping batch {batch_number}] "
      f"mapped {mapped_count}/{len(batch)}",
      flush=True,
    )
  return pd.DataFrame(rows)


def atomic_write_frame(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def write_json(payload: dict[str, object], path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  temporary.write_text(
    json.dumps(payload, indent=2, sort_keys=True),
    encoding="utf-8",
  )
  temporary.replace(path)


def cached_string_result(
  cache_tsv: Path,
  cache_json: Path,
  expected: dict[str, object],
) -> tuple[pd.DataFrame | None, dict[str, object] | None]:
  if not cache_json.exists():
    return None, None
  try:
    metadata = json.loads(cache_json.read_text(encoding="utf-8"))
  except Exception:
    return None, None
  comparable_keys = (
    "api_root",
    "string_id",
    "required_score",
    "limit",
    "network_type",
  )
  if any(metadata.get(key) != expected.get(key) for key in comparable_keys):
    return None, None
  status = str(metadata.get("status") or "")
  if status == "ok":
    if not cache_tsv.exists():
      return None, None
    try:
      return pd.read_csv(cache_tsv, sep="\t"), metadata
    except Exception:
      return None, None
  if status == "no_interactions":
    return pd.DataFrame(columns=STRING_REQUIRED_COLUMNS), metadata
  return None, None


def download_string(
  genes: list[str],
  output: Path,
  mapping_output: Path,
  status_output: Path,
  cache_dir: Path,
  required_score: int,
  limit: int,
  caller_identity: str,
  sleep_seconds: float,
  retries: int,
  map_chunk_size: int,
  refresh: bool,
  strict_requests: bool,
) -> tuple[pd.DataFrame, dict[str, object]]:
  """Download STRING partners with stable-version mapping and per-gene resume cache."""
  output.parent.mkdir(parents=True, exist_ok=True)
  cache_dir.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.10.2 STRING acquisition",
    "Accept": "text/tab-separated-values, application/json",
  })
  api_root, version_metadata = discover_string_api_root(session, retries)
  print(f"STRING API root: {api_root}", flush=True)

  mapping = map_string_ids(
    session,
    genes,
    api_root=api_root,
    caller_identity=caller_identity,
    retries=retries,
    chunk_size=map_chunk_size,
  )
  atomic_write_frame(mapping, mapping_output)

  frames: list[pd.DataFrame] = []
  status_rows: list[dict[str, object]] = []
  request_failures = 0
  for index, record in enumerate(mapping.to_dict("records"), start=1):
    gene = canonical_gene_name(record.get("query_gene"))
    string_id = str(record.get("string_id") or "").strip()
    mapping_status = str(record.get("mapping_status") or "")
    if mapping_status != "mapped" or not string_id:
      if mapping_status == "mapping_request_failed":
        request_failures += 1
      status_rows.append({
        "query_gene": gene,
        "string_id": string_id,
        "preferred_name": record.get("preferred_name"),
        "status": mapping_status or "unmapped",
        "edge_rows": 0,
        "from_cache": False,
        "detail": record.get("detail"),
        "api_root": api_root,
      })
      print(
        f"[STRING {index}/{len(genes)}] {gene}: "
        f"{mapping_status or 'unmapped'}",
        flush=True,
      )
      continue

    cache_tsv = cache_dir / f"{safe_name(gene)}.tsv"
    cache_json = cache_dir / f"{safe_name(gene)}.json"
    expected = {
      "api_root": api_root,
      "string_id": string_id,
      "required_score": required_score,
      "limit": limit,
      "network_type": "functional",
    }
    cached_frame: pd.DataFrame | None = None
    cached_metadata: dict[str, object] | None = None
    if not refresh:
      cached_frame, cached_metadata = cached_string_result(
        cache_tsv,
        cache_json,
        expected,
      )
    if cached_frame is not None and cached_metadata is not None:
      if not cached_frame.empty:
        frames.append(cached_frame)
      status_rows.append({
        "query_gene": gene,
        "string_id": string_id,
        "preferred_name": record.get("preferred_name"),
        "status": cached_metadata.get("status"),
        "edge_rows": int(cached_metadata.get("edge_rows") or 0),
        "from_cache": True,
        "detail": cached_metadata.get("detail", ""),
        "api_root": api_root,
      })
      print(
        f"[STRING {index}/{len(genes)}] {gene}: "
        f"cache {cached_metadata.get('status')} "
        f"({cached_metadata.get('edge_rows', 0)} rows)",
        flush=True,
      )
      continue

    status = "request_failed"
    detail = ""
    edge_rows = 0
    frame = pd.DataFrame(columns=STRING_REQUIRED_COLUMNS)
    try:
      response = request_with_retries(
        session,
        "POST",
        f"{api_root}/tsv/interaction_partners",
        retries=retries,
        data={
          "identifiers": string_id,
          "species": 9606,
          "required_score": required_score,
          "limit": limit,
          "network_type": "functional",
          "caller_identity": caller_identity,
        },
      )
      text = response.text.strip()
      if text:
        frame = pd.read_csv(io.StringIO(text), sep="\t")
        required = {"preferredName_A", "preferredName_B", "score"}
        missing = sorted(required - set(frame.columns))
        if missing:
          raise ValueError(
            f"Unexpected STRING partner schema; missing {missing}; "
            f"columns={frame.columns.tolist()[:20]}"
          )
        frame["query_gene"] = gene
        frame["query_string_id"] = string_id
        edge_rows = len(frame)
        status = "ok"
        atomic_write_frame(frame, cache_tsv)
        frames.append(frame)
      else:
        status = "no_interactions"
        cache_tsv.unlink(missing_ok=True)
    except Exception as exc:
      request_failures += 1
      detail = str(exc)

    cache_metadata = {
      **expected,
      "query_gene": gene,
      "preferred_name": record.get("preferred_name"),
      "status": status,
      "edge_rows": edge_rows,
      "detail": detail,
    }
    write_json(cache_metadata, cache_json)
    status_rows.append({
      "query_gene": gene,
      "string_id": string_id,
      "preferred_name": record.get("preferred_name"),
      "status": status,
      "edge_rows": edge_rows,
      "from_cache": False,
      "detail": detail,
      "api_root": api_root,
    })
    print(
      f"[STRING {index}/{len(genes)}] {gene}: {status} ({edge_rows} rows)",
      flush=True,
    )
    if sleep_seconds:
      time.sleep(sleep_seconds)

  result = (
    pd.concat(frames, ignore_index=True).drop_duplicates()
    if frames
    else pd.DataFrame(columns=STRING_REQUIRED_COLUMNS)
  )
  for column in STRING_REQUIRED_COLUMNS:
    if column not in result:
      result[column] = pd.Series(dtype=object)
  atomic_write_frame(result, output)
  status_frame = pd.DataFrame(status_rows)
  atomic_write_frame(status_frame, status_output)

  status_counts = (
    status_frame["status"].value_counts(dropna=False).to_dict()
    if not status_frame.empty
    else {}
  )
  summary = {
    **version_metadata,
    "mapping_endpoint": f"{api_root}/tsv/get_string_ids",
    "interaction_partners_endpoint": f"{api_root}/tsv/interaction_partners",
    "candidate_gene_count": len(genes),
    "mapped_gene_count": int(mapping["mapping_status"].eq("mapped").sum()),
    "unmapped_gene_count": int(mapping["mapping_status"].eq("unmapped").sum()),
    "mapping_request_failure_count": int(
      mapping["mapping_status"].eq("mapping_request_failed").sum()
    ),
    "request_failure_count": request_failures,
    "edge_rows": len(result),
    "status_counts": status_counts,
    "mapping_output": str(mapping_output),
    "status_output": str(status_output),
    "cache_dir": str(cache_dir),
  }
  if strict_requests and request_failures:
    raise RuntimeError(
      "STRING acquisition had persistent request failures after writing resume "
      f"caches and status tables; failures={request_failures}; inspect {status_output}"
    )
  return result, summary


def download_dorothea(output: Path) -> pd.DataFrame:
  """Download high/medium-confidence human DoRothEA interactions.

  OmniPath deployments have accepted both ``organisms`` and ``organism`` in
  different client generations. The downloader tries the plural documented web
  service parameter first, validates the returned schema, and falls back to the
  singular spelling without silently accepting an HTML error page.
  """
  output.parent.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  common = {
    "datasets": "dorothea",
    "genesymbols": 1,
    "format": "tsv",
    "dorothea_levels": "A,B,C",
    "fields": "sources,references,dorothea_level",
  }
  errors: list[str] = []
  for organism_key in ("organisms", "organism"):
    params = {**common, organism_key: 9606}
    try:
      response = request_with_retries(
        session,
        "GET",
        OMNIPATH_ENDPOINT,
        params=params,
      )
      text = response.text.strip()
      frame = pd.read_csv(io.StringIO(text), sep="\t")
      source_present = any(
        column in frame.columns for column in ("source_genesymbol", "source", "tf")
      )
      target_present = any(
        column in frame.columns for column in ("target_genesymbol", "target", "gene")
      )
      if frame.empty or not source_present or not target_present:
        raise ValueError(
          f"unexpected OmniPath schema: rows={len(frame)}, columns={frame.columns.tolist()[:12]}"
        )
      output.write_text(text + "\n", encoding="utf-8")
      return frame
    except Exception as exc:
      errors.append(f"{organism_key}: {exc}")
  raise RuntimeError("DoRothEA acquisition failed; " + " | ".join(errors))


def download_hpa(output: Path) -> pd.DataFrame:
  output.parent.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  response = request_with_retries(session, "GET", HPA_URL)
  with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
    members = [
      name for name in archive.namelist()
      if name.lower().endswith((".tsv", ".txt"))
    ]
    if not members:
      raise RuntimeError("HPA archive does not contain a TSV file")
    content = archive.read(members[0])
  output.write_bytes(content)
  return pd.read_csv(output, sep="\t")


def download_uniprot(
  genes: list[str],
  output: Path,
  chunk_size: int = 40,
) -> pd.DataFrame:
  output.parent.mkdir(parents=True, exist_ok=True)
  session = requests.Session()
  frames: list[pd.DataFrame] = []
  fields = ",".join([
    "accession",
    "id",
    "gene_primary",
    "protein_name",
    "ec",
    "cc_function",
    "cc_catalytic_activity",
    "cc_cofactor",
    "cc_subcellular_location",
    "xref_pdb",
  ])
  for start in range(0, len(genes), chunk_size):
    chunk = genes[start:start + chunk_size]
    gene_query = " OR ".join(f"gene_exact:{gene}" for gene in chunk)
    response = request_with_retries(
      session,
      "GET",
      UNIPROT_ENDPOINT,
      params={
        "query": f"({gene_query}) AND organism_id:9606 AND reviewed:true",
        "format": "tsv",
        "fields": fields,
        "size": 500,
      },
    )
    text = response.text.strip()
    if text:
      frames.append(pd.read_csv(io.StringIO(text), sep="\t"))
    print(f"[UniProt {min(start + chunk_size, len(genes))}/{len(genes)}]", flush=True)
  result = (
    pd.concat(frames, ignore_index=True).drop_duplicates()
    if frames
    else pd.DataFrame()
  )
  result.to_csv(output, sep="\t", index=False)
  return result


def load_or_download(
  path: Path,
  downloader: Callable[[], pd.DataFrame],
  skip_download: bool,
) -> pd.DataFrame:
  if skip_download:
    if not path.exists():
      raise FileNotFoundError(
        f"--skip-download requested but file is absent: {path}"
      )
    return pd.read_csv(path, sep="\t")
  return downloader()


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--raw-dir",
    default="data/raw/human_functional_evidence",
  )
  parser.add_argument(
    "--output",
    default="data/processed/expanded_pair_functional_evidence.tsv",
  )
  parser.add_argument("--string-required-score", type=int, default=700)
  parser.add_argument("--string-limit", type=int, default=100)
  parser.add_argument("--string-sleep", type=float, default=1.0)
  parser.add_argument("--string-retries", type=int, default=7)
  parser.add_argument("--string-map-chunk-size", type=int, default=200)
  parser.add_argument("--caller-identity", default="rses-onco")
  parser.add_argument("--refresh-string", action="store_true")
  parser.add_argument("--strict-string-requests", action="store_true")
  parser.add_argument("--skip-download", action="store_true")
  args = parser.parse_args()

  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t")
  required = {"pair_id", "lost_gene", "target_gene"}
  missing = sorted(required - set(candidates.columns))
  if missing:
    raise ValueError(f"Candidate universe lacks columns: {missing}")
  genes = sorted({
    canonical_gene_name(value)
    for column in ("lost_gene", "target_gene")
    for value in candidates[column].dropna()
    if canonical_gene_name(value)
  })

  raw_dir = resolve_path(args.raw_dir)
  raw_dir.mkdir(parents=True, exist_ok=True)
  string_path = raw_dir / "string_interaction_partners.tsv"
  string_mapping_path = raw_dir / "string_id_mapping.tsv"
  string_status_path = raw_dir / "string_acquisition_status.tsv"
  string_cache_dir = raw_dir / "string_cache"
  dorothea_path = raw_dir / "omnipath_dorothea.tsv"
  hpa_path = raw_dir / "hpa_subcellular_location.tsv"
  uniprot_path = raw_dir / "uniprot_reviewed_annotations.tsv"
  metadata_path = raw_dir / "source_metadata.json"

  if args.skip_download:
    if not string_path.exists():
      raise FileNotFoundError(
        f"--skip-download requested but file is absent: {string_path}"
      )
    string_edges = pd.read_csv(string_path, sep="\t")
    string_summary: dict[str, object] = {
      "discovery_status": "skip_download",
      "edge_rows": len(string_edges),
    }
  else:
    string_edges, string_summary = download_string(
      genes,
      string_path,
      string_mapping_path,
      string_status_path,
      string_cache_dir,
      args.string_required_score,
      args.string_limit,
      args.caller_identity,
      args.string_sleep,
      args.string_retries,
      args.string_map_chunk_size,
      args.refresh_string,
      args.strict_string_requests,
    )

  dorothea = load_or_download(
    dorothea_path,
    lambda: download_dorothea(dorothea_path),
    args.skip_download,
  )
  hpa = load_or_download(
    hpa_path,
    lambda: download_hpa(hpa_path),
    args.skip_download,
  )
  uniprot = load_or_download(
    uniprot_path,
    lambda: download_uniprot(genes, uniprot_path),
    args.skip_download,
  )

  rows: list[dict[str, object]] = []
  for record in candidates.to_dict("records"):
    gene_a = canonical_gene_name(record["lost_gene"])
    gene_b = canonical_gene_name(record["target_gene"])
    string_metrics = string_pair_metrics(string_edges, gene_a, gene_b)
    regulatory_metrics = regulatory_pair_metrics(dorothea, gene_a, gene_b)
    localization_metrics = localization_pair_metrics(hpa, gene_a, gene_b)

    structural_divergence = None
    cluster_a = record.get("lost_structural_cluster")
    cluster_b = record.get("target_structural_cluster")
    if pd.notna(cluster_a) and pd.notna(cluster_b):
      structural_divergence = float(str(cluster_a) != str(cluster_b))

    rows.append({
      "pair_id": record["pair_id"],
      "lost_gene": gene_a,
      "target_gene": gene_b,
      "source_class": record.get("source_class"),
      "component_localization": localization_metrics.divergence,
      "localization_shared": localization_metrics.shared_count,
      "localization_exclusive_lost": localization_metrics.exclusive_a_count,
      "localization_exclusive_target": localization_metrics.exclusive_b_count,
      "component_biochemical_structural": structural_divergence,
      "component_interaction_network": string_metrics.divergence,
      "string_direct_score": direct_string_score(string_edges, gene_a, gene_b),
      "string_neighbor_jaccard": string_metrics.jaccard,
      "string_shared_neighbors": string_metrics.shared_count,
      "component_regulatory_network": regulatory_metrics.divergence,
      "regulator_jaccard": regulatory_metrics.jaccard,
      "shared_regulators": regulatory_metrics.shared_count,
    })

  result = pd.DataFrame(rows)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  result.to_csv(output, sep="\t", index=False)
  metadata_path.write_text(
    json.dumps(
      {
        "string": string_summary,
        "string_species": 9606,
        "string_required_score": args.string_required_score,
        "string_limit_per_query": args.string_limit,
        "string_mapping_output": str(string_mapping_path),
        "string_status_output": str(string_status_path),
        "omnipath_endpoint": OMNIPATH_ENDPOINT,
        "omnipath_dataset": "dorothea",
        "omnipath_confidence_levels": ["A", "B", "C"],
        "hpa_url": HPA_URL,
        "uniprot_endpoint": UNIPROT_ENDPOINT,
        "candidate_gene_count": len(genes),
        "string_edge_rows": len(string_edges),
        "dorothea_rows": len(dorothea),
        "hpa_rows": len(hpa),
        "uniprot_rows": len(uniprot),
      },
      indent=2,
      sort_keys=True,
    ),
    encoding="utf-8",
  )
  print(f"Candidate genes: {len(genes):,}")
  print(f"STRING edges: {len(string_edges):,}")
  print(f"DoRothEA interactions: {len(dorothea):,}")
  print(f"HPA localization records: {len(hpa):,}")
  print(f"UniProt reviewed records: {len(uniprot):,}")
  print(f"Wrote {len(result):,} pair evidence rows to {output}")
  print(f"Wrote {metadata_path}")


if __name__ == "__main__":
  main()
