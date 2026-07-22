from __future__ import annotations

import hashlib
import io
import json
import random
import time
from pathlib import Path
from typing import Any, Callable

import pandas as pd
import requests

from .utils import canonical_gene_name

OMNIPATH_ENDPOINT = "https://omnipathdb.org/interactions"
DOROTHEA_FALLBACK_COMMIT = "833165d3c790ced3a3e3852899e93412c63f0f44"
DOROTHEA_FALLBACK_URL = (
  "https://raw.githubusercontent.com/saezlab/dorothea-py/"
  f"{DOROTHEA_FALLBACK_COMMIT}/dorothea/data/dorothea_hs.pkl"
)
DOROTHEA_LEVELS = ("A", "B", "C")
RETRYABLE_STATUS = {408, 425, 429, 500, 502, 503, 504}


class DorotheaRequestFailure(RuntimeError):
  def __init__(self, message: str, status_code: int | None = None):
    super().__init__(message)
    self.status_code = status_code


def _retry_delay(response: requests.Response | None, attempt: int) -> float:
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
  retries: int = 4,
  timeout: int = 240,
  sleep: Callable[[float], None] = time.sleep,
  **kwargs: object,
) -> requests.Response:
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
      if status_code in RETRYABLE_STATUS:
        raise requests.HTTPError(
          f"retryable HTTP {status_code}", response=response,
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
      delay = _retry_delay(response, attempt)
      print(
        f"DoRothEA request attempt {attempt}/{retries} failed for {url}: "
        f"{exc}; retrying in {delay:.1f}s",
        flush=True,
      )
      sleep(delay)
  raise DorotheaRequestFailure(
    f"Request failed after {attempts_used}/{retries} attempts: "
    f"{method} {url}: {error}",
    status_code=status_code,
  )


def _resolve_column(frame: pd.DataFrame, names: tuple[str, ...]) -> str | None:
  return next((name for name in names if name in frame.columns), None)


def normalize_dorothea_frame(
  frame: pd.DataFrame,
  *,
  acquisition_source: str,
  levels: tuple[str, ...] = DOROTHEA_LEVELS,
) -> pd.DataFrame:
  if not isinstance(frame, pd.DataFrame) or frame.empty:
    raise ValueError("DoRothEA frame is empty")

  source_column = _resolve_column(
    frame,
    ("source_genesymbol", "tf", "source"),
  )
  target_column = _resolve_column(
    frame,
    ("target_genesymbol", "target", "gene"),
  )
  level_column = _resolve_column(
    frame,
    ("dorothea_level", "confidence", "level"),
  )
  if source_column is None or target_column is None:
    raise ValueError(
      "DoRothEA table lacks TF/target gene-symbol columns: "
      f"{frame.columns.tolist()[:20]}"
    )
  if level_column is None:
    raise ValueError(
      "DoRothEA table lacks a confidence-level column: "
      f"{frame.columns.tolist()[:20]}"
    )

  result = frame.copy()
  result["source_genesymbol"] = result[source_column].map(canonical_gene_name)
  result["target_genesymbol"] = result[target_column].map(canonical_gene_name)
  result["dorothea_level"] = (
    result[level_column].astype(str).str.strip().str.upper()
  )
  result = result.loc[result["dorothea_level"].isin(levels)].copy()
  result = result.loc[
    result["source_genesymbol"].ne("")
    & result["target_genesymbol"].ne("")
  ].copy()

  if "mor" not in result.columns:
    if "is_stimulation" in result.columns or "is_inhibition" in result.columns:
      stimulation = pd.to_numeric(
        result.get("is_stimulation", 0), errors="coerce"
      ).fillna(0)
      inhibition = pd.to_numeric(
        result.get("is_inhibition", 0), errors="coerce"
      ).fillna(0)
      result["mor"] = stimulation - inhibition
    else:
      result["mor"] = pd.NA
  result["acquisition_source"] = acquisition_source

  result = result.drop_duplicates(
    ["source_genesymbol", "target_genesymbol", "dorothea_level"],
    keep="first",
  )
  result = result.sort_values(
    ["source_genesymbol", "target_genesymbol", "dorothea_level"]
  ).reset_index(drop=True)
  if result.empty:
    raise ValueError(
      f"DoRothEA table has no A/B/C interactions after normalization from {acquisition_source}"
    )
  return result


def _atomic_write_frame(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def _atomic_write_bytes(content: bytes, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  temporary.write_bytes(content)
  temporary.replace(path)


def _atomic_write_json(payload: dict[str, object], path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  temporary.write_text(
    json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
  )
  temporary.replace(path)


def _read_fallback_pickle(content: bytes) -> pd.DataFrame:
  value: Any = pd.read_pickle(io.BytesIO(content))
  if not isinstance(value, pd.DataFrame):
    raise TypeError(
      f"Pinned DoRothEA pickle contains {type(value).__name__}, expected DataFrame"
    )
  return value


def acquire_dorothea(
  output: Path,
  *,
  metadata_output: Path | None = None,
  fallback_pickle: Path | None = None,
  retries: int = 4,
  refresh: bool = False,
  allow_fallback: bool = True,
  session: requests.Session | None = None,
  sleep: Callable[[float], None] = time.sleep,
) -> tuple[pd.DataFrame, dict[str, object]]:
  output = Path(output)
  metadata_output = metadata_output or output.with_suffix(".metadata.json")
  fallback_pickle = fallback_pickle or output.parent / "dorothea_hs_official.pkl"

  if output.exists() and not refresh:
    cached = normalize_dorothea_frame(
      pd.read_csv(output, sep="\t"),
      acquisition_source="cached_normalized_tsv",
    )
    summary = {
      "status": "cache",
      "source": "cached_normalized_tsv",
      "output": str(output),
      "row_count": len(cached),
      "confidence_levels": list(DOROTHEA_LEVELS),
    }
    _atomic_write_json(summary, metadata_output)
    print(f"DoRothEA: cache ({len(cached):,} A/B/C interactions)", flush=True)
    return cached, summary

  client = session or requests.Session()
  client.headers.update({
    "User-Agent": "RSES-Onco/0.10.3 DoRothEA acquisition",
    "Accept": "text/tab-separated-values, application/octet-stream",
  })
  common = {
    "datasets": "dorothea",
    "genesymbols": 1,
    "format": "tsv",
    "dorothea_levels": ",".join(DOROTHEA_LEVELS),
    "fields": "sources,references,dorothea_level",
  }
  errors: list[str] = []
  server_failure = False
  for organism_key in ("organisms", "organism"):
    params = {**common, organism_key: 9606}
    try:
      response = request_with_retries(
        client,
        "GET",
        OMNIPATH_ENDPOINT,
        params=params,
        retries=retries,
        sleep=sleep,
      )
      frame = normalize_dorothea_frame(
        pd.read_csv(io.StringIO(response.text), sep="\t"),
        acquisition_source="omnipath_webservice",
      )
      _atomic_write_frame(frame, output)
      summary = {
        "status": "ok",
        "source": "omnipath_webservice",
        "endpoint": OMNIPATH_ENDPOINT,
        "organism_parameter": organism_key,
        "row_count": len(frame),
        "confidence_levels": list(DOROTHEA_LEVELS),
        "fallback_used": False,
      }
      _atomic_write_json(summary, metadata_output)
      print(f"DoRothEA: OmniPath ({len(frame):,} A/B/C interactions)", flush=True)
      return frame, summary
    except DorotheaRequestFailure as exc:
      errors.append(f"{organism_key}: {exc}")
      if exc.status_code in {500, 502, 503, 504}:
        server_failure = True
        break
    except Exception as exc:
      errors.append(f"{organism_key}: {exc}")

  if not allow_fallback:
    raise RuntimeError(
      "DoRothEA OmniPath acquisition failed and fallback is disabled; "
      + " | ".join(errors)
    )

  fallback_content: bytes | None = None
  fallback_origin = ""
  if fallback_pickle.exists() and not refresh:
    fallback_content = fallback_pickle.read_bytes()
    fallback_origin = "cached_pinned_official_pickle"
  else:
    try:
      response = request_with_retries(
        client,
        "GET",
        DOROTHEA_FALLBACK_URL,
        retries=retries,
        timeout=300,
        sleep=sleep,
      )
      fallback_content = response.content
      _atomic_write_bytes(fallback_content, fallback_pickle)
      fallback_origin = "pinned_saezlab_dorothea_py"
    except Exception as exc:
      errors.append(f"official_fallback: {exc}")

  if fallback_content is None:
    raise RuntimeError(
      "DoRothEA acquisition failed for OmniPath and the pinned official fallback; "
      + " | ".join(errors)
    )

  frame = normalize_dorothea_frame(
    _read_fallback_pickle(fallback_content),
    acquisition_source=fallback_origin,
  )
  _atomic_write_frame(frame, output)
  observed_sha256 = hashlib.sha256(fallback_content).hexdigest()
  summary = {
    "status": "fallback",
    "source": fallback_origin,
    "omnipath_endpoint": OMNIPATH_ENDPOINT,
    "omnipath_errors": errors,
    "omnipath_server_failure": server_failure,
    "fallback_repository": "saezlab/dorothea-py",
    "fallback_commit": DOROTHEA_FALLBACK_COMMIT,
    "fallback_url": DOROTHEA_FALLBACK_URL,
    "fallback_pickle": str(fallback_pickle),
    "fallback_sha256": observed_sha256,
    "row_count": len(frame),
    "confidence_levels": list(DOROTHEA_LEVELS),
    "fallback_used": True,
  }
  _atomic_write_json(summary, metadata_output)
  print(
    f"DoRothEA: official pinned fallback ({len(frame):,} A/B/C interactions)",
    flush=True,
  )
  return frame, summary
