#!/usr/bin/env python3
"""Download DoRothEA from OmniPath's official static-table backup.

The primary OmniPath REST endpoint occasionally returns HTTP 5xx responses. The
OmnipathR project publishes an official no-TLS static backup specifically for
server or client connectivity failures. This script downloads the human
DoRothEA interaction table, validates the schema, filters confidence levels,
and writes a reusable local TSV for the expanded RSES-Onco pipeline.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import random
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATIC_URL = (
  "http://no-tls.static.omnipathdb.org/resources/"
  "interactions_dorothea_9606.tsv.gz"
)
DEFAULT_OUTPUT = (
  "data/raw/human_functional_evidence/"
  "dorothea_official_static_9606.tsv"
)
DEFAULT_STATUS_OUTPUT = (
  "data/raw/human_functional_evidence/"
  "dorothea_official_static_9606_status.json"
)
RETRYABLE_STATUS = {408, 425, 429, 500, 502, 503, 504}
SOURCE_COLUMNS = ("source_genesymbol", "source", "tf")
TARGET_COLUMNS = ("target_genesymbol", "target", "gene")


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_write_frame(frame: pd.DataFrame, output: Path) -> None:
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)


def atomic_write_json(payload: dict[str, object], output: Path) -> None:
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  temporary.write_text(
    json.dumps(payload, indent=2, sort_keys=True),
    encoding="utf-8",
  )
  temporary.replace(output)


def normalize_levels(value: str) -> tuple[str, ...]:
  levels = tuple(
    dict.fromkeys(
      item.strip().upper()
      for item in value.split(",")
      if item.strip()
    )
  )
  if not levels:
    raise ValueError("At least one DoRothEA confidence level is required")
  invalid = sorted(set(levels) - {"A", "B", "C", "D", "E"})
  if invalid:
    raise ValueError(f"Unsupported DoRothEA confidence levels: {invalid}")
  return levels


def validate_and_filter(
  frame: pd.DataFrame,
  levels: tuple[str, ...],
) -> pd.DataFrame:
  source_present = any(column in frame.columns for column in SOURCE_COLUMNS)
  target_present = any(column in frame.columns for column in TARGET_COLUMNS)
  required_missing = []
  if not source_present:
    required_missing.append(f"one of {SOURCE_COLUMNS}")
  if not target_present:
    required_missing.append(f"one of {TARGET_COLUMNS}")
  if "dorothea_level" not in frame.columns:
    required_missing.append("dorothea_level")
  if frame.empty or required_missing:
    raise ValueError(
      "Unexpected official static DoRothEA schema: "
      f"rows={len(frame)}, missing={required_missing}, "
      f"columns={frame.columns.tolist()[:30]}"
    )

  pattern = f"[{re.escape(''.join(levels))}]"
  mask = (
    frame["dorothea_level"]
      .fillna("")
      .astype(str)
      .str.upper()
      .str.contains(pattern, regex=True)
  )
  filtered = frame.loc[mask].copy()
  if filtered.empty:
    raise ValueError(
      "Official static DoRothEA table contains no records for levels "
      f"{','.join(levels)}"
    )
  filtered = filtered.drop_duplicates().reset_index(drop=True)
  return filtered


def parse_static_archive(
  content: bytes,
  levels: tuple[str, ...],
) -> pd.DataFrame:
  if not content:
    raise ValueError("Official static DoRothEA response is empty")
  if content.lstrip().lower().startswith((b"<html", b"<!doctype")):
    raise ValueError("Official static DoRothEA response is HTML, not gzip TSV")
  try:
    payload = gzip.decompress(content)
  except (OSError, EOFError) as exc:
    raise ValueError(
      "Official static DoRothEA response is not a valid gzip archive"
    ) from exc
  frame = pd.read_csv(io.BytesIO(payload), sep="\t", low_memory=False)
  return validate_and_filter(frame, levels)


def read_valid_cache(
  output: Path,
  levels: tuple[str, ...],
) -> pd.DataFrame:
  frame = pd.read_csv(output, sep="\t", low_memory=False)
  return validate_and_filter(frame, levels)


def retry_delay(attempt: int) -> float:
  return min(30.0, (2.0 ** attempt) + random.uniform(0.0, 0.75))


def download_static_archive(
  url: str,
  *,
  retries: int,
  timeout: int,
) -> bytes:
  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.10.4 official DoRothEA static fallback",
    "Accept": "application/gzip, application/octet-stream, text/tab-separated-values",
  })
  last_error: Exception | None = None
  attempts_used = 0
  for attempt in range(1, retries + 1):
    attempts_used = attempt
    response: requests.Response | None = None
    try:
      response = session.get(url, timeout=timeout)
      if response.status_code in RETRYABLE_STATUS:
        raise requests.HTTPError(
          f"retryable HTTP {response.status_code}",
          response=response,
        )
      response.raise_for_status()
      return response.content
    except Exception as exc:
      last_error = exc
      status = (
        exc.response.status_code
        if isinstance(exc, requests.HTTPError) and exc.response is not None
        else None
      )
      retryable = status in RETRYABLE_STATUS or not isinstance(
        exc,
        requests.HTTPError,
      )
      if attempt >= retries or not retryable:
        break
      delay = retry_delay(attempt)
      print(
        f"DoRothEA static attempt {attempt}/{retries} failed: {exc}; "
        f"retrying in {delay:.1f}s",
        flush=True,
      )
      time.sleep(delay)
  raise RuntimeError(
    "Official static DoRothEA request failed after "
    f"{attempts_used}/{retries} attempts: {url}: {last_error}"
  )


def sha256_bytes(content: bytes) -> str:
  return hashlib.sha256(content).hexdigest()


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--url", default=DEFAULT_STATIC_URL)
  parser.add_argument("--output", default=DEFAULT_OUTPUT)
  parser.add_argument("--status-output", default=DEFAULT_STATUS_OUTPUT)
  parser.add_argument("--levels", default="A,B,C")
  parser.add_argument("--retries", type=int, default=4)
  parser.add_argument("--timeout", type=int, default=120)
  parser.add_argument("--refresh", action="store_true")
  args = parser.parse_args()

  output = resolve_path(args.output)
  status_output = resolve_path(args.status_output)
  levels = normalize_levels(args.levels)

  if output.exists() and not args.refresh:
    try:
      frame = read_valid_cache(output, levels)
      summary = {
        "status": "cache",
        "available": True,
        "rows": len(frame),
        "source_url": args.url,
        "output": str(output),
        "levels": list(levels),
      }
      atomic_write_json(summary, status_output)
      print(
        f"DoRothEA official static: reused validated cache "
        f"({len(frame):,} rows): {output}",
        flush=True,
      )
      return
    except Exception as exc:
      print(
        f"DoRothEA official static cache invalid; refreshing: {exc}",
        flush=True,
      )

  print(
    f"DoRothEA official static: downloading {args.url}",
    flush=True,
  )
  content = download_static_archive(
    args.url,
    retries=args.retries,
    timeout=args.timeout,
  )
  frame = parse_static_archive(content, levels)
  atomic_write_frame(frame, output)
  summary = {
    "status": "downloaded_static_official",
    "available": True,
    "rows": len(frame),
    "source_url": args.url,
    "output": str(output),
    "levels": list(levels),
    "archive_sha256": sha256_bytes(content),
    "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
    "scientific_interpretation": (
      "Official OmniPath static DoRothEA backup used because the primary REST "
      "service was unavailable or intentionally bypassed."
    ),
  }
  atomic_write_json(summary, status_output)
  print(
    f"DoRothEA official static: downloaded {len(frame):,} rows",
    flush=True,
  )
  print(f"DoRothEA official static: wrote {output}", flush=True)


if __name__ == "__main__":
  main()
