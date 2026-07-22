#!/usr/bin/env python3
"""Download and validate the current Human Protein Atlas subcellular table.

The Human Protein Atlas moved the archive from
``/download/subcellular_location.tsv.zip`` to
``/download/tsv/subcellular_location.tsv.zip``. This utility writes the
uncompressed TSV to the cache path consumed by the resilient functional-evidence
pipeline, allowing a failed run to resume without repeating STRING acquisition.
"""
from __future__ import annotations

import argparse
import io
import json
import random
import time
import zipfile
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
HPA_URLS = (
  "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip",
  "https://www.proteinatlas.org/download/subcellular_location.tsv.zip",
)
RETRYABLE_STATUS = {408, 425, 429, 500, 502, 503, 504}
REQUIRED_COLUMNS = {
  "Gene",
  "Gene name",
  "Reliability",
  "Enhanced",
  "Supported",
  "Approved",
  "Uncertain",
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def retry_delay(attempt: int) -> float:
  return min(30.0, (2.0 ** attempt) + random.uniform(0.0, 0.75))


def request_archive(
  session: requests.Session,
  url: str,
  retries: int,
  timeout: int,
) -> requests.Response:
  last_error: Exception | None = None
  for attempt in range(1, retries + 1):
    try:
      response = session.get(url, timeout=timeout)
      if response.status_code in RETRYABLE_STATUS:
        raise requests.HTTPError(
          f"retryable HTTP {response.status_code}",
          response=response,
        )
      response.raise_for_status()
      return response
    except Exception as exc:
      last_error = exc
      status = (
        exc.response.status_code
        if isinstance(exc, requests.HTTPError) and exc.response is not None
        else None
      )
      retryable = status in RETRYABLE_STATUS or not isinstance(exc, requests.HTTPError)
      if attempt >= retries or not retryable:
        break
      delay = retry_delay(attempt)
      print(
        f"HPA attempt {attempt}/{retries} failed for {url}: {exc}; "
        f"retrying in {delay:.1f}s",
        flush=True,
      )
      time.sleep(delay)
  raise RuntimeError(f"HPA request failed for {url}: {last_error}")


def parse_hpa_archive(content: bytes) -> tuple[pd.DataFrame, str]:
  try:
    archive = zipfile.ZipFile(io.BytesIO(content))
  except zipfile.BadZipFile as exc:
    raise ValueError("HPA response is not a valid ZIP archive") from exc
  with archive:
    members = [
      name
      for name in archive.namelist()
      if name.lower().endswith((".tsv", ".txt"))
    ]
    if not members:
      raise ValueError("HPA archive contains no TSV or TXT member")
    member = members[0]
    payload = archive.read(member)
  frame = pd.read_csv(io.BytesIO(payload), sep="\t")
  missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
  if frame.empty or missing:
    raise ValueError(
      "Unexpected HPA subcellular schema: "
      f"rows={len(frame)}, missing={missing}, columns={frame.columns.tolist()[:20]}"
    )
  return frame, member


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


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--output",
    default="data/raw/human_functional_evidence/hpa_subcellular_location.tsv",
  )
  parser.add_argument(
    "--status-output",
    default="data/raw/human_functional_evidence/hpa_subcellular_status.json",
  )
  parser.add_argument("--retries", type=int, default=5)
  parser.add_argument("--timeout", type=int, default=240)
  parser.add_argument("--refresh", action="store_true")
  args = parser.parse_args()

  output = resolve_path(args.output)
  status_output = resolve_path(args.status_output)

  if output.exists() and not args.refresh:
    frame = pd.read_csv(output, sep="\t")
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    if not frame.empty and not missing:
      print(f"HPA: reused validated cache ({len(frame):,} rows): {output}")
      atomic_write_json(
        {
          "status": "cache",
          "available": True,
          "rows": len(frame),
          "output": str(output),
        },
        status_output,
      )
      return

  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.10.3 HPA subcellular acquisition",
    "Accept": "application/zip, application/octet-stream",
  })
  failures: list[dict[str, str]] = []
  for url in HPA_URLS:
    try:
      response = request_archive(session, url, args.retries, args.timeout)
      frame, member = parse_hpa_archive(response.content)
      atomic_write_frame(frame, output)
      atomic_write_json(
        {
          "status": "downloaded",
          "available": True,
          "rows": len(frame),
          "source_url": url,
          "archive_member": member,
          "output": str(output),
          "failures": failures,
        },
        status_output,
      )
      print(f"HPA: downloaded {len(frame):,} rows from {url}")
      print(f"HPA: wrote {output}")
      return
    except Exception as exc:
      failures.append({"url": url, "detail": str(exc)})
      print(f"HPA URL failed: {url}: {exc}", flush=True)

  atomic_write_json(
    {
      "status": "unavailable",
      "available": False,
      "rows": 0,
      "output": str(output),
      "failures": failures,
    },
    status_output,
  )
  raise SystemExit(
    "Human Protein Atlas subcellular acquisition failed for all known URLs; "
    f"inspect {status_output}"
  )


if __name__ == "__main__":
  main()
