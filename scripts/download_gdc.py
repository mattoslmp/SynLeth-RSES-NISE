#!/usr/bin/env python3
"""Query, download and validate open GDC ASCAT3 primary-tumor gene-level CN files."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import requests

PROJECTS = ["TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-LUAD", "TCGA-LUSC"]
FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"
ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def md5sum(path: Path) -> str:
  digest = hashlib.md5()
  with path.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(block)
  return digest.hexdigest()


def query_project(project: str, workflow: str) -> list[dict]:
  filters = {
    "op": "and",
    "content": [
      {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
      {"op": "in", "content": {"field": "data_category", "value": ["Copy Number Variation"]}},
      {"op": "in", "content": {"field": "data_type", "value": ["Gene Level Copy Number"]}},
      {"op": "in", "content": {"field": "analysis.workflow_type", "value": [workflow]}},
      {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
      {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ],
  }
  params = {
    "filters": json.dumps(filters),
    "fields": (
      "file_id,file_name,data_type,data_format,analysis.workflow_type,file_size,md5sum,"
      "cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type"
    ),
    "format": "JSON",
    "size": "10000",
  }
  response = requests.get(FILES_ENDPOINT, params=params, timeout=180)
  response.raise_for_status()
  return response.json()["data"]["hits"]


def validate_local(path: Path, hit: dict) -> tuple[bool, str]:
  if not path.exists():
    return False, "missing"
  expected_size = int(hit.get("file_size") or -1)
  if expected_size >= 0 and path.stat().st_size != expected_size:
    return False, f"size mismatch {path.stat().st_size} != {expected_size}"
  expected_md5 = str(hit.get("md5sum") or "").lower()
  if expected_md5:
    observed = md5sum(path).lower()
    if observed != expected_md5:
      return False, f"md5 mismatch {observed} != {expected_md5}"
  return True, "valid"


def download_one(hit: dict, target: Path, retries: int) -> None:
  valid, message = validate_local(target, hit)
  if valid:
    print(f"SKIP valid: {target}", flush=True)
    return
  if target.exists():
    print(f"REPLACE invalid ({message}): {target}", flush=True)
  part = target.with_suffix(target.suffix + ".part")
  for attempt in range(1, retries + 1):
    try:
      with requests.get(
        f"{DATA_ENDPOINT}/{hit['file_id']}",
        stream=True,
        timeout=(60, 600),
      ) as response:
        response.raise_for_status()
        with part.open("wb") as handle:
          for chunk in response.iter_content(1024 * 1024):
            if chunk:
              handle.write(chunk)
      part.replace(target)
      valid, message = validate_local(target, hit)
      if not valid:
        raise OSError(message)
      return
    except Exception:
      part.unlink(missing_ok=True)
      if attempt >= retries:
        raise
      wait = 2 ** attempt
      print(f"Retry {attempt}/{retries} after {wait}s: {target.name}", flush=True)
      time.sleep(wait)


def query_manifest(projects: list[str], workflow: str) -> dict[str, list[dict]]:
  all_hits: dict[str, list[dict]] = {}
  for project in projects:
    hits = query_project(project, workflow)
    all_hits[project] = hits
    print(
      f"{project}: {len(hits)} open {workflow} primary-tumor gene-level files",
      flush=True,
    )
  return all_hits


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--output-dir", default="data/raw/gdc")
  parser.add_argument("--manifest-only", action="store_true")
  parser.add_argument("--validate-only", action="store_true")
  parser.add_argument(
    "--use-existing-manifest",
    action="store_true",
    help="Do not requery GDC; download from the existing JSON manifest.",
  )
  parser.add_argument("--manifest", default=None)
  parser.add_argument("--projects", nargs="*", default=PROJECTS)
  parser.add_argument("--workflow", default="ASCAT3")
  parser.add_argument("--retries", type=int, default=3)
  args = parser.parse_args()

  if args.manifest_only and args.validate_only:
    parser.error("--manifest-only and --validate-only cannot be combined")

  outdir = resolve_path(args.output_dir)
  outdir.mkdir(parents=True, exist_ok=True)
  manifest = (
    resolve_path(args.manifest)
    if args.manifest
    else outdir / "gdc_gene_level_copy_number_manifest.json"
  )

  if args.validate_only or args.use_existing_manifest:
    if not manifest.exists():
      raise SystemExit(f"Manifest not found: {manifest}")
    all_hits = json.loads(manifest.read_text(encoding="utf-8"))
    print(f"Using existing manifest: {manifest}", flush=True)
  else:
    all_hits = query_manifest(args.projects, args.workflow)
    manifest.write_text(json.dumps(all_hits, indent=2), encoding="utf-8")
    print(f"Wrote {manifest}", flush=True)
    if args.manifest_only:
      return

  failures: list[str] = []
  for project, hits in all_hits.items():
    project_dir = outdir / project
    project_dir.mkdir(exist_ok=True)
    for index, hit in enumerate(hits, start=1):
      target = project_dir / hit["file_name"]
      if args.validate_only:
        valid, message = validate_local(target, hit)
        print(
          f"{'OK' if valid else 'FAIL'} [{project} {index}/{len(hits)}] "
          f"{target.name}: {message}",
          flush=True,
        )
        if not valid:
          failures.append(str(target))
      else:
        print(f"[{project} {index}/{len(hits)}] {target.name}", flush=True)
        download_one(hit, target, args.retries)
  if failures:
    raise SystemExit(f"Validation failed for {len(failures)} files")
  print("All selected GDC files passed size and MD5 validation.", flush=True)


if __name__ == "__main__":
  main()
