#!/usr/bin/env python3
"""Query GDC for public gene-level copy-number files for selected TCGA projects.

The GDC API response schema can evolve. This acquisition script records the
complete JSON manifest before downloading and fails loudly if no public files
match. Restricted files require a GDC token and are intentionally not assumed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import requests

PROJECTS = ["TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-LUAD", "TCGA-LUSC"]
FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"
ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def query_project(project: str) -> list[dict]:
  filters = {
    "op": "and",
    "content": [
      {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
      {"op": "in", "content": {"field": "data_category", "value": ["Copy Number Variation"]}},
      {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ],
  }
  params = {
    "filters": json.dumps(filters),
    "fields": "file_id,file_name,data_type,analysis.workflow_type,file_size,md5sum,cases.submitter_id",
    "format": "JSON",
    "size": "10000",
  }
  response = requests.get(FILES_ENDPOINT, params=params, timeout=180)
  response.raise_for_status()
  hits = response.json()["data"]["hits"]
  return [h for h in hits if "gene" in str(h.get("data_type", "")).casefold() and "copy" in str(h.get("data_type", "")).casefold()]


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--output-dir", default="data/raw/gdc")
  parser.add_argument("--manifest-only", action="store_true")
  parser.add_argument("--projects", nargs="*", default=PROJECTS)
  args = parser.parse_args()
  outdir = resolve_path(args.output_dir)
  outdir.mkdir(parents=True, exist_ok=True)
  all_hits = {}
  for project in args.projects:
    hits = query_project(project)
    all_hits[project] = hits
    print(f"{project}: {len(hits)} public gene-level copy-number files")
  manifest = outdir / "gdc_gene_level_copy_number_manifest.json"
  manifest.write_text(json.dumps(all_hits, indent=2), encoding="utf-8")
  if args.manifest_only:
    return
  for project, hits in all_hits.items():
    project_dir = outdir / project
    project_dir.mkdir(exist_ok=True)
    for hit in hits:
      file_id = hit["file_id"]
      target = project_dir / hit["file_name"]
      if target.exists():
        continue
      with requests.get(f"{DATA_ENDPOINT}/{file_id}", stream=True, timeout=600) as response:
        response.raise_for_status()
        with target.open("wb") as handle:
          for chunk in response.iter_content(1024 * 1024):
            handle.write(chunk)


if __name__ == "__main__":
  main()
