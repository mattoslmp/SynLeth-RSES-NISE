#!/usr/bin/env python3
"""Acquire open GDC/TCGA methylation beta-value files and official probe annotations.

The script queries the GDC API for primary-tumor `Methylation Beta Value` files in
TCGA-COAD, TCGA-READ, TCGA-STAD, TCGA-LUAD and TCGA-LUSC. It writes a traceable
manifest, downloads files resumably with size/MD5 validation, and retrieves the
official GDC GENCODE-v36 HM27, HM450 and EPIC probe-annotation manifests.

Examples
--------
Create a manifest only::

  python scripts/download_gdc_methylation.py --stage manifest

Download/validate all open files and annotation manifests::

  python scripts/download_gdc_methylation.py --stage all
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from typing import Any

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
API = "https://api.gdc.cancer.gov"
PROJECTS = (
  "TCGA-COAD",
  "TCGA-READ",
  "TCGA-STAD",
  "TCGA-LUAD",
  "TCGA-LUSC",
)
ANNOTATIONS = {
  "EPIC": {
    "uuid": "5ce8ae8f-3386-4d12-9035-152742aa07e0",
    "file_name": "EPIC.hg38.manifest.gencode.v36.tsv.gz",
    "md5": "071d925096dce531739cfb955605217b",
  },
  "HM27": {
    "uuid": "e5182c42-bdc6-433e-9b4a-7b7c6696ce89",
    "file_name": "HM27.hg38.manifest.gencode.v36.tsv.gz",
    "md5": "9d4e032a9bd13127ffb9782f66450fd6",
  },
  "HM450": {
    "uuid": "021a2330-951d-474f-af24-1acd77e7664f",
    "file_name": "HM450.hg38.manifest.gencode.v36.tsv.gz",
    "md5": "e163fc110043abb5a7ef623816383bb9",
  },
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def md5sum(path: Path) -> str:
  digest = hashlib.md5()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def flatten_hit(hit: dict[str, Any]) -> dict[str, Any]:
  cases = hit.get("cases") or []
  case = cases[0] if cases else {}
  samples = case.get("samples") or []
  primary = next(
    (
      sample
      for sample in samples
      if str(sample.get("sample_type", "")).casefold() == "primary tumor"
    ),
    samples[0] if samples else {},
  )
  project = case.get("project") or {}
  analysis = hit.get("analysis") or {}
  return {
    "file_id": hit.get("file_id"),
    "file_name": hit.get("file_name"),
    "md5sum": hit.get("md5sum"),
    "file_size": hit.get("file_size"),
    "data_category": hit.get("data_category"),
    "data_type": hit.get("data_type"),
    "data_format": hit.get("data_format"),
    "experimental_strategy": hit.get("experimental_strategy"),
    "platform": hit.get("platform"),
    "workflow_type": analysis.get("workflow_type"),
    "project_id": project.get("project_id"),
    "case_id": case.get("case_id"),
    "case_submitter_id": case.get("submitter_id"),
    "sample_id": primary.get("sample_id"),
    "sample_submitter_id": primary.get("submitter_id"),
    "sample_type": primary.get("sample_type"),
    "access": hit.get("access"),
  }


def query_manifest(session: requests.Session) -> pd.DataFrame:
  filters = {
    "op": "and",
    "content": [
      {
        "op": "in",
        "content": {
          "field": "cases.project.project_id",
          "value": list(PROJECTS),
        },
      },
      {
        "op": "=",
        "content": {
          "field": "data_category",
          "value": "DNA Methylation",
        },
      },
      {
        "op": "=",
        "content": {
          "field": "data_type",
          "value": "Methylation Beta Value",
        },
      },
      {
        "op": "=",
        "content": {
          "field": "access",
          "value": "open",
        },
      },
      {
        "op": "in",
        "content": {
          "field": "cases.samples.sample_type",
          "value": ["Primary Tumor"],
        },
      },
    ],
  }
  fields = ",".join([
    "file_id",
    "file_name",
    "md5sum",
    "file_size",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "platform",
    "access",
    "analysis.workflow_type",
    "cases.case_id",
    "cases.submitter_id",
    "cases.project.project_id",
    "cases.samples.sample_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
  ])
  response = session.get(
    f"{API}/files",
    params={
      "filters": json.dumps(filters, separators=(",", ":")),
      "fields": fields,
      "format": "JSON",
      "size": "10000",
    },
    timeout=120,
  )
  response.raise_for_status()
  hits = response.json().get("data", {}).get("hits", [])
  frame = pd.DataFrame(flatten_hit(hit) for hit in hits)
  if frame.empty:
    raise RuntimeError("The GDC API returned no open primary-tumor methylation beta files")
  return frame.sort_values(
    ["project_id", "sample_submitter_id", "platform", "file_id"],
    na_position="last",
  ).reset_index(drop=True)


def download_one(
  session: requests.Session,
  file_id: str,
  destination: Path,
  expected_size: int | None,
  expected_md5: str | None,
  retries: int,
) -> dict[str, Any]:
  destination.parent.mkdir(parents=True, exist_ok=True)
  if destination.exists():
    size_ok = expected_size is None or destination.stat().st_size == expected_size
    md5_ok = expected_md5 is None or md5sum(destination) == expected_md5
    if size_ok and md5_ok:
      return {"status": "reused", "path": str(destination), "file_id": file_id}
  temporary = destination.with_suffix(destination.suffix + ".part")
  last_error: Exception | None = None
  for attempt in range(1, retries + 1):
    try:
      with session.get(f"{API}/data/{file_id}", stream=True, timeout=300) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
          for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
              handle.write(chunk)
      if expected_size is not None and temporary.stat().st_size != expected_size:
        raise RuntimeError(
          f"size mismatch for {file_id}: {temporary.stat().st_size} != {expected_size}"
        )
      if expected_md5 is not None and md5sum(temporary) != expected_md5:
        raise RuntimeError(f"MD5 mismatch for {file_id}")
      temporary.replace(destination)
      return {"status": "downloaded", "path": str(destination), "file_id": file_id}
    except Exception as exc:  # noqa: BLE001 - preserve remote error context
      last_error = exc
      temporary.unlink(missing_ok=True)
      if attempt < retries:
        time.sleep(min(30, 2**attempt))
  raise RuntimeError(f"Failed to download {file_id}: {last_error}")


def download_annotations(
  session: requests.Session,
  annotation_dir: Path,
  retries: int,
) -> list[dict[str, Any]]:
  records = []
  for platform, metadata in ANNOTATIONS.items():
    destination = annotation_dir / str(metadata["file_name"])
    result = download_one(
      session,
      str(metadata["uuid"]),
      destination,
      None,
      str(metadata["md5"]),
      retries,
    )
    result["platform"] = platform
    result["source_role"] = "GDC_GENCODE_v36_methylation_probe_annotation"
    records.append(result)
  return records


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--stage",
    choices=("manifest", "download", "validate", "all"),
    default="manifest",
  )
  parser.add_argument("--output-dir", default="data/raw/methylation")
  parser.add_argument("--manifest", default=None)
  parser.add_argument("--workers", type=int, default=4)
  parser.add_argument("--retries", type=int, default=4)
  args = parser.parse_args()

  output_dir = resolve_path(args.output_dir)
  output_dir.mkdir(parents=True, exist_ok=True)
  manifest_path = (
    resolve_path(args.manifest)
    if args.manifest
    else output_dir / "gdc_methylation_manifest.tsv"
  )
  status_path = output_dir / "gdc_methylation_download_status.tsv"
  metadata_path = output_dir / "gdc_methylation_source_metadata.json"
  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.11.1 methylation acquisition",
    "X-Tool": "RSES-Onco",
  })

  if args.stage in {"manifest", "all"}:
    manifest = query_manifest(session)
    manifest["local_path"] = manifest.apply(
      lambda row: str(
        output_dir
        / "files"
        / str(row["project_id"])
        / str(row["file_id"])
        / str(row["file_name"])
      ),
      axis=1,
    )
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(f"Wrote GDC methylation manifest: {manifest_path} ({len(manifest):,} files)")
  else:
    if not manifest_path.exists() or manifest_path.stat().st_size == 0:
      raise FileNotFoundError(f"Methylation manifest is absent: {manifest_path}")
    manifest = pd.read_csv(manifest_path, sep="\t", low_memory=False)

  results: list[dict[str, Any]] = []
  if args.stage in {"download", "all"}:
    results.extend(download_annotations(session, output_dir / "annotations", args.retries))
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
      futures = []
      for row in manifest.to_dict("records"):
        futures.append(executor.submit(
          download_one,
          requests.Session(),
          str(row["file_id"]),
          resolve_path(str(row["local_path"])),
          int(row["file_size"]) if pd.notna(row.get("file_size")) else None,
          str(row["md5sum"]) if pd.notna(row.get("md5sum")) else None,
          args.retries,
        ))
      for future in as_completed(futures):
        result = future.result()
        results.append(result)
        print(f"[{result['status']}] {result['file_id']}", flush=True)
    pd.DataFrame(results).to_csv(status_path, sep="\t", index=False)

  if args.stage in {"validate", "all"}:
    validation = []
    for row in manifest.to_dict("records"):
      path = resolve_path(str(row["local_path"]))
      expected_size = int(row["file_size"]) if pd.notna(row.get("file_size")) else None
      expected_md5 = str(row["md5sum"]) if pd.notna(row.get("md5sum")) else None
      status = "ok"
      reason = ""
      if not path.exists():
        status = "missing"
        reason = "file_not_downloaded"
      elif expected_size is not None and path.stat().st_size != expected_size:
        status = "invalid"
        reason = "file_size_mismatch"
      elif expected_md5 is not None and md5sum(path) != expected_md5:
        status = "invalid"
        reason = "md5_mismatch"
      validation.append({
        "file_id": row["file_id"],
        "path": str(path),
        "status": status,
        "reason": reason,
      })
    validation_frame = pd.DataFrame(validation)
    validation_frame.to_csv(status_path, sep="\t", index=False)
    if not validation_frame["status"].eq("ok").all():
      raise RuntimeError("One or more GDC methylation files failed validation")

  metadata = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "source": "NCI Genomic Data Commons",
    "data_category": "DNA Methylation",
    "data_type": "Methylation Beta Value",
    "workflow": "GDC SeSAMe Methylation Array Harmonization",
    "projects": list(PROJECTS),
    "sample_type": "Primary Tumor",
    "manifest": str(manifest_path),
    "file_count": int(len(manifest)),
    "annotation_files": ANNOTATIONS,
    "interpretation_boundary": (
      "Promoter beta values are epigenetic context and are not proof of gene silencing."
    ),
  }
  metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
  print(f"Wrote source metadata: {metadata_path}")


if __name__ == "__main__":
  main()
