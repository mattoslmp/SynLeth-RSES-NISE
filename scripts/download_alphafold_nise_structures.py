#!/usr/bin/env python3
"""Download AlphaFold DB structures for all curated human NISE proteins.

The downloader uses the AlphaFold DB prediction API, preserves model/fragments,
writes through temporary files, validates non-empty downloads, calculates SHA-256
checksums and records all metadata required for reproducibility.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
from typing import Any

import pandas as pd
import requests

from rses_onco.structural import file_sha256, mean_plddt_from_pdb, plddt_class, write_json

ROOT = Path(__file__).resolve().parents[1]
AF_API = "https://alphafold.ebi.ac.uk/api/prediction"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def request(session: requests.Session, url: str, retries: int = 4) -> requests.Response:
  error: Exception | None = None
  for attempt in range(1, retries + 1):
    try:
      response = session.get(url, timeout=240)
      response.raise_for_status()
      return response
    except Exception as exc:
      error = exc
      if attempt == retries:
        break
      time.sleep(min(20, 2 ** attempt))
  raise RuntimeError(f"Request failed after {retries} attempts: {url}: {error}")


def download_file(
  session: requests.Session,
  url: str | None,
  destination: Path,
  refresh: bool,
) -> tuple[str, int, str] | None:
  if not url:
    return None
  if destination.exists() and destination.stat().st_size > 100 and not refresh:
    return str(destination), destination.stat().st_size, file_sha256(destination)
  destination.parent.mkdir(parents=True, exist_ok=True)
  part = destination.with_suffix(destination.suffix + ".part")
  response = request(session, url)
  part.write_bytes(response.content)
  if part.stat().st_size <= 100:
    part.unlink(missing_ok=True)
    raise RuntimeError(f"Downloaded file is unexpectedly small: {url}")
  part.replace(destination)
  return str(destination), destination.stat().st_size, file_sha256(destination)


def choose_records(payload: Any) -> list[dict[str, Any]]:
  if isinstance(payload, dict):
    payload = [payload]
  if not isinstance(payload, list):
    return []
  records = [record for record in payload if isinstance(record, dict)]
  return sorted(
    records,
    key=lambda record: (
      int(record.get("uniprotStart") or 1),
      int(record.get("uniprotEnd") or 10**9),
      str(record.get("entryId") or ""),
    ),
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--proteins",
    default="data/curated/human_nise_bonafide_2017.tsv",
  )
  parser.add_argument(
    "--output-dir",
    default="data/raw/structures/alphafold",
  )
  parser.add_argument(
    "--manifest",
    default="data/processed/structures/alphafold_structure_manifest.tsv",
  )
  parser.add_argument("--refresh", action="store_true")
  parser.add_argument("--sleep", type=float, default=0.10)
  args = parser.parse_args()

  proteins = pd.read_csv(resolve_path(args.proteins), sep="\t")
  required = {
    "group_id", "ec_number", "activity", "gene_symbol",
    "uniprot_accession", "structural_cluster",
  }
  missing = sorted(required - set(proteins.columns))
  if missing:
    raise ValueError(f"Protein table lacks columns: {missing}")
  proteins = proteins.drop_duplicates("uniprot_accession").sort_values(
    ["group_id", "structural_cluster", "gene_symbol"]
  )
  output_dir = resolve_path(args.output_dir)
  metadata_dir = output_dir / "metadata"
  session = requests.Session()
  session.headers.update({
    "User-Agent": "RSES-Onco/0.10 AlphaFold structure downloader",
    "Accept": "application/json",
  })

  rows: list[dict[str, Any]] = []
  failures: list[dict[str, str]] = []
  for index, protein in enumerate(proteins.to_dict("records"), start=1):
    accession = str(protein["uniprot_accession"])
    gene = str(protein["gene_symbol"])
    try:
      metadata_response = request(session, f"{AF_API}/{accession}")
      metadata = metadata_response.json()
      write_json(metadata_dir / f"{accession}.json", metadata)
      records = choose_records(metadata)
      if not records:
        raise RuntimeError("AlphaFold DB returned no model records")
      for record_index, record in enumerate(records, start=1):
        entry_id = str(record.get("entryId") or f"AF-{accession}-F{record_index}")
        safe_entry = entry_id.replace("/", "_")
        pdb_path = output_dir / "pdb" / f"{safe_entry}.pdb"
        cif_path = output_dir / "cif" / f"{safe_entry}.cif"
        pae_path = output_dir / "pae" / f"{safe_entry}.json"
        pdb_info = download_file(session, record.get("pdbUrl"), pdb_path, args.refresh)
        cif_info = download_file(session, record.get("cifUrl"), cif_path, args.refresh)
        pae_info = download_file(session, record.get("paeDocUrl"), pae_path, args.refresh)
        mean_plddt = mean_plddt_from_pdb(pdb_path) if pdb_info else None
        rows.append({
          **protein,
          "entry_id": entry_id,
          "fragment_index": record_index,
          "uniprot_start": record.get("uniprotStart"),
          "uniprot_end": record.get("uniprotEnd"),
          "model_created_date": record.get("modelCreatedDate"),
          "latest_version": record.get("latestVersion"),
          "pdb_path": pdb_info[0] if pdb_info else None,
          "pdb_size_bytes": pdb_info[1] if pdb_info else None,
          "pdb_sha256": pdb_info[2] if pdb_info else None,
          "cif_path": cif_info[0] if cif_info else None,
          "cif_size_bytes": cif_info[1] if cif_info else None,
          "cif_sha256": cif_info[2] if cif_info else None,
          "pae_path": pae_info[0] if pae_info else None,
          "pae_size_bytes": pae_info[1] if pae_info else None,
          "pae_sha256": pae_info[2] if pae_info else None,
          "mean_plddt": mean_plddt,
          "confidence_class": plddt_class(mean_plddt),
          "metadata_path": str(metadata_dir / f"{accession}.json"),
          "source_api": AF_API,
          "status": "ok",
        })
      print(f"[AlphaFold {index}/{len(proteins)}] {gene} {accession}: {len(records)} model record(s)", flush=True)
    except Exception as exc:
      failures.append({
        "gene_symbol": gene,
        "uniprot_accession": accession,
        "status": "failed",
        "message": str(exc),
      })
      print(f"[AlphaFold {index}/{len(proteins)}] FAILED {gene} {accession}: {exc}", flush=True)
    if args.sleep:
      time.sleep(args.sleep)

  manifest = resolve_path(args.manifest)
  manifest.parent.mkdir(parents=True, exist_ok=True)
  frame = pd.DataFrame(rows)
  if not frame.empty:
    frame = frame.sort_values(["group_id", "gene_symbol", "fragment_index"])
  frame.to_csv(manifest, sep="\t", index=False)
  pd.DataFrame(failures).to_csv(
    manifest.with_name("alphafold_structure_failures.tsv"), sep="\t", index=False
  )
  summary = {
    "curated_proteins": int(len(proteins)),
    "successful_proteins": int(frame["uniprot_accession"].nunique()) if not frame.empty else 0,
    "downloaded_model_records": int(len(frame)),
    "failed_proteins": int(len(failures)),
    "api": AF_API,
  }
  write_json(manifest.with_suffix(".summary.json"), summary)
  print(json.dumps(summary, indent=2, sort_keys=True))
  print(f"Wrote {manifest}")


if __name__ == "__main__":
  main()
