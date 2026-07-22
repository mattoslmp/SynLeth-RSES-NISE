#!/usr/bin/env python3
"""Download canonical human promoter coordinates and sequences from Ensembl REST.

Promoter sequence is defined relative to the canonical transcript TSS as 2,000 bp
upstream and 500 bp downstream on the transcript strand. Coordinates and sequences
are cached and source failures remain explicit; they are never converted into
negative regulatory evidence.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from typing import Any

import pandas as pd
import requests

from rses_onco.utils import canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
SERVER = "https://rest.ensembl.org"


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def atomic_text(text: str, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  temporary.write_text(text, encoding="utf-8")
  temporary.replace(path)


def request_json(
  session: requests.Session,
  method: str,
  url: str,
  *,
  payload: dict[str, Any] | None = None,
  retries: int = 5,
) -> Any:
  delay = 1.0
  for attempt in range(retries):
    try:
      response = session.request(
        method,
        url,
        json=payload,
        headers={
          "Content-Type": "application/json",
          "Accept": "application/json",
        },
        timeout=90,
      )
      if response.ok:
        return response.json()
      if (
        response.status_code not in {408, 425, 429}
        and response.status_code < 500
      ):
        response.raise_for_status()
    except requests.RequestException:
      if attempt + 1 >= retries:
        raise
    time.sleep(delay)
    delay = min(delay * 2, 30)
  raise RuntimeError(
    f"Ensembl request failed after {retries} attempts: {url}"
  )


def request_sequence(
  session: requests.Session,
  region: str,
  retries: int = 5,
) -> str:
  url = f"{SERVER}/sequence/region/homo_sapiens/{region}"
  delay = 1.0
  for attempt in range(retries):
    try:
      response = session.get(
        url,
        headers={
          "Content-Type": "text/plain",
          "Accept": "text/plain",
        },
        timeout=90,
      )
      if response.ok:
        return response.text.strip().upper()
      if (
        response.status_code not in {408, 425, 429}
        and response.status_code < 500
      ):
        response.raise_for_status()
    except requests.RequestException:
      if attempt + 1 >= retries:
        raise
    time.sleep(delay)
    delay = min(delay * 2, 30)
  raise RuntimeError(f"Ensembl sequence request failed: {region}")


def promoter_coordinates(
  record: dict[str, Any],
  upstream: int,
  downstream: int,
) -> dict[str, Any]:
  transcripts = record.get("Transcript") or []
  canonical_id = str(
    record.get("canonical_transcript") or ""
  ).split(".")[0]
  canonical = next(
    (
      item
      for item in transcripts
      if str(item.get("id")) == canonical_id
      or int(item.get("is_canonical") or 0) == 1
    ),
    None,
  )
  if canonical is None:
    canonical = next(
      (
        item
        for item in transcripts
        if item.get("biotype") == "protein_coding"
      ),
      None,
    )
  if canonical is None:
    raise ValueError("no_canonical_or_protein_coding_transcript")
  strand = int(canonical["strand"])
  chromosome = str(canonical["seq_region_name"])
  tss = int(
    canonical["start"]
    if strand == 1
    else canonical["end"]
  )
  if strand == 1:
    start = max(1, tss - upstream)
    end = tss + downstream
  else:
    start = max(1, tss - downstream)
    end = tss + upstream
  return {
    "ensembl_gene_id": record.get("id"),
    "canonical_transcript_id": canonical.get("id"),
    "assembly": (
      canonical.get("assembly_name")
      or record.get("assembly_name")
    ),
    "chromosome": chromosome,
    "strand": strand,
    "tss": tss,
    "promoter_start": start,
    "promoter_end": end,
    "region": f"{chromosome}:{start}..{end}:{strand}",
  }


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--output",
    default="data/raw/regulatory/ensembl_promoters.tsv",
  )
  parser.add_argument(
    "--fasta",
    default="data/raw/regulatory/ensembl_promoters.fa",
  )
  parser.add_argument(
    "--cache-dir",
    default="data/raw/regulatory/ensembl_promoter_cache",
  )
  parser.add_argument(
    "--status-output",
    default="data/raw/regulatory/ensembl_promoters_status.json",
  )
  parser.add_argument("--upstream", type=int, default=2000)
  parser.add_argument("--downstream", type=int, default=500)
  args = parser.parse_args()

  candidates = pd.read_csv(
    resolve_path(args.candidates),
    sep="\t",
    low_memory=False,
  )
  genes = sorted({
    canonical_gene_name(value)
    for column in ("lost_gene", "target_gene")
    if column in candidates
    for value in candidates[column].dropna()
    if canonical_gene_name(value)
  })
  output = resolve_path(args.output)
  fasta = resolve_path(args.fasta)
  cache_dir = resolve_path(args.cache_dir)
  status_output = resolve_path(args.status_output)
  cache_dir.mkdir(parents=True, exist_ok=True)

  session = requests.Session()
  lookup: dict[str, Any] = {}
  for offset in range(0, len(genes), 200):
    batch = genes[offset:offset + 200]
    uncached = []
    for gene in batch:
      cache = cache_dir / f"{gene}.lookup.json"
      if cache.exists() and cache.stat().st_size:
        lookup[gene] = json.loads(
          cache.read_text(encoding="utf-8")
        )
      else:
        uncached.append(gene)
    if uncached:
      response = request_json(
        session,
        "POST",
        f"{SERVER}/lookup/symbol/homo_sapiens?expand=1",
        payload={"symbols": uncached},
      )
      for gene in uncached:
        value = response.get(gene)
        lookup[gene] = value
        atomic_text(
          json.dumps(value, indent=2, sort_keys=True),
          cache_dir / f"{gene}.lookup.json",
        )
      time.sleep(0.2)

  rows: list[dict[str, Any]] = []
  fasta_lines: list[str] = []
  for gene in genes:
    record = lookup.get(gene)
    base = {
      "gene": gene,
      "source": "Ensembl REST",
      "source_url": (
        f"{SERVER}/lookup/symbol/homo_sapiens/{gene}?expand=1"
      ),
      "accessed_at_utc": datetime.now(timezone.utc).isoformat(),
      "promoter_upstream_bp": args.upstream,
      "promoter_downstream_bp": args.downstream,
    }
    if not isinstance(record, dict):
      rows.append({**base, "status": "symbol_unresolved"})
      continue
    try:
      coordinates = promoter_coordinates(
        record,
        args.upstream,
        args.downstream,
      )
      sequence_cache = cache_dir / f"{gene}.promoter.fa"
      if sequence_cache.exists() and sequence_cache.stat().st_size:
        sequence = "".join(
          line.strip()
          for line in sequence_cache.read_text(
            encoding="utf-8"
          ).splitlines()
          if not line.startswith(">")
        ).upper()
      else:
        sequence = request_sequence(
          session,
          coordinates["region"],
        )
        atomic_text(
          f">{gene}|{coordinates['region']}\n{sequence}\n",
          sequence_cache,
        )
        time.sleep(0.1)
      digest = hashlib.sha256(
        sequence.encode("utf-8")
      ).hexdigest()
      rows.append({
        **base,
        **coordinates,
        "sequence_length": len(sequence),
        "sequence_sha256": digest,
        "status": "available",
      })
      fasta_lines.extend([
        f">{gene}|{coordinates['region']}",
        sequence,
      ])
    except Exception as error:
      rows.append({
        **base,
        "status": "technical_failure",
        "failure_reason": str(error),
      })

  frame = pd.DataFrame(rows)
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)
  atomic_text("\n".join(fasta_lines) + "\n", fasta)
  available = frame["status"].eq("available")
  status = {
    "source": "Ensembl REST",
    "lookup_endpoint": (
      f"{SERVER}/lookup/symbol/homo_sapiens"
    ),
    "sequence_endpoint": (
      f"{SERVER}/sequence/region/homo_sapiens"
    ),
    "candidate_gene_count": len(genes),
    "available_promoters": int(available.sum()),
    "failed_or_unresolved": int((~available).sum()),
    "output": str(output),
    "fasta": str(fasta),
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
  }
  atomic_text(
    json.dumps(status, indent=2, sort_keys=True),
    status_output,
  )
  print(frame["status"].value_counts(dropna=False).to_string())
  print(f"Wrote {output}")
  print(f"Wrote {fasta}")


if __name__ == "__main__":
  main()
