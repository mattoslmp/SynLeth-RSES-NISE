#!/usr/bin/env python3
"""Scan Ensembl promoter sequences with JASPAR CORE motifs using MEME FIMO.

The result is predicted motif occurrence only. It is never labelled as direct TF
binding, promoter occupancy, transcriptional activation or repression.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def file_sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--motifs", required=True)
  parser.add_argument("--promoters", required=True)
  parser.add_argument(
    "--output",
    default="data/processed/regulatory/jaspar_promoter_motif_hits.tsv",
  )
  parser.add_argument(
    "--summary-output",
    default="data/processed/regulatory/jaspar_promoter_tf_summary.tsv",
  )
  parser.add_argument(
    "--status-output",
    default="data/processed/regulatory/jaspar_promoter_scan_status.json",
  )
  parser.add_argument("--fimo", default="fimo")
  parser.add_argument("--threshold", type=float, default=1e-4)
  parser.add_argument("--q-value-threshold", type=float, default=0.05)
  parser.add_argument("--refresh", action="store_true")
  args = parser.parse_args()

  motifs = resolve_path(args.motifs)
  promoters = resolve_path(args.promoters)
  output = resolve_path(args.output)
  summary_output = resolve_path(args.summary_output)
  status_output = resolve_path(args.status_output)
  for path in (motifs, promoters):
    if not path.exists() or path.stat().st_size == 0:
      raise FileNotFoundError(path)
  executable = shutil.which(args.fimo)
  if executable is None:
    raise RuntimeError(
      "FIMO is required for promoter motif scanning. Install the MEME suite in "
      "the rses-onco environment, for example `conda install -c bioconda meme`."
    )

  if (
    not args.refresh
    and output.exists()
    and output.stat().st_size
    and summary_output.exists()
    and summary_output.stat().st_size
  ):
    print(f"Promoter motif scan: cache ({output})")
    return

  with tempfile.TemporaryDirectory(prefix="rses_fimo_") as directory:
    run_dir = Path(directory) / "fimo"
    command = [
      executable,
      "--oc",
      str(run_dir),
      "--thresh",
      str(args.threshold),
      "--max-stored-scores",
      "5000000",
      str(motifs),
      str(promoters),
    ]
    subprocess.run(command, cwd=ROOT, check=True)
    fimo_path = run_dir / "fimo.tsv"
    if not fimo_path.exists():
      raise FileNotFoundError(f"FIMO did not create {fimo_path}")
    frame = pd.read_csv(fimo_path, sep="\t", comment="#", low_memory=False)

  if frame.empty:
    hits = pd.DataFrame(columns=[
      "gene", "motif_id", "transcription_factor", "start", "stop",
      "strand", "score", "p_value", "q_value", "matched_sequence",
      "evidence_type", "direct_binding_claim",
    ])
  else:
    rename = {
      "motif_id": "motif_id",
      "motif_alt_id": "transcription_factor",
      "p-value": "p_value",
      "q-value": "q_value",
      "matched_sequence": "matched_sequence",
    }
    frame = frame.rename(columns=rename)
    frame["gene"] = (
      frame["sequence_name"].astype(str).str.split("|", regex=False).str[0]
    )
    frame["transcription_factor"] = (
      frame["transcription_factor"]
        .fillna("")
        .astype(str)
        .str.split("::", regex=False)
        .str[0]
        .str.strip()
        .str.upper()
    )
    frame["q_value"] = pd.to_numeric(frame["q_value"], errors="coerce")
    hits = frame.loc[
      frame["q_value"].notna()
      & frame["q_value"].le(args.q_value_threshold)
      & frame["transcription_factor"].ne("")
    ].copy()
    hits["evidence_type"] = "JASPAR_promoter_motif_prediction"
    hits["direct_binding_claim"] = False
    ordered = [
      "gene", "motif_id", "transcription_factor", "start", "stop",
      "strand", "score", "p_value", "q_value", "matched_sequence",
      "evidence_type", "direct_binding_claim",
    ]
    hits = hits[[column for column in ordered if column in hits]]
  atomic_tsv(hits, output)

  if hits.empty:
    summary = pd.DataFrame(columns=[
      "gene", "transcription_factor", "motif_count", "best_p_value",
      "best_q_value", "evidence_type", "direct_binding_claim",
    ])
  else:
    summary = (
      hits.groupby(["gene", "transcription_factor"], as_index=False)
        .agg(
          motif_count=("motif_id", "size"),
          best_p_value=("p_value", "min"),
          best_q_value=("q_value", "min"),
        )
    )
    summary["evidence_type"] = "JASPAR_promoter_motif_prediction"
    summary["direct_binding_claim"] = False
  atomic_tsv(summary, summary_output)

  status = {
    "source": "JASPAR 2026 CORE vertebrates + Ensembl canonical promoters",
    "method": "MEME FIMO",
    "fimo_executable": executable,
    "threshold": args.threshold,
    "q_value_threshold": args.q_value_threshold,
    "hit_rows": len(hits),
    "gene_tf_pairs": len(summary),
    "motif_sha256": file_sha256(motifs),
    "promoter_fasta_sha256": file_sha256(promoters),
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "interpretation": (
      "Motif occurrence is predicted cis-regulatory support and is not direct "
      "TF binding, promoter occupancy or causal regulation."
    ),
  }
  status_output.parent.mkdir(parents=True, exist_ok=True)
  temporary = status_output.with_suffix(status_output.suffix + ".tmp")
  temporary.write_text(json.dumps(status, indent=2, sort_keys=True), encoding="utf-8")
  temporary.replace(status_output)
  print(f"Promoter motif hits: {len(hits):,}")
  print(f"Unique promoter gene-TF pairs: {len(summary):,}")
  print(f"Wrote {output}")
  print(f"Wrote {summary_output}")


if __name__ == "__main__":
  main()
