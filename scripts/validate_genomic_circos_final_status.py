#!/usr/bin/env python3
"""Validate finalized genomic Circos counts and SHA-256 provenance."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--circos-dir", default="data/processed/circos")
  args = parser.parse_args()
  directory = resolve(args.circos_dir)
  status_path = directory / "genomic_circos_status.json"
  provenance_path = directory / "genomic_circos_source_provenance.tsv"
  if not status_path.exists() or status_path.stat().st_size == 0:
    raise FileNotFoundError(status_path)
  if not provenance_path.exists() or provenance_path.stat().st_size == 0:
    raise FileNotFoundError(provenance_path)

  status = json.loads(status_path.read_text(encoding="utf-8"))
  provenance = pd.read_csv(
    provenance_path,
    sep="\t",
    low_memory=False,
  )
  errors: list[str] = []
  expected_values = {
    "version": "genomic-circos-v0.11.2",
    "assembly": "GRCh38",
    "coordinate_missing": 0,
    "duplicate_pair_chords": 0,
    "tracks": 35,
    "panel_a_tracks": 14,
    "panel_b_tracks": 21,
  }
  for key, expected in expected_values.items():
    observed = status.get(key)
    if observed != expected:
      errors.append(f"{key}={observed!r}, expected {expected!r}")
  if int(status.get("candidate_pair_chords", 0)) <= 0:
    errors.append("candidate_pair_chords must be positive")
  if int(status.get("genes", 0)) <= 0:
    errors.append("genes must be positive")
  if int(status.get("ring_rows", 0)) != int(status.get("genes", 0)) * 35:
    errors.append("ring_rows must equal genes multiplied by 35 tracks")
  if int(status.get("expression_gene_cancer_rows", 0)) != int(
    status.get("genes", 0)
  ) * 3:
    errors.append(
      "expression_gene_cancer_rows must equal genes multiplied by 3 cancers"
    )
  required_roles = {
    "final_ranking",
    "candidate_universe",
    "ensembl_canonical_coordinates",
    "depmap_expression",
    "depmap_model_metadata",
    "wgcna_pair_metrics",
    "gene_coordinates",
    "pair_links",
    "ring_values",
    "track_definitions",
    "expression_summary",
    "expression_model_values",
  }
  observed_roles = set(provenance["role"].astype(str))
  if required_roles != observed_roles:
    errors.append(
      "provenance roles mismatch; missing="
      f"{sorted(required_roles - observed_roles)}, extra="
      f"{sorted(observed_roles - required_roles)}"
    )
  if not provenance["sha256"].astype(str).str.fullmatch(
    r"[0-9a-f]{64}"
  ).all():
    errors.append("one or more provenance SHA-256 values are invalid")
  if not provenance["bytes"].gt(0).all():
    errors.append("one or more provenance files have non-positive byte size")

  if errors:
    raise SystemExit(
      "Genomic Circos final-status validation failed:\n"
      + "\n".join(f"- {error}" for error in errors)
    )
  print("Genomic Circos final-status and provenance validation passed.")
  print(json.dumps(status, indent=2))


if __name__ == "__main__":
  main()
