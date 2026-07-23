#!/usr/bin/env python3
"""Finalize genomic Circos status and source provenance after all enrichment stages."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def require(path: Path) -> None:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(path)


def sha256(path: Path) -> str:
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
  parser.add_argument("--output-dir", default="data/processed/circos")
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--promoters",
    default="data/raw/regulatory/ensembl_promoters.tsv",
  )
  parser.add_argument(
    "--expression",
    default=(
      "data/raw/depmap/"
      "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    ),
  )
  parser.add_argument("--models", default="data/raw/depmap/Model.csv")
  parser.add_argument(
    "--wgcna",
    default=(
      "data/processed/regulatory/wgcna/"
      "wgcna_pair_metrics_all_cancers.tsv"
    ),
  )
  args = parser.parse_args()

  output_dir = resolve(args.output_dir)
  generated_paths = {
    "gene_coordinates": output_dir / "genomic_circos_gene_coordinates.tsv",
    "pair_links": output_dir / "genomic_circos_pair_links.tsv",
    "ring_values": output_dir / "genomic_circos_ring_values.tsv",
    "track_definitions": output_dir / "genomic_circos_track_definitions.tsv",
    "expression_summary": output_dir / "genomic_circos_expression_summary.tsv",
    "expression_model_values": output_dir / "genomic_circos_expression_model_values.tsv",
  }
  source_paths = {
    "final_ranking": resolve(args.ranking),
    "candidate_universe": resolve(args.candidates),
    "ensembl_canonical_coordinates": resolve(args.promoters),
    "depmap_expression": resolve(args.expression),
    "depmap_model_metadata": resolve(args.models),
    "wgcna_pair_metrics": resolve(args.wgcna),
  }
  for path in [*generated_paths.values(), *source_paths.values()]:
    require(path)

  coordinates = pd.read_csv(
    generated_paths["gene_coordinates"], sep="\t", low_memory=False
  )
  links = pd.read_csv(
    generated_paths["pair_links"], sep="\t", low_memory=False
  )
  rings = pd.read_csv(
    generated_paths["ring_values"], sep="\t", low_memory=False
  )
  tracks = pd.read_csv(
    generated_paths["track_definitions"], sep="\t", low_memory=False
  )
  expression_summary = pd.read_csv(
    generated_paths["expression_summary"], sep="\t", low_memory=False
  )
  expression_model = pd.read_csv(
    generated_paths["expression_model_values"], sep="\t", low_memory=False
  )

  if len(tracks) != 35:
    raise RuntimeError(f"Expected 35 final Circos tracks; observed {len(tracks)}")
  if links["pair_id"].astype(str).duplicated().any():
    raise RuntimeError("Final Circos link table contains duplicate pair IDs")
  if coordinates["gene"].astype(str).duplicated().any():
    raise RuntimeError("Final Circos coordinate table contains duplicate genes")

  provenance_rows = []
  for role, path in {**source_paths, **generated_paths}.items():
    provenance_rows.append({
      "role": role,
      "path": str(path),
      "bytes": int(path.stat().st_size),
      "sha256": sha256(path),
      "source_or_generated": (
        "source" if role in source_paths else "generated"
      ),
      "used_by": (
        "scripts/build_genomic_circos_inputs.py;"
        "scripts/enrich_genomic_circos_internal_layers.py;"
        "scripts/complete_genomic_circos_expression_summary.py;"
        "scripts/complete_genomic_circos_links.py;"
        "scripts/make_genomic_circos_figure_resilient.py"
      ),
      "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
    })
  provenance_path = output_dir / "genomic_circos_source_provenance.tsv"
  atomic_tsv(pd.DataFrame(provenance_rows), provenance_path)

  measurement = (
    expression_model.get(
      "is_measurement",
      pd.Series(True, index=expression_model.index),
    )
    .astype(str)
    .str.casefold()
    .isin({"true", "1", "yes"})
  )
  status = {
    "version": "genomic-circos-v0.11.2",
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "assembly": "GRCh38",
    "genes": int(coordinates["gene"].nunique()),
    "coordinate_missing": int(
      (~coordinates["coordinate_status"].astype(str).eq("available")).sum()
    ),
    "candidate_pair_chords": int(links["pair_id"].nunique()),
    "duplicate_pair_chords": int(
      links["pair_id"].astype(str).duplicated().sum()
    ),
    "score_missing_chords": int(
      links["link_status"].astype(str).eq("score_missing").sum()
    ),
    "nise_chords": int(
      links["pair_class"].astype(str).str.contains("NISE").sum()
    ),
    "homologous_paralog_chords": int(
      links["pair_class"].astype(str).eq("homologous_paralog").sum()
    ),
    "tracks": int(len(tracks)),
    "panel_a_tracks": int(tracks["panel"].astype(str).eq("A").sum()),
    "panel_b_tracks": int(tracks["panel"].astype(str).eq("B").sum()),
    "ring_rows": int(len(rings)),
    "observed_ring_rows": int(
      rings["evidence_status"].astype(str).eq("observed").sum()
    ),
    "missing_or_noneligible_ring_rows": int(
      rings["evidence_status"]
      .astype(str)
      .eq("missing_or_not_eligible")
      .sum()
    ),
    "expression_gene_cancer_rows": int(len(expression_summary)),
    "observed_expression_measurements": int(measurement.sum()),
    "unavailable_expression_sentinels": int((~measurement).sum()),
    "missing_data_rule": (
      "missing and non-eligible evidence remain NA; no zero imputation"
    ),
    "link_rule": (
      "exactly one chord per simple NISE/paralog candidate pair"
    ),
    "outputs": {
      **{key: str(path) for key, path in generated_paths.items()},
      "source_provenance": str(provenance_path),
    },
  }
  status_path = output_dir / "genomic_circos_status.json"
  status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
  print(json.dumps(status, indent=2))


if __name__ == "__main__":
  main()
