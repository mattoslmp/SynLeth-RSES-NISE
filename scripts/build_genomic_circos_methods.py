#!/usr/bin/env python3
"""Generate supplementary methods for the genomic Circos and its score rings."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--tracks",
    default=(
      "data/processed/circos/"
      "genomic_circos_track_definitions.tsv"
    ),
  )
  parser.add_argument(
    "--output",
    default=(
      "article_outputs/manuscript_assets/"
      "supplementary_methods/GENOMIC_CIRCOS_METHODS.md"
    ),
  )
  args = parser.parse_args()

  tracks_path = resolve(args.tracks)
  if not tracks_path.exists() or tracks_path.stat().st_size == 0:
    raise FileNotFoundError(tracks_path)
  tracks = pd.read_csv(tracks_path, sep="\t", low_memory=False)
  lines = [
    "# Genomic Circos representation of NISE and paralog hypotheses",
    "",
    (
      "The genomic Circos figure includes every coordinate-complete "
      "simple-gene hypothesis classified as a non-homologous isofunctional "
      "enzyme (NISE), homologous paralog or both. Canonical human genomic "
      "coordinates are derived from the Ensembl canonical-transcript lookup "
      "used by the promoter workflow and are displayed on the GRCh38 "
      "chromosome ideogram. Every included gene is shown as a genomic tick. "
      "Red chords denote NISE relationships and black chords denote "
      "homologous-paralog relationships. Link width and transparency are "
      "proportional to the maximum cancer-specific coverage-adjusted "
      "RSES-Onco score for the pair."
    ),
    "",
    (
      "For every gene, each ring value is the maximum observed value across "
      "all associated pair × cancer rows. The exact median, minimum, number "
      "of observed rows and number of eligible rows are retained in "
      "Supplementary Table S47. Missing or non-eligible evidence remains NA "
      "and is rendered as a hollow marker rather than numeric zero."
    ),
    "",
    "## Ring definitions",
    "",
    (
      "| Ring | Panel | Domain family | Parent domain | Source column | "
      "Aggregation |"
    ),
    "|---|---|---|---|---|---|",
  ]
  for row in tracks.sort_values([
    "panel",
    "ring_order",
  ]).to_dict("records"):
    lines.append(
      f"| {row['track_id']} — {row['track_label']} | "
      f"{row['panel']} | {row['domain_family']} | "
      f"{row['parent_domain']} | `{row['source_column']}` | "
      f"{row['aggregation']} |"
    )
  lines.extend([
    "",
    "## Expression data transparency",
    "",
    (
      "All model-level DepMap expression values used for Circos genes are "
      "exported as Supplementary Table S50, retaining ModelID, cancer "
      "context, gene, log2(TPM+1) value and source file. Supplementary "
      "Table S49 provides the corresponding per-gene/per-cancer summary. "
      "These tables are not reconstructed from the figure and are generated "
      "directly from the same expression matrix used by the analysis."
    ),
    "",
    "## Reproducibility",
    "",
    (
      "Supplementary Tables S45–S52 contain genomic coordinates, pair "
      "links, ring values, ring definitions, expression summaries, "
      "model-level expression values, the complete script catalogue and "
      "source provenance. The exact combined source TSV used to render "
      "Figure S70 is copied byte-for-byte to the figure-data directory and "
      "registered in the publication manifest."
    ),
  ])
  output = resolve(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  output.write_text("\n".join(lines) + "\n", encoding="utf-8")
  print(f"Wrote genomic Circos supplementary methods: {output}")


if __name__ == "__main__":
  main()
