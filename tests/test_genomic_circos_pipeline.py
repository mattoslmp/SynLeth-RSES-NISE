from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import pandas as pd
import yaml

from rses_onco.circos import GenomeLayout, normalize_chromosome

ROOT = Path(__file__).resolve().parents[1]


def test_genome_layout_maps_grch38_positions() -> None:
  layout = GenomeLayout()
  assert normalize_chromosome("chr1") == "1"
  assert normalize_chromosome("23") == "X"
  start = layout.theta("1", 1)
  end = layout.theta("1", layout.lengths["1"])
  assert start != end
  x, y = layout.point(start, 1.0)
  assert abs((x * x + y * y) - 1.0) < 1e-8


def test_genomic_circos_figure_generates_strict_triplet(
  tmp_path: Path,
) -> None:
  processed = tmp_path / "data/processed/circos"
  processed.mkdir(parents=True)
  article = tmp_path / "article_outputs"
  config = tmp_path / "article_assets.yaml"
  config.write_text(
    yaml.safe_dump({
      "supplementary_figures": [{
        "id": "Figure_S70",
        "file": "Figure_S70_genomic_circos_rses_onco",
        "title": "Genomic Circos of NISE and paralog hypotheses",
        "caption": "Synthetic strict-layout Circos test.",
      }],
    }),
    encoding="utf-8",
  )

  coordinates = pd.DataFrame([
    {
      "gene": "GENEA",
      "gene_class": "NISE",
      "chromosome": "1",
      "genomic_position": 10_000_000,
    },
    {
      "gene": "GENEB",
      "gene_class": "homologous_paralog",
      "chromosome": "2",
      "genomic_position": 20_000_000,
    },
    {
      "gene": "GENEC",
      "gene_class": "NISE_and_paralog",
      "chromosome": "X",
      "genomic_position": 30_000_000,
    },
  ])
  links = pd.DataFrame([
    {
      "pair_id": "P1",
      "lost_gene": "GENEA",
      "target_gene": "GENEB",
      "pair_class": "NISE",
      "lost_chromosome": "1",
      "lost_position": 10_000_000,
      "target_chromosome": "2",
      "target_position": 20_000_000,
      "maximum_coverage_adjusted_rses": 0.8,
      "link_width": 2.0,
      "link_alpha": 0.4,
      "link_color": "#C62828",
    },
    {
      "pair_id": "P2",
      "lost_gene": "GENEB",
      "target_gene": "GENEC",
      "pair_class": "homologous_paralog",
      "lost_chromosome": "2",
      "lost_position": 20_000_000,
      "target_chromosome": "X",
      "target_position": 30_000_000,
      "maximum_coverage_adjusted_rses": 0.6,
      "link_width": 1.5,
      "link_alpha": 0.3,
      "link_color": "#111111",
    },
  ])
  tracks = pd.DataFrame([
    {
      "track_id": "A01",
      "track_label": "Coverage-adjusted RSES",
      "source_column": "coverage_adjusted_rses",
      "panel": "A",
      "domain_family": "RSES-Onco",
      "parent_domain": "global",
      "ring_order": 1,
      "aggregation": "maximum",
    },
    {
      "track_id": "A02",
      "track_label": "Evidence coverage",
      "source_column": "evidence_coverage",
      "panel": "A",
      "domain_family": "RSES-Onco",
      "parent_domain": "coverage",
      "ring_order": 2,
      "aggregation": "maximum",
    },
    {
      "track_id": "B01",
      "track_label": "Expression context",
      "source_column": "microniche_expression_context",
      "panel": "B",
      "domain_family": "Functional microniche",
      "parent_domain": "expression_context",
      "ring_order": 1,
      "aggregation": "maximum",
    },
    {
      "track_id": "B02",
      "track_label": "Promoter methylation",
      "source_column": "component_promoter_methylation_context",
      "panel": "B",
      "domain_family": "Regulatory network",
      "parent_domain": "promoter_methylation",
      "ring_order": 2,
      "aggregation": "maximum",
    },
  ])
  ring_rows = []
  for gene in coordinates["gene"]:
    for index, track in tracks.iterrows():
      ring_rows.append({
        "gene": gene,
        **track.to_dict(),
        "value": 0.2 + 0.1 * index,
      })
  rings = pd.DataFrame(ring_rows)
  paths = {
    "coordinates": processed / "genomic_circos_gene_coordinates.tsv",
    "links": processed / "genomic_circos_pair_links.tsv",
    "rings": processed / "genomic_circos_ring_values.tsv",
    "tracks": processed / "genomic_circos_track_definitions.tsv",
  }
  coordinates.to_csv(paths["coordinates"], sep="\t", index=False)
  links.to_csv(paths["links"], sep="\t", index=False)
  rings.to_csv(paths["rings"], sep="\t", index=False)
  tracks.to_csv(paths["tracks"], sep="\t", index=False)

  subprocess.run([
    sys.executable,
    "scripts/make_genomic_circos_figure.py",
    "--config",
    str(config),
    "--coordinates",
    str(paths["coordinates"]),
    "--links",
    str(paths["links"]),
    "--ring-values",
    str(paths["rings"]),
    "--tracks",
    str(paths["tracks"]),
    "--output-root",
    str(article),
    "--strict-layout",
  ], cwd=ROOT, check=True)

  base = (
    article
    / "figures/supplementary/Figure_S70_genomic_circos_rses_onco"
  )
  for extension in ("png", "pdf", "svg", "layout_audit.json"):
    path = base.with_suffix(f".{extension}")
    assert path.exists() and path.stat().st_size > 0
  manifest = pd.read_csv(
    article / "manifests/genomic_circos_figure_manifest.tsv",
    sep="\t",
  )
  assert manifest.iloc[0]["figure_id"] == "Figure_S70"
  assert manifest.iloc[0]["layout_status"] == "pass"
