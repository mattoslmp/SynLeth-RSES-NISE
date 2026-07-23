from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import pandas as pd
import yaml

from scripts.build_genomic_circos_inputs import TRACKS
from scripts.enrich_genomic_circos_internal_layers import (
  ADDITIONAL_TRACKS,
)

ROOT = Path(__file__).resolve().parents[1]


def test_all_circos_score_rings_pass_strict_layout(
  tmp_path: Path,
) -> None:
  data = tmp_path / "circos"
  data.mkdir()
  article = tmp_path / "article"
  config = tmp_path / "config.yaml"
  config.write_text(
    yaml.safe_dump({
      "supplementary_figures": [{
        "id": "Figure_S70",
        "file": "Figure_S70_genomic_circos_rses_onco",
        "title": "Circos",
        "caption": "All-ring layout regression test.",
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
  links = pd.DataFrame([{
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
  }])
  tracks = pd.DataFrame(
    [*TRACKS, *ADDITIONAL_TRACKS],
    columns=[
      "track_id",
      "track_label",
      "source_column",
      "panel",
      "domain_family",
      "parent_domain",
      "ring_order",
    ],
  ).sort_values(["panel", "ring_order"])
  assert len(tracks) == 35
  tracks["aggregation"] = "maximum"
  rings = pd.DataFrame([
    {
      "gene": gene,
      **track,
      "value": 0.5,
    }
    for gene in coordinates["gene"]
    for track in tracks.to_dict("records")
  ])
  paths = {
    "coordinates": data / "coordinates.tsv",
    "links": data / "links.tsv",
    "rings": data / "rings.tsv",
    "tracks": data / "tracks.tsv",
  }
  coordinates.to_csv(paths["coordinates"], sep="\t", index=False)
  links.to_csv(paths["links"], sep="\t", index=False)
  rings.to_csv(paths["rings"], sep="\t", index=False)
  tracks.to_csv(paths["tracks"], sep="\t", index=False)

  subprocess.run([
    sys.executable,
    "scripts/make_genomic_circos_figure_resilient.py",
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

  audit = (
    article
    / "figures/supplementary/"
    "Figure_S70_genomic_circos_rses_onco.layout_audit.json"
  )
  assert audit.exists()
  assert '"status": "pass"' in audit.read_text(encoding="utf-8")
  assert '"warnings": []' in audit.read_text(encoding="utf-8")
