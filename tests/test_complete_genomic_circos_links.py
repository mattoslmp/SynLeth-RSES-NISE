from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_missing_score_pair_still_receives_a_circos_chord(
  tmp_path: Path,
) -> None:
  candidates = tmp_path / "candidates.tsv"
  coordinates = tmp_path / "coordinates.tsv"
  ranking = tmp_path / "ranking.tsv"
  links = tmp_path / "links.tsv"

  pd.DataFrame([
    {
      "pair_id": "P1",
      "lost_gene": "GENEA",
      "target_gene": "GENEB",
      "source_class": "NISE",
      "relation_type": "nise",
      "ensembl_homology_type": "",
    },
    {
      "pair_id": "P2",
      "lost_gene": "GENEC",
      "target_gene": "GENED",
      "source_class": "PARALOG",
      "relation_type": "paralog",
      "ensembl_homology_type": "within_species_paralog",
    },
  ]).to_csv(candidates, sep="\t", index=False)
  pd.DataFrame([
    {
      "gene": "GENEA",
      "chromosome": "1",
      "genomic_position": 10,
    },
    {
      "gene": "GENEB",
      "chromosome": "2",
      "genomic_position": 20,
    },
    {
      "gene": "GENEC",
      "chromosome": "3",
      "genomic_position": 30,
    },
    {
      "gene": "GENED",
      "chromosome": "4",
      "genomic_position": 40,
    },
  ]).to_csv(coordinates, sep="\t", index=False)
  pd.DataFrame([
    {
      "pair_id": "P1",
      "cancer": "colon",
      "coverage_adjusted_rses": 0.8,
      "evidence_coverage": 0.9,
    },
  ]).to_csv(ranking, sep="\t", index=False)
  pd.DataFrame([
    {
      "pair_id": "P1",
      "lost_gene": "GENEA",
      "target_gene": "GENEB",
      "pair_class": "NISE",
      "lost_chromosome": "1",
      "lost_position": 10,
      "target_chromosome": "2",
      "target_position": 20,
      "cancers": "colon",
      "maximum_coverage_adjusted_rses": 0.8,
      "median_coverage_adjusted_rses": 0.8,
      "maximum_evidence_coverage": 0.9,
      "link_width": 2.45,
      "link_alpha": 0.412,
      "link_color": "#C62828",
      "link_status": "available",
    },
  ]).to_csv(links, sep="\t", index=False)

  subprocess.run([
    sys.executable,
    "scripts/complete_genomic_circos_links.py",
    "--candidates",
    str(candidates),
    "--coordinates",
    str(coordinates),
    "--ranking",
    str(ranking),
    "--links",
    str(links),
  ], cwd=ROOT, check=True)

  result = pd.read_csv(links, sep="\t", low_memory=False)
  assert set(result["pair_id"]) == {"P1", "P2"}
  p2 = result.set_index("pair_id").loc["P2"]
  assert p2["link_status"] == "score_missing"
  assert p2["link_color"] == "#111111"
  assert pd.isna(p2["maximum_coverage_adjusted_rses"])
