from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  frame.to_csv(path, sep="\t", index=False)


def test_finalize_and_validate_circos_status(
  tmp_path: Path,
) -> None:
  output = tmp_path / "circos"
  coordinates = pd.DataFrame([
    {
      "gene": "A",
      "coordinate_status": "available",
    },
    {
      "gene": "B",
      "coordinate_status": "available",
    },
  ])
  links = pd.DataFrame([
    {
      "pair_id": "P1",
      "pair_class": "NISE",
      "link_status": "available",
    },
  ])
  tracks = pd.DataFrame([
    {
      "track_id": f"A{index:02d}" if index < 14 else f"B{index - 13:02d}",
      "panel": "A" if index < 14 else "B",
    }
    for index in range(35)
  ])
  rings = pd.DataFrame([
    {
      "gene": gene,
      "track_id": track,
      "evidence_status": "observed",
    }
    for gene in coordinates["gene"]
    for track in tracks["track_id"]
  ])
  expression_summary = pd.DataFrame([
    {"cancer": cancer, "gene": gene}
    for cancer in ("colon", "stomach", "lung")
    for gene in coordinates["gene"]
  ])
  expression_model = pd.DataFrame([
    {
      "ModelID": "ACH-1",
      "gene": "A",
      "is_measurement": True,
    },
    {
      "ModelID": "UNAVAILABLE::colon::B",
      "gene": "B",
      "is_measurement": False,
    },
  ])
  generated = {
    "genomic_circos_gene_coordinates.tsv": coordinates,
    "genomic_circos_pair_links.tsv": links,
    "genomic_circos_ring_values.tsv": rings,
    "genomic_circos_track_definitions.tsv": tracks,
    "genomic_circos_expression_summary.tsv": expression_summary,
    "genomic_circos_expression_model_values.tsv": expression_model,
  }
  for name, frame in generated.items():
    write_tsv(output / name, frame)

  source_paths = {}
  for name in (
    "ranking.tsv",
    "candidates.tsv",
    "promoters.tsv",
    "expression.csv",
    "models.csv",
    "wgcna.tsv",
  ):
    path = tmp_path / name
    path.write_text("column\nvalue\n", encoding="utf-8")
    source_paths[name] = path

  subprocess.run([
    sys.executable,
    "scripts/finalize_genomic_circos_status.py",
    "--output-dir",
    str(output),
    "--ranking",
    str(source_paths["ranking.tsv"]),
    "--candidates",
    str(source_paths["candidates.tsv"]),
    "--promoters",
    str(source_paths["promoters.tsv"]),
    "--expression",
    str(source_paths["expression.csv"]),
    "--models",
    str(source_paths["models.csv"]),
    "--wgcna",
    str(source_paths["wgcna.tsv"]),
  ], cwd=ROOT, check=True)
  subprocess.run([
    sys.executable,
    "scripts/validate_genomic_circos_final_status.py",
    "--circos-dir",
    str(output),
  ], cwd=ROOT, check=True)

  status = json.loads(
    (output / "genomic_circos_status.json").read_text(
      encoding="utf-8"
    )
  )
  assert status["tracks"] == 35
  assert status["panel_a_tracks"] == 14
  assert status["panel_b_tracks"] == 21
  assert status["unavailable_expression_sentinels"] == 1
  provenance = pd.read_csv(
    output / "genomic_circos_source_provenance.tsv",
    sep="\t",
  )
  assert "wgcna_pair_metrics" in set(provenance["role"])
