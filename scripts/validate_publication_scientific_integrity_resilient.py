#!/usr/bin/env python3
"""Run scientific-integrity validation with registered-output freshness checks."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.validate_publication_scientific_integrity as target


def validate_run_freshness(article_root: Path, marker: Path | None) -> None:
  if marker is None:
    return
  if not marker.exists():
    raise FileNotFoundError(f"Run marker not found: {marker}")
  threshold = marker.stat().st_mtime
  mandatory: list[Path] = []

  figure_manifest = article_root / "manifests/figure_manifest.tsv"
  table_manifest = article_root / "manifests/table_manifest.tsv"
  mandatory.extend([figure_manifest, table_manifest])
  if figure_manifest.exists():
    figures = pd.read_csv(figure_manifest, sep="\t")
    for record in figures.to_dict("records"):
      base = Path(str(record["base_path"]))
      if not base.is_absolute():
        base = ROOT / base
      source = Path(str(record["source_data_path"]))
      if not source.is_absolute():
        source = ROOT / source
      mandatory.extend([source, base.with_suffix(".png"), base.with_suffix(".pdf"), base.with_suffix(".svg"), base.with_suffix(".layout_audit.json")])
  if table_manifest.exists():
    tables = pd.read_csv(table_manifest, sep="\t")
    for value in tables["path"].astype(str):
      path = Path(value)
      mandatory.append(path if path.is_absolute() else ROOT / path)

  mandatory.extend([
    article_root / "tables/qc/candidate_domain_evidence_audit.tsv",
    article_root / "tables/qc/evidence_overlap_registry.tsv",
    article_root / "tables/score_components/rses_onco_score_decomposition.tsv",
    article_root / "tables/robustness/leave_one_domain_out.tsv",
    article_root / "tables/figure_data/figure_source_data_inventory.tsv",
    article_root / "tables/supporting_evidence/supporting_evidence_manifest.tsv",
  ])
  stale = []
  for path in mandatory:
    if not path.exists() or path.stat().st_size == 0:
      stale.append(f"missing_or_empty:{path}")
    elif path.stat().st_mtime + 1e-6 < threshold:
      stale.append(f"stale:{path}")
  if stale:
    raise RuntimeError(
      "Mandatory registered publication outputs were not regenerated in this assets-only run:\n"
      + "\n".join(stale[:150])
    )


target.validate_run_freshness = validate_run_freshness


if __name__ == "__main__":
  target.main()
