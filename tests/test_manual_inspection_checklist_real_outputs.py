from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_manual_inspection_checklist_when_real_outputs_exist() -> None:
  figure_manifest = ROOT / "article_outputs/manifests/figure_manifest.tsv"
  if not figure_manifest.exists():
    pytest.skip("Real publication figures have not been generated")
  subprocess.run(
    [
      sys.executable,
      "-u",
      "scripts/create_manual_visual_inspection_checklist.py",
      "--article-root",
      "article_outputs",
    ],
    cwd=ROOT,
    check=True,
  )
  checklist_path = (
    ROOT
    / "article_outputs/review_records/"
      "MANUAL_VISUAL_INSPECTION_CHECKLIST.tsv"
  )
  checklist = pd.read_csv(checklist_path, sep="\t", low_memory=False)
  figures = pd.read_csv(figure_manifest, sep="\t", low_memory=False)
  if len(figures) != 80:
    pytest.skip("Existing real outputs predate the v0.11.0 publication contract")
  assert len(checklist) == len(figures) == 80
  assert set(checklist["figure_id"].astype(str)) == set(
    figures["figure_id"].astype(str)
  )
  assert checklist["manual_review_status"].eq(
    "pending_manual_review"
  ).all()
