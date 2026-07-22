from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_extended_supporting_evidence_when_real_outputs_exist() -> None:
  model_manifest = (
    ROOT
    / "article_outputs/tables/supporting_evidence/"
      "model_level_supporting_evidence_manifest.tsv"
  )
  raw_manifest = (
    ROOT
    / "article_outputs/tables/supporting_evidence/networks/raw_sources/"
      "raw_functional_evidence_manifest.tsv"
  )
  if not model_manifest.exists() and not raw_manifest.exists():
    pytest.skip("Real-data extended supporting evidence has not been generated")
  if not model_manifest.exists() or not raw_manifest.exists():
    pytest.fail(
      "Only part of the mandatory extended supporting-evidence package exists"
    )
  subprocess.run(
    [
      sys.executable,
      "-u",
      "scripts/validate_extended_supporting_evidence.py",
      "--article-root",
      "article_outputs",
    ],
    cwd=ROOT,
    check=True,
  )
