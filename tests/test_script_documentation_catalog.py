from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_script_catalog_covers_every_pipeline_source(
  tmp_path: Path,
) -> None:
  markdown = tmp_path / "SCRIPT_CATALOG.md"
  manifest = tmp_path / "script_manifest.tsv"
  processed = tmp_path / "pipeline_script_catalog.tsv"
  subprocess.run([
    sys.executable,
    "scripts/build_script_documentation.py",
    "--output-md",
    str(markdown),
    "--output-tsv",
    str(manifest),
    "--processed-output",
    str(processed),
  ], cwd=ROOT, check=True)
  frame = pd.read_csv(manifest, sep="\t", low_memory=False)
  expected = {
    path.relative_to(ROOT).as_posix()
    for directory in (ROOT / "scripts", ROOT / "src/rses_onco")
    for path in directory.rglob("*")
    if path.is_file()
    and path.suffix in {".py", ".sh", ".R", ".r"}
    and "__pycache__" not in path.parts
  }
  assert set(frame["script_path"]) == expected
  assert frame["purpose"].astype(str).str.len().gt(0).all()
  assert frame["sha256"].astype(str).str.fullmatch(
    r"[0-9a-f]{64}"
  ).all()
  assert markdown.exists() and markdown.stat().st_size > 0
  assert processed.exists() and processed.stat().st_size > 0
