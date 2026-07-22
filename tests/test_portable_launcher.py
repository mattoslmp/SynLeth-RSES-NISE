from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_portable_launcher_resolves_script_namespace() -> None:
  result = subprocess.run(
    [
      "bash",
      str(ROOT / "scripts/run_expanded_pipeline_portable.sh"),
      "--check-python-imports",
    ],
    cwd=ROOT,
    check=False,
    capture_output=True,
    text=True,
  )
  assert result.returncode == 0, result.stderr
  assert "RSES-Onco script namespace imports: OK" in result.stdout
