#!/usr/bin/env python3
"""Run functional-evidence acquisition with cached/pinned DoRothEA fallback.

The underlying publication logic remains in download_human_functional_evidence.py.
This wrapper replaces only its DoRothEA acquisition function, preserving STRING,
HPA, UniProt and pair-level calculations exactly as implemented there.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rses_onco.dorothea import acquire_dorothea  # noqa: E402
from scripts import download_human_functional_evidence as base  # noqa: E402


DOROTHEA_SUMMARY: dict[str, object] = {}


def env_flag(name: str, default: bool = False) -> bool:
  value = os.environ.get(name)
  if value is None:
    return default
  return value.strip().casefold() in {"1", "true", "yes", "on"}


def argument_value(name: str, default: str) -> str:
  try:
    index = sys.argv.index(name)
  except ValueError:
    return default
  if index + 1 >= len(sys.argv):
    return default
  return sys.argv[index + 1]


def resilient_dorothea(output: Path):
  retries = int(os.environ.get("RSES_DOROTHEA_RETRIES", "4"))
  refresh = env_flag("RSES_REFRESH_DOROTHEA", False)
  allow_fallback = not env_flag("RSES_DISABLE_DOROTHEA_FALLBACK", False)
  frame, summary = acquire_dorothea(
    Path(output),
    metadata_output=Path(output).parent / "dorothea_acquisition_metadata.json",
    fallback_pickle=Path(output).parent / "dorothea_hs_official.pkl",
    retries=retries,
    refresh=refresh,
    allow_fallback=allow_fallback,
  )
  DOROTHEA_SUMMARY.clear()
  DOROTHEA_SUMMARY.update(summary)
  return frame


def update_combined_metadata() -> None:
  raw_dir = Path(argument_value(
    "--raw-dir", "data/raw/human_functional_evidence"
  ))
  if not raw_dir.is_absolute():
    raw_dir = ROOT / raw_dir
  metadata_path = raw_dir / "source_metadata.json"
  if not metadata_path.exists():
    return
  payload = json.loads(metadata_path.read_text(encoding="utf-8"))
  payload["dorothea_acquisition"] = DOROTHEA_SUMMARY
  temporary = metadata_path.with_suffix(metadata_path.suffix + ".tmp")
  temporary.write_text(
    json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
  )
  temporary.replace(metadata_path)


def main() -> None:
  base.download_dorothea = resilient_dorothea
  base.main()
  update_combined_metadata()


if __name__ == "__main__":
  main()
