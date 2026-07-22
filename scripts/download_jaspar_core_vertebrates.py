#!/usr/bin/env python3
"""Download the official JASPAR 2026 CORE vertebrate non-redundant MEME file."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_URL = (
  "https://jaspar.elixir.no/download/data/2026/CORE/"
  "JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
)


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--url", default=DEFAULT_URL)
  parser.add_argument(
    "--output",
    default="data/raw/regulatory/JASPAR2026_CORE_vertebrates_non-redundant.meme",
  )
  parser.add_argument(
    "--status-output",
    default="data/raw/regulatory/jaspar2026_core_status.json",
  )
  parser.add_argument("--refresh", action="store_true")
  args = parser.parse_args()

  output = resolve_path(args.output)
  status_output = resolve_path(args.status_output)
  output.parent.mkdir(parents=True, exist_ok=True)
  source_status = "cache"
  if args.refresh or not output.exists() or output.stat().st_size == 0:
    response = requests.get(args.url, timeout=180)
    response.raise_for_status()
    text = response.text
    if "MEME version" not in text or "MOTIF" not in text:
      raise ValueError("Downloaded JASPAR file is not a MEME motif file")
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(output)
    source_status = "downloaded"
  content = output.read_text(encoding="utf-8", errors="replace")
  motif_count = sum(line.startswith("MOTIF ") for line in content.splitlines())
  if motif_count == 0:
    raise ValueError("JASPAR MEME file contains no motifs")
  status = {
    "source": "JASPAR 2026 CORE vertebrates non-redundant",
    "url": args.url,
    "status": source_status,
    "motif_count": motif_count,
    "sha256": sha256(output),
    "output": str(output),
    "accessed_at_utc": datetime.now(timezone.utc).isoformat(),
    "interpretation": (
      "JASPAR matrices support predicted promoter motif occurrence; motif occurrence "
      "is not direct TF binding or transcriptional causality."
    ),
  }
  temporary = status_output.with_suffix(status_output.suffix + ".tmp")
  temporary.write_text(json.dumps(status, indent=2, sort_keys=True), encoding="utf-8")
  temporary.replace(status_output)
  print(f"JASPAR motifs: {motif_count:,} ({source_status})")
  print(f"Wrote {output}")


if __name__ == "__main__":
  main()
