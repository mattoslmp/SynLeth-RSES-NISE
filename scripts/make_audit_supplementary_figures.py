#!/usr/bin/env python3
"""Generate supplementary evidence-audit and RSES robustness figures S33-S38."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from rses_onco.publication import write_figure_manifest, write_legends_markdown
from scripts.publication_audit_figures import (
  figure_s33,
  figure_s34,
  figure_s37,
  figure_s38,
)
from scripts.publication_robustness_figures import figure_s35, figure_s36

ROOT = Path(__file__).resolve().parents[1]
AUDIT_FIGURE_BUILDERS = [
  figure_s33,
  figure_s34,
  figure_s35,
  figure_s36,
  figure_s37,
  figure_s38,
]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument(
    "--strict-layout",
    action=argparse.BooleanOptionalAction,
    default=True,
  )
  args = parser.parse_args()

  output_root = resolve_path(args.output_root)
  records = [
    builder(output_root, args.strict_layout)
    for builder in AUDIT_FIGURE_BUILDERS
  ]
  manifest_dir = output_root / "manifests"
  legend_dir = output_root / "manuscript_assets"
  write_figure_manifest(
    records,
    manifest_dir / "audit_supplementary_figure_manifest.tsv",
  )
  write_legends_markdown(
    records,
    legend_dir / "audit_supplementary_figure_legends.md",
  )
  frame = pd.DataFrame([asdict(record) for record in records])
  print(
    frame[["figure_id", "layout_status", "base_path"]].to_string(
      index=False
    )
  )
  if len(records) != 6 or not frame["layout_status"].eq("pass").all():
    raise RuntimeError(
      "Audit supplementary figure generation did not produce six passing figures"
    )
  print(
    "Wrote audit supplementary figures to "
    f"{output_root / 'figures/supplementary'}"
  )


if __name__ == "__main__":
  main()
