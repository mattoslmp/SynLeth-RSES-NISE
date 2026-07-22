#!/usr/bin/env python3
"""Generate and validate every main, supplementary, structural and audit figure."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def run(command: list[str]) -> None:
  print("+", " ".join(command), flush=True)
  subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--config", default="config/article_assets.yaml")
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--candidates", default="data/processed/expanded_candidate_universe.tsv")
  parser.add_argument("--discovery", default="results/expanded_26Q1/discovery/all_target_dependency_screen.tsv")
  parser.add_argument("--pharmacology", default="results/expanded_26Q1/pharmacology/pharmacology_ranked_hypotheses.tsv")
  parser.add_argument("--proteins", default="data/curated/human_nise_bonafide_2017.tsv")
  parser.add_argument("--structure-manifest", default="data/processed/structures/alphafold_structure_manifest.tsv")
  parser.add_argument("--render-manifest", default="data/processed/structures/nise_structure_render_manifest.tsv")
  parser.add_argument("--structural-annotations", default="data/processed/structures/nise_structural_residue_annotations.tsv")
  parser.add_argument("--structural-coverage", default="data/processed/structures/nise_structural_annotation_coverage.tsv")
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--strict-layout", action=argparse.BooleanOptionalAction, default=True)
  parser.add_argument("--top-n", type=int, default=15)
  args = parser.parse_args()

  strict_flag = "--strict-layout" if args.strict_layout else "--no-strict-layout"
  common = [
    "--config", args.config,
    "--ranking", args.ranking,
    "--candidates", args.candidates,
    "--discovery", args.discovery,
    "--output-root", args.output_root,
    strict_flag,
  ]
  run([
    sys.executable, "-u", "scripts/make_main_figures_resilient.py", *common,
    "--pharmacology", args.pharmacology, "--top-n", str(args.top_n),
  ])
  run([
    sys.executable, "-u", "scripts/make_supplementary_figures_resilient.py", *common,
  ])
  run([
    sys.executable, "-u", "scripts/make_nise_structure_figures.py",
    "--ranking", args.ranking,
    "--proteins", args.proteins,
    "--structure-manifest", args.structure_manifest,
    "--render-manifest", args.render_manifest,
    "--annotations", args.structural_annotations,
    "--coverage", args.structural_coverage,
    "--output-root", args.output_root,
    strict_flag,
  ])
  run([
    sys.executable, "-u", "scripts/make_audit_supplementary_figures.py",
    "--output-root", args.output_root,
    strict_flag,
  ])

  output_root = resolve_path(args.output_root)
  manifest_names = [
    "main_figure_manifest.tsv",
    "supplementary_figure_manifest.tsv",
    "structural_main_figure_manifest.tsv",
    "structural_supplementary_figure_manifest.tsv",
    "audit_supplementary_figure_manifest.tsv",
  ]
  manifests = []
  for name in manifest_names:
    path = output_root / "manifests" / name
    if not path.exists() or path.stat().st_size == 0:
      raise FileNotFoundError(f"Missing or empty figure manifest: {path}")
    manifests.append(pd.read_csv(path, sep="\t"))
  combined = pd.concat(manifests, ignore_index=True)
  expected = 46
  if len(combined) != expected:
    raise RuntimeError(f"Expected {expected} registered figures; observed {len(combined)}")
  if combined["figure_id"].duplicated().any():
    duplicated = combined.loc[combined["figure_id"].duplicated(keep=False), "figure_id"]
    raise RuntimeError("Duplicated figure identifiers: " + ", ".join(sorted(set(duplicated))))
  if args.strict_layout and not combined["layout_status"].eq("pass").all():
    failed = combined.loc[~combined["layout_status"].eq("pass"), ["figure_id", "layout_warnings"]]
    raise RuntimeError("Figure layout audit did not pass:\n" + failed.to_string(index=False))
  combined.to_csv(output_root / "manifests/figure_manifest.tsv", sep="\t", index=False)

  missing = []
  for record in combined.to_dict("records"):
    base = Path(record["base_path"])
    if not base.is_absolute():
      base = ROOT / base
    source_path = Path(str(record["source_data_path"]))
    if not source_path.is_absolute():
      source_path = ROOT / source_path
    if not source_path.exists() or source_path.stat().st_size == 0:
      missing.append(str(source_path))
    for extension in ("png", "pdf", "svg"):
      path = base.with_suffix(f".{extension}")
      if not path.exists() or path.stat().st_size == 0:
        missing.append(str(path))
  if missing:
    raise RuntimeError("Missing or empty figure files/source data:\n" + "\n".join(missing))

  legend_names = [
    "main_figure_legends.md",
    "supplementary_figure_legends.md",
    "structural_main_figure_legends.md",
    "structural_supplementary_figure_legends.md",
    "audit_supplementary_figure_legends.md",
  ]
  legends = []
  for name in legend_names:
    path = output_root / "manuscript_assets" / name
    if path.exists():
      legends.append(path.read_text(encoding="utf-8"))
  (output_root / "manuscript_assets/all_figure_legends.md").write_text(
    "\n\n".join(legends), encoding="utf-8"
  )
  print(f"Validated {len(combined)} figures and {len(combined) * 3} exported image files.")
  print(f"Combined manifest: {output_root / 'manifests/figure_manifest.tsv'}")


if __name__ == "__main__":
  main()
