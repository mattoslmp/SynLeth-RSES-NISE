#!/usr/bin/env python3
"""Create publication-asset inventories, provenance and SHA-256 checksums."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys

import matplotlib
import numpy
import pandas
import scipy

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def git_value(arguments: list[str]) -> str | None:
  try:
    return subprocess.check_output(
      ["git", *arguments], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
    ).strip()
  except Exception:
    return None


def inventory_files(article_root: Path) -> list[Path]:
  excluded_names = {
    "SHA256SUMS.txt",
    "publication_file_inventory.tsv",
    "publication_provenance.json",
  }
  return sorted(
    path
    for path in article_root.rglob("*")
    if path.is_file() and path.name not in excluded_names
  )


def classify(path: Path, article_root: Path) -> str:
  relative = path.relative_to(article_root)
  parts = relative.parts
  if parts[:2] == ("figures", "main"):
    return "main_figure"
  if parts[:2] == ("figures", "supplementary"):
    return "supplementary_figure"
  table_categories = {
    ("tables", "main"): "main_table",
    ("tables", "supplementary"): "supplementary_table",
    ("tables", "qc"): "quality_control_table",
    ("tables", "score_components"): "score_component_table",
    ("tables", "robustness"): "robustness_table",
    ("tables", "figure_data"): "exact_figure_source_table",
    ("tables", "supporting_evidence"): "supporting_evidence_table",
  }
  for prefix, category in table_categories.items():
    if parts[:2] == prefix:
      return category
  if parts and parts[0] == "source_data":
    return "source_data"
  if parts and parts[0] == "workbooks":
    return "workbook"
  if parts and parts[0] == "manifests":
    return "manifest"
  if parts and parts[0] == "manuscript_assets":
    return "manuscript_asset"
  if parts and parts[0] == "structure_atlas":
    return "structure_render"
  return "other"


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument(
    "--input",
    action="append",
    default=[],
    help="Important analysis input to fingerprint; repeatable.",
  )
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  manifest_dir = article_root / "manifests"
  manifest_dir.mkdir(parents=True, exist_ok=True)
  files = inventory_files(article_root)
  inventory_rows = []
  for path in files:
    inventory_rows.append({
      "path": str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path),
      "relative_to_article_root": str(path.relative_to(article_root)),
      "category": classify(path, article_root),
      "size_bytes": path.stat().st_size,
      "sha256": sha256(path),
      "modified_at_epoch": path.stat().st_mtime,
    })
  inventory = pandas.DataFrame(inventory_rows)
  inventory_path = manifest_dir / "publication_file_inventory.tsv"
  inventory.to_csv(inventory_path, sep="\t", index=False)

  input_rows = []
  for value in args.input:
    path = resolve_path(value)
    input_rows.append({
      "path": str(path.relative_to(ROOT) if path.exists() and path.is_relative_to(ROOT) else path),
      "exists": path.exists(),
      "size_bytes": path.stat().st_size if path.exists() else None,
      "sha256": sha256(path) if path.exists() and path.is_file() else None,
    })
  pandas.DataFrame(input_rows).to_csv(
    manifest_dir / "analysis_input_fingerprints.tsv", sep="\t", index=False
  )

  category_counts = (
    inventory.groupby("category").size().astype(int).to_dict()
    if not inventory.empty else {}
  )
  provenance = {
    "repository_root": str(ROOT),
    "git_commit": git_value(["rev-parse", "HEAD"]),
    "git_branch": git_value(["rev-parse", "--abbrev-ref", "HEAD"]),
    "git_status_porcelain": git_value(["status", "--porcelain"]),
    "python": sys.version,
    "platform": platform.platform(),
    "packages": {
      "numpy": numpy.__version__,
      "pandas": pandas.__version__,
      "scipy": scipy.__version__,
      "matplotlib": matplotlib.__version__,
    },
    "article_root": str(article_root),
    "asset_count": len(files),
    "category_counts": category_counts,
    "main_figure_files": int(
      inventory.loc[
        inventory["category"].eq("main_figure")
        & inventory["path"].str.endswith((".png", ".pdf", ".svg")),
      ].shape[0]
    ) if not inventory.empty else 0,
    "supplementary_figure_files": int(
      inventory.loc[
        inventory["category"].eq("supplementary_figure")
        & inventory["path"].str.endswith((".png", ".pdf", ".svg")),
      ].shape[0]
    ) if not inventory.empty else 0,
  }
  (manifest_dir / "publication_provenance.json").write_text(
    json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
  )

  checksum_files = inventory_files(article_root)
  checksum_rows = [
    f"{sha256(path)}  {path.relative_to(article_root)}"
    for path in checksum_files
  ]
  (manifest_dir / "SHA256SUMS.txt").write_text(
    "\n".join(checksum_rows) + ("\n" if checksum_rows else ""),
    encoding="utf-8",
  )

  readme = f"""# RSES-Onco publication assets

This directory was generated entirely by repository scripts.

- Main figures: `figures/main/`
- Supplementary figures: `figures/supplementary/`
- Main tables: `tables/main/`
- Supplementary tables: `tables/supplementary/`
- Candidate-domain coverage and missingness audits: `tables/qc/`
- Score decomposition: `tables/score_components/`
- Robustness analyses: `tables/robustness/`
- Exact source table used by each figure: `tables/figure_data/`
- Expression, network, phenotype, tumor-event, structural and pharmacology support: `tables/supporting_evidence/`
- Original figure/table source data: `source_data/`
- Figure legends, score formulas and reproduction methods: `manuscript_assets/`
- Workbooks: `workbooks/`
- Provenance, manifests and SHA-256 checksums: `manifests/`

Automated layout audits are stored beside every figure. Scientific-integrity validation
checks missingness semantics, formula reproduction, overlap control, FDR boundaries,
biologically valid event frequencies and exact figure-source correspondence. These
checks do not replace manual inspection of every rendered page at 100% zoom.

Pharmacology outputs are computational experimental-priority hypotheses, not medical
advice, clinical efficacy evidence or claims of cure.

Git commit: `{provenance.get('git_commit')}`
Assets inventoried: `{len(files)}`
"""
  (article_root / "README.md").write_text(readme, encoding="utf-8")
  print(f"Inventoried {len(files):,} publication assets")
  print(f"Wrote {inventory_path}")
  print(f"Wrote {manifest_dir / 'SHA256SUMS.txt'}")


if __name__ == "__main__":
  main()
