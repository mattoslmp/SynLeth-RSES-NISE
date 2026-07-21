#!/usr/bin/env python3
"""Render publication-grade whole-structure and functional-site views with PyMOL.

All residue labels are placed in the final scripted figure layout rather than on
top of the molecular render, preventing label/structure overlap. Important known
residues are shown as colored sticks and CA spheres. AlphaFold confidence is read
from the PDB B-factor field and displayed using the standard confidence bands.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any

import pandas as pd

from rses_onco.structural import build_pymol_script, file_sha256, write_json

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def find_pymol(requested: str) -> str:
  executable = shutil.which(requested)
  if executable:
    return executable
  for candidate in ("pymol", "pymol-open-source"):
    executable = shutil.which(candidate)
    if executable:
      return executable
  raise FileNotFoundError(
    "PyMOL was not found. Install with: conda install -c conda-forge pymol-open-source"
  )


def run_pymol(executable: str, script_path: Path, log_path: Path) -> None:
  log_path.parent.mkdir(parents=True, exist_ok=True)
  process = subprocess.run(
    [executable, "-cq", str(script_path)],
    cwd=ROOT,
    text=True,
    capture_output=True,
  )
  log_path.write_text(
    (process.stdout or "") + "\n--- STDERR ---\n" + (process.stderr or ""),
    encoding="utf-8",
  )
  if process.returncode != 0:
    raise RuntimeError(
      f"PyMOL failed with code {process.returncode}; see {log_path}"
    )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--structure-manifest",
    default="data/processed/structures/alphafold_structure_manifest.tsv",
  )
  parser.add_argument(
    "--annotations",
    default="data/processed/structures/nise_structural_residue_annotations.tsv",
  )
  parser.add_argument(
    "--output-dir",
    default="article_outputs/structure_atlas/individual",
  )
  parser.add_argument(
    "--render-manifest",
    default="data/processed/structures/nise_structure_render_manifest.tsv",
  )
  parser.add_argument("--pymol", default="pymol")
  parser.add_argument("--width", type=int, default=2800)
  parser.add_argument("--height", type=int, default=2200)
  parser.add_argument("--refresh", action="store_true")
  parser.add_argument("--strict", action="store_true")
  args = parser.parse_args()

  structures = pd.read_csv(resolve_path(args.structure_manifest), sep="\t")
  annotation_path = resolve_path(args.annotations)
  if annotation_path.exists() and annotation_path.stat().st_size > 0:
    annotations = pd.read_csv(annotation_path, sep="\t")
  else:
    annotations = pd.DataFrame()
  structures = structures.loc[structures["status"].astype(str).eq("ok")].copy()
  if structures.empty:
    raise SystemExit("No validated AlphaFold structures were available")

  pymol = find_pymol(args.pymol)
  output_dir = resolve_path(args.output_dir)
  script_dir = output_dir / "pymol_scripts"
  log_dir = output_dir / "logs"
  output_dir.mkdir(parents=True, exist_ok=True)
  rows: list[dict[str, Any]] = []
  failures: list[dict[str, str]] = []

  # Use the first/complete record for each accession; fragments remain in the raw manifest.
  structures = structures.sort_values(
    ["uniprot_accession", "uniprot_start", "uniprot_end", "fragment_index"],
    na_position="last",
  ).drop_duplicates("uniprot_accession", keep="first")

  for index, record in enumerate(structures.to_dict("records"), start=1):
    gene = str(record["gene_symbol"])
    accession = str(record["uniprot_accession"])
    pdb_path = resolve_path(str(record["pdb_path"]))
    subset = (
      annotations.loc[annotations["uniprot_accession"].astype(str).eq(accession)].copy()
      if not annotations.empty and "uniprot_accession" in annotations
      else pd.DataFrame()
    )
    gene_dir = output_dir / str(record["group_id"]) / gene
    gene_dir.mkdir(parents=True, exist_ok=True)
    views = {
      "whole": gene_dir / f"{gene}_{accession}_whole.png",
      "site": gene_dir / f"{gene}_{accession}_site.png",
    }
    for view, png_path in views.items():
      script_path = script_dir / f"{gene}_{accession}_{view}.pml"
      log_path = log_dir / f"{gene}_{accession}_{view}.log"
      try:
        if not png_path.exists() or png_path.stat().st_size < 1000 or args.refresh:
          script = build_pymol_script(
            pdb_path,
            png_path,
            subset,
            view=view,
            width=args.width,
            height=args.height,
          )
          script_path.parent.mkdir(parents=True, exist_ok=True)
          script_path.write_text(script, encoding="utf-8")
          run_pymol(pymol, script_path, log_path)
        if not png_path.exists() or png_path.stat().st_size < 1000:
          raise RuntimeError("Rendered image is missing or unexpectedly small")
        rows.append({
          **record,
          "view": view,
          "render_path": str(png_path),
          "render_size_bytes": png_path.stat().st_size,
          "render_sha256": file_sha256(png_path),
          "pymol_script": str(script_path),
          "pymol_log": str(log_path),
          "annotated_residue_count": int(subset["residue_number"].nunique()) if not subset.empty else 0,
          "annotation_types": ";".join(sorted(set(subset["annotation_type"].astype(str)))) if not subset.empty else "",
          "status": "ok",
        })
      except Exception as exc:
        failure = {
          "gene_symbol": gene,
          "uniprot_accession": accession,
          "view": view,
          "status": "failed",
          "message": str(exc),
        }
        failures.append(failure)
        print(f"[Render {index}/{len(structures)}] FAILED {gene} {view}: {exc}", flush=True)
        if args.strict:
          raise
    print(
      f"[Render {index}/{len(structures)}] {gene} {accession}: "
      f"{int(subset['residue_number'].nunique()) if not subset.empty else 0} annotated residues",
      flush=True,
    )

  manifest = resolve_path(args.render_manifest)
  manifest.parent.mkdir(parents=True, exist_ok=True)
  frame = pd.DataFrame(rows)
  if not frame.empty:
    frame = frame.sort_values(["group_id", "gene_symbol", "view"])
  frame.to_csv(manifest, sep="\t", index=False)
  pd.DataFrame(failures).to_csv(
    manifest.with_name("nise_structure_render_failures.tsv"), sep="\t", index=False
  )
  summary = {
    "proteins_requested": int(structures["uniprot_accession"].nunique()),
    "successful_proteins": int(frame["uniprot_accession"].nunique()) if not frame.empty else 0,
    "successful_renders": int(len(frame)),
    "failed_renders": int(len(failures)),
    "pymol_executable": pymol,
    "image_dimensions": [args.width, args.height],
    "dpi": 600,
  }
  write_json(manifest.with_suffix(".summary.json"), summary)
  print(json.dumps(summary, indent=2, sort_keys=True))
  print(f"Wrote {manifest}")


if __name__ == "__main__":
  main()
