#!/usr/bin/env python3
"""Validate generated DOCX/PDF publication documents and page separation."""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str | Path) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def require(path: Path, minimum: int = 1) -> None:
  if not path.exists() or not path.is_file() or path.stat().st_size < minimum:
    raise FileNotFoundError(f"Missing, invalid or too-small publication document: {path}")


def pdf_pages(pdf: Path, pdfinfo: str) -> int:
  result = subprocess.check_output([pdfinfo, str(pdf)], text=True)
  match = re.search(r"^Pages:\s+(\d+)", result, flags=re.MULTILINE)
  if not match:
    raise RuntimeError(f"Could not determine page count for {pdf}")
  return int(match.group(1))


def figure_pages(pdf: Path, identifiers: list[str], pdftotext: str, pdfinfo: str) -> dict[str, int]:
  pages = pdf_pages(pdf, pdfinfo)
  found: dict[str, int] = {}
  for page in range(1, pages + 1):
    text = subprocess.check_output([
      pdftotext, "-f", str(page), "-l", str(page), "-layout", str(pdf), "-"
    ], text=True, errors="replace")
    for identifier in identifiers:
      if re.search(rf"\b{re.escape(identifier)}\.", text) and identifier not in found:
        found[identifier] = page
  return found


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument("--document-dir", default="article_outputs/documents")
  parser.add_argument("--require-page-renders", action=argparse.BooleanOptionalAction, default=False)
  args = parser.parse_args()

  article_root = resolve_path(args.article_root)
  document_dir = resolve_path(args.document_dir)
  manifest = pd.read_csv(article_root / "manifests/figure_manifest.tsv", sep="\t", low_memory=False)
  article_docx = document_dir / "RSES_Onco_manuscript.docx"
  supplement_docx = document_dir / "RSES_Onco_supplementary_material.docx"
  article_pdf = document_dir / "RSES_Onco_manuscript.pdf"
  supplement_pdf = document_dir / "RSES_Onco_supplementary_material.pdf"
  for path in (article_docx, supplement_docx, article_pdf, supplement_pdf):
    require(path, minimum=1000)

  pdftotext = shutil.which("pdftotext")
  pdfinfo = shutil.which("pdfinfo")
  if not pdftotext or not pdfinfo:
    raise RuntimeError("pdftotext and pdfinfo are required for document validation")

  main_ids = sorted(
    manifest.loc[manifest["category"].astype(str).eq("main"), "figure_id"].astype(str),
    key=lambda value: int(re.search(r"\d+", value).group()),
  )
  supplementary_ids = sorted(
    manifest.loc[manifest["category"].astype(str).eq("supplementary"), "figure_id"].astype(str),
    key=lambda value: int(re.search(r"\d+", value).group()),
  )
  main_pages = figure_pages(article_pdf, main_ids, pdftotext, pdfinfo)
  supplement_pages = figure_pages(supplement_pdf, supplementary_ids, pdftotext, pdfinfo)
  missing_main = sorted(set(main_ids) - set(main_pages))
  missing_supplementary = sorted(set(supplementary_ids) - set(supplement_pages))
  if missing_main or missing_supplementary:
    raise RuntimeError(
      f"Missing figure captions in rendered documents; main={missing_main}, supplementary={missing_supplementary}"
    )
  if supplement_pages["Figure_S68"] == supplement_pages["Figure_S69"]:
    raise RuntimeError("Figures S68 and S69 were rendered on the same page")

  page_rows = [
    {"document": "main", "figure_id": key, "page": value}
    for key, value in main_pages.items()
  ] + [
    {"document": "supplementary", "figure_id": key, "page": value}
    for key, value in supplement_pages.items()
  ]
  pd.DataFrame(page_rows).to_csv(
    document_dir / "document_figure_page_map.tsv", sep="\t", index=False
  )

  if args.require_page_renders:
    for stem in ("RSES_Onco_manuscript", "RSES_Onco_supplementary_material"):
      page_dir = document_dir / f"{stem}_pages"
      if not page_dir.is_dir() or not list(page_dir.glob("page-*.png")):
        raise RuntimeError(f"Missing rendered page images for 100% inspection: {page_dir}")

  print(f"Validated {len(main_pages)} main and {len(supplement_pages)} supplementary figure captions.")
  print(f"Figure S68 page: {supplement_pages['Figure_S68']}; Figure S69 page: {supplement_pages['Figure_S69']}")


if __name__ == "__main__":
  main()
