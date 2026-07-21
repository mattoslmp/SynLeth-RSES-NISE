#!/usr/bin/env python3
"""Build the complete organized main/supplementary article workbook."""
from __future__ import annotations

import argparse
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def safe_sheet_name(name: str, used: set[str]) -> str:
  base = name[:31]
  candidate = base
  suffix = 1
  while candidate in used:
    trailer = f"_{suffix}"
    candidate = base[: 31 - len(trailer)] + trailer
    suffix += 1
  used.add(candidate)
  return candidate


def format_workbook(path: Path) -> None:
  workbook = load_workbook(path)
  header_fill = PatternFill("solid", fgColor="17365D")
  header_font = Font(color="FFFFFF", bold=True, size=10)
  title_fill = PatternFill("solid", fgColor="DDEBF7")
  for worksheet in workbook.worksheets:
    worksheet.freeze_panes = "A2"
    worksheet.auto_filter.ref = worksheet.dimensions
    for cell in worksheet[1]:
      cell.fill = header_fill
      cell.font = header_font
      cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    worksheet.row_dimensions[1].height = 36
    for column_index, column_cells in enumerate(worksheet.columns, start=1):
      maximum = 0
      for cell in list(column_cells)[:300]:
        value = "" if cell.value is None else str(cell.value)
        maximum = max(maximum, max((len(line) for line in value.splitlines()), default=0))
        if cell.row > 1:
          cell.alignment = Alignment(vertical="top", wrap_text=True)
      worksheet.column_dimensions[get_column_letter(column_index)].width = min(max(maximum + 2, 11), 55)
  if "Contents" in workbook.sheetnames:
    contents = workbook["Contents"]
    for cell in contents[1]:
      cell.fill = title_fill
      cell.font = Font(bold=True, color="17365D")
  workbook.save(path)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument(
    "--output",
    default="article_outputs/workbooks/RSES_Onco_Article_Tables_and_Evidence.xlsx",
  )
  args = parser.parse_args()

  article_root = resolve_path(args.article_root)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  main_tables = sorted((article_root / "tables" / "main").glob("*.tsv"))
  supplementary_tables = sorted(
    (article_root / "tables" / "supplementary").glob("*.tsv")
  )
  manifests = sorted((article_root / "manifests").glob("*.tsv"))
  if not main_tables or not supplementary_tables:
    raise RuntimeError(
      "Article tables are absent. Run scripts/export_article_tables.py first."
    )

  contents_rows = []
  used: set[str] = set()
  with pd.ExcelWriter(output, engine="openpyxl") as writer:
    for category, paths in (
      ("Main", main_tables),
      ("Supplementary", supplementary_tables),
      ("Manifest", manifests),
    ):
      for path in paths:
        frame = pd.read_csv(path, sep="\t")
        sheet = safe_sheet_name(path.stem, used)
        frame.to_excel(writer, sheet_name=sheet, index=False)
        contents_rows.append({
          "category": category,
          "sheet": sheet,
          "source_file": str(path),
          "rows": len(frame),
          "columns": len(frame.columns),
        })
    contents = pd.DataFrame(contents_rows)
    contents.to_excel(writer, sheet_name="Contents", index=False)
  format_workbook(output)
  print(f"Workbook sheets: {len(contents_rows) + 1}")
  print(f"Wrote {output}")


if __name__ == "__main__":
  main()
