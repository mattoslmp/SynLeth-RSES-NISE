#!/usr/bin/env python3
"""Build a multi-sheet Excel workbook from benchmark or expanded RSES-Onco outputs."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists():
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t")
  except EmptyDataError:
    return pd.DataFrame()


def read_first(parent: Path, names: list[str]) -> pd.DataFrame:
  for name in names:
    frame = read_optional(parent / name)
    if not frame.empty:
      return frame
  return pd.DataFrame()


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--input",
    default="results/empirical_26Q1/full/empirical_rses_onco_by_cancer.tsv",
  )
  parser.add_argument(
    "--output",
    default="supplementary/RSES_Onco_Empirical_26Q1.xlsx",
  )
  args = parser.parse_args()

  input_path = resolve_path(args.input)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  parent = input_path.parent

  tables = {
    "RSES_Onco_Ranking": pd.read_csv(input_path, sep="\t"),
    "Dependency_Contrasts": read_first(parent, [
      "dependency_contrasts.tsv",
      "expanded_dependency_contrasts.tsv",
    ]),
    "Expression_Compensation": read_first(parent, [
      "expression_contrasts.tsv",
      "expanded_expression_compensation.tsv",
    ]),
    "Expression_Context": read_optional(parent / "expanded_expression_context_profiles.tsv"),
    "CRISPR_Phenotype_Profiles": read_optional(parent / "expanded_crispr_phenotype_profiles.tsv"),
    "Skipped_Biomarkers": read_first(parent, [
      "skipped_complex_biomarkers.tsv",
      "expanded_skipped_complex_biomarkers.tsv",
    ]),
    "GDC_Matrix_QC": read_optional(parent / "gdc_matrix_qc.tsv"),
    "TCGA_Event_Summary": read_optional(parent / "tcga_gene_event_summary.tsv"),
    "Candidate_Universe": read_optional(ROOT / "data/processed/expanded_candidate_universe.tsv"),
    "Class_Member_Inventory": read_optional(ROOT / "data/processed/expanded_class_member_inventory.tsv"),
    "Functional_Evidence": read_optional(ROOT / "data/processed/expanded_pair_functional_evidence.tsv"),
  }

  with pd.ExcelWriter(output, engine="openpyxl") as writer:
    for sheet_name, frame in tables.items():
      frame.to_excel(writer, sheet_name=sheet_name[:31], index=False)
      worksheet = writer.book[sheet_name[:31]]
      worksheet.freeze_panes = "A2"
      worksheet.auto_filter.ref = worksheet.dimensions
      for column_cells in worksheet.columns:
        maximum = max(
          len(str(cell.value)) if cell.value is not None else 0
          for cell in column_cells[:200]
        )
        worksheet.column_dimensions[column_cells[0].column_letter].width = min(
          max(maximum + 2, 11), 55
        )

  print(f"Wrote {output}", flush=True)


if __name__ == "__main__":
  main()
