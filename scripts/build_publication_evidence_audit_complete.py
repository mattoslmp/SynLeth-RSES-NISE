#!/usr/bin/env python3
"""Run the evidence audit and finalize complete coverage and overlap views."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.build_publication_evidence_audit as target
from rses_onco.audit import coverage_summary


def argument_value(name: str, default: str) -> str:
  prefix = f"--{name}="
  arguments = sys.argv[1:]
  for index, argument in enumerate(arguments):
    if argument.startswith(prefix):
      return argument.split("=", 1)[1]
    if argument == f"--{name}" and index + 1 < len(arguments):
      return arguments[index + 1]
  return default


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t", low_memory=False)
  except EmptyDataError:
    return pd.DataFrame()


def source_table(path: Path, label: str) -> pd.DataFrame:
  frame = read_optional(path)
  if frame.empty:
    return frame
  frame = frame.copy()
  if "source" in frame.columns:
    frame["source_record_identifier"] = frame["source"]
  frame["source"] = label
  return frame


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def main() -> None:
  target.main()
  output_root = resolve_path(argument_value("output-root", "article_outputs"))
  audit_path = output_root / "tables/qc/candidate_domain_evidence_audit.tsv"
  audit = pd.read_csv(audit_path, sep="\t", low_memory=False)

  complete = coverage_summary(
    audit,
    [
      "domain_family", "domain", "domain_label", "evidence_source",
      "cancer", "mechanistic_class",
    ],
  )
  coverage_path = (
    output_root
    / "tables/supplementary/Table_S19_coverage_domain_source_cancer_class.tsv"
  )
  atomic_tsv(complete, coverage_path)
  if complete.empty or not coverage_path.exists() or coverage_path.stat().st_size == 0:
    raise RuntimeError("Complete domain/source/cancer/class coverage table is empty")

  raw_dir = ROOT / "data/raw/human_functional_evidence"
  functional_path = resolve_path(
    argument_value(
      "functional-evidence",
      "data/processed/expanded_pair_functional_evidence.tsv",
    )
  )
  pharmacology_path = resolve_path(
    argument_value(
      "pharmacology-evidence",
      "data/processed/pharmacology/pharmacology_evidence_long.tsv",
    )
  )
  named_tables = [
    ("functional_pair_evidence", read_optional(functional_path), functional_path),
    ("pharmacology_evidence", read_optional(pharmacology_path), pharmacology_path),
    (
      "STRING",
      source_table(raw_dir / "string_interaction_partners.tsv", "STRING"),
      raw_dir / "string_interaction_partners.tsv",
    ),
    (
      "DoRothEA_OmniPath",
      source_table(raw_dir / "omnipath_dorothea.tsv", "DoRothEA/OmniPath"),
      raw_dir / "omnipath_dorothea.tsv",
    ),
    (
      "Human_Protein_Atlas",
      source_table(raw_dir / "hpa_subcellular_location.tsv", "Human Protein Atlas"),
      raw_dir / "hpa_subcellular_location.tsv",
    ),
    (
      "UniProtKB",
      source_table(raw_dir / "uniprot_reviewed_annotations.tsv", "UniProtKB"),
      raw_dir / "uniprot_reviewed_annotations.tsv",
    ),
  ]
  named_tables = [item for item in named_tables if not item[1].empty]
  overlap, overlap_summary = target.build_overlap_registry(named_tables)
  overlap_path = output_root / "tables/qc/evidence_overlap_registry.tsv"
  overlap_summary_path = output_root / "tables/qc/evidence_overlap_summary.tsv"
  supplementary_overlap_path = (
    output_root / "tables/supplementary/Table_S21_evidence_overlap_and_roles.tsv"
  )
  atomic_tsv(overlap, overlap_path)
  atomic_tsv(overlap_summary, overlap_summary_path)
  atomic_tsv(overlap, supplementary_overlap_path)
  if overlap.empty:
    raise RuntimeError("Evidence-overlap registry is empty after adding raw sources")

  print(
    f"Finalized complete coverage view: {coverage_path} "
    f"({len(complete):,} rows)"
  )
  print(
    f"Finalized raw-source overlap registry: {overlap_path} "
    f"({len(overlap):,} rows)"
  )


if __name__ == "__main__":
  main()
