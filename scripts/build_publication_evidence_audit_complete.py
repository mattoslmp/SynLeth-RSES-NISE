#!/usr/bin/env python3
"""Run the evidence audit and finalize complete domain/source/cancer/class views."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.build_publication_evidence_audit as target
from rses_onco.audit import coverage_summary


def argument_value(name: str, default: str) -> str:
  prefix = f"--{name}="
  for index, argument in enumerate(sys.argv[1:]):
    if argument.startswith(prefix):
      return argument.split("=", 1)[1]
    if argument == f"--{name}" and index + 2 <= len(sys.argv[1:]):
      return sys.argv[1:][index + 1]
  return default


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


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
  path = output_root / "tables/supplementary/Table_S19_coverage_domain_source_cancer_class.tsv"
  temporary = path.with_suffix(path.suffix + ".tmp")
  complete.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  if complete.empty or not path.exists() or path.stat().st_size == 0:
    raise RuntimeError("Complete domain/source/cancer/class coverage table is empty")
  print(f"Finalized complete coverage view: {path} ({len(complete):,} rows)")


if __name__ == "__main__":
  main()
