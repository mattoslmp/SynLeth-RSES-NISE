#!/usr/bin/env python3
"""Validate GDC promoter-methylation evidence and RSES-Onco integration."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--evidence",
    default="data/processed/methylation/pair_promoter_methylation_evidence.tsv",
  )
  parser.add_argument("--ranking", required=True)
  args = parser.parse_args()

  evidence = pd.read_csv(resolve_path(args.evidence), sep="\t", low_memory=False)
  ranking = pd.read_csv(resolve_path(args.ranking), sep="\t", low_memory=False)
  errors: list[str] = []

  required_evidence = {
    "pair_id",
    "cancer",
    "evidence_status",
    "promoter_methylation_context_score",
    "methylation_source",
    "direct_gene_silencing_claim",
  }
  missing = sorted(required_evidence - set(evidence.columns))
  if missing:
    errors.append("missing_evidence_columns:" + ",".join(missing))
  if not missing:
    observed = evidence["evidence_status"].astype(str).eq("observed")
    scores = pd.to_numeric(
      evidence.loc[observed, "promoter_methylation_context_score"],
      errors="coerce",
    )
    if scores.isna().any() or not scores.between(0, 1, inclusive="both").all():
      errors.append("methylation_scores_outside_0_1_or_missing")
    direct_claim = evidence["direct_gene_silencing_claim"].astype(str).str.casefold()
    if direct_claim.isin({"true", "1", "yes"}).any():
      errors.append("unsupported_direct_gene_silencing_claim")
    duplicates = evidence.duplicated(["pair_id", "cancer"])
    if duplicates.any():
      errors.append("duplicate_pair_cancer_methylation_rows")

  required_ranking = {
    "pair_id",
    "cancer",
    "component_expression_compensation",
    "expression_compensation_expression_only",
    "promoter_methylation_context_score",
    "expression_methylation_subcoverage",
    "score_version",
    "expression_regulatory_semantics_version",
  }
  missing_ranking = sorted(required_ranking - set(ranking.columns))
  if missing_ranking:
    errors.append("missing_ranking_columns:" + ",".join(missing_ranking))
  else:
    if set(ranking["score_version"].dropna().astype(str)) != {
      "RSES-Onco-expanded-v0.11.1"
    }:
      errors.append("unexpected_score_version")
    if set(
      ranking["expression_regulatory_semantics_version"]
      .dropna()
      .astype(str)
    ) != {"eligibility-aware-wgcna-regulatory-methylation-v4"}:
      errors.append("unexpected_expression_regulatory_semantics")
    subcoverage = pd.to_numeric(
      ranking["expression_methylation_subcoverage"],
      errors="coerce",
    )
    finite = subcoverage[np.isfinite(subcoverage)]
    if not finite.between(0, 1, inclusive="both").all():
      errors.append("expression_methylation_subcoverage_outside_0_1")

  if errors:
    raise SystemExit(
      "Methylation validation failed:\n"
      + "\n".join(f"- {error}" for error in errors)
    )
  print("Methylation evidence and score integration validation passed.")
  print(evidence["evidence_status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
  main()
