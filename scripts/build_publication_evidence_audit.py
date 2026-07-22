#!/usr/bin/env python3
"""Build publication-grade coverage, missingness and score-decomposition audits.

Inputs
------
- cancer-specific RSES-Onco ranking;
- functional evidence table;
- pharmacology evidence table;
- acquisition metadata/status JSON files.

Outputs
-------
Open TSV tables under ``article_outputs/tables/qc`` and
``article_outputs/tables/score_components``. Missing evidence is retained as
missing, non-eligible domains are explicit, and technical source failures are
never interpreted as biological negatives.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

from rses_onco.audit import (
  build_candidate_domain_audit,
  coverage_summary,
  evidence_category_table,
  missingness_summary,
  score_decomposition,
)

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists():
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t", low_memory=False)
  except EmptyDataError:
    return pd.DataFrame()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> Path:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  return path


def atomic_json(payload: object, path: Path) -> Path:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
  temporary.replace(path)
  return path


def text_value(value: object) -> str:
  if value is None:
    return ""
  try:
    if pd.isna(value):
      return ""
  except (TypeError, ValueError):
    pass
  return str(value).strip()


def normalize_doi(value: object) -> str:
  text = text_value(value).casefold()
  text = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", text)
  return text.strip()


def evidence_identifier(record: dict[str, object], row_number: int) -> tuple[str, str]:
  for column in (
    "evidence_id", "interaction_id", "assay_chembl_id", "document_chembl_id",
    "publication_id", "pmid", "pmids", "doi", "primary_doi", "supporting_doi",
  ):
    value = text_value(record.get(column))
    if value:
      return value, column
  payload = "|".join(
    text_value(record.get(column))
    for column in (
      "source", "target_gene", "drug_id", "drug_name", "interaction_type",
      "mechanism_of_action", "references", "raw_record",
    )
  )
  if payload.strip("|"):
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20], "content_hash"
  return f"row_{row_number}", "row_number"


def publication_key(record: dict[str, object]) -> str:
  doi = next(
    (normalize_doi(record.get(column)) for column in ("doi", "primary_doi", "supporting_doi") if normalize_doi(record.get(column))),
    "",
  )
  if doi:
    return f"doi:{doi}"
  for column in ("pmid", "pmids", "references", "document_chembl_id"):
    value = text_value(record.get(column))
    if value:
      return f"{column}:{value}"
  return "not_recorded"


def evidence_domain(source: str, record: dict[str, object]) -> str:
  source_cf = source.casefold()
  if any(token in source_cf for token in ("string", "bioplex", "ppi")):
    return "interaction_network"
  if any(token in source_cf for token in ("dorothea", "omnipath", "regul")):
    return "regulatory_network"
  if any(token in source_cf for token in ("chembl", "opentarget", "dgidb", "mychem", "pharos", "civic")):
    return "pharmacology"
  if any(token in source_cf for token in ("uniprot", "pdb", "alphafold", "mcsa")):
    return "biochemical_structural"
  if any(token in source_cf for token in ("hpa", "protein atlas")):
    return "localization"
  if text_value(record.get("domain")):
    return text_value(record.get("domain"))
  return "unclassified"


def evidence_role(source: str, domain: str) -> str:
  if domain in {"interaction_network", "regulatory_network", "biochemical_structural", "localization"}:
    return "score"
  if domain == "pharmacology":
    return "prioritization"
  return "interpretation"


def build_overlap_registry(
  named_tables: list[tuple[str, pd.DataFrame, Path]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
  rows: list[dict[str, object]] = []
  for table_name, frame, source_path in named_tables:
    for index, record in enumerate(frame.to_dict("records")):
      source = text_value(record.get("source")) or table_name
      identifier, identifier_method = evidence_identifier(record, index)
      pub_key = publication_key(record)
      domain = evidence_domain(source, record)
      role = evidence_role(source, domain)
      rows.append({
        "evidence_id": identifier,
        "evidence_id_method": identifier_method,
        "original_publication_or_dataset": pub_key,
        "aggregating_source": source,
        "source_table": str(source_path),
        "domain": domain,
        "evidence_role": role,
        "target_gene": record.get("target_gene"),
        "pair_id": record.get("pair_id"),
        "drug_id": record.get("drug_id"),
        "drug_name": record.get("drug_name"),
        "raw_source_record": record.get("raw_record"),
      })
  registry = pd.DataFrame(rows)
  if registry.empty:
    return registry, pd.DataFrame(columns=[
      "original_publication_or_dataset", "domain", "aggregator_count",
      "evidence_rows", "overlap_class", "score_weight_multiplier", "justification",
    ])

  grouping_key = registry["original_publication_or_dataset"].astype(str)
  fallback_key = (
    registry["domain"].astype(str) + "|"
    + registry["target_gene"].astype(str) + "|"
    + registry["drug_id"].astype(str) + "|"
    + registry["evidence_id"].astype(str)
  )
  registry["deduplication_key"] = np.where(
    grouping_key.ne("not_recorded"), grouping_key, fallback_key,
  )
  counts = registry.groupby(["deduplication_key", "domain"], dropna=False).agg(
    aggregators=("aggregating_source", "nunique"),
    rows=("evidence_id", "size"),
  ).reset_index()
  registry = registry.merge(counts, on=["deduplication_key", "domain"], how="left")
  registry["overlap_class"] = np.where(
    registry["aggregators"] > 1,
    "overlapping_aggregators",
    np.where(registry["rows"] > 1, "duplicate_records_within_source", "independent_record"),
  )
  registry["score_weight_multiplier"] = np.where(
    registry["overlap_class"].eq("independent_record"), 1.0,
    1.0 / registry["rows"].clip(lower=1),
  )
  registry["use_in_score"] = registry["evidence_role"].eq("score")
  registry["independent_evidence_unit"] = registry["deduplication_key"]
  registry["justification"] = np.where(
    registry["overlap_class"].eq("independent_record"),
    "Single traceable evidence unit.",
    "Repeated representations share one total evidence weight; additional rows are interpretative confirmation.",
  )

  summary = registry.groupby(
    ["deduplication_key", "original_publication_or_dataset", "domain", "overlap_class"],
    dropna=False,
  ).agg(
    aggregator_count=("aggregating_source", "nunique"),
    evidence_rows=("evidence_id", "size"),
    assigned_total_weight=("score_weight_multiplier", "sum"),
    roles=("evidence_role", lambda values: ";".join(sorted(set(map(str, values))))),
    sources=("aggregating_source", lambda values: ";".join(sorted(set(map(str, values))))),
  ).reset_index()
  summary["assigned_total_weight"] = summary["assigned_total_weight"].clip(upper=1.0)
  summary["justification"] = np.where(
    summary["overlap_class"].eq("independent_record"),
    "Independent evidence unit.",
    "Overlapping rows are retained for traceability but receive no more than one combined evidence unit.",
  )
  return registry, summary


def validate_audit(audit: pd.DataFrame, decomposition: pd.DataFrame) -> list[str]:
  failures: list[str] = []
  required = {
    "candidate_id", "cancer", "domain", "eligible", "evidence_state",
    "component_original", "component_normalized", "domain_weight",
    "final_score_contribution", "absence_reason",
  }
  missing = sorted(required - set(audit.columns))
  if missing:
    failures.append(f"audit missing columns: {missing}")
  allowed = {
    "observed_evidence", "negative_evidence", "neutral_evidence", "missing",
    "not_eligible", "technical_failure", "insufficient_sample",
  }
  observed_states = set(audit.get("evidence_state", pd.Series(dtype=str)).dropna().astype(str))
  if not observed_states.issubset(allowed):
    failures.append(f"unexpected evidence states: {sorted(observed_states - allowed)}")
  absent = audit["component_normalized"].isna()
  if audit.loc[absent, "final_score_contribution"].notna().any():
    failures.append("missing component received a score contribution")
  not_eligible = audit["evidence_state"].eq("not_eligible")
  if audit.loc[not_eligible, "component_normalized"].notna().any():
    failures.append("non-eligible domain has a numeric component")
  if decomposition.empty:
    failures.append("score decomposition is empty")
  else:
    for pipeline, recomputed, label in (
      ("pipeline_observed_score", "recomputed_observed_score", "observed score"),
      ("pipeline_coverage", "recomputed_coverage", "coverage"),
      ("pipeline_adjusted_score", "recomputed_adjusted_score", "adjusted score"),
    ):
      valid = decomposition[[pipeline, recomputed]].dropna()
      if not valid.empty and not np.allclose(valid[pipeline], valid[recomputed], atol=1e-10, rtol=1e-8):
        failures.append(f"pipeline and recomputed {label} disagree")
  return failures


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--functional-evidence", default="data/processed/expanded_pair_functional_evidence.tsv")
  parser.add_argument("--pharmacology-evidence", default="data/processed/pharmacology/pharmacology_evidence_long.tsv")
  parser.add_argument("--source-metadata-root", default="data/raw")
  parser.add_argument("--output-root", default="article_outputs")
  args = parser.parse_args()

  ranking_path = resolve_path(args.ranking)
  if not ranking_path.exists() or ranking_path.stat().st_size == 0:
    raise FileNotFoundError(f"Missing or empty ranking: {ranking_path}")
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  if ranking.empty:
    raise ValueError("The ranking contains no rows")

  functional_path = resolve_path(args.functional_evidence)
  pharmacology_path = resolve_path(args.pharmacology_evidence)
  functional = read_optional(functional_path)
  pharmacology = read_optional(pharmacology_path)
  output_root = resolve_path(args.output_root)
  qc_dir = output_root / "tables" / "qc"
  score_dir = output_root / "tables" / "score_components"
  supplementary_dir = output_root / "tables" / "supplementary"

  audit = build_candidate_domain_audit(ranking, resolve_path(args.source_metadata_root))
  decomposition = score_decomposition(audit)
  categories = evidence_category_table(ranking)

  summaries = {
    "coverage_by_domain": coverage_summary(audit, ["domain_family", "domain", "domain_label"]),
    "coverage_by_source": coverage_summary(audit, ["evidence_source"]),
    "coverage_by_cancer": coverage_summary(audit, ["domain_family", "cancer"]),
    "coverage_by_mechanistic_class": coverage_summary(audit, ["domain_family", "mechanistic_class"]),
    "coverage_by_domain_cancer_class": coverage_summary(
      audit, ["domain_family", "domain", "domain_label", "cancer", "mechanistic_class"]
    ),
    "missingness_reasons": missingness_summary(audit),
  }

  overlap, overlap_summary = build_overlap_registry([
    ("functional_evidence", functional, functional_path),
    ("pharmacology_evidence", pharmacology, pharmacology_path),
  ])

  outputs = [
    atomic_tsv(audit, qc_dir / "candidate_domain_evidence_audit.tsv"),
    atomic_tsv(decomposition, score_dir / "rses_onco_score_decomposition.tsv"),
    atomic_tsv(categories, qc_dir / "evidence_category_assignments.tsv"),
    atomic_tsv(overlap, qc_dir / "evidence_overlap_registry.tsv"),
    atomic_tsv(overlap_summary, qc_dir / "evidence_overlap_summary.tsv"),
  ]
  for name, frame in summaries.items():
    outputs.append(atomic_tsv(frame, qc_dir / f"{name}.tsv"))

  complete_coverage = summaries["coverage_by_domain_cancer_class"]
  outputs.append(atomic_tsv(
    complete_coverage,
    supplementary_dir / "Table_S19_coverage_domain_source_cancer_class.tsv",
  ))
  outputs.append(atomic_tsv(
    decomposition,
    supplementary_dir / "Table_S20_score_decomposition_and_coverage.tsv",
  ))
  outputs.append(atomic_tsv(
    overlap,
    supplementary_dir / "Table_S21_evidence_overlap_and_roles.tsv",
  ))
  outputs.append(atomic_tsv(
    audit,
    supplementary_dir / "Table_S22_candidate_domain_missingness_audit.tsv",
  ))

  failures = validate_audit(audit, decomposition)
  report = {
    "status": "failed" if failures else "passed",
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "ranking": str(ranking_path),
    "ranking_rows": len(ranking),
    "audit_rows": len(audit),
    "score_decomposition_rows": len(decomposition),
    "functional_evidence_rows": len(functional),
    "pharmacology_evidence_rows": len(pharmacology),
    "outputs": [str(path) for path in outputs],
    "failures": failures,
    "scientific_rules": {
      "missing_is_zero": False,
      "noneligible_enters_denominator": False,
      "overlap_total_weight_capped_at_one": True,
      "technical_failure_is_negative_evidence": False,
    },
  }
  atomic_json(report, qc_dir / "evidence_audit_validation.json")

  print(f"Candidate-domain audit rows: {len(audit):,}")
  print(f"Score decomposition rows: {len(decomposition):,}")
  print(f"Evidence-overlap registry rows: {len(overlap):,}")
  for path in outputs:
    if not path.exists() or path.stat().st_size == 0:
      raise RuntimeError(f"Mandatory audit output is missing or empty: {path}")
    print(f"Wrote {path}")
  if failures:
    raise RuntimeError("Scientific evidence audit failed: " + "; ".join(failures))
  print("Publication evidence audit passed.")


if __name__ == "__main__":
  main()
