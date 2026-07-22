#!/usr/bin/env python3
"""Validate scientific integrity, missingness semantics and figure-source traceability."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INVALID_DISPLAY = {"", "nan", "none", "na", "null", "<na>"}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def read_tsv(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(f"Missing or empty required TSV: {path}")
  return pd.read_csv(path, sep="\t", low_memory=False)


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
  missing = sorted(columns - set(frame.columns))
  if missing:
    raise ValueError(f"{label} missing columns: {missing}")


def invalid_text(series: pd.Series) -> pd.Series:
  return series.astype(str).str.strip().str.casefold().isin(INVALID_DISPLAY)


def validate_audit(article_root: Path) -> None:
  frame = read_tsv(article_root / "tables/qc/candidate_domain_evidence_audit.tsv")
  require_columns(frame, {
    "candidate_id", "cancer", "domain", "eligible", "evidence_state",
    "component_normalized", "final_score_contribution", "absence_reason",
  }, "candidate-domain audit")
  allowed = {
    "observed_evidence", "negative_evidence", "neutral_evidence", "missing",
    "not_eligible", "technical_failure", "insufficient_sample",
  }
  unexpected = set(frame["evidence_state"].dropna().astype(str)) - allowed
  if unexpected:
    raise ValueError(f"Unexpected evidence states: {sorted(unexpected)}")
  missing_component = frame["component_normalized"].isna()
  if frame.loc[missing_component, "final_score_contribution"].notna().any():
    raise ValueError("A missing component received a numeric final contribution")
  if frame.loc[frame["evidence_state"].eq("not_eligible"), "component_normalized"].notna().any():
    raise ValueError("A non-eligible domain received a numeric component")
  observed_zero = frame["component_normalized"].eq(0)
  if not frame.loc[observed_zero, "evidence_state"].eq("negative_evidence").all():
    raise ValueError("Observed zero components are not classified as negative evidence")
  absent = frame["evidence_state"].isin({"missing", "technical_failure", "insufficient_sample"})
  if frame.loc[absent, "absence_reason"].fillna("").astype(str).str.strip().eq("").any():
    raise ValueError("Unavailable evidence lacks an explicit absence reason")


def validate_scores(article_root: Path) -> None:
  frame = read_tsv(article_root / "tables/score_components/rses_onco_score_decomposition.tsv")
  pairs = [
    ("pipeline_observed_score", "recomputed_observed_score"),
    ("pipeline_coverage", "recomputed_coverage"),
    ("pipeline_adjusted_score", "recomputed_adjusted_score"),
  ]
  require_columns(frame, {item for pair in pairs for item in pair}, "score decomposition")
  for pipeline, recomputed in pairs:
    valid = frame[[pipeline, recomputed]].dropna()
    if not valid.empty and not np.allclose(valid[pipeline], valid[recomputed], atol=1e-10, rtol=1e-8):
      raise ValueError(f"Score formula mismatch: {pipeline} versus {recomputed}")
  coverage = pd.to_numeric(frame["pipeline_coverage"], errors="coerce")
  if ((coverage < 0) | (coverage > 1)).any():
    raise ValueError("Evidence coverage outside [0,1]")


def validate_overlap(article_root: Path) -> None:
  frame = read_tsv(article_root / "tables/qc/evidence_overlap_summary.tsv")
  require_columns(frame, {"overlap_class", "assigned_total_weight", "roles"}, "evidence overlap summary")
  if (pd.to_numeric(frame["assigned_total_weight"], errors="coerce") > 1.0 + 1e-12).any():
    raise ValueError("Overlapping evidence exceeds one combined evidence unit")


def validate_figures(article_root: Path) -> None:
  manifest = read_tsv(article_root / "manifests/figure_manifest.tsv")
  require_columns(manifest, {"figure_id", "category", "base_path", "source_data_path", "layout_status"}, "figure manifest")
  if len(manifest) != 46:
    raise ValueError(f"Expected 46 registered figures; observed {len(manifest)}")
  if manifest["figure_id"].duplicated().any():
    raise ValueError("Duplicated figure identifiers")
  if not manifest["layout_status"].eq("pass").all():
    failed = manifest.loc[~manifest["layout_status"].eq("pass"), ["figure_id", "layout_warnings"]]
    raise ValueError("Layout audit failures:\n" + failed.to_string(index=False))
  for record in manifest.to_dict("records"):
    base = Path(str(record["base_path"]))
    if not base.is_absolute():
      base = ROOT / base
    source = Path(str(record["source_data_path"]))
    if not source.is_absolute():
      source = ROOT / source
    for path in [source, *[base.with_suffix(f".{extension}") for extension in ("png", "pdf", "svg")]]:
      if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty figure asset: {path}")

  main_source = article_root / "source_data/figures/main"
  fig3 = read_tsv(main_source / "Figure_3_cancer_specific_ranking_source_data.tsv")
  require_columns(fig3, {
    "display_pair_label", "rses_onco", "coverage_adjusted_rses", "evidence_coverage",
    "eligible_domain_count", "evidence_domain_count", "source_class", "statistical_status", "support_level",
  }, "Figure 3 source data")
  if invalid_text(fig3["display_pair_label"]).any():
    raise ValueError("Figure 3 contains unresolved/invalid display labels")

  fig4 = read_tsv(main_source / "Figure_4_tcga_depmap_integration_source_data.tsv")
  if "tcga_homdel_frequency" in fig4:
    values = pd.to_numeric(fig4["tcga_homdel_frequency"], errors="coerce").dropna()
    if ((values < 0) | (values > 1)).any():
      raise ValueError("Figure 4 contains biologically invalid event frequencies")

  fig5 = read_tsv(main_source / "Figure_5_functional_microniches_source_data.tsv")
  require_columns(fig5, {"domain", "display_status", "display_value", "eligible", "absence_reason", "evidence_source"}, "Figure 5 source data")
  unavailable = fig5["display_status"].isin({"missing", "technical_failure", "insufficient_sample", "not_eligible"})
  if fig5.loc[unavailable, "display_value"].notna().any():
    raise ValueError("Figure 5 unavailable/non-eligible cell has a numeric displayed value")

  fig6 = read_tsv(main_source / "Figure_6_class_discoveries_source_data.tsv")
  discovery = fig6.loc[fig6.get("panel", pd.Series(index=fig6.index, dtype=object)).astype(str).eq("B_all_target_discovery")]
  if not discovery.empty:
    q_col = next((column for column in ("q_value_bh_within_loss_cancer", "q_value_bh") if column in discovery), None)
    if q_col is None or not (pd.to_numeric(discovery[q_col], errors="coerce") < 0.05).all():
      raise ValueError("Figure 6 discovery panel contains a non-FDR-supported row")
    if not (pd.to_numeric(discovery["delta_effect"], errors="coerce") < 0).all():
      raise ValueError("Figure 6 discovery panel contains a non-supportive effect direction")

  fig7 = read_tsv(main_source / "Figure_7_pharmacology_source_data.tsv")
  require_columns(fig7, {"compound_display", "compound_resolution", "display_pair_label"}, "Figure 7 source data")
  if invalid_text(fig7["compound_display"]).any() or fig7["compound_resolution"].eq("target_only").any():
    raise ValueError("Figure 7 main table contains unresolved/target-only compound evidence")


def validate_tables(article_root: Path) -> None:
  manifest = read_tsv(article_root / "manifests/table_manifest.tsv")
  require_columns(manifest, {"table_id", "category", "path", "rows", "status"}, "table manifest")
  if int(manifest["category"].eq("main").sum()) != 4:
    raise ValueError("Expected four main tables")
  if int(manifest["category"].eq("supplementary").sum()) != 25:
    raise ValueError("Expected 25 supplementary tables")
  for record in manifest.to_dict("records"):
    path = Path(str(record["path"]))
    if not path.is_absolute():
      path = ROOT / path
    if not path.exists() or path.stat().st_size == 0:
      raise FileNotFoundError(f"Missing or empty registered table: {path}")


def validate_figure_catalog(article_root: Path) -> None:
  inventory = read_tsv(article_root / "tables/figure_data/figure_source_data_inventory.tsv")
  figures = read_tsv(article_root / "manifests/figure_manifest.tsv")
  if len(inventory) != len(figures):
    raise ValueError("Not every figure has an exact source-data catalogue entry")
  if set(inventory["figure_id"].astype(str)) != set(figures["figure_id"].astype(str)):
    raise ValueError("Figure source-data inventory and manifest identifiers disagree")


def validate_run_freshness(article_root: Path, marker: Path | None) -> None:
  if marker is None:
    return
  if not marker.exists():
    raise FileNotFoundError(f"Run marker not found: {marker}")
  threshold = marker.stat().st_mtime
  mandatory: list[Path] = []
  figure_manifest = article_root / "manifests/figure_manifest.tsv"
  table_manifest = article_root / "manifests/table_manifest.tsv"
  mandatory.extend([figure_manifest, table_manifest])
  if figure_manifest.exists():
    figures = pd.read_csv(figure_manifest, sep="\t")
    for record in figures.to_dict("records"):
      base = Path(str(record["base_path"]))
      if not base.is_absolute():
        base = ROOT / base
      source = Path(str(record["source_data_path"]))
      if not source.is_absolute():
        source = ROOT / source
      mandatory.extend([
        source, base.with_suffix(".png"), base.with_suffix(".pdf"),
        base.with_suffix(".svg"), base.with_suffix(".layout_audit.json"),
      ])
  if table_manifest.exists():
    tables = pd.read_csv(table_manifest, sep="\t")
    for value in tables["path"].astype(str):
      path = Path(value)
      mandatory.append(path if path.is_absolute() else ROOT / path)
  mandatory.extend([
    article_root / "tables/qc/candidate_domain_evidence_audit.tsv",
    article_root / "tables/qc/evidence_overlap_registry.tsv",
    article_root / "tables/score_components/rses_onco_score_decomposition.tsv",
    article_root / "tables/robustness/leave_one_domain_out.tsv",
    article_root / "tables/figure_data/figure_source_data_inventory.tsv",
    article_root / "tables/supporting_evidence/supporting_evidence_manifest.tsv",
  ])
  stale = []
  for path in mandatory:
    if not path.exists() or path.stat().st_size == 0:
      stale.append(f"missing_or_empty:{path}")
    elif path.stat().st_mtime + 1e-6 < threshold:
      stale.append(f"stale:{path}")
  if stale:
    raise RuntimeError(
      "Mandatory registered publication outputs were not regenerated in this assets-only run:\n"
      + "\n".join(stale[:150])
    )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument("--run-marker", default=None)
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  marker = resolve_path(args.run_marker) if args.run_marker else None
  validate_audit(article_root)
  validate_scores(article_root)
  validate_overlap(article_root)
  validate_figures(article_root)
  validate_tables(article_root)
  validate_figure_catalog(article_root)
  validate_run_freshness(article_root, marker)
  report = {
    "status": "passed",
    "missing_values_preserved": True,
    "overlap_weight_capped": True,
    "figure_source_tables_complete": True,
    "tcga_frequency_domain_valid": True,
    "figure_5_missingness_explicit": True,
    "fdr_semantics_checked": True,
    "clinical_efficacy_claimed": False,
  }
  report_path = article_root / "manifests/scientific_integrity_validation.json"
  report_path.parent.mkdir(parents=True, exist_ok=True)
  report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
  print("Publication scientific-integrity validation passed.")
  print(f"Report: {report_path}")


if __name__ == "__main__":
  main()
