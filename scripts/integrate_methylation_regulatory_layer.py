#!/usr/bin/env python3
"""Integrate promoter methylation into the existing regulatory microniche domain.

Methylation does not receive a new independent RSES-Onco domain weight. Instead,
it shares the existing regulatory-network domain with DoRothEA TF associations,
TF-expression consistency and JASPAR/FIMO promoter-motif support. Missing
methylation remains missing and lowers regulatory subcoverage; it is never
converted to biological zero.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from rses_onco.depmap import normalize_model_id_column
from rses_onco.expanded import coverage_aware_score
from rses_onco.methylation import (
  build_methylation_pair_metrics,
  read_promoter_methylation,
)
from rses_onco.utils import canonical_gene_name


REGULATORY_METHYLATION_SUBWEIGHTS = {
  "tf_association_divergence": 0.32,
  "tf_expression_profile_divergence": 0.28,
  "promoter_motif_divergence": 0.20,
  "promoter_methylation_context": 0.20,
}


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if np.isfinite(result) else None


def read_copy_number(path: Path) -> pd.DataFrame:
  frame = normalize_model_id_column(
    pd.read_csv(path, low_memory=False),
    "copy_number",
  )
  metadata = {
    "ModelID", "ProfileID", "PROFILEID", "is_default_entry",
    "IsDefaultEntryForModel", "IsDefaultEntryForMC", "ModelConditionID",
    "SequencingID",
  }
  frame = frame.rename(columns={
    column: canonical_gene_name(column)
    for column in frame.columns
    if column not in metadata
  })
  if frame["ModelID"].duplicated().any():
    for column in ("IsDefaultEntryForModel", "is_default_entry"):
      if column in frame:
        flag = frame[column].astype(str).str.casefold().isin({"1", "true", "yes"})
        selected = frame.loc[flag].copy()
        if not selected.empty and not selected["ModelID"].duplicated().any():
          frame = selected
          break
    else:
      frame = frame.groupby("ModelID", as_index=False).median(numeric_only=True)
  return frame


def missing_metrics(
  pairs: pd.DataFrame,
  cancers: tuple[str, ...],
  reason: str,
) -> pd.DataFrame:
  rows = []
  for cancer in cancers:
    for record in pairs.to_dict("records"):
      rows.append({
        "pair_id": record["pair_id"],
        "cancer": cancer,
        "lost_gene": record["lost_gene"],
        "target_gene": record["target_gene"],
        "methylation_profile_n_models": 0,
        "methylation_pair_spearman_rho": None,
        "methylation_pair_median_absolute_difference": None,
        "methylation_pair_profile_divergence": None,
        "methylation_n_loss": 0,
        "methylation_n_intact": 0,
        "methylation_target_median_loss": None,
        "methylation_target_median_intact": None,
        "methylation_target_delta_loss_minus_intact": None,
        "methylation_target_hypomethylation_support": None,
        "methylation_p_value": None,
        "methylation_q_value_bh": None,
        "methylation_q_value_bh_within_cancer": None,
        "methylation_raw": None,
        "methylation_coverage": 0.0,
        "component_promoter_methylation_context": None,
        "methylation_observed_subcomponents": 0,
        "methylation_lost_promoter_feature_count": 0,
        "methylation_target_promoter_feature_count": 0,
        "methylation_absence_reason": reason,
        "methylation_evidence_type": (
          "CCLE_RRBS_weighted_1kb_upstream_TSS_promoter_methylation"
        ),
        "methylation_interpretation": "eligible_missing_evidence_not_negative_evidence",
      })
  return pd.DataFrame(rows)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--methylation", default="")
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument("--candidates", required=True)
  parser.add_argument("--input", required=True)
  parser.add_argument("--output", required=True)
  parser.add_argument(
    "--metrics-output",
    default="data/processed/regulatory/promoter_methylation_pair_metrics.tsv",
  )
  parser.add_argument(
    "--status-output",
    default="data/processed/regulatory/promoter_methylation_status.json",
  )
  parser.add_argument("--loss-threshold", type=float, default=0.30)
  parser.add_argument("--min-group-size", type=int, default=3)
  parser.add_argument("--difference-saturation", type=float, default=0.25)
  args = parser.parse_args()

  evidence = pd.read_csv(resolve_path(args.input), sep="\t", low_memory=False)
  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t", low_memory=False)
  pairs = candidates[["pair_id", "lost_gene", "target_gene"]].copy()
  pairs["lost_gene"] = pairs["lost_gene"].map(canonical_gene_name)
  pairs["target_gene"] = pairs["target_gene"].map(canonical_gene_name)
  pairs = pairs.loc[
    pairs["pair_id"].notna()
    & pairs["lost_gene"].ne("")
    & pairs["target_gene"].ne("")
  ].drop_duplicates("pair_id")
  models = normalize_model_id_column(
    pd.read_csv(resolve_path(args.models), low_memory=False),
    "Model.csv",
  )
  copy_number = read_copy_number(resolve_path(args.copy_number))
  cancers = tuple(sorted(set(evidence["cancer"].dropna().astype(str))))
  if not cancers:
    cancers = ("colon", "stomach", "lung")

  methylation_path = resolve_path(args.methylation) if args.methylation else None
  source_available = bool(
    methylation_path is not None
    and methylation_path.exists()
    and methylation_path.stat().st_size > 0
  )
  diagnostics: dict[str, object]
  if source_available:
    methylation = read_promoter_methylation(methylation_path, models)
    metrics = build_methylation_pair_metrics(
      methylation,
      copy_number,
      models,
      pairs,
      cancers=cancers,
      loss_threshold=args.loss_threshold,
      min_group_size=args.min_group_size,
      difference_saturation=args.difference_saturation,
    )
    diagnostics = dict(methylation.diagnostics)
    diagnostics["source_status"] = "available"
  else:
    metrics = missing_metrics(
      pairs,
      cancers,
      "methylation_source_not_provided_or_file_absent",
    )
    diagnostics = {
      "source_path": str(methylation_path or ""),
      "source_status": "eligible_source_missing",
      "missing_data_rule": "preserved_as_NA_and_lowers_regulatory_subcoverage",
    }

  metrics_output = resolve_path(args.metrics_output)
  metrics_output.parent.mkdir(parents=True, exist_ok=True)
  metrics.to_csv(metrics_output, sep="\t", index=False)

  methylation_columns = [
    column
    for column in metrics.columns
    if column not in {"lost_gene", "target_gene"}
  ]
  duplicated = [
    column
    for column in methylation_columns
    if column in evidence.columns and column not in {"pair_id", "cancer"}
  ]
  evidence = evidence.drop(columns=duplicated, errors="ignore")
  enriched = evidence.merge(
    metrics[methylation_columns],
    on=["pair_id", "cancer"],
    how="left",
  )

  rows = []
  for record in enriched.to_dict("records"):
    values = {
      "tf_association_divergence": numeric(
        record.get("regulatory_tf_association_divergence")
      ),
      "tf_expression_profile_divergence": numeric(
        record.get("regulatory_tf_expression_profile_divergence")
      ),
      "promoter_motif_divergence": numeric(
        record.get("regulatory_promoter_motif_divergence")
      ),
      "promoter_methylation_context": numeric(
        record.get("component_promoter_methylation_context")
      ),
    }
    score = coverage_aware_score(values, REGULATORY_METHYLATION_SUBWEIGHTS)
    updated = dict(record)
    updated.update({
      "regulatory_methylation_context": values["promoter_methylation_context"],
      "regulatory_network_raw": score.observed_score,
      "regulatory_network_coverage": score.coverage,
      "component_regulatory_network_composite": score.adjusted_score,
      "component_regulatory_network": score.adjusted_score,
      "regulatory_observed_subcomponents": score.n_domains,
      "regulatory_network_method": (
        "DoRothEA_TF_target_plus_TF_expression_consistency_plus_"
        "JASPAR_promoter_motif_prediction_plus_CCLE_promoter_methylation"
      ),
      "regulatory_layer_version": "wgcna-promoter-methylation-regulatory-v3",
      "methylation_source_status": diagnostics["source_status"],
      "methylation_direct_causal_claim": False,
    })
    rows.append(updated)
  result = pd.DataFrame(rows)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  result.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)

  status = {
    "version": "wgcna-promoter-methylation-regulatory-v3",
    "source_status": diagnostics["source_status"],
    "source_path": diagnostics.get("source_path", ""),
    "regulatory_subweights": REGULATORY_METHYLATION_SUBWEIGHTS,
    "methylation_subweights": {
      "pair_profile_divergence": 0.50,
      "conditional_target_hypomethylation": 0.50,
    },
    "loss_threshold": args.loss_threshold,
    "min_group_size": args.min_group_size,
    "difference_saturation": args.difference_saturation,
    "direct_causal_methylation_claim": False,
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    **diagnostics,
  }
  status_path = resolve_path(args.status_output)
  status_path.parent.mkdir(parents=True, exist_ok=True)
  status_path.write_text(
    json.dumps(status, indent=2, sort_keys=True),
    encoding="utf-8",
  )
  print(f"Methylation metrics: {metrics_output} ({len(metrics):,} rows)")
  print(
    "Integrated methylation-aware regulatory evidence: "
    f"{output} ({len(result):,} rows)"
  )
  print(f"Methylation source status: {diagnostics['source_status']}")


if __name__ == "__main__":
  main()
