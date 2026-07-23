#!/usr/bin/env python3
"""Integrate TCGA/GDC promoter methylation into the regulatory microniche domain."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
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
  return result if math.isfinite(result) else None


def explicit_bool(value: object) -> bool:
  if isinstance(value, bool):
    return value
  return str(value).strip().casefold() in {
    "1", "true", "yes", "eligible", "available"
  }


def eligibility_aware_subscore(
  components: dict[str, float | None],
  weights: dict[str, float],
  eligible: set[str],
) -> dict[str, float | int]:
  eligible_weight = float(sum(weights[name] for name in eligible))
  numerator = 0.0
  observed_weight = 0.0
  observed = 0
  for name, weight in weights.items():
    if name not in eligible:
      continue
    value = numeric(components.get(name))
    if value is None:
      continue
    clipped = float(np.clip(value, 0.0, 1.0))
    numerator += weight * clipped
    observed_weight += weight
    observed += 1
  raw = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / eligible_weight if eligible_weight else float("nan")
  adjusted = raw * coverage if math.isfinite(raw) and math.isfinite(coverage) else float("nan")
  return {
    "raw": raw,
    "coverage": coverage,
    "adjusted": adjusted,
    "observed_subcomponents": observed,
    "eligible_subcomponents": len(eligible),
    "observed_weight": observed_weight,
    "eligible_weight": eligible_weight,
  }


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--regulatory-evidence",
    default=(
      "data/processed/regulatory/"
      "expanded_pair_functional_evidence_by_cancer.tsv"
    ),
  )
  parser.add_argument(
    "--methylation-metrics",
    default=(
      "data/processed/epigenetics/methylation/"
      "tcga_nise_methylation_pair_metrics.tsv"
    ),
  )
  parser.add_argument(
    "--output",
    default=(
      "data/processed/regulatory/"
      "expanded_pair_functional_evidence_by_cancer.tsv"
    ),
  )
  parser.add_argument(
    "--status-output",
    default=(
      "data/processed/regulatory/"
      "methylation_regulatory_integration_status.json"
    ),
  )
  args = parser.parse_args()

  evidence = pd.read_csv(
    resolve_path(args.regulatory_evidence),
    sep="\t",
    low_memory=False,
  )
  metrics_path = resolve_path(args.methylation_metrics)
  methylation = (
    pd.read_csv(metrics_path, sep="\t", low_memory=False)
    if metrics_path.exists() and metrics_path.stat().st_size
    else pd.DataFrame(columns=["pair_id", "cancer"])
  )
  required = {"pair_id", "cancer"}
  if not required.issubset(evidence.columns):
    raise ValueError("Regulatory evidence requires pair_id and cancer")
  if not required.issubset(methylation.columns):
    raise ValueError("Methylation metrics require pair_id and cancer")

  methylation_columns = [
    column for column in methylation.columns
    if column not in evidence.columns or column in required
  ]
  merged = evidence.merge(
    methylation[methylation_columns],
    on=["pair_id", "cancer"],
    how="left",
  )

  results: list[dict[str, object]] = []
  for record in merged.to_dict("records"):
    methylation_eligible = explicit_bool(
      record.get("methylation_source_eligible", False)
    )
    eligible = {
      "tf_association_divergence",
      "tf_expression_profile_divergence",
      "promoter_motif_divergence",
    }
    if methylation_eligible:
      eligible.add("promoter_methylation_context")
    components = {
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
    score = eligibility_aware_subscore(
      components,
      REGULATORY_METHYLATION_SUBWEIGHTS,
      eligible,
    )
    updated = dict(record)
    updated.update({
      "regulatory_promoter_methylation_context": components[
        "promoter_methylation_context"
      ],
      "regulatory_methylation_eligible": methylation_eligible,
      "regulatory_network_raw": score["raw"],
      "regulatory_network_coverage": score["coverage"],
      "component_regulatory_network_composite": score["adjusted"],
      "component_regulatory_network": score["adjusted"],
      "regulatory_observed_subcomponents": score[
        "observed_subcomponents"
      ],
      "regulatory_eligible_subcomponents": score[
        "eligible_subcomponents"
      ],
      "regulatory_observed_weight": score["observed_weight"],
      "regulatory_eligible_weight": score["eligible_weight"],
      "regulatory_network_method": (
        "DoRothEA_TF_target_plus_TF_expression_consistency_plus_"
        "JASPAR_promoter_motif_prediction_plus_"
        "TCGA_GDC_promoter_methylation_context"
      ),
      "regulatory_layer_version": (
        "wgcna-promoter-methylation-regulatory-v3"
      ),
      "methylation_direct_silencing_claim": False,
    })
    results.append(updated)

  output_frame = pd.DataFrame(results)
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  output_frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)

  status = {
    "version": "wgcna-promoter-methylation-regulatory-v3",
    "regulatory_subweights": REGULATORY_METHYLATION_SUBWEIGHTS,
    "top_level_regulatory_weight_changed": False,
    "methylation_source": "NCI_GDC_via_UCSC_Xena_GDC_Hub",
    "repbase_used": False,
    "methylation_direct_silencing_claim": False,
    "output_rows": len(output_frame),
    "methylation_eligible_rows": int(
      output_frame["regulatory_methylation_eligible"]
      .fillna(False)
      .sum()
    ),
    "methylation_observed_rows": int(
      output_frame["regulatory_promoter_methylation_context"]
      .notna()
      .sum()
    ),
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
  }
  status_path = resolve_path(args.status_output)
  status_path.parent.mkdir(parents=True, exist_ok=True)
  status_path.write_text(
    json.dumps(status, indent=2, sort_keys=True),
    encoding="utf-8",
  )
  print(
    "Integrated promoter methylation into regulatory evidence: "
    f"{output} ({len(output_frame):,} rows)"
  )


if __name__ == "__main__":
  main()
