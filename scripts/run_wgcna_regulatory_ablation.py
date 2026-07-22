#!/usr/bin/env python3
"""Ablate WGCNA, pairwise expression, TF and promoter regulatory subcomponents."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  FUNCTIONAL_MICRONICHE_WEIGHTS,
  coverage_aware_score,
  expanded_onco_score,
  functional_microniche_score,
)

ROOT = Path(__file__).resolve().parents[1]
EXPRESSION_WEIGHTS = {
  "pairwise": 0.50,
  "wgcna": 0.50,
}
REGULATORY_WEIGHTS = {
  "dorothea": 0.40,
  "tf_expression": 0.35,
  "promoter": 0.25,
}
SCENARIOS = (
  "baseline",
  "without_wgcna",
  "without_pairwise_expression",
  "without_dorothea_regulator_sets",
  "without_tf_expression_consistency",
  "without_promoter_motifs",
  "without_regulatory_domain",
)


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def numeric(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if np.isfinite(result) else None


def explicit_bool(value: object) -> bool:
  if isinstance(value, bool):
    return value
  return str(value).strip().casefold() in {
    "1",
    "true",
    "yes",
    "eligible",
  }


def subscore(
  values: dict[str, float | None],
  weights: dict[str, float],
) -> float | None:
  result = coverage_aware_score(values, weights)
  return (
    result.adjusted_score
    if np.isfinite(result.adjusted_score)
    else None
  )


def score_scenario(
  record: dict[str, object],
  scenario: str,
) -> float | None:
  expression_values = {
    "pairwise": numeric(record.get("pairwise_expression_context")),
    "wgcna": numeric(record.get("wgcna_expression_network")),
  }
  regulatory_values = {
    "dorothea": numeric(
      record.get("regulatory_tf_association_divergence")
    ),
    "tf_expression": numeric(
      record.get("regulatory_tf_expression_profile_divergence")
    ),
    "promoter": numeric(
      record.get("regulatory_promoter_motif_divergence")
    ),
  }
  if scenario == "without_wgcna":
    expression_values["wgcna"] = None
  elif scenario == "without_pairwise_expression":
    expression_values["pairwise"] = None
  elif scenario == "without_dorothea_regulator_sets":
    regulatory_values["dorothea"] = None
  elif scenario == "without_tf_expression_consistency":
    regulatory_values["tf_expression"] = None
  elif scenario == "without_promoter_motifs":
    regulatory_values["promoter"] = None

  microniche_components = {
    domain: numeric(record.get(f"microniche_{domain}"))
    for domain in FUNCTIONAL_MICRONICHE_WEIGHTS
  }
  microniche_components["expression_context"] = subscore(
    expression_values,
    EXPRESSION_WEIGHTS,
  )
  microniche_components["regulatory_network"] = (
    None
    if scenario == "without_regulatory_domain"
    else subscore(regulatory_values, REGULATORY_WEIGHTS)
  )
  eligible_microniche = {
    domain
    for domain in FUNCTIONAL_MICRONICHE_WEIGHTS
    if explicit_bool(record.get(f"eligible_microniche_{domain}", True))
  }
  if scenario == "without_regulatory_domain":
    eligible_microniche.discard("regulatory_network")
  microniche = functional_microniche_score(
    microniche_components,
    eligible_domains=eligible_microniche,
  )

  onco_components = {
    domain: numeric(record.get(f"component_{domain}"))
    for domain in EXPANDED_ONCO_WEIGHTS
  }
  onco_components["functional_microniche"] = (
    microniche.adjusted_score
    if np.isfinite(microniche.adjusted_score)
    else None
  )
  eligible_onco = {
    domain
    for domain in EXPANDED_ONCO_WEIGHTS
    if explicit_bool(record.get(f"eligible_component_{domain}", True))
  }
  onco = expanded_onco_score(
    onco_components,
    eligible_domains=eligible_onco,
  )
  return (
    onco.adjusted_score
    if np.isfinite(onco.adjusted_score)
    else None
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument("--article-root", default="article_outputs")
  parser.add_argument("--top-k", type=int, default=20)
  args = parser.parse_args()

  ranking = pd.read_csv(
    resolve_path(args.ranking),
    sep="\t",
    low_memory=False,
  )
  required = {
    "cancer",
    "pair_id",
    "coverage_adjusted_rses",
    "pairwise_expression_context",
    "wgcna_expression_network",
    "regulatory_tf_association_divergence",
    "regulatory_tf_expression_profile_divergence",
    "regulatory_promoter_motif_divergence",
  }
  missing = sorted(required - set(ranking.columns))
  if missing:
    raise ValueError(
      f"Ranking lacks WGCNA/regulatory ablation fields: {missing}"
    )

  scenario_rows = []
  summary_rows = []
  baseline = ranking[[
    "cancer",
    "pair_id",
    "coverage_adjusted_rses",
  ]].copy()
  for scenario in SCENARIOS:
    values = [
      score_scenario(record, scenario)
      for record in ranking.to_dict("records")
    ]
    frame = ranking[[
      column
      for column in (
        "cancer",
        "pair_id",
        "source_class",
        "score_comparability_group",
      )
      if column in ranking
    ]].copy()
    frame["scenario"] = scenario
    frame["scenario_adjusted_score"] = values
    frame["scenario_rank"] = frame.groupby(
      [
        column
        for column in (
          "cancer",
          "score_comparability_group",
        )
        if column in frame
      ],
      dropna=False,
    )["scenario_adjusted_score"].rank(
      method="min",
      ascending=False,
    )
    scenario_rows.append(frame)

    merged = baseline.merge(
      frame[["cancer", "pair_id", "scenario_adjusted_score"]],
      on=["cancer", "pair_id"],
      how="inner",
    )
    for cancer, group in merged.groupby("cancer", dropna=False):
      valid = group[[
        "coverage_adjusted_rses",
        "scenario_adjusted_score",
      ]].dropna()
      rho = (
        float(
          spearmanr(
            valid.iloc[:, 0],
            valid.iloc[:, 1],
          ).statistic
        )
        if len(valid) >= 3
        else np.nan
      )
      base_top = set(
        group.sort_values(
          "coverage_adjusted_rses",
          ascending=False,
        ).head(args.top_k)["pair_id"].astype(str)
      )
      scenario_top = set(
        group.sort_values(
          "scenario_adjusted_score",
          ascending=False,
        ).head(args.top_k)["pair_id"].astype(str)
      )
      union = base_top | scenario_top
      summary_rows.append({
        "scenario": scenario,
        "cancer": cancer,
        "n_pairs": len(valid),
        "spearman_rho": rho,
        "top_k": args.top_k,
        "top_k_overlap": len(base_top & scenario_top),
        "top_k_jaccard": (
          len(base_top & scenario_top) / len(union)
          if union
          else np.nan
        ),
      })

  root = resolve_path(args.article_root)
  output_dir = root / "tables/robustness"
  output_dir.mkdir(parents=True, exist_ok=True)
  scenarios = pd.concat(scenario_rows, ignore_index=True)
  summary = pd.DataFrame(summary_rows)
  scenarios.to_csv(
    output_dir / "wgcna_regulatory_ablation_scores.tsv",
    sep="\t",
    index=False,
  )
  summary.to_csv(
    output_dir / "wgcna_regulatory_ablation_summary.tsv",
    sep="\t",
    index=False,
  )
  print("WGCNA/regulatory ablation analysis passed.")
  print(summary.to_string(index=False))


if __name__ == "__main__":
  main()
