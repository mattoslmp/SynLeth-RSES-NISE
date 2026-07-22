#!/usr/bin/env python3
"""Run deterministic RSES-Onco robustness and missing-data sensitivity analyses.

This script does not create new biological evidence. It recomputes rankings from the
already observed component matrix under leave-one-domain-out and controlled weight
perturbations. Missing values remain missing. Non-eligible domains do not enter the
observed denominator.
"""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from rses_onco.expanded import EXPANDED_ONCO_WEIGHTS

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


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


def numeric_series(frame: pd.DataFrame, column: str) -> pd.Series:
  if column not in frame:
    return pd.Series(np.nan, index=frame.index, dtype=float)
  return pd.to_numeric(frame[column], errors="coerce").clip(lower=0, upper=1)


def recompute(frame: pd.DataFrame, weights: dict[str, float]) -> pd.DataFrame:
  total_weight = float(sum(weights.values()))
  numerator = pd.Series(0.0, index=frame.index)
  observed_weight = pd.Series(0.0, index=frame.index)
  observed_domains = pd.Series(0, index=frame.index, dtype=int)
  for domain, weight in weights.items():
    values = numeric_series(frame, f"component_{domain}")
    observed = values.notna()
    numerator.loc[observed] += float(weight) * values.loc[observed]
    observed_weight.loc[observed] += float(weight)
    observed_domains.loc[observed] += 1
  raw = numerator / observed_weight.replace(0, np.nan)
  coverage = observed_weight / total_weight if total_weight else np.nan
  adjusted = raw * coverage
  return pd.DataFrame({
    "recomputed_raw_score": raw,
    "recomputed_coverage": coverage,
    "recomputed_adjusted_score": adjusted,
    "recomputed_observed_domains": observed_domains,
  })


def rank_within_cancer(frame: pd.DataFrame, score: str) -> pd.Series:
  return frame.groupby("cancer", dropna=False)[score].rank(method="min", ascending=False)


def top_k_set(frame: pd.DataFrame, score: str, cancer: str, k: int) -> set[str]:
  subset = frame.loc[frame["cancer"].astype(str).eq(cancer)].sort_values(score, ascending=False).head(k)
  return set(subset["pair_id"].astype(str))


def rank_correlation(base: pd.DataFrame, alternative: pd.DataFrame, score: str) -> pd.DataFrame:
  rows = []
  merged = base[["cancer", "pair_id", "coverage_adjusted_rses"]].merge(
    alternative[["cancer", "pair_id", score]],
    on=["cancer", "pair_id"],
    how="inner",
  )
  for cancer, group in merged.groupby("cancer", dropna=False):
    valid = group[["coverage_adjusted_rses", score]].dropna()
    rho = float(spearmanr(valid.iloc[:, 0], valid.iloc[:, 1]).statistic) if len(valid) >= 3 else np.nan
    rows.append({"cancer": cancer, "n_pairs": len(valid), "spearman_rho": rho})
  return pd.DataFrame(rows)


def controlled_weight_scenarios(delta: float) -> list[tuple[str, dict[str, float]]]:
  scenarios = [("baseline", dict(EXPANDED_ONCO_WEIGHTS))]
  for domain in EXPANDED_ONCO_WEIGHTS:
    for direction, multiplier in (("down", 1.0 - delta), ("up", 1.0 + delta)):
      weights = dict(EXPANDED_ONCO_WEIGHTS)
      weights[domain] *= multiplier
      total = sum(weights.values())
      weights = {key: value / total for key, value in weights.items()}
      scenarios.append((f"{domain}_{direction}_{int(delta * 100)}pct", weights))
  return scenarios


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--ranking", default="results/expanded_26Q1/full/expanded_rses_onco.tsv")
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--top-k", type=int, default=20)
  parser.add_argument("--weight-delta", type=float, default=0.20)
  args = parser.parse_args()

  ranking_path = resolve_path(args.ranking)
  if not ranking_path.exists() or ranking_path.stat().st_size == 0:
    raise FileNotFoundError(f"Missing or empty ranking: {ranking_path}")
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  required = {"cancer", "pair_id", "coverage_adjusted_rses"}
  missing = sorted(required - set(ranking.columns))
  if missing:
    raise ValueError(f"Ranking missing required columns: {missing}")

  output_root = resolve_path(args.output_root)
  robustness_dir = output_root / "tables" / "robustness"
  supplementary_dir = output_root / "tables" / "supplementary"

  baseline = recompute(ranking, dict(EXPANDED_ONCO_WEIGHTS))
  baseline_rows = pd.concat([
    ranking[[column for column in (
      "cancer", "pair_id", "source_class", "lost_feature", "analysis_lost_gene",
      "lost_gene", "analysis_target_gene", "target_gene", "rses_onco",
      "evidence_coverage", "coverage_adjusted_rses",
    ) if column in ranking]].reset_index(drop=True),
    baseline.reset_index(drop=True),
  ], axis=1)
  baseline_rows["pipeline_rank"] = rank_within_cancer(baseline_rows, "coverage_adjusted_rses")
  baseline_rows["recomputed_rank"] = rank_within_cancer(baseline_rows, "recomputed_adjusted_score")

  valid = baseline_rows[["coverage_adjusted_rses", "recomputed_adjusted_score"]].dropna()
  if not valid.empty and not np.allclose(
    valid["coverage_adjusted_rses"], valid["recomputed_adjusted_score"], atol=1e-10, rtol=1e-8,
  ):
    raise RuntimeError("Baseline recomputation does not reproduce pipeline adjusted scores")

  leave_one_rows = []
  stability_rows = []
  correlations = []
  cancers = sorted(ranking["cancer"].dropna().astype(str).unique())
  for excluded in EXPANDED_ONCO_WEIGHTS:
    weights = {key: value for key, value in EXPANDED_ONCO_WEIGHTS.items() if key != excluded}
    scores = recompute(ranking, weights)
    alternative = pd.concat([
      ranking[["cancer", "pair_id", "source_class"]].reset_index(drop=True),
      scores.reset_index(drop=True),
    ], axis=1)
    alternative["excluded_domain"] = excluded
    alternative["alternative_rank"] = rank_within_cancer(alternative, "recomputed_adjusted_score")
    leave_one_rows.append(alternative)
    corr = rank_correlation(ranking, alternative, "recomputed_adjusted_score")
    corr["scenario"] = f"leave_out_{excluded}"
    correlations.append(corr)
    for cancer in cancers:
      baseline_set = top_k_set(ranking, "coverage_adjusted_rses", cancer, args.top_k)
      alternative_set = top_k_set(alternative, "recomputed_adjusted_score", cancer, args.top_k)
      union = baseline_set | alternative_set
      stability_rows.append({
        "scenario": f"leave_out_{excluded}",
        "cancer": cancer,
        "top_k": args.top_k,
        "overlap_n": len(baseline_set & alternative_set),
        "jaccard": len(baseline_set & alternative_set) / len(union) if union else np.nan,
        "baseline_top_pairs": ";".join(sorted(baseline_set)),
        "alternative_top_pairs": ";".join(sorted(alternative_set)),
      })

  weight_rows = []
  for scenario, weights in controlled_weight_scenarios(args.weight_delta):
    scores = recompute(ranking, weights)
    alternative = pd.concat([
      ranking[["cancer", "pair_id", "source_class"]].reset_index(drop=True),
      scores.reset_index(drop=True),
    ], axis=1)
    alternative["scenario"] = scenario
    alternative["weights"] = ";".join(f"{key}={value:.8f}" for key, value in weights.items())
    alternative["alternative_rank"] = rank_within_cancer(alternative, "recomputed_adjusted_score")
    weight_rows.append(alternative)
    corr = rank_correlation(ranking, alternative, "recomputed_adjusted_score")
    corr["scenario"] = f"weights_{scenario}"
    correlations.append(corr)
    for cancer in cancers:
      baseline_set = top_k_set(ranking, "coverage_adjusted_rses", cancer, args.top_k)
      alternative_set = top_k_set(alternative, "recomputed_adjusted_score", cancer, args.top_k)
      union = baseline_set | alternative_set
      stability_rows.append({
        "scenario": f"weights_{scenario}",
        "cancer": cancer,
        "top_k": args.top_k,
        "overlap_n": len(baseline_set & alternative_set),
        "jaccard": len(baseline_set & alternative_set) / len(union) if union else np.nan,
        "baseline_top_pairs": ";".join(sorted(baseline_set)),
        "alternative_top_pairs": ";".join(sorted(alternative_set)),
      })

  leave_one = pd.concat(leave_one_rows, ignore_index=True)
  weight_sensitivity = pd.concat(weight_rows, ignore_index=True)
  correlations_frame = pd.concat(correlations, ignore_index=True)
  stability = pd.DataFrame(stability_rows)

  missingness_rows = []
  for record in ranking.to_dict("records"):
    observed_values = [
      pd.to_numeric(pd.Series([record.get(f"component_{domain}")]), errors="coerce").iloc[0]
      for domain in EXPANDED_ONCO_WEIGHTS
    ]
    n_missing = int(sum(pd.isna(value) for value in observed_values))
    missingness_rows.append({
      "cancer": record.get("cancer"),
      "pair_id": record.get("pair_id"),
      "source_class": record.get("source_class"),
      "observed_domains": len(EXPANDED_ONCO_WEIGHTS) - n_missing,
      "missing_domains": n_missing,
      "pipeline_adjusted_score": record.get("coverage_adjusted_rses"),
      "pipeline_coverage": record.get("evidence_coverage"),
      "missing_data_interpretation": (
        "Missing components omitted from the observed denominator and penalized through coverage; no zero imputation."
      ),
    })
  missingness_sensitivity = pd.DataFrame(missingness_rows)

  outputs = [
    atomic_tsv(baseline_rows, robustness_dir / "raw_vs_coverage_adjusted.tsv"),
    atomic_tsv(leave_one, robustness_dir / "leave_one_domain_out.tsv"),
    atomic_tsv(weight_sensitivity, robustness_dir / "controlled_weight_perturbation.tsv"),
    atomic_tsv(correlations_frame, robustness_dir / "ranking_correlations.tsv"),
    atomic_tsv(stability, robustness_dir / "top_k_stability.tsv"),
    atomic_tsv(missingness_sensitivity, robustness_dir / "missing_data_sensitivity.tsv"),
    atomic_tsv(stability, supplementary_dir / "Table_S23_rses_ranking_stability.tsv"),
    atomic_tsv(leave_one, supplementary_dir / "Table_S24_leave_one_domain_out.tsv"),
    atomic_tsv(weight_sensitivity, supplementary_dir / "Table_S25_weight_sensitivity.tsv"),
  ]
  summary = {
    "ranking_rows": len(ranking),
    "leave_one_domain_out_rows": len(leave_one),
    "weight_sensitivity_rows": len(weight_sensitivity),
    "correlation_rows": len(correlations_frame),
    "stability_rows": len(stability),
    "top_k": args.top_k,
    "weight_delta": args.weight_delta,
    "rules": {
      "zero_imputation": False,
      "missing_domains_in_observed_denominator": False,
      "deterministic": True,
      "clinical_claims": False,
    },
    "outputs": [str(path) for path in outputs],
  }
  atomic_json(summary, robustness_dir / "robustness_summary.json")
  for path in outputs:
    if not path.exists() or path.stat().st_size == 0:
      raise RuntimeError(f"Missing or empty robustness output: {path}")
    print(f"Wrote {path}")
  print("RSES-Onco robustness analyses passed.")


if __name__ == "__main__":
  main()
