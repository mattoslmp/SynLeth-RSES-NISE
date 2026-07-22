#!/usr/bin/env python3
"""Run deterministic eligibility-aware RSES-Onco robustness analyses.

The script recomputes rankings from observed components under leave-one-domain-out
and controlled weight perturbations. Missing values remain missing. Non-eligible
domains enter neither the observed nor eligible denominator. Rank stability is
reported within cancer and score-comparability groups so composite-event hypotheses
are not silently compared as if they had the same evaluable domain universe as
gene-pair hypotheses.
"""
from __future__ import annotations

import argparse
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
  temporary.write_text(
    json.dumps(payload, indent=2, sort_keys=True),
    encoding="utf-8",
  )
  temporary.replace(path)
  return path


def numeric_series(frame: pd.DataFrame, column: str) -> pd.Series:
  if column not in frame:
    return pd.Series(np.nan, index=frame.index, dtype=float)
  return pd.to_numeric(frame[column], errors="coerce").clip(lower=0, upper=1)


def boolean_series(frame: pd.DataFrame, column: str) -> pd.Series:
  if column not in frame:
    return pd.Series(True, index=frame.index, dtype=bool)
  values = frame[column]
  if values.dtype == bool:
    return values.fillna(False)
  text = values.fillna("").astype(str).str.strip().str.casefold()
  return text.isin({"1", "true", "yes", "eligible"})


def recompute(frame: pd.DataFrame, weights: dict[str, float]) -> pd.DataFrame:
  numerator = pd.Series(0.0, index=frame.index, dtype=float)
  observed_weight = pd.Series(0.0, index=frame.index, dtype=float)
  eligible_weight = pd.Series(0.0, index=frame.index, dtype=float)
  observed_domains = pd.Series(0, index=frame.index, dtype=int)
  eligible_domains = pd.Series(0, index=frame.index, dtype=int)

  for domain, weight in weights.items():
    values = numeric_series(frame, f"component_{domain}")
    eligible = boolean_series(frame, f"eligible_component_{domain}")
    observed = eligible & values.notna()
    eligible_weight.loc[eligible] += float(weight)
    eligible_domains.loc[eligible] += 1
    numerator.loc[observed] += float(weight) * values.loc[observed]
    observed_weight.loc[observed] += float(weight)
    observed_domains.loc[observed] += 1

  raw = numerator / observed_weight.replace(0, np.nan)
  coverage = observed_weight / eligible_weight.replace(0, np.nan)
  adjusted = numerator / eligible_weight.replace(0, np.nan)
  return pd.DataFrame({
    "recomputed_weighted_numerator": numerator,
    "recomputed_observed_weight": observed_weight,
    "recomputed_eligible_weight": eligible_weight,
    "recomputed_raw_score": raw,
    "recomputed_coverage": coverage,
    "recomputed_adjusted_score": adjusted,
    "recomputed_observed_domains": observed_domains,
    "recomputed_eligible_domains": eligible_domains,
  })


def context_columns(frame: pd.DataFrame) -> list[str]:
  columns = ["cancer"]
  if "score_comparability_group" in frame:
    columns.append("score_comparability_group")
  return columns


def rank_within_context(frame: pd.DataFrame, score: str) -> pd.Series:
  return frame.groupby(context_columns(frame), dropna=False)[score].rank(
    method="min",
    ascending=False,
  )


def top_k_set(
  frame: pd.DataFrame,
  score: str,
  cancer: str,
  comparability_group: str,
  k: int,
) -> set[str]:
  mask = frame["cancer"].astype(str).eq(cancer)
  if "score_comparability_group" in frame:
    mask &= frame["score_comparability_group"].astype(str).eq(
      comparability_group
    )
  subset = frame.loc[mask].sort_values(score, ascending=False).head(k)
  return set(subset["pair_id"].astype(str))


def rank_correlation(
  base: pd.DataFrame,
  alternative: pd.DataFrame,
  score: str,
  group_columns: list[str],
  scope: str,
) -> pd.DataFrame:
  rows = []
  merge_columns = ["cancer", "pair_id"]
  if "score_comparability_group" in base and "score_comparability_group" in alternative:
    merge_columns.append("score_comparability_group")
  carry = list(dict.fromkeys([*merge_columns, *group_columns]))
  merged = base[[*carry, "coverage_adjusted_rses"]].merge(
    alternative[[*merge_columns, score]],
    on=merge_columns,
    how="inner",
  )
  for values, group in merged.groupby(group_columns, dropna=False):
    if not isinstance(values, tuple):
      values = (values,)
    valid = group[["coverage_adjusted_rses", score]].dropna()
    rho = (
      float(spearmanr(valid.iloc[:, 0], valid.iloc[:, 1]).statistic)
      if len(valid) >= 3
      else np.nan
    )
    row = dict(zip(group_columns, values))
    row.update({
      "comparison_scope": scope,
      "n_pairs": len(valid),
      "spearman_rho": rho,
    })
    rows.append(row)
  return pd.DataFrame(rows)


def controlled_weight_scenarios(
  delta: float,
) -> list[tuple[str, dict[str, float]]]:
  scenarios = [("baseline", dict(EXPANDED_ONCO_WEIGHTS))]
  for domain in EXPANDED_ONCO_WEIGHTS:
    for direction, multiplier in (
      ("down", 1.0 - delta),
      ("up", 1.0 + delta),
    ):
      weights = dict(EXPANDED_ONCO_WEIGHTS)
      weights[domain] *= multiplier
      total = sum(weights.values())
      weights = {key: value / total for key, value in weights.items()}
      scenarios.append(
        (f"{domain}_{direction}_{int(delta * 100)}pct", weights)
      )
  return scenarios


def base_columns(ranking: pd.DataFrame) -> list[str]:
  return [
    column
    for column in (
      "cancer",
      "pair_id",
      "source_class",
      "score_comparability_group",
      "hypothesis_type",
      "lost_feature",
      "analysis_lost_gene",
      "lost_gene",
      "analysis_target_gene",
      "target_gene",
      "rses_onco",
      "evidence_coverage",
      "coverage_adjusted_rses",
      "n_domains",
      "eligible_domains",
      *[f"eligible_component_{domain}" for domain in EXPANDED_ONCO_WEIGHTS],
    )
    if column in ranking
  ]


def add_correlations(
  correlations: list[pd.DataFrame],
  ranking: pd.DataFrame,
  alternative: pd.DataFrame,
  scenario: str,
) -> None:
  comparability_columns = ["cancer"]
  if "score_comparability_group" in ranking:
    comparability_columns.append("score_comparability_group")
  frame = rank_correlation(
    ranking,
    alternative,
    "recomputed_adjusted_score",
    comparability_columns,
    "cancer_comparability_group",
  )
  frame["scenario"] = scenario
  correlations.append(frame)

  if "source_class" in ranking:
    class_columns = [*comparability_columns, "source_class"]
    by_class = rank_correlation(
      ranking,
      alternative,
      "recomputed_adjusted_score",
      class_columns,
      "cancer_comparability_group_mechanistic_class",
    )
    by_class["scenario"] = scenario
    correlations.append(by_class)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--top-k", type=int, default=20)
  parser.add_argument("--weight-delta", type=float, default=0.20)
  args = parser.parse_args()

  ranking_path = resolve_path(args.ranking)
  if not ranking_path.exists() or ranking_path.stat().st_size == 0:
    raise FileNotFoundError(f"Missing or empty ranking: {ranking_path}")
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  required = {
    "cancer",
    "pair_id",
    "coverage_adjusted_rses",
    "score_comparability_group",
    "scoring_semantics_version",
    *[f"eligible_component_{domain}" for domain in EXPANDED_ONCO_WEIGHTS],
  }
  missing = sorted(required - set(ranking.columns))
  if missing:
    raise ValueError(
      "Ranking lacks eligibility-aware fields; rerun the expanded scoring stage. "
      f"Missing: {missing}"
    )
  versions = set(
    ranking["scoring_semantics_version"].dropna().astype(str)
  )
  if versions != {"eligibility-aware-v1"}:
    raise ValueError(
      "Unsupported scoring semantics for robustness analyses: "
      f"{sorted(versions)}"
    )

  output_root = resolve_path(args.output_root)
  robustness_dir = output_root / "tables" / "robustness"
  supplementary_dir = output_root / "tables" / "supplementary"

  baseline = recompute(ranking, dict(EXPANDED_ONCO_WEIGHTS))
  baseline_rows = pd.concat(
    [
      ranking[base_columns(ranking)].reset_index(drop=True),
      baseline.reset_index(drop=True),
    ],
    axis=1,
  )
  baseline_rows["pipeline_rank_within_comparability_group"] = (
    rank_within_context(baseline_rows, "coverage_adjusted_rses")
  )
  baseline_rows["recomputed_rank_within_comparability_group"] = (
    rank_within_context(baseline_rows, "recomputed_adjusted_score")
  )

  valid = baseline_rows[
    ["coverage_adjusted_rses", "recomputed_adjusted_score"]
  ].dropna()
  if not valid.empty and not np.allclose(
    valid["coverage_adjusted_rses"],
    valid["recomputed_adjusted_score"],
    atol=1e-10,
    rtol=1e-8,
  ):
    differences = (
      valid["coverage_adjusted_rses"]
      - valid["recomputed_adjusted_score"]
    ).abs()
    raise RuntimeError(
      "Baseline recomputation does not reproduce eligibility-aware pipeline "
      f"adjusted scores; maximum absolute difference={differences.max()}"
    )

  leave_one_rows: list[pd.DataFrame] = []
  weight_rows: list[pd.DataFrame] = []
  stability_rows: list[dict[str, object]] = []
  correlations: list[pd.DataFrame] = []
  contexts = (
    ranking[["cancer", "score_comparability_group"]]
      .drop_duplicates()
      .astype(str)
      .itertuples(index=False, name=None)
  )
  context_values = list(contexts)

  for excluded in EXPANDED_ONCO_WEIGHTS:
    weights = {
      key: value
      for key, value in EXPANDED_ONCO_WEIGHTS.items()
      if key != excluded
    }
    scores = recompute(ranking, weights)
    alternative = pd.concat(
      [ranking[base_columns(ranking)].reset_index(drop=True), scores],
      axis=1,
    )
    alternative["excluded_domain"] = excluded
    alternative["alternative_rank_within_comparability_group"] = (
      rank_within_context(alternative, "recomputed_adjusted_score")
    )
    leave_one_rows.append(alternative)
    scenario = f"leave_out_{excluded}"
    add_correlations(correlations, ranking, alternative, scenario)
    for cancer, comparability_group in context_values:
      baseline_set = top_k_set(
        ranking,
        "coverage_adjusted_rses",
        cancer,
        comparability_group,
        args.top_k,
      )
      alternative_set = top_k_set(
        alternative,
        "recomputed_adjusted_score",
        cancer,
        comparability_group,
        args.top_k,
      )
      union = baseline_set | alternative_set
      stability_rows.append({
        "scenario": scenario,
        "cancer": cancer,
        "score_comparability_group": comparability_group,
        "top_k": args.top_k,
        "overlap_n": len(baseline_set & alternative_set),
        "jaccard": (
          len(baseline_set & alternative_set) / len(union)
          if union
          else np.nan
        ),
        "baseline_top_pairs": ";".join(sorted(baseline_set)),
        "alternative_top_pairs": ";".join(sorted(alternative_set)),
      })

  for scenario_name, weights in controlled_weight_scenarios(
    args.weight_delta
  ):
    scores = recompute(ranking, weights)
    alternative = pd.concat(
      [ranking[base_columns(ranking)].reset_index(drop=True), scores],
      axis=1,
    )
    alternative["scenario"] = scenario_name
    alternative["weights"] = ";".join(
      f"{key}={value:.8f}" for key, value in weights.items()
    )
    alternative["alternative_rank_within_comparability_group"] = (
      rank_within_context(alternative, "recomputed_adjusted_score")
    )
    weight_rows.append(alternative)
    scenario = f"weights_{scenario_name}"
    add_correlations(correlations, ranking, alternative, scenario)
    for cancer, comparability_group in context_values:
      baseline_set = top_k_set(
        ranking,
        "coverage_adjusted_rses",
        cancer,
        comparability_group,
        args.top_k,
      )
      alternative_set = top_k_set(
        alternative,
        "recomputed_adjusted_score",
        cancer,
        comparability_group,
        args.top_k,
      )
      union = baseline_set | alternative_set
      stability_rows.append({
        "scenario": scenario,
        "cancer": cancer,
        "score_comparability_group": comparability_group,
        "top_k": args.top_k,
        "overlap_n": len(baseline_set & alternative_set),
        "jaccard": (
          len(baseline_set & alternative_set) / len(union)
          if union
          else np.nan
        ),
        "baseline_top_pairs": ";".join(sorted(baseline_set)),
        "alternative_top_pairs": ";".join(sorted(alternative_set)),
      })

  leave_one = pd.concat(leave_one_rows, ignore_index=True)
  weight_sensitivity = pd.concat(weight_rows, ignore_index=True)
  correlations_frame = pd.concat(correlations, ignore_index=True, sort=False)
  stability = pd.DataFrame(stability_rows)

  missingness_rows = []
  for record in ranking.to_dict("records"):
    eligible_domains = [
      domain
      for domain in EXPANDED_ONCO_WEIGHTS
      if boolean_series(
        pd.DataFrame([record]),
        f"eligible_component_{domain}",
      ).iloc[0]
    ]
    observed_domains = [
      domain
      for domain in eligible_domains
      if pd.notna(
        pd.to_numeric(
          pd.Series([record.get(f"component_{domain}")]),
          errors="coerce",
        ).iloc[0]
      )
    ]
    missingness_rows.append({
      "cancer": record.get("cancer"),
      "pair_id": record.get("pair_id"),
      "source_class": record.get("source_class"),
      "score_comparability_group": record.get(
        "score_comparability_group"
      ),
      "eligible_domains": len(eligible_domains),
      "observed_domains": len(observed_domains),
      "missing_eligible_domains": (
        len(eligible_domains) - len(observed_domains)
      ),
      "noneligible_domains": (
        len(EXPANDED_ONCO_WEIGHTS) - len(eligible_domains)
      ),
      "pipeline_adjusted_score": record.get("coverage_adjusted_rses"),
      "pipeline_coverage": record.get("evidence_coverage"),
      "missing_data_interpretation": (
        "Missing eligible components are omitted from the observed denominator "
        "and penalized through coverage; non-eligible domains enter neither "
        "denominator and comparisons are stratified by comparability group."
      ),
    })
  missingness_sensitivity = pd.DataFrame(missingness_rows)

  outputs = [
    atomic_tsv(
      baseline_rows,
      robustness_dir / "raw_vs_coverage_adjusted.tsv",
    ),
    atomic_tsv(
      leave_one,
      robustness_dir / "leave_one_domain_out.tsv",
    ),
    atomic_tsv(
      weight_sensitivity,
      robustness_dir / "controlled_weight_perturbation.tsv",
    ),
    atomic_tsv(
      correlations_frame,
      robustness_dir / "ranking_correlations.tsv",
    ),
    atomic_tsv(stability, robustness_dir / "top_k_stability.tsv"),
    atomic_tsv(
      missingness_sensitivity,
      robustness_dir / "missing_data_sensitivity.tsv",
    ),
    atomic_tsv(
      stability,
      supplementary_dir / "Table_S23_rses_ranking_stability.tsv",
    ),
    atomic_tsv(
      leave_one,
      supplementary_dir / "Table_S24_leave_one_domain_out.tsv",
    ),
    atomic_tsv(
      weight_sensitivity,
      supplementary_dir / "Table_S25_weight_sensitivity.tsv",
    ),
  ]
  summary = {
    "ranking_rows": len(ranking),
    "leave_one_domain_out_rows": len(leave_one),
    "weight_sensitivity_rows": len(weight_sensitivity),
    "correlation_rows": len(correlations_frame),
    "stability_rows": len(stability),
    "top_k": args.top_k,
    "weight_delta": args.weight_delta,
    "comparison_contexts": len(context_values),
    "rules": {
      "zero_imputation": False,
      "missing_domains_in_observed_denominator": False,
      "noneligible_domains_in_eligible_denominator": False,
      "cross_comparability_group_ranking": False,
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
  print("Eligibility-aware RSES-Onco robustness analyses passed.")


if __name__ == "__main__":
  main()
