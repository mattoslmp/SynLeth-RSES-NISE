from __future__ import annotations

from pathlib import Path
from typing import Mapping

import pandas as pd

from .score import rses_onco_score


def score_candidate_table(
  candidates: pd.DataFrame,
  empirical_components: Mapping[str, Mapping[str, float | None]] | None = None,
) -> pd.DataFrame:
  """Score candidate pairs with optional empirical component overrides."""
  empirical_components = empirical_components or {}
  rows = []
  for record in candidates.to_dict("records"):
    pair_id = str(record["pair_id"])
    validation = sum(float(record[k]) for k in (
      "genetic_screen", "isogenic_validation", "in_vivo", "clinical_tractability"
    )) / 4.0
    components = {
      "tumor_event": float(record["lineage_relevance"]),
      "dependency": float(record["genetic_screen"]),
      "selectivity": float(record["isogenic_validation"]),
      "expression_compensation": None,
      "functional_relation": float(record["relation_confidence"]),
      "validation_tractability": validation,
    }
    components.update(empirical_components.get(pair_id, {}))
    scored = rses_onco_score(components)
    out = dict(record)
    out.update({f"component_{k}": v for k, v in components.items()})
    out.update({
      "rses_onco": scored.observed_score,
      "evidence_coverage": scored.coverage,
      "coverage_adjusted_rses": scored.adjusted_score,
      "n_domains": scored.n_domains,
      "priority_class": scored.interpretation,
    })
    rows.append(out)
  result = pd.DataFrame(rows)
  return result.sort_values(["coverage_adjusted_rses", "rses_onco"], ascending=False)


def load_reference_candidates(path: str | Path) -> pd.DataFrame:
  return pd.read_csv(path, sep="\t")
