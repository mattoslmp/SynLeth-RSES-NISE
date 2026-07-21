from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .depmap import cancer_model_ids
from .utils import canonical_gene_name, clamp01


FUNCTIONAL_MICRONICHE_WEIGHTS = {
  "expression_context": 0.20,
  "localization": 0.15,
  "biochemical_structural": 0.15,
  "genetic_phenotype": 0.20,
  "interaction_network": 0.15,
  "regulatory_network": 0.15,
}

EXPANDED_ONCO_WEIGHTS = {
  "tumor_event": 0.16,
  "dependency": 0.22,
  "selectivity": 0.14,
  "expression_compensation": 0.08,
  "functional_relation": 0.06,
  "functional_microniche": 0.16,
  "validation_tractability": 0.18,
}


@dataclass(frozen=True)
class CoverageAwareResult:
  observed_score: float
  coverage: float
  adjusted_score: float
  n_domains: int


@dataclass(frozen=True)
class PairProfileMetrics:
  gene_a: str
  gene_b: str
  cancer: str
  n_models: int
  spearman_rho: float | None
  median_absolute_difference: float | None
  set_jaccard: float | None
  divergence: float | None


def coverage_aware_score(
  components: Mapping[str, float | None],
  weights: Mapping[str, float],
) -> CoverageAwareResult:
  total_weight = float(sum(weights.values()))
  numerator = 0.0
  observed_weight = 0.0
  n_domains = 0
  for domain, weight in weights.items():
    value = clamp01(components.get(domain))
    if value is None:
      continue
    numerator += float(weight) * float(value)
    observed_weight += float(weight)
    n_domains += 1
  observed = numerator / observed_weight if observed_weight else float("nan")
  coverage = observed_weight / total_weight if total_weight else float("nan")
  adjusted = observed * coverage if np.isfinite(observed) else float("nan")
  return CoverageAwareResult(observed, coverage, adjusted, n_domains)


def functional_microniche_score(
  components: Mapping[str, float | None],
) -> CoverageAwareResult:
  return coverage_aware_score(components, FUNCTIONAL_MICRONICHE_WEIGHTS)


def expanded_onco_score(
  components: Mapping[str, float | None],
) -> CoverageAwareResult:
  return coverage_aware_score(components, EXPANDED_ONCO_WEIGHTS)


def _normalized_pair_table(
  frame: pd.DataFrame,
  models: pd.DataFrame,
  gene_a: str,
  gene_b: str,
  cancer: str,
) -> pd.DataFrame:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  if gene_a not in frame or gene_b not in frame:
    return pd.DataFrame(columns=["ModelID", gene_a, gene_b])
  selected = set(cancer_model_ids(models, cancer))
  if not selected:
    return pd.DataFrame(columns=["ModelID", gene_a, gene_b])
  table = frame[["ModelID", gene_a, gene_b]].copy()
  table = table.loc[table["ModelID"].astype(str).isin(selected)]
  table[gene_a] = pd.to_numeric(table[gene_a], errors="coerce")
  table[gene_b] = pd.to_numeric(table[gene_b], errors="coerce")
  return table.dropna(subset=[gene_a, gene_b])


def expression_profile_metrics(
  expression: pd.DataFrame,
  models: pd.DataFrame,
  gene_a: str,
  gene_b: str,
  cancer: str,
  separation_saturation: float = 2.0,
) -> PairProfileMetrics:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  table = _normalized_pair_table(expression, models, gene_a, gene_b, cancer)
  if len(table) < 3:
    return PairProfileMetrics(gene_a, gene_b, cancer, len(table), None, None, None, None)
  rho = float(spearmanr(table[gene_a], table[gene_b], nan_policy="omit").statistic)
  if not np.isfinite(rho):
    rho = None
  median_difference = float(np.median(np.abs(table[gene_a] - table[gene_b])))
  correlation_divergence = None if rho is None else (1.0 - rho) / 2.0
  separation = float(np.clip(median_difference / separation_saturation, 0, 1))
  available = [value for value in (correlation_divergence, separation) if value is not None]
  divergence = float(np.mean(available)) if available else None
  return PairProfileMetrics(
    gene_a, gene_b, cancer, len(table), rho, median_difference, None, divergence,
  )


def phenotype_profile_metrics(
  gene_effect: pd.DataFrame,
  models: pd.DataFrame,
  gene_a: str,
  gene_b: str,
  cancer: str,
  dependency_threshold: float = -0.5,
  separation_saturation: float = 1.0,
) -> PairProfileMetrics:
  gene_a = canonical_gene_name(gene_a)
  gene_b = canonical_gene_name(gene_b)
  table = _normalized_pair_table(gene_effect, models, gene_a, gene_b, cancer)
  if len(table) < 3:
    return PairProfileMetrics(gene_a, gene_b, cancer, len(table), None, None, None, None)
  rho = float(spearmanr(table[gene_a], table[gene_b], nan_policy="omit").statistic)
  if not np.isfinite(rho):
    rho = None
  median_difference = float(np.median(np.abs(table[gene_a] - table[gene_b])))
  set_a = set(table.loc[table[gene_a] < dependency_threshold, "ModelID"].astype(str))
  set_b = set(table.loc[table[gene_b] < dependency_threshold, "ModelID"].astype(str))
  union = set_a | set_b
  jaccard = float(len(set_a & set_b) / len(union)) if union else None
  correlation_divergence = None if rho is None else (1.0 - rho) / 2.0
  separation = float(np.clip(median_difference / separation_saturation, 0, 1))
  set_divergence = None if jaccard is None else 1.0 - jaccard
  available = [
    value for value in (correlation_divergence, separation, set_divergence)
    if value is not None
  ]
  divergence = float(np.mean(available)) if available else None
  return PairProfileMetrics(
    gene_a, gene_b, cancer, len(table), rho, median_difference, jaccard, divergence,
  )


def build_directed_nise_candidates(nise_pairs: pd.DataFrame) -> pd.DataFrame:
  required = {
    "group_id", "ec_number", "activity", "gene_a", "uniprot_a", "cluster_a",
    "gene_b", "uniprot_b", "cluster_b",
  }
  missing = sorted(required - set(nise_pairs.columns))
  if missing:
    raise ValueError(f"NISE pair table lacks columns: {missing}")

  rows: list[dict[str, object]] = []
  for record in nise_pairs.to_dict("records"):
    directions = [
      (
        record["gene_a"], record["uniprot_a"], record["cluster_a"],
        record["gene_b"], record["uniprot_b"], record["cluster_b"],
      ),
      (
        record["gene_b"], record["uniprot_b"], record["cluster_b"],
        record["gene_a"], record["uniprot_a"], record["cluster_a"],
      ),
    ]
    for lost, lost_uniprot, lost_cluster, target, target_uniprot, target_cluster in directions:
      lost = canonical_gene_name(lost)
      target = canonical_gene_name(target)
      pair_id = f"NISE_{record['group_id']}_{lost}_TO_{target}"
      rows.append({
        "pair_id": pair_id,
        "lost_feature": f"{lost} loss",
        "lost_gene": lost,
        "target_gene": target,
        "source_class": "NISE",
        "relation_type": "human_NISE",
        "mechanism": (
          f"Loss of {lost} may expose context-specific dependence on the "
          f"non-homologous isofunctional alternative {target} for {record['activity']}."
        ),
        "group_id": record["group_id"],
        "ec_number": record["ec_number"],
        "activity": record["activity"],
        "lost_uniprot": lost_uniprot,
        "target_uniprot": target_uniprot,
        "lost_structural_cluster": lost_cluster,
        "target_structural_cluster": target_cluster,
        "colon": 1,
        "stomach": 1,
        "lung": 1,
        "relation_confidence": 1.0,
        "genetic_screen": 0.0,
        "isogenic_validation": 0.0,
        "in_vivo": 0.0,
        "clinical_tractability": 0.25,
        "lineage_relevance": 0.0,
        "evidence_stage": "systematic discovery",
        "primary_doi": "10.1093/gbe/evx119",
        "supporting_doi": "10.1186/s13104-026-07742-5",
        "status": "unvalidated NISE direction",
      })
  result = pd.DataFrame(rows)
  return result.drop_duplicates(["lost_gene", "target_gene", "source_class"])


def merge_candidate_sources(
  nise_candidates: pd.DataFrame,
  benchmark_candidates: pd.DataFrame,
  additional_sources: Iterable[pd.DataFrame] = (),
) -> pd.DataFrame:
  frames = [nise_candidates.copy(), benchmark_candidates.copy(), *[x.copy() for x in additional_sources]]
  all_columns = sorted(set().union(*(set(frame.columns) for frame in frames)))
  normalized = []
  for frame in frames:
    for column in all_columns:
      if column not in frame:
        frame[column] = np.nan
    if "lost_gene" not in frame or frame["lost_gene"].isna().all():
      frame["lost_gene"] = frame["lost_feature"].astype(str).str.extract(
        r"^([A-Za-z0-9-]+)", expand=False
      )
    frame["lost_gene"] = frame["lost_gene"].map(canonical_gene_name)
    frame["target_gene"] = frame["target_gene"].map(canonical_gene_name)
    if "source_class" not in frame or frame["source_class"].isna().all():
      frame["source_class"] = frame.get("relation_type", "curated")
    normalized.append(frame[all_columns])
  combined = pd.concat(normalized, ignore_index=True)
  combined = combined.sort_values(
    ["lost_gene", "target_gene", "source_class", "pair_id"], na_position="last"
  )
  return combined.drop_duplicates(["lost_gene", "target_gene", "source_class"], keep="first")


def class_member_inventory(candidates: pd.DataFrame) -> pd.DataFrame:
  rows: list[dict[str, object]] = []
  for record in candidates.to_dict("records"):
    source_class = str(record.get("source_class") or record.get("relation_type") or "unclassified")
    for role, column in (("lost_or_biomarker", "lost_gene"), ("target", "target_gene")):
      gene = canonical_gene_name(record.get(column))
      if not gene:
        continue
      rows.append({
        "source_class": source_class,
        "gene": gene,
        "role": role,
        "pair_id": record.get("pair_id"),
      })
  return pd.DataFrame(rows).drop_duplicates().sort_values(["source_class", "gene", "role"])


def load_optional_table(path: str | Path | None) -> pd.DataFrame | None:
  if path is None:
    return None
  path = Path(path)
  if not path.exists():
    raise FileNotFoundError(path)
  separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
  return pd.read_csv(path, sep=separator)
