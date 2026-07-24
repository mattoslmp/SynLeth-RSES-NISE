from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


CANCER_ALIASES = {
  "colon": {"colon", "colorectal", "large intestine"},
  "stomach": {"stomach", "gastric"},
  "lung": {"lung", "non-small cell lung", "small cell lung"},
}

DIRECT_SCORE_SOURCE_KEYS = {
  "crispr_dependency",
  "absolute_cn",
  "loh",
  "damaging_mutations",
  "mutation_table",
  "proteomics_gygi",
  "proteomics_sanger",
  "proteomics_olink",
  "proteomics_rppa",
  "proteomics_rppa500",
  "rnai_demeter2",
}

EXPLORATORY_SOURCE_KEYS = {
  "hotspot_mutations",
  "fusions",
  "metabolomics",
  "ssgsea",
  "mirna",
  "global_chromatin",
  "molecular_subtypes",
  "subtype_matrix",
  "omics_signatures",
  "metmap_125",
  "metmap_500",
  "metmap_penetrance",
}


@dataclass(frozen=True)
class SourceSpec:
  key: str
  filename: str
  role: str
  value_direction: str = "higher"
  required: bool = False


@dataclass(frozen=True)
class Contrast:
  n_loss: int
  n_intact: int
  median_loss: float
  median_intact: float
  delta: float
  standardized_delta: float
  p_value: float | None
  score: float | None


def clamp01(value: object) -> float | None:
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  if not np.isfinite(result):
    return None
  return float(np.clip(result, 0.0, 1.0))


def sha256_file(path: str | Path) -> str:
  digest = hashlib.sha256()
  with Path(path).open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def read_table(path: str | Path, nrows: int | None = None) -> pd.DataFrame:
  path = Path(path)
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(path)
  suffix = path.suffix.lower()
  if suffix in {".tsv", ".txt"}:
    return pd.read_csv(path, sep="\t", low_memory=False, nrows=nrows)
  return pd.read_csv(path, low_memory=False, nrows=nrows)


def normalize_text(value: object) -> str:
  if value is None or pd.isna(value):
    return ""
  return re.sub(r"\s+", " ", str(value).strip())


def canonical_gene(value: object) -> str:
  text = normalize_text(value)
  if not text:
    return ""
  text = re.sub(r"\s*\([^)]*\)\s*$", "", text)
  text = re.sub(r"\s*\[[^]]*\]\s*$", "", text)
  text = text.split("|")[0].strip()
  return text.upper()


def canonical_feature(value: object) -> str:
  text = normalize_text(value)
  text = re.sub(r"^Unnamed:\s*\d+$", "", text, flags=re.I)
  return text


def find_first_column(columns: Iterable[object], candidates: Iterable[str]) -> str | None:
  original = [str(column) for column in columns]
  exact = {column: column for column in original}
  folded = {column.casefold(): column for column in original}
  for candidate in candidates:
    if candidate in exact:
      return exact[candidate]
    if candidate.casefold() in folded:
      return folded[candidate.casefold()]
  return None


def build_model_lookup(models: pd.DataFrame) -> dict[str, str]:
  if "ModelID" not in models.columns:
    raise ValueError("Model.csv must contain ModelID")
  alias_columns = [
    column
    for column in (
      "ModelID",
      "CCLEName",
      "CellLineName",
      "StrippedCellLineName",
      "ModelConditionID",
      "COSMICID",
      "COSMIC_ID",
      "OncotreeLineage",
    )
    if column in models.columns
  ]
  lookup: dict[str, str] = {}
  for record in models[alias_columns].to_dict("records"):
    model_id = normalize_text(record.get("ModelID"))
    if not model_id:
      continue
    lookup[model_id.casefold()] = model_id
    for column in alias_columns:
      value = normalize_text(record.get(column))
      if value:
        lookup.setdefault(value.casefold(), model_id)
  return lookup


def map_model_id(value: object, lookup: Mapping[str, str]) -> str | None:
  text = normalize_text(value)
  if not text:
    return None
  if text.upper().startswith("ACH-"):
    return text.upper()
  return lookup.get(text.casefold())


def _recognized_model_fraction(values: Iterable[object], lookup: Mapping[str, str]) -> float:
  values = list(values)
  if not values:
    return 0.0
  recognized = sum(map_model_id(value, lookup) is not None for value in values)
  return recognized / len(values)


def read_model_feature_matrix(
  path: str | Path,
  models: pd.DataFrame,
  *,
  gene_features: bool = False,
) -> pd.DataFrame:
  """Read a custom-download matrix and return ModelID-indexed numeric values.

  Both model-by-feature and feature-by-model layouts are supported. Duplicate
  feature columns are collapsed by median after canonicalization.
  """
  frame = read_table(path)
  if frame.empty:
    return pd.DataFrame()
  lookup = build_model_lookup(models)
  id_candidates = (
    "ModelID",
    "model_id",
    "DepMap_ID",
    "depmap_id",
    "Broad_ID",
    "CCLEName",
    "CellLineName",
    "cell_line",
    "Unnamed: 0",
  )
  id_column = find_first_column(frame.columns, id_candidates)
  if id_column is None:
    id_column = str(frame.columns[0])

  row_fraction = _recognized_model_fraction(frame[id_column].head(500), lookup)
  column_fraction = _recognized_model_fraction(frame.columns[1:501], lookup)

  if column_fraction > row_fraction and column_fraction >= 0.25:
    feature_column = id_column
    transposed = frame.set_index(feature_column).T
    transposed.index = [map_model_id(value, lookup) for value in transposed.index]
    transposed = transposed.loc[pd.notna(transposed.index)]
    matrix = transposed
  else:
    model_ids = frame[id_column].map(lambda value: map_model_id(value, lookup))
    matrix = frame.drop(columns=[id_column]).copy()
    matrix.index = model_ids
    matrix = matrix.loc[pd.notna(matrix.index)]

  if matrix.empty:
    return matrix
  matrix.index = matrix.index.astype(str)
  matrix = matrix.apply(pd.to_numeric, errors="coerce")
  matrix = matrix.loc[:, matrix.notna().any(axis=0)]
  if gene_features:
    names = [canonical_gene(column) for column in matrix.columns]
  else:
    names = [canonical_feature(column) for column in matrix.columns]
  matrix.columns = names
  matrix = matrix.loc[:, [bool(column) for column in matrix.columns]]
  if matrix.columns.duplicated().any():
    matrix = matrix.T.groupby(level=0).median().T
  if matrix.index.duplicated().any():
    matrix = matrix.groupby(level=0).median()
  return matrix.sort_index()


def cancer_model_ids(models: pd.DataFrame, cancer: str) -> set[str]:
  aliases = CANCER_ALIASES.get(str(cancer).casefold(), {str(cancer).casefold()})
  lineage_columns = [
    column
    for column in (
      "OncotreeLineage",
      "OncotreePrimaryDisease",
      "lineage",
      "Lineage",
      "PrimaryDisease",
      "ModelType",
    )
    if column in models.columns
  ]
  if not lineage_columns:
    return set(models["ModelID"].astype(str))
  mask = pd.Series(False, index=models.index)
  for column in lineage_columns:
    normalized = models[column].astype(str).str.casefold()
    for alias in aliases:
      mask |= normalized.str.contains(re.escape(alias), na=False)
  return set(models.loc[mask, "ModelID"].astype(str))


def robust_scale(values: pd.Series) -> pd.Series:
  numeric = pd.to_numeric(values, errors="coerce")
  median = numeric.median()
  mad = (numeric - median).abs().median()
  if not np.isfinite(mad) or mad <= 1e-12:
    std = numeric.std(ddof=0)
    scale = std if np.isfinite(std) and std > 1e-12 else 1.0
  else:
    scale = 1.4826 * mad
  return (numeric - median) / scale


def contrast_score(
  values: pd.Series,
  loss_mask: pd.Series,
  *,
  lower_is_supportive: bool,
  min_group_size: int = 3,
  saturation: float = 1.0,
) -> Contrast | None:
  aligned = pd.DataFrame({
    "value": pd.to_numeric(values, errors="coerce"),
    "loss": loss_mask.astype(bool),
  }).dropna(subset=["value"])
  loss = aligned.loc[aligned["loss"], "value"]
  intact = aligned.loc[~aligned["loss"], "value"]
  if len(loss) < min_group_size or len(intact) < min_group_size:
    return None
  median_loss = float(loss.median())
  median_intact = float(intact.median())
  delta = median_loss - median_intact
  pooled = pd.concat([loss, intact], ignore_index=True)
  mad = float((pooled - pooled.median()).abs().median())
  scale = 1.4826 * mad if np.isfinite(mad) and mad > 1e-12 else float(pooled.std(ddof=0))
  if not np.isfinite(scale) or scale <= 1e-12:
    scale = 1.0
  standardized = float(delta / scale)
  directed = -standardized if lower_is_supportive else standardized
  score = float(np.clip(max(directed, 0.0) / max(saturation, 1e-8), 0.0, 1.0))
  try:
    p_value = float(mannwhitneyu(loss, intact, alternative="two-sided").pvalue)
  except ValueError:
    p_value = None
  return Contrast(
    n_loss=len(loss),
    n_intact=len(intact),
    median_loss=median_loss,
    median_intact=median_intact,
    delta=delta,
    standardized_delta=standardized,
    p_value=p_value,
    score=score,
  )


def _binary_component(matrix: pd.DataFrame | None, gene: str) -> pd.Series | None:
  if matrix is None or matrix.empty or gene not in matrix.columns:
    return None
  values = pd.to_numeric(matrix[gene], errors="coerce")
  unique = set(values.dropna().unique())
  if unique and unique.issubset({0, 1, 0.0, 1.0, False, True}):
    return values.astype(float)
  return (values > 0).astype(float).where(values.notna())


def build_functional_loss_table(
  candidates: pd.DataFrame,
  models: pd.DataFrame,
  *,
  relative_cn: pd.DataFrame | None = None,
  absolute_cn: pd.DataFrame | None = None,
  loh: pd.DataFrame | None = None,
  damaging: pd.DataFrame | None = None,
  hotspot: pd.DataFrame | None = None,
  fusions: pd.DataFrame | None = None,
  loss_threshold: float = 0.30,
) -> pd.DataFrame:
  genes = sorted({
    canonical_gene(value)
    for column in ("lost_gene", "target_gene")
    if column in candidates.columns
    for value in candidates[column]
    if canonical_gene(value)
  })
  model_ids = models["ModelID"].astype(str)
  rows: list[pd.DataFrame] = []
  for gene in genes:
    table = pd.DataFrame(index=model_ids)
    table.index.name = "ModelID"
    table["gene"] = gene

    if relative_cn is not None and gene in relative_cn.columns:
      values = pd.to_numeric(relative_cn[gene], errors="coerce").reindex(model_ids)
      table["relative_cn"] = values
      relative_component = (values < loss_threshold).astype(float).where(values.notna())
    else:
      table["relative_cn"] = np.nan
      relative_component = pd.Series(np.nan, index=model_ids)

    if absolute_cn is not None and gene in absolute_cn.columns:
      values = pd.to_numeric(absolute_cn[gene], errors="coerce").reindex(model_ids)
      table["absolute_cn"] = values
      absolute_component = pd.Series(
        np.select(
          [values <= 0.5, values <= 1.5],
          [1.0, 0.45],
          default=0.0,
        ),
        index=model_ids,
      ).where(values.notna())
    else:
      table["absolute_cn"] = np.nan
      absolute_component = pd.Series(np.nan, index=model_ids)

    component_map: dict[str, pd.Series] = {
      "relative_cn_loss": relative_component,
      "absolute_cn_loss": absolute_component,
    }
    for name, matrix in (
      ("loh", loh),
      ("damaging_mutation", damaging),
      ("hotspot_mutation", hotspot),
      ("fusion", fusions),
    ):
      component = _binary_component(matrix, gene)
      if component is None:
        component = pd.Series(np.nan, index=model_ids)
      else:
        component = component.reindex(model_ids)
      table[name] = component
      component_map[name] = component

    combined = pd.DataFrame(component_map, index=model_ids)
    biallelic = (
      (combined["absolute_cn_loss"] >= 0.95)
      | ((combined["loh"] >= 0.5) & (combined["damaging_mutation"] >= 0.5))
    )
    probable = (
      (combined["relative_cn_loss"] >= 0.5)
      | (combined["absolute_cn_loss"] >= 0.4)
      | (combined["damaging_mutation"] >= 0.5)
    )
    score = pd.concat([
      combined["relative_cn_loss"].fillna(0.0),
      combined["absolute_cn_loss"].fillna(0.0),
      (combined["loh"].fillna(0.0) * 0.45),
      (combined["damaging_mutation"].fillna(0.0) * 0.65),
      (combined["hotspot_mutation"].fillna(0.0) * 0.10),
      (combined["fusion"].fillna(0.0) * 0.15),
    ], axis=1).max(axis=1)
    score = score.where(combined.notna().any(axis=1))
    score.loc[biallelic] = 1.0
    score.loc[~biallelic & probable] = np.maximum(
      score.loc[~biallelic & probable].fillna(0.0),
      0.60,
    )
    table["functional_loss_score"] = score
    table["functional_loss_state"] = np.select(
      [score >= 0.85, score >= 0.60, score >= 0.30, score < 0.30],
      ["biallelic_or_homdel", "probable_functional_loss", "partial_loss", "intact"],
      default="missing",
    )
    table["event_source_count"] = combined.notna().sum(axis=1)
    rows.append(table.reset_index())
  if not rows:
    return pd.DataFrame()
  return pd.concat(rows, ignore_index=True)


def _matrix_for_gene(matrix: pd.DataFrame | None, gene: str) -> pd.Series | None:
  if matrix is None or matrix.empty or gene not in matrix.columns:
    return None
  result = pd.to_numeric(matrix[gene], errors="coerce")
  result.index = result.index.astype(str)
  return result


def _global_shift(
  matrix: pd.DataFrame | None,
  model_ids: set[str],
  loss_models: set[str],
  *,
  max_features: int = 500,
  min_group_size: int = 3,
) -> tuple[float | None, int, int, int]:
  if matrix is None or matrix.empty:
    return None, 0, 0, 0
  selected = matrix.loc[matrix.index.astype(str).isin(model_ids)].copy()
  if selected.empty:
    return None, 0, 0, 0
  loss_mask = selected.index.astype(str).isin(loss_models)
  n_loss = int(loss_mask.sum())
  n_intact = int((~loss_mask).sum())
  if n_loss < min_group_size or n_intact < min_group_size:
    return None, n_loss, n_intact, 0
  variances = selected.var(axis=0, skipna=True).sort_values(ascending=False)
  columns = list(variances.head(max_features).index)
  shifts = []
  for column in columns:
    contrast = contrast_score(
      selected[column],
      pd.Series(loss_mask, index=selected.index),
      lower_is_supportive=False,
      min_group_size=min_group_size,
      saturation=1.5,
    )
    if contrast is not None:
      shifts.append(abs(contrast.standardized_delta))
  if not shifts:
    return None, n_loss, n_intact, 0
  score = float(np.clip(np.median(shifts) / 1.5, 0.0, 1.0))
  return score, n_loss, n_intact, len(shifts)


def build_pair_evidence(
  ranking: pd.DataFrame,
  models: pd.DataFrame,
  loss_table: pd.DataFrame,
  *,
  dependency_probability: pd.DataFrame | None = None,
  proteomics: Mapping[str, pd.DataFrame] | None = None,
  rnai: pd.DataFrame | None = None,
  metabolomics: pd.DataFrame | None = None,
  covariates: Mapping[str, pd.DataFrame] | None = None,
  min_group_size: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
  proteomics = dict(proteomics or {})
  covariates = dict(covariates or {})
  loss_index = loss_table.set_index(["ModelID", "gene"], drop=False)
  rows: list[dict[str, object]] = []
  protein_rows: list[dict[str, object]] = []
  covariate_rows: list[dict[str, object]] = []
  shift_cache: dict[tuple[str, str, str], tuple[float | None, int, int, int]] = {}

  key_columns = [column for column in ("pair_id", "cancer", "lost_gene", "target_gene") if column in ranking.columns]
  hypotheses = ranking[key_columns].drop_duplicates()
  for record in hypotheses.to_dict("records"):
    pair_id = str(record.get("pair_id") or "")
    cancer = str(record.get("cancer") or "")
    lost_gene = canonical_gene(record.get("lost_gene"))
    target_gene = canonical_gene(record.get("target_gene"))
    cancer_models = cancer_model_ids(models, cancer)
    if not pair_id or not lost_gene or not target_gene or not cancer_models:
      continue

    try:
      gene_loss = loss_index.xs(lost_gene, level="gene").copy()
    except KeyError:
      gene_loss = pd.DataFrame()
    if gene_loss.empty:
      rows.append({
        "pair_id": pair_id,
        "cancer": cancer,
        "lost_gene": lost_gene,
        "target_gene": target_gene,
        "evidence_status": "missing_loss_state",
      })
      continue
    gene_loss = gene_loss.loc[gene_loss["ModelID"].astype(str).isin(cancer_models)]
    observed_loss = pd.to_numeric(gene_loss["functional_loss_score"], errors="coerce")
    evaluable = observed_loss.notna()
    loss_mask = observed_loss >= 0.60
    loss_models = set(gene_loss.loc[loss_mask & evaluable, "ModelID"].astype(str))
    intact_models = set(gene_loss.loc[(observed_loss < 0.30) & evaluable, "ModelID"].astype(str))
    event_frequency = (
      len(loss_models) / (len(loss_models) + len(intact_models))
      if len(loss_models) + len(intact_models) > 0
      else None
    )
    event_score = (
      float(np.clip(event_frequency / 0.10, 0.0, 1.0))
      if event_frequency is not None
      else None
    )

    row: dict[str, object] = {
      "pair_id": pair_id,
      "cancer": cancer,
      "lost_gene": lost_gene,
      "target_gene": target_gene,
      "integrated_loss_n": len(loss_models),
      "integrated_intact_n": len(intact_models),
      "integrated_loss_frequency": event_frequency,
      "integrated_functional_loss_support": event_score,
      "integrated_event_source_count_median": float(gene_loss["event_source_count"].median()),
      "evidence_status": "ok" if len(loss_models) >= min_group_size and len(intact_models) >= min_group_size else "insufficient_groups",
    }

    dep_values = _matrix_for_gene(dependency_probability, target_gene)
    dep_contrast = None
    if dep_values is not None:
      selected_ids = sorted((loss_models | intact_models) & set(dep_values.dropna().index.astype(str)))
      if selected_ids:
        dep_contrast = contrast_score(
          dep_values.reindex(selected_ids),
          pd.Series([model in loss_models for model in selected_ids], index=selected_ids),
          lower_is_supportive=False,
          min_group_size=min_group_size,
          saturation=1.0,
        )
    if dep_contrast is not None:
      loss_probability = float(np.clip(dep_contrast.median_loss, 0.0, 1.0))
      row.update({
        "dependency_probability_n_loss": dep_contrast.n_loss,
        "dependency_probability_n_intact": dep_contrast.n_intact,
        "dependency_probability_median_loss": dep_contrast.median_loss,
        "dependency_probability_median_intact": dep_contrast.median_intact,
        "dependency_probability_delta": dep_contrast.delta,
        "dependency_probability_p_value": dep_contrast.p_value,
        "dependency_probability_support": float(np.mean([
          loss_probability,
          dep_contrast.score or 0.0,
        ])),
      })

    protein_scores = []
    for source, matrix in proteomics.items():
      values = _matrix_for_gene(matrix, target_gene)
      source_record: dict[str, object] = {
        "pair_id": pair_id,
        "cancer": cancer,
        "lost_gene": lost_gene,
        "target_gene": target_gene,
        "source": source,
      }
      if values is None:
        source_record["status"] = "target_not_measured"
        protein_rows.append(source_record)
        continue
      selected_ids = sorted((loss_models | intact_models) & set(values.dropna().index.astype(str)))
      contrast = contrast_score(
        values.reindex(selected_ids),
        pd.Series([model in loss_models for model in selected_ids], index=selected_ids),
        lower_is_supportive=False,
        min_group_size=min_group_size,
        saturation=1.0,
      )
      if contrast is None:
        source_record["status"] = "insufficient_groups"
      else:
        source_record.update({
          "status": "ok",
          "n_loss": contrast.n_loss,
          "n_intact": contrast.n_intact,
          "median_loss": contrast.median_loss,
          "median_intact": contrast.median_intact,
          "delta": contrast.delta,
          "standardized_delta": contrast.standardized_delta,
          "p_value": contrast.p_value,
          "protein_compensation_support": contrast.score,
        })
        if contrast.score is not None:
          protein_scores.append(contrast.score)
      protein_rows.append(source_record)
    if protein_scores:
      row["protein_compensation_support"] = float(np.median(protein_scores))
      row["protein_source_coverage"] = len(protein_scores) / max(len(proteomics), 1)
      row["protein_source_count"] = len(protein_scores)
      row["protein_source_concordance"] = float(np.mean(np.asarray(protein_scores) > 0.0))

    rnai_values = _matrix_for_gene(rnai, target_gene)
    rnai_contrast = None
    if rnai_values is not None:
      selected_ids = sorted((loss_models | intact_models) & set(rnai_values.dropna().index.astype(str)))
      rnai_contrast = contrast_score(
        rnai_values.reindex(selected_ids),
        pd.Series([model in loss_models for model in selected_ids], index=selected_ids),
        lower_is_supportive=True,
        min_group_size=min_group_size,
        saturation=1.0,
      )
    if rnai_contrast is not None:
      row.update({
        "rnai_n_loss": rnai_contrast.n_loss,
        "rnai_n_intact": rnai_contrast.n_intact,
        "rnai_delta": rnai_contrast.delta,
        "rnai_standardized_delta": rnai_contrast.standardized_delta,
        "rnai_p_value": rnai_contrast.p_value,
        "rnai_orthogonal_support": rnai_contrast.score,
      })

    if metabolomics is not None:
      cache_key = (lost_gene, cancer, "metabolomics")
      if cache_key not in shift_cache:
        shift_cache[cache_key] = _global_shift(
          metabolomics,
          cancer_models,
          loss_models,
          min_group_size=min_group_size,
        )
      shift, n_loss, n_intact, n_features = shift_cache[cache_key]
      row.update({
        "metabolomic_state_shift_exploratory": shift,
        "metabolomic_n_loss": n_loss,
        "metabolomic_n_intact": n_intact,
        "metabolomic_feature_count": n_features,
        "metabolomic_scored_in_primary_rses": False,
      })

    for source, matrix in covariates.items():
      cache_key = (lost_gene, cancer, source)
      if cache_key not in shift_cache:
        shift_cache[cache_key] = _global_shift(
          matrix,
          cancer_models,
          loss_models,
          max_features=300,
          min_group_size=min_group_size,
        )
      shift, n_loss, n_intact, n_features = shift_cache[cache_key]
      covariate_rows.append({
        "pair_id": pair_id,
        "cancer": cancer,
        "lost_gene": lost_gene,
        "target_gene": target_gene,
        "source": source,
        "global_shift_exploratory": shift,
        "n_loss": n_loss,
        "n_intact": n_intact,
        "feature_count": n_features,
        "scored_in_primary_rses": False,
      })

    scored_layers = [
      row.get("integrated_functional_loss_support"),
      row.get("dependency_probability_support"),
      row.get("protein_compensation_support"),
      row.get("rnai_orthogonal_support"),
    ]
    row["extended_scored_layer_count"] = sum(clamp01(value) is not None for value in scored_layers)
    row["extended_scored_layer_coverage"] = row["extended_scored_layer_count"] / 4.0
    rows.append(row)

  return (
    pd.DataFrame(rows),
    pd.DataFrame(protein_rows),
    pd.DataFrame(covariate_rows),
  )


def coverage_consensus(
  components: Mapping[str, float | None],
  weights: Mapping[str, float],
) -> tuple[float | None, float]:
  numerator = 0.0
  observed_weight = 0.0
  total_weight = float(sum(weights.values()))
  for key, weight in weights.items():
    value = clamp01(components.get(key))
    if value is None:
      continue
    numerator += weight * value
    observed_weight += weight
  if observed_weight <= 0:
    return None, 0.0
  observed = numerator / observed_weight
  coverage = observed_weight / total_weight if total_weight else 0.0
  return float(observed * coverage), float(coverage)


def write_status_json(path: str | Path, payload: Mapping[str, object]) -> Path:
  path = Path(path)
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(json.dumps(dict(payload), indent=2, sort_keys=True), encoding="utf-8")
  return path


def read_long_event_matrix(
  path: str | Path,
  models: pd.DataFrame,
  *,
  event: str,
) -> pd.DataFrame:
  """Convert a long mutation/fusion table to a binary ModelID × gene matrix."""
  frame = read_table(path)
  if frame.empty:
    return pd.DataFrame()
  lookup = build_model_lookup(models)
  model_column = find_first_column(
    frame.columns,
    (
      "ModelID",
      "DepMap_ID",
      "depmap_id",
      "model_id",
      "CCLEName",
      "Tumor_Sample_Barcode",
      "cell_line",
    ),
  )
  gene_column = find_first_column(
    frame.columns,
    (
      "HugoSymbol",
      "Hugo_Symbol",
      "gene",
      "Gene",
      "GeneSymbol",
      "gene_symbol",
      "Symbol",
    ),
  )
  if model_column is None or gene_column is None:
    return read_model_feature_matrix(path, models, gene_features=True)
  standardized = pd.DataFrame({
    "ModelID": frame[model_column].map(lambda value: map_model_id(value, lookup)),
    "gene": frame[gene_column].map(canonical_gene),
  })
  if event == "mutation":
    consequence_column = find_first_column(
      frame.columns,
      (
        "Variant_Classification",
        "Consequence",
        "consequence",
        "Protein_Change",
        "VariantType",
        "Variant_Type",
      ),
    )
    if consequence_column is not None:
      consequences = frame[consequence_column].astype(str).str.casefold()
      clear_lof_terms = (
        "frame_shift",
        "frameshift",
        "nonsense",
        "stop_gained",
        "stop_lost",
        "splice_acceptor",
        "splice_donor",
        "splice_site",
        "start_lost",
        "translation_start_site",
      )
      clear_lof = pd.Series(False, index=frame.index)
      for term in clear_lof_terms:
        clear_lof |= consequences.str.contains(term, na=False)

      annotation_columns = [
        column
        for column in (
          "isDamaging",
          "damaging",
          "LikelyLoF",
          "likely_lof",
          "PolyPhen",
          "SIFT",
          "VEP_Impact",
          "IMPACT",
          "Oncogenic",
        )
        if column in frame.columns
      ]
      annotated_damaging = pd.Series(False, index=frame.index)
      for column in annotation_columns:
        values = frame[column].astype(str).str.casefold()
        annotated_damaging |= values.str.contains(
          r"(^|[^a-z])(true|1|yes|damaging|deleterious|probably_damaging|high)([^a-z]|$)",
          regex=True,
          na=False,
        )
      missense_or_inframe = consequences.str.contains(
        r"missense|in_frame|inframe", regex=True, na=False
      )
      mask = clear_lof | (missense_or_inframe & annotated_damaging)
      standardized = standardized.loc[mask]
  elif event == "fusion":
    partner_column = find_first_column(
      frame.columns,
      ("FusionName", "fusion", "Fusion", "PartnerGene", "partner_gene"),
    )
    if partner_column is not None:
      standardized["partner"] = frame.loc[standardized.index, partner_column].astype(str)
  standardized = standardized.dropna(subset=["ModelID"])
  standardized = standardized.loc[standardized["gene"].astype(bool)]
  if standardized.empty:
    return pd.DataFrame()
  standardized["value"] = 1.0
  matrix = standardized.pivot_table(
    index="ModelID",
    columns="gene",
    values="value",
    aggfunc="max",
    fill_value=0.0,
  )
  return matrix.astype(float)


def standardize_matrix_long(
  matrix: pd.DataFrame,
  *,
  source: str,
  metric: str,
  source_file: str | Path,
) -> pd.DataFrame:
  if matrix.empty:
    return pd.DataFrame(columns=[
      "source",
      "model_id",
      "feature",
      "response_value",
      "response_metric",
      "source_file",
    ])
  long = matrix.reset_index(names="model_id").melt(
    id_vars=["model_id"],
    var_name="feature",
    value_name="response_value",
  )
  long["response_value"] = pd.to_numeric(long["response_value"], errors="coerce")
  long = long.dropna(subset=["response_value"])
  long.insert(0, "source", source)
  long["response_metric"] = metric
  long["source_file"] = str(source_file)
  return long


def parse_combination_feature(feature: object) -> tuple[str, str]:
  text = normalize_text(feature)
  for separator in ("::", " + ", "|", "__", ";", " / "):
    if separator in text:
      left, right = text.split(separator, 1)
      return left.strip(), right.strip()
  match = re.match(r"^(.+?)\s*\+\s*(.+)$", text)
  if match:
    return match.group(1).strip(), match.group(2).strip()
  return "", text


def build_gdsc_combination_table(
  matrices: Mapping[str, pd.DataFrame],
  source_files: Mapping[str, str | Path],
) -> pd.DataFrame:
  """Align GDSC combination matrices without confusing AUC and viability."""
  renamed = {
    "anchor_viability": "anchor_viability",
    "combination_auc": "combination_auc",
    "combination_viability": "combination_viability",
    "library_auc": "library_auc",
    "library_viability": "library_viability",
  }
  pieces = []
  for key, output_column in renamed.items():
    matrix = matrices.get(key)
    if matrix is None or matrix.empty:
      continue
    long = matrix.reset_index(names="model_id").melt(
      id_vars=["model_id"],
      var_name="combination_feature",
      value_name=output_column,
    )
    long[output_column] = pd.to_numeric(long[output_column], errors="coerce")
    pieces.append(long)
  if not pieces:
    return pd.DataFrame()
  combined = pieces[0]
  for piece in pieces[1:]:
    combined = combined.merge(
      piece,
      on=["model_id", "combination_feature"],
      how="outer",
      validate="one_to_one",
    )
  parsed = combined["combination_feature"].map(parse_combination_feature)
  combined["anchor_drug"] = [value[0] for value in parsed]
  combined["library_drug"] = [value[1] for value in parsed]
  if {"combination_auc", "library_auc"}.issubset(combined.columns):
    combined["auc_sensitization"] = (
      combined["library_auc"] - combined["combination_auc"]
    )
  if {
    "anchor_viability",
    "library_viability",
    "combination_viability",
  }.issubset(combined.columns):
    anchor = pd.to_numeric(combined["anchor_viability"], errors="coerce")
    library = pd.to_numeric(combined["library_viability"], errors="coerce")
    combination = pd.to_numeric(combined["combination_viability"], errors="coerce")
    scale = 100.0 if max(anchor.max(skipna=True), library.max(skipna=True), 1.0) > 2.0 else 1.0
    anchor = anchor / scale
    library = library / scale
    combination = combination / scale
    combined["bliss_expected_viability"] = anchor * library
    combined["bliss_excess_effect"] = combined["bliss_expected_viability"] - combination
  combined["source_files"] = ";".join(
    str(source_files[key])
    for key in renamed
    if key in source_files
  )
  return combined


def source_inventory(
  specs: Iterable[SourceSpec],
  data_dir: str | Path,
) -> pd.DataFrame:
  data_dir = Path(data_dir)
  rows = []
  for spec in specs:
    path = data_dir / spec.filename
    row = {
      "source_key": spec.key,
      "filename": spec.filename,
      "role": spec.role,
      "required": spec.required,
      "path": str(path),
      "exists": path.exists(),
      "size_bytes": path.stat().st_size if path.exists() else 0,
      "sha256": sha256_file(path) if path.exists() and path.is_file() else "",
      "direct_score_layer": spec.key in DIRECT_SCORE_SOURCE_KEYS,
      "exploratory_or_validation_layer": spec.key in EXPLORATORY_SOURCE_KEYS,
    }
    rows.append(row)
  return pd.DataFrame(rows)


def build_global_context_evidence(
  ranking: pd.DataFrame,
  models: pd.DataFrame,
  loss_table: pd.DataFrame,
  matrix: pd.DataFrame,
  *,
  source: str,
  min_group_size: int = 3,
  max_features: int = 300,
) -> pd.DataFrame:
  """Build non-causal context shifts for robustness, never primary scoring."""
  if matrix.empty or loss_table.empty:
    return pd.DataFrame()
  loss_index = loss_table.set_index(["ModelID", "gene"], drop=False)
  rows = []
  cache: dict[tuple[str, str], tuple[float | None, int, int, int]] = {}
  columns = [
    column
    for column in ("pair_id", "cancer", "lost_gene", "target_gene")
    if column in ranking.columns
  ]
  for record in ranking[columns].drop_duplicates().to_dict("records"):
    pair_id = str(record.get("pair_id") or "")
    cancer = str(record.get("cancer") or "")
    lost_gene = canonical_gene(record.get("lost_gene"))
    target_gene = canonical_gene(record.get("target_gene"))
    if not pair_id or not lost_gene:
      continue
    key = (lost_gene, cancer)
    if key not in cache:
      try:
        gene_loss = loss_index.xs(lost_gene, level="gene").copy()
      except KeyError:
        gene_loss = pd.DataFrame()
      if gene_loss.empty:
        cache[key] = (None, 0, 0, 0)
      else:
        cancer_models = cancer_model_ids(models, cancer)
        gene_loss = gene_loss.loc[
          gene_loss["ModelID"].astype(str).isin(cancer_models)
        ]
        score = pd.to_numeric(
          gene_loss["functional_loss_score"], errors="coerce"
        )
        loss_models = set(
          gene_loss.loc[score >= 0.60, "ModelID"].astype(str)
        )
        cache[key] = _global_shift(
          matrix,
          cancer_models,
          loss_models,
          max_features=max_features,
          min_group_size=min_group_size,
        )
    shift, n_loss, n_intact, n_features = cache[key]
    rows.append({
      "pair_id": pair_id,
      "cancer": cancer,
      "lost_gene": lost_gene,
      "target_gene": target_gene,
      "source": source,
      "global_shift_exploratory": shift,
      "n_loss": n_loss,
      "n_intact": n_intact,
      "feature_count": n_features,
      "scored_in_primary_rses": False,
      "interpretation": (
        "multivariate_context_shift_not_mechanistic_gene_evidence"
      ),
    })
  return pd.DataFrame(rows)


def bh_adjust(values: Iterable[object]) -> np.ndarray:
  series = pd.to_numeric(pd.Series(list(values)), errors="coerce")
  result = np.full(len(series), np.nan, dtype=float)
  mask = series.notna().to_numpy()
  observed = series[mask].to_numpy(dtype=float)
  if observed.size == 0:
    return result
  order = np.argsort(observed)
  ranked = observed[order]
  adjusted = ranked * observed.size / np.arange(1, observed.size + 1)
  adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
  adjusted = np.clip(adjusted, 0.0, 1.0)
  restored = np.empty_like(adjusted)
  restored[order] = adjusted
  result[np.flatnonzero(mask)] = restored
  return result
