from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .utils import canonical_gene_name


MODEL_ID_ALIASES = (
  "ModelID",
  "DepMap_ID",
  "DepMapID",
  "DepMap ID",
  "model_id",
)
MODEL_ID_PATTERN = re.compile(r"^ACH-[A-Za-z0-9]+$", re.I)


CANCER_CONTEXT_LABELS = {
  "colon": "Colorectal",
  "colorectal": "Colorectal",
  "bowel": "Colorectal",
  "stomach": "Gastric",
  "gastric": "Gastric",
  "esophagus/stomach": "Gastric",
  "lung": "Lung",
}


def _normalized_text(frame: pd.DataFrame, column: str) -> pd.Series:
  if column not in frame:
    return pd.Series("", index=frame.index, dtype="object")
  return frame[column].fillna("").astype(str).str.strip()


def cancer_model_mask(models: pd.DataFrame, cancer: str) -> pd.Series:
  """Return a 26Q1-compatible mask for the requested cancer cohort.

  DepMap/OncoTree does not label colorectal models as ``Colorectal`` or
  gastric models as ``Gastric`` in ``OncotreeLineage``. In current releases,
  colorectal adenocarcinoma belongs to lineage ``Bowel`` and stomach
  adenocarcinoma belongs to lineage ``Esophagus/Stomach``. The gastric mask
  additionally requires a stomach/gastric subtype or STAD code so esophageal
  models are not silently included.
  """
  key = str(cancer).strip().casefold()
  lineage = _normalized_text(models, "OncotreeLineage").str.casefold()
  primary = _normalized_text(models, "OncotreePrimaryDisease").str.casefold()
  subtype = _normalized_text(models, "OncotreeSubtype")
  subtype_cf = subtype.str.casefold()
  code = _normalized_text(models, "OncotreeCode").str.upper()
  model_type = _normalized_text(models, "DepmapModelType").str.upper()

  if key in {"colon", "colorectal", "bowel"}:
    primary_match = primary.eq("colorectal adenocarcinoma") | primary.str.contains(
      r"\bcolorectal\b", regex=True, na=False
    )
    fallback = lineage.eq("bowel") & subtype_cf.str.contains(
      r"\b(?:colon|rectal|colorectal)\b", regex=True, na=False
    )
    legacy = lineage.isin({"colorectal", "colon", "rectal"})
    return primary_match | fallback | legacy

  if key in {"stomach", "gastric", "esophagus/stomach"}:
    lineage_match = lineage.eq("esophagus/stomach")
    subtype_match = subtype_cf.str.contains(
      r"\b(?:stomach|gastric)\b", regex=True, na=False
    )
    primary_match = primary.str.contains(
      r"\b(?:stomach|gastric)\b", regex=True, na=False
    )
    code_match = code.eq("STAD") | model_type.eq("STAD")
    legacy = lineage.isin({"gastric", "stomach"})
    return (lineage_match & (subtype_match | primary_match | code_match)) | legacy

  if key == "lung":
    return lineage.eq("lung")

  # Backward-compatible exact lineage matching for other user-defined groups.
  return lineage.eq(key)


def cancer_model_ids(models: pd.DataFrame, cancer: str) -> pd.Series:
  """Return normalized ModelIDs belonging to a cancer-specific cohort."""
  if "ModelID" not in models:
    models = normalize_model_id_column(models, "Model.csv")
  return models.loc[cancer_model_mask(models, cancer), "ModelID"].astype(str).str.strip()


@dataclass(frozen=True)
class DependencyContrast:
  lost_gene: str
  target_gene: str
  lineage: str
  n_loss: int
  n_intact: int
  median_effect_loss: float
  median_effect_intact: float
  delta_effect: float
  p_value: float


@dataclass(frozen=True)
class ExpressionContrast:
  lost_gene: str
  target_gene: str
  lineage: str
  n_loss: int
  n_intact: int
  median_expression_loss: float
  median_expression_intact: float
  delta_expression: float
  p_value: float


def _looks_like_model_ids(values: pd.Series) -> bool:
  observed = values.dropna().astype(str).str.strip()
  if observed.empty:
    return False
  matches = observed.str.fullmatch(MODEL_ID_PATTERN)
  return bool(matches.mean() >= 0.80)


def detect_model_id_column(
  path: str | Path,
  dataset_name: str = "DepMap dataset",
) -> tuple[str, bool]:
  """Return the raw CSV column that contains ModelID.

  DepMap release matrices have historically represented the row identifier in
  three ways: a literal ``ModelID`` column, a legacy ``DepMap_ID`` column, or
  an unnamed CSV index column whose values are ACH-* model identifiers. This
  detector accepts all three without rewriting the large source files.
  """
  path = Path(path)
  header = pd.read_csv(path, nrows=0)
  columns = [str(column) for column in header.columns]

  for alias in MODEL_ID_ALIASES:
    if alias in columns:
      return alias, alias != "ModelID"

  if not columns:
    raise ValueError(f"{dataset_name} has no readable CSV columns: {path}")

  first = columns[0]
  first_is_index_like = first.strip() == "" or first.startswith("Unnamed:")
  if first_is_index_like:
    sample = pd.read_csv(path, usecols=[first], nrows=50)[first]
    if _looks_like_model_ids(sample):
      return first, True

  raise ValueError(
    f"{dataset_name} does not expose a usable ModelID column. "
    f"First columns: {columns[:6]}. File: {path}"
  )


def normalize_model_id_column(
  frame: pd.DataFrame,
  dataset_name: str = "DepMap dataset",
) -> pd.DataFrame:
  """Normalize ModelID/DepMap_ID/unnamed ACH-* index to ``ModelID``."""
  columns = [str(column) for column in frame.columns]
  for alias in MODEL_ID_ALIASES:
    if alias in columns:
      if alias == "ModelID":
        result = frame
      else:
        result = frame.rename(columns={alias: "ModelID"})
      result["ModelID"] = result["ModelID"].astype(str).str.strip()
      return result

  if not columns:
    raise ValueError(f"{dataset_name} has no columns")

  first = columns[0]
  first_is_index_like = first.strip() == "" or first.startswith("Unnamed:")
  if first_is_index_like and _looks_like_model_ids(frame.iloc[:, 0]):
    result = frame.rename(columns={frame.columns[0]: "ModelID"})
    result["ModelID"] = result["ModelID"].astype(str).str.strip()
    return result

  raise ValueError(
    f"{dataset_name} requires ModelID. First columns: {columns[:6]}"
  )


def read_model_ids(
  path: str | Path,
  dataset_name: str = "DepMap dataset",
) -> pd.Series:
  """Read only the model identifier column and normalize it to strings."""
  raw_column, _ = detect_model_id_column(path, dataset_name)
  values = pd.read_csv(path, usecols=[raw_column])[raw_column]
  return values.dropna().astype(str).str.strip()


def _rename_gene_columns(frame: pd.DataFrame) -> pd.DataFrame:
  metadata = {
    "ModelID", "ProfileID", "PROFILEID", "is_default_entry",
    "IsDefaultEntryForModel", "IsDefaultEntryForMC",
    "ModelConditionID", "SequencingID",
  }
  renamed = {c: canonical_gene_name(c) for c in frame.columns if c not in metadata}
  return frame.rename(columns=renamed)


def _default_model_rows(frame: pd.DataFrame, name: str) -> pd.DataFrame:
  if "ModelID" not in frame:
    raise ValueError(f"{name} requires a ModelID column")
  if not frame["ModelID"].duplicated().any():
    return frame
  for column in ("IsDefaultEntryForModel", "is_default_entry"):
    if column in frame:
      flag = frame[column].astype(str).str.casefold().isin({"1", "true", "yes"})
      selected = frame.loc[flag].copy()
      if not selected.empty and not selected["ModelID"].duplicated().any():
        return selected
  duplicates = int(frame["ModelID"].duplicated().sum())
  raise ValueError(
    f"{name} contains {duplicates} duplicate ModelID rows and no unique default-model selection could be made"
  )


def _read_matrix(path: str | Path, name: str) -> pd.DataFrame:
  frame = pd.read_csv(path)
  frame = normalize_model_id_column(frame, name)
  return _rename_gene_columns(frame)


def read_depmap_inputs(
  gene_effect_path: str | Path,
  copy_number_path: str | Path,
  model_path: str | Path,
  expression_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
  """Read DepMap matrices using normalized ModelID as the stable join key."""
  effect = _default_model_rows(_read_matrix(gene_effect_path, "gene_effect"), "gene_effect")
  copy_number = _default_model_rows(_read_matrix(copy_number_path, "copy_number"), "copy_number")
  models = normalize_model_id_column(pd.read_csv(model_path), "Model.csv")
  expression = None
  if expression_path:
    expression = _default_model_rows(_read_matrix(expression_path, "expression"), "expression")
  if models["ModelID"].duplicated().any():
    raise ValueError("Model.csv contains duplicate ModelID values")
  return effect, copy_number, models, expression


def _lineage_join(
  copy_number: pd.DataFrame,
  models: pd.DataFrame,
  lost_gene: str,
  lineage: str,
  lineage_column: str,
) -> pd.DataFrame:
  # ``lineage`` accepts cancer keys (colon, stomach, lung) as well as exact
  # lineage labels for backward compatibility. ``lineage_column`` is retained
  # in the public signature but current standard cohorts use the full metadata
  # rules in ``cancer_model_mask``.
  if lineage_column not in models:
    raise ValueError(f"Model.csv lacks lineage column: {lineage_column}")
  selected_ids = cancer_model_ids(models, lineage)
  metadata = pd.DataFrame({"ModelID": selected_ids.drop_duplicates()})
  return metadata.merge(copy_number[["ModelID", lost_gene]], on="ModelID", how="inner")


def dependency_contrast(
  effect: pd.DataFrame,
  copy_number: pd.DataFrame,
  models: pd.DataFrame,
  lost_gene: str,
  target_gene: str,
  lineage: str,
  lineage_column: str = "OncotreeLineage",
  loss_threshold: float = 0.30,
  min_group_size: int = 3,
) -> DependencyContrast | None:
  """Compare target CRISPR effect in copy-number-loss versus intact models."""
  lost_gene = canonical_gene_name(lost_gene)
  target_gene = canonical_gene_name(target_gene)
  if lost_gene not in copy_number or target_gene not in effect:
    return None
  joined = (
    _lineage_join(copy_number, models, lost_gene, lineage, lineage_column)
      .merge(effect[["ModelID", target_gene]], on="ModelID", how="inner")
      .dropna(subset=[lost_gene, target_gene])
  )
  loss = joined.loc[joined[lost_gene] < loss_threshold, target_gene].astype(float)
  intact = joined.loc[joined[lost_gene] >= loss_threshold, target_gene].astype(float)
  if len(loss) < min_group_size or len(intact) < min_group_size:
    return None
  p_value = float(mannwhitneyu(loss, intact, alternative="less").pvalue)
  median_loss = float(np.median(loss))
  median_intact = float(np.median(intact))
  return DependencyContrast(
    lost_gene=lost_gene,
    target_gene=target_gene,
    lineage=lineage,
    n_loss=len(loss),
    n_intact=len(intact),
    median_effect_loss=median_loss,
    median_effect_intact=median_intact,
    delta_effect=median_loss - median_intact,
    p_value=p_value,
  )


def expression_contrast(
  expression: pd.DataFrame,
  copy_number: pd.DataFrame,
  models: pd.DataFrame,
  lost_gene: str,
  target_gene: str,
  lineage: str,
  lineage_column: str = "OncotreeLineage",
  loss_threshold: float = 0.30,
  min_group_size: int = 3,
) -> ExpressionContrast | None:
  """Compare target log2(TPM+1) expression in loss versus intact models."""
  lost_gene = canonical_gene_name(lost_gene)
  target_gene = canonical_gene_name(target_gene)
  if lost_gene not in copy_number or target_gene not in expression:
    return None
  joined = (
    _lineage_join(copy_number, models, lost_gene, lineage, lineage_column)
      .merge(expression[["ModelID", target_gene]], on="ModelID", how="inner")
      .dropna(subset=[lost_gene, target_gene])
  )
  loss = joined.loc[joined[lost_gene] < loss_threshold, target_gene].astype(float)
  intact = joined.loc[joined[lost_gene] >= loss_threshold, target_gene].astype(float)
  if len(loss) < min_group_size or len(intact) < min_group_size:
    return None
  p_value = float(mannwhitneyu(loss, intact, alternative="greater").pvalue)
  median_loss = float(np.median(loss))
  median_intact = float(np.median(intact))
  return ExpressionContrast(
    lost_gene=lost_gene,
    target_gene=target_gene,
    lineage=lineage,
    n_loss=len(loss),
    n_intact=len(intact),
    median_expression_loss=median_loss,
    median_expression_intact=median_intact,
    delta_expression=median_loss - median_intact,
    p_value=p_value,
  )


def dependency_components(contrast: DependencyContrast) -> dict[str, float]:
  dependency = np.clip((-contrast.median_effect_loss - 0.25) / 0.75, 0, 1)
  selectivity = np.clip((-contrast.delta_effect) / 0.75, 0, 1)
  return {"dependency": float(dependency), "selectivity": float(selectivity)}


def expression_component(contrast: ExpressionContrast, saturation_log2: float = 2.0) -> float:
  """Normalize target induction; 2 log2(TPM+1) units saturates by default."""
  return float(np.clip(contrast.delta_expression / saturation_log2, 0, 1))
