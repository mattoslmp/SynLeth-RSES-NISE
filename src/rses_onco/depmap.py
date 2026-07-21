from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .utils import canonical_gene_name


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


def _rename_gene_columns(frame: pd.DataFrame) -> pd.DataFrame:
  renamed = {c: canonical_gene_name(c) for c in frame.columns if c != "ModelID"}
  return frame.rename(columns=renamed)


def read_depmap_inputs(
  gene_effect_path: str | Path,
  copy_number_path: str | Path,
  model_path: str | Path,
  expression_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
  """Read DepMap matrices using ModelID as the stable join key."""
  effect = _rename_gene_columns(pd.read_csv(gene_effect_path))
  copy_number = _rename_gene_columns(pd.read_csv(copy_number_path))
  models = pd.read_csv(model_path)
  expression = _rename_gene_columns(pd.read_csv(expression_path)) if expression_path else None
  for name, frame in {"gene_effect": effect, "copy_number": copy_number}.items():
    if "ModelID" not in frame:
      raise ValueError(f"{name} requires a ModelID column")
  if "ModelID" not in models:
    raise ValueError("Model.csv requires a ModelID column")
  return effect, copy_number, models, expression


def dependency_contrast(
  effect: pd.DataFrame,
  copy_number: pd.DataFrame,
  models: pd.DataFrame,
  lost_gene: str,
  target_gene: str,
  lineage: str,
  lineage_column: str = "OncotreeLineage",
  loss_threshold: float = -1.5,
  min_group_size: int = 3,
) -> DependencyContrast | None:
  """Compare target CRISPR gene effect in loss versus intact cell lines."""
  lost_gene = canonical_gene_name(lost_gene)
  target_gene = canonical_gene_name(target_gene)
  if lost_gene not in copy_number or target_gene not in effect:
    return None
  metadata = models[["ModelID", lineage_column]].copy()
  metadata[lineage_column] = metadata[lineage_column].fillna("").astype(str)
  ids = metadata.loc[metadata[lineage_column].str.casefold() == lineage.casefold(), "ModelID"]
  joined = (
    metadata.loc[metadata.ModelID.isin(ids)]
      .merge(copy_number[["ModelID", lost_gene]], on="ModelID", how="inner")
      .merge(effect[["ModelID", target_gene]], on="ModelID", how="inner")
      .dropna(subset=[lost_gene, target_gene])
  )
  loss = joined.loc[joined[lost_gene] <= loss_threshold, target_gene].astype(float)
  intact = joined.loc[joined[lost_gene] > loss_threshold, target_gene].astype(float)
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


def dependency_components(contrast: DependencyContrast) -> dict[str, float]:
  """Normalize a dependency contrast to interpretable score components."""
  dependency = np.clip((-contrast.median_effect_loss - 0.25) / 0.75, 0, 1)
  selectivity = np.clip((-contrast.delta_effect) / 0.75, 0, 1)
  return {"dependency": float(dependency), "selectivity": float(selectivity)}
