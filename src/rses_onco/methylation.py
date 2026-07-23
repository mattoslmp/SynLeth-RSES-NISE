from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from .depmap import cancer_model_ids, normalize_model_id_column
from .utils import bh_adjust, canonical_gene_name


METHYLATION_SUBWEIGHTS = {
  "pair_profile_divergence": 0.50,
  "conditional_target_hypomethylation": 0.50,
}

MODEL_LABEL_ALIASES = (
  "ModelID", "DepMap_ID", "DepMapID", "DepMap ID", "model_id",
  "CCLEName", "CCLE_Name", "CCLE name", "cell_line", "CellLine",
  "Cell Line",
)
GENE_LABEL_ALIASES = (
  "gene", "Gene", "gene_symbol", "GeneSymbol", "symbol", "HugoSymbol",
  "feature", "Feature", "locus", "Locus",
)
VALUE_ALIASES = (
  "methylation", "Methylation", "beta_value", "Beta Value", "value", "Value",
)


@dataclass(frozen=True)
class MethylationMatrix:
  matrix: pd.DataFrame
  promoter_feature_counts: dict[str, int]
  diagnostics: dict[str, object]


def _read_delimited(path: Path) -> pd.DataFrame:
  suffixes = "".join(path.suffixes).casefold()
  if suffixes.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
  if suffixes.endswith((".csv", ".csv.gz")):
    return pd.read_csv(path, compression="infer", low_memory=False)
  return pd.read_csv(
    path,
    sep=None,
    engine="python",
    compression="infer",
    low_memory=False,
  )


def _first_present(columns: Iterable[object], aliases: Iterable[str]) -> object | None:
  lookup = {str(column).strip().casefold(): column for column in columns}
  for alias in aliases:
    match = lookup.get(alias.strip().casefold())
    if match is not None:
      return match
  return None


def promoter_gene_name(value: object) -> str:
  """Extract the leading gene symbol from a promoter-feature label."""
  text = str(value).strip()
  if not text or text.casefold() in {"nan", "none"}:
    return ""
  text = re.sub(r"^[^A-Za-z0-9]+", "", text)
  token = re.split(r"[\s(|;,:]+", text, maxsplit=1)[0]
  return canonical_gene_name(token)


def _model_alias_map(models: pd.DataFrame) -> dict[str, str]:
  models = normalize_model_id_column(models.copy(), "Model.csv")
  aliases = {
    "ModelID", "CCLEName", "CCLE_Name", "StrippedCellLineName",
    "CellLineName", "SangerModelID", "COSMICID",
  }
  mapping: dict[str, str] = {}
  for record in models.to_dict("records"):
    model_id = str(record.get("ModelID") or "").strip()
    if not model_id:
      continue
    mapping[model_id.casefold()] = model_id
    for column in aliases:
      value = str(record.get(column) or "").strip()
      if value:
        mapping[value.casefold()] = model_id
  return mapping


def _collapse_feature_columns(
  frame: pd.DataFrame,
  id_column: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
  groups: dict[str, list[object]] = {}
  for column in frame.columns:
    if column == id_column:
      continue
    gene = promoter_gene_name(column)
    if gene:
      groups.setdefault(gene, []).append(column)
  output = pd.DataFrame({"ModelID": frame[id_column].astype(str).str.strip()})
  counts: dict[str, int] = {}
  for gene, columns in sorted(groups.items()):
    values = frame[columns].apply(pd.to_numeric, errors="coerce")
    output[gene] = values.median(axis=1, skipna=True)
    counts[gene] = len(columns)
  return output, counts


def _normalize_model_rows(
  frame: pd.DataFrame,
  models: pd.DataFrame,
  id_column: object,
) -> tuple[pd.DataFrame, dict[str, int]]:
  alias_map = _model_alias_map(models)
  raw = frame[id_column].astype(str).str.strip()
  working = frame.copy()
  working["__ModelID__"] = raw.map(
    lambda value: alias_map.get(value.casefold(), value)
  )
  working = working.loc[
    working["__ModelID__"].astype(str).str.match(r"^ACH-", na=False)
  ]
  output, counts = _collapse_feature_columns(working, "__ModelID__")
  output = output.groupby("ModelID", as_index=False).median(numeric_only=True)
  return output, counts


def _normalize_feature_rows(
  frame: pd.DataFrame,
  models: pd.DataFrame,
  feature_column: object,
) -> tuple[pd.DataFrame, dict[str, int]]:
  alias_map = _model_alias_map(models)
  model_columns: list[object] = []
  model_ids: list[str] = []
  for column in frame.columns:
    if column == feature_column:
      continue
    mapped = alias_map.get(str(column).strip().casefold())
    if mapped:
      model_columns.append(column)
      model_ids.append(mapped)
  if not model_columns:
    raise ValueError(
      "Could not map methylation columns to DepMap ModelID values; "
      "provide a ModelID-indexed custom download or matching Model.csv aliases"
    )
  rows = frame[[feature_column, *model_columns]].copy()
  rows["__gene__"] = rows[feature_column].map(promoter_gene_name)
  rows = rows.loc[rows["__gene__"].ne("")]
  counts = rows.groupby("__gene__").size().astype(int).to_dict()
  numeric = rows[model_columns].apply(pd.to_numeric, errors="coerce")
  numeric["__gene__"] = rows["__gene__"].to_numpy()
  collapsed = numeric.groupby("__gene__").median(numeric_only=True)
  collapsed.columns = model_ids
  output = collapsed.T
  output.index.name = "ModelID"
  output = output.groupby(level=0).median(numeric_only=True).reset_index()
  return output, counts


def read_promoter_methylation(
  path: str | Path,
  models: pd.DataFrame,
) -> MethylationMatrix:
  """Read DepMap/CCLE 1-kb-upstream-TSS methylation in common layouts.

  Model-by-feature, long ModelID/gene/value, and historical feature-by-cell-line
  layouts are supported. Multiple promoter/TSS features assigned to the same gene
  are collapsed by the median within each model. Missing values remain missing.
  """
  path = Path(path)
  if not path.exists() or path.stat().st_size == 0:
    raise FileNotFoundError(path)
  frame = _read_delimited(path)
  if frame.empty:
    raise ValueError(f"Methylation file is empty: {path}")

  model_column = _first_present(frame.columns, MODEL_LABEL_ALIASES)
  gene_column = _first_present(frame.columns, GENE_LABEL_ALIASES)
  value_column = _first_present(frame.columns, VALUE_ALIASES)

  if model_column is not None and gene_column is not None and value_column is not None:
    alias_map = _model_alias_map(models)
    long = frame[[model_column, gene_column, value_column]].copy()
    long["ModelID"] = long[model_column].astype(str).str.strip().map(
      lambda value: alias_map.get(value.casefold(), value)
    )
    long["gene"] = long[gene_column].map(promoter_gene_name)
    long["methylation"] = pd.to_numeric(long[value_column], errors="coerce")
    long = long.loc[
      long["ModelID"].astype(str).str.match(r"^ACH-", na=False)
      & long["gene"].ne("")
    ]
    counts = long.groupby("gene").size().astype(int).to_dict()
    output = long.pivot_table(
      index="ModelID",
      columns="gene",
      values="methylation",
      aggfunc="median",
    ).reset_index()
    layout = "long_model_gene_value"
  elif model_column is not None:
    output, counts = _normalize_model_rows(frame, models, model_column)
    layout = "model_rows_feature_columns"
  else:
    feature_column = gene_column if gene_column is not None else frame.columns[0]
    output, counts = _normalize_feature_rows(frame, models, feature_column)
    layout = "feature_rows_model_columns"

  gene_columns = [column for column in output.columns if column != "ModelID"]
  for column in gene_columns:
    output[column] = pd.to_numeric(output[column], errors="coerce")
  values = output[gene_columns].to_numpy(dtype=float, copy=False)
  finite = values[np.isfinite(values)]
  if finite.size:
    minimum = float(finite.min())
    maximum = float(finite.max())
    if minimum < -1e-8 or maximum > 1 + 1e-8:
      raise ValueError(
        "Promoter methylation values must be beta-like ratios in [0,1]; "
        f"observed range {minimum:.6g} to {maximum:.6g}"
      )
  diagnostics = {
    "source_path": str(path),
    "layout": layout,
    "models": int(output["ModelID"].nunique()),
    "genes": len(gene_columns),
    "observed_values": int(output[gene_columns].notna().sum().sum()),
    "multiple_promoter_gene_count": int(sum(value > 1 for value in counts.values())),
    "aggregation": "median_across_multiple_1kb_upstream_TSS_features_per_gene",
    "missing_data_rule": "preserved_as_NA",
    "methylation_scale": "promoter_methylation_ratio_0_to_1",
  }
  return MethylationMatrix(output, counts, diagnostics)


def _coverage_adjusted_subscore(
  values: dict[str, float | None],
  weights: dict[str, float],
) -> tuple[float | None, float, float | None, int]:
  observed = [
    name
    for name, value in values.items()
    if value is not None and np.isfinite(value)
  ]
  observed_weight = sum(weights[name] for name in observed)
  total_weight = sum(weights.values())
  if not observed_weight or not total_weight:
    return None, 0.0, None, 0
  raw = sum(
    weights[name] * float(np.clip(values[name], 0, 1))
    for name in observed
  ) / observed_weight
  coverage = observed_weight / total_weight
  return float(raw), float(coverage), float(raw * coverage), len(observed)


def build_methylation_pair_metrics(
  methylation: MethylationMatrix,
  copy_number: pd.DataFrame,
  models: pd.DataFrame,
  pairs: pd.DataFrame,
  cancers: Iterable[str] = ("colon", "stomach", "lung"),
  loss_threshold: float = 0.30,
  min_group_size: int = 3,
  difference_saturation: float = 0.25,
) -> pd.DataFrame:
  """Build cancer-specific epigenetic microniche metrics.

  Pair-profile divergence combines Spearman divergence and median absolute beta
  difference. Conditional support measures lower target-promoter methylation in
  lost-gene-loss models relative to intact models. It is an association and not
  proof of causal epigenetic silencing or derepression.
  """
  matrix = methylation.matrix.copy()
  copy_number = normalize_model_id_column(copy_number.copy(), "copy_number")
  rows: list[dict[str, object]] = []
  for cancer in cancers:
    selected = set(cancer_model_ids(models, cancer).astype(str))
    cancer_methylation = matrix.loc[matrix["ModelID"].astype(str).isin(selected)]
    cancer_copy = copy_number.loc[copy_number["ModelID"].astype(str).isin(selected)]
    for record in pairs.to_dict("records"):
      pair_id = str(record.get("pair_id") or "")
      lost = canonical_gene_name(record.get("lost_gene"))
      target = canonical_gene_name(record.get("target_gene"))
      profile_divergence: float | None = None
      profile_rho: float | None = None
      profile_mad: float | None = None
      profile_n = 0
      profile_reason = ""

      if lost not in cancer_methylation.columns or target not in cancer_methylation.columns:
        profile_reason = "lost_or_target_gene_unmapped_in_methylation"
      else:
        profile = cancer_methylation[[lost, target]].apply(
          pd.to_numeric,
          errors="coerce",
        ).dropna()
        profile_n = len(profile)
        if profile_n < 3:
          profile_reason = "insufficient_complete_methylation_profiles"
        else:
          rho = spearmanr(
            profile[lost],
            profile[target],
            nan_policy="omit",
          ).statistic
          if np.isfinite(rho):
            profile_rho = float(rho)
          profile_mad = float(np.median(np.abs(profile[lost] - profile[target])))
          available = []
          if profile_rho is not None:
            available.append((1.0 - profile_rho) / 2.0)
          available.append(
            float(np.clip(profile_mad / difference_saturation, 0, 1))
          )
          profile_divergence = float(np.mean(available)) if available else None

      conditional_support: float | None = None
      n_loss = 0
      n_intact = 0
      median_loss: float | None = None
      median_intact: float | None = None
      delta: float | None = None
      p_value: float | None = None
      conditional_reason = ""
      if target not in cancer_methylation.columns:
        conditional_reason = "target_gene_unmapped_in_methylation"
      elif lost not in cancer_copy.columns:
        conditional_reason = "lost_gene_unmapped_in_copy_number"
      else:
        table = cancer_copy[["ModelID", lost]].merge(
          cancer_methylation[["ModelID", target]],
          on="ModelID",
          how="inner",
        )
        table[lost] = pd.to_numeric(table[lost], errors="coerce")
        table[target] = pd.to_numeric(table[target], errors="coerce")
        table = table.dropna(subset=[lost, target])
        loss_values = table.loc[table[lost] < loss_threshold, target]
        intact_values = table.loc[table[lost] >= loss_threshold, target]
        n_loss = len(loss_values)
        n_intact = len(intact_values)
        if n_loss < min_group_size or n_intact < min_group_size:
          conditional_reason = "insufficient_loss_or_intact_group_size"
        else:
          median_loss = float(loss_values.median())
          median_intact = float(intact_values.median())
          delta = median_loss - median_intact
          conditional_support = float(
            np.clip((-delta) / difference_saturation, 0, 1)
          )
          p_value = float(
            mannwhitneyu(
              loss_values,
              intact_values,
              alternative="two-sided",
            ).pvalue
          )

      raw, coverage, adjusted, observed = _coverage_adjusted_subscore(
        {
          "pair_profile_divergence": profile_divergence,
          "conditional_target_hypomethylation": conditional_support,
        },
        METHYLATION_SUBWEIGHTS,
      )
      reasons = [reason for reason in (profile_reason, conditional_reason) if reason]
      rows.append({
        "pair_id": pair_id,
        "cancer": cancer,
        "lost_gene": lost,
        "target_gene": target,
        "methylation_profile_n_models": profile_n,
        "methylation_pair_spearman_rho": profile_rho,
        "methylation_pair_median_absolute_difference": profile_mad,
        "methylation_pair_profile_divergence": profile_divergence,
        "methylation_n_loss": n_loss,
        "methylation_n_intact": n_intact,
        "methylation_target_median_loss": median_loss,
        "methylation_target_median_intact": median_intact,
        "methylation_target_delta_loss_minus_intact": delta,
        "methylation_target_hypomethylation_support": conditional_support,
        "methylation_p_value": p_value,
        "methylation_raw": raw,
        "methylation_coverage": coverage,
        "component_promoter_methylation_context": adjusted,
        "methylation_observed_subcomponents": observed,
        "methylation_lost_promoter_feature_count": (
          methylation.promoter_feature_counts.get(lost, 0)
        ),
        "methylation_target_promoter_feature_count": (
          methylation.promoter_feature_counts.get(target, 0)
        ),
        "methylation_absence_reason": ";".join(dict.fromkeys(reasons)),
        "methylation_evidence_type": (
          "CCLE_RRBS_weighted_1kb_upstream_TSS_promoter_methylation"
        ),
        "methylation_interpretation": (
          "profile_divergence_plus_conditional_target_hypomethylation;"
          "association_not_causal_silencing_proof"
        ),
      })
  result = pd.DataFrame(rows)
  if not result.empty:
    result["methylation_q_value_bh"] = bh_adjust(result["methylation_p_value"])
    result["methylation_q_value_bh_within_cancer"] = (
      result.groupby("cancer", group_keys=False)["methylation_p_value"]
        .transform(lambda values: bh_adjust(values))
    )
  return result
