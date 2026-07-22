#!/usr/bin/env python3
"""Build model-level expression, coexpression, compensation and CRISPR support tables.

Only gene-pair hypotheses with real DepMap measurements are evaluated. Composite
biomarkers are retained in an explicit non-eligible status table and are never
converted into single-gene events. Correlation is reported descriptively and is not
interpreted as compensation; compensation is operationally tested as target
expression in loss versus intact models.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from rses_onco.depmap import cancer_model_ids, read_depmap_inputs
from rses_onco.utils import bh_adjust, canonical_gene_name

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SCORING_SEMANTICS = "eligibility-aware-v1"


@dataclass(frozen=True)
class OutputRecord:
  evidence_family: str
  path: str
  rows: int
  columns: int
  sha256: str
  status: str
  interpretation_boundary: str


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> Path:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)
  return path


def safe_spearman(first: pd.Series, second: pd.Series) -> tuple[int, float | None, float | None]:
  table = pd.DataFrame({"first": first, "second": second}).apply(
    pd.to_numeric,
    errors="coerce",
  ).dropna()
  if len(table) < 3 or table["first"].nunique() < 2 or table["second"].nunique() < 2:
    return len(table), None, None
  result = spearmanr(table["first"], table["second"], nan_policy="omit")
  rho = float(result.statistic) if np.isfinite(result.statistic) else None
  p_value = float(result.pvalue) if np.isfinite(result.pvalue) else None
  return len(table), rho, p_value


def model_metadata_columns(models: pd.DataFrame) -> list[str]:
  preferred = (
    "ModelID", "CellLineName", "ModelType", "OncotreeLineage",
    "OncotreePrimaryDisease", "OncotreeSubtype", "Sex", "Age",
    "PrimaryOrMetastasis", "SangerModelID", "RRID",
  )
  return [column for column in preferred if column in models]


def pair_table(
  *,
  cancer: str,
  pair_id: str,
  lost_gene: str,
  target_gene: str,
  effect: pd.DataFrame,
  copy_number: pd.DataFrame,
  expression: pd.DataFrame,
  models: pd.DataFrame,
  loss_threshold: float,
) -> tuple[pd.DataFrame, str]:
  required_mapping = {
    "lost_copy_number": lost_gene in copy_number,
    "target_copy_number": target_gene in copy_number,
    "lost_expression": lost_gene in expression,
    "target_expression": target_gene in expression,
    "lost_gene_effect": lost_gene in effect,
    "target_gene_effect": target_gene in effect,
  }
  if not required_mapping["lost_copy_number"]:
    return pd.DataFrame(), "lost_gene_unmapped_in_copy_number"
  selected = set(cancer_model_ids(models, cancer).astype(str))
  if not selected:
    return pd.DataFrame(), "no_compatible_cancer_models"

  result = copy_number[["ModelID", lost_gene]].rename(
    columns={lost_gene: "lost_copy_number"}
  )
  if target_gene in copy_number:
    target_cn = copy_number[["ModelID", target_gene]].rename(
      columns={target_gene: "target_copy_number"}
    )
    result = result.merge(target_cn, on="ModelID", how="left")
  else:
    result["target_copy_number"] = np.nan

  for matrix, gene, name in (
    (expression, lost_gene, "lost_expression"),
    (expression, target_gene, "target_expression"),
    (effect, lost_gene, "lost_gene_effect"),
    (effect, target_gene, "target_gene_effect"),
  ):
    if gene in matrix:
      one = matrix[["ModelID", gene]].rename(columns={gene: name})
      result = result.merge(one, on="ModelID", how="left")
    else:
      result[name] = np.nan

  metadata = models[model_metadata_columns(models)].drop_duplicates("ModelID")
  result = result.merge(metadata, on="ModelID", how="left")
  result = result.loc[result["ModelID"].astype(str).isin(selected)].copy()
  for column in (
    "lost_copy_number", "target_copy_number", "lost_expression",
    "target_expression", "lost_gene_effect", "target_gene_effect",
  ):
    result[column] = pd.to_numeric(result[column], errors="coerce")
  result["loss_group"] = np.where(
    result["lost_copy_number"].isna(),
    "copy_number_unavailable",
    np.where(
      result["lost_copy_number"] < loss_threshold,
      "loss",
      "intact",
    ),
  )
  result["cancer"] = cancer
  result["pair_id"] = pair_id
  result["lost_gene"] = lost_gene
  result["target_gene"] = target_gene
  result["loss_threshold"] = loss_threshold
  result["mapping_lost_expression"] = required_mapping["lost_expression"]
  result["mapping_target_expression"] = required_mapping["target_expression"]
  result["mapping_lost_gene_effect"] = required_mapping["lost_gene_effect"]
  result["mapping_target_gene_effect"] = required_mapping["target_gene_effect"]
  if result.empty:
    return result, "no_models_after_cancer_and_identifier_join"
  return result, "available"


def correlation_rows(model_level: pd.DataFrame) -> list[dict[str, object]]:
  rows: list[dict[str, object]] = []
  for (cancer, pair_id, lost_gene, target_gene), group in model_level.groupby(
    ["cancer", "pair_id", "lost_gene", "target_gene"],
    dropna=False,
  ):
    for measurement, first, second in (
      ("expression", "lost_expression", "target_expression"),
      ("crispr_gene_effect", "lost_gene_effect", "target_gene_effect"),
    ):
      for stratum in ("all", "loss", "intact"):
        subset = group if stratum == "all" else group.loc[group["loss_group"].eq(stratum)]
        n, rho, p_value = safe_spearman(subset[first], subset[second])
        rows.append({
          "cancer": cancer,
          "pair_id": pair_id,
          "lost_gene": lost_gene,
          "target_gene": target_gene,
          "measurement": measurement,
          "event_stratum": stratum,
          "n_models": n,
          "spearman_rho": rho,
          "p_value": p_value,
          "analysis_status": "available" if rho is not None else "insufficient_or_constant_values",
          "interpretation_boundary": (
            "Correlation is descriptive and is not automatic evidence of transcriptional compensation or synthetic lethality."
          ),
        })
  return rows


def contrast_rows(model_level: pd.DataFrame) -> list[dict[str, object]]:
  rows: list[dict[str, object]] = []
  for (cancer, pair_id, lost_gene, target_gene), group in model_level.groupby(
    ["cancer", "pair_id", "lost_gene", "target_gene"],
    dropna=False,
  ):
    loss = group.loc[group["loss_group"].eq("loss")]
    intact = group.loc[group["loss_group"].eq("intact")]
    for analysis, column, alternative, supportive_direction in (
      (
        "transcriptional_compensation",
        "target_expression",
        "greater",
        "higher_target_expression_in_loss_group",
      ),
      (
        "conditional_dependency",
        "target_gene_effect",
        "less",
        "more_negative_target_gene_effect_in_loss_group",
      ),
    ):
      loss_values = pd.to_numeric(loss[column], errors="coerce").dropna()
      intact_values = pd.to_numeric(intact[column], errors="coerce").dropna()
      status = "available"
      p_value = None
      if len(loss_values) < 3 or len(intact_values) < 3:
        status = "insufficient_group_size"
      else:
        p_value = float(
          mannwhitneyu(
            loss_values,
            intact_values,
            alternative=alternative,
          ).pvalue
        )
      median_loss = float(np.median(loss_values)) if len(loss_values) else None
      median_intact = float(np.median(intact_values)) if len(intact_values) else None
      delta = (
        median_loss - median_intact
        if median_loss is not None and median_intact is not None
        else None
      )
      rows.append({
        "cancer": cancer,
        "pair_id": pair_id,
        "lost_gene": lost_gene,
        "target_gene": target_gene,
        "analysis": analysis,
        "measurement": column,
        "n_loss": len(loss_values),
        "n_intact": len(intact_values),
        "median_loss": median_loss,
        "median_intact": median_intact,
        "delta_loss_minus_intact": delta,
        "p_value": p_value,
        "analysis_status": status,
        "supportive_direction": supportive_direction,
        "interpretation_boundary": (
          "The comparison is operational and context-specific; significance does not establish clinical efficacy."
        ),
      })
  return rows


def add_fdr(frame: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
  if frame.empty or "p_value" not in frame:
    return frame
  frame = frame.copy()
  frame["q_value_bh"] = np.nan
  for _, indices in frame.groupby(group_columns, dropna=False).groups.items():
    index = list(indices)
    values = pd.to_numeric(frame.loc[index, "p_value"], errors="coerce")
    valid = values.notna()
    if valid.any():
      frame.loc[np.array(index)[valid.to_numpy()], "q_value_bh"] = bh_adjust(
        values.loc[valid]
      )
  return frame


def output_record(
  family: str,
  path: Path,
  frame: pd.DataFrame,
  boundary: str,
) -> OutputRecord:
  return OutputRecord(
    evidence_family=family,
    path=str(path),
    rows=len(frame),
    columns=len(frame.columns),
    sha256=sha256(path),
    status="available" if not frame.empty else "empty_no_eligible_records",
    interpretation_boundary=boundary,
  )


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--gene-effect", required=True)
  parser.add_argument("--copy-number", required=True)
  parser.add_argument("--models", required=True)
  parser.add_argument("--expression", required=True)
  parser.add_argument(
    "--ranking",
    default="results/expanded_26Q1/full/expanded_rses_onco.tsv",
  )
  parser.add_argument("--output-root", default="article_outputs")
  parser.add_argument("--loss-threshold", type=float, default=0.30)
  args = parser.parse_args()

  ranking_path = resolve_path(args.ranking)
  ranking = pd.read_csv(ranking_path, sep="\t", low_memory=False)
  required = {
    "cancer", "pair_id", "analysis_lost_gene", "analysis_target_gene",
    "hypothesis_type", "scoring_semantics_version",
  }
  missing = sorted(required - set(ranking.columns))
  if missing:
    raise ValueError(
      "Ranking lacks model-level eligibility fields; rerun expanded scoring. "
      f"Missing: {missing}"
    )
  versions = set(ranking["scoring_semantics_version"].dropna().astype(str))
  if versions != {EXPECTED_SCORING_SEMANTICS}:
    raise ValueError(
      f"Expected {EXPECTED_SCORING_SEMANTICS}; observed {sorted(versions)}"
    )

  effect, copy_number, models, expression = read_depmap_inputs(
    resolve_path(args.gene_effect),
    resolve_path(args.copy_number),
    resolve_path(args.models),
    resolve_path(args.expression),
  )
  if expression is None:
    raise ValueError("Expression matrix is required")

  contexts = ranking[[
    "cancer", "pair_id", "analysis_lost_gene", "analysis_target_gene",
    "hypothesis_type", "source_class", "score_comparability_group",
  ]].drop_duplicates()
  model_frames: list[pd.DataFrame] = []
  status_rows: list[dict[str, object]] = []
  for record in contexts.to_dict("records"):
    cancer = str(record["cancer"])
    pair_id = str(record["pair_id"])
    lost_gene = canonical_gene_name(record.get("analysis_lost_gene"))
    target_gene = canonical_gene_name(record.get("analysis_target_gene"))
    if record.get("hypothesis_type") != "gene" or not lost_gene or not target_gene:
      status_rows.append({
        **record,
        "analysis_status": "not_eligible",
        "absence_reason": "composite_event_not_executable_as_single_gene_model_level_analysis",
      })
      continue
    table, status = pair_table(
      cancer=cancer,
      pair_id=pair_id,
      lost_gene=lost_gene,
      target_gene=target_gene,
      effect=effect,
      copy_number=copy_number,
      expression=expression,
      models=models,
      loss_threshold=args.loss_threshold,
    )
    status_rows.append({
      **record,
      "analysis_status": status,
      "absence_reason": "" if status == "available" else status,
      "model_rows": len(table),
    })
    if not table.empty:
      table["source_class"] = record.get("source_class")
      table["score_comparability_group"] = record.get(
        "score_comparability_group"
      )
      model_frames.append(table)

  model_level = (
    pd.concat(model_frames, ignore_index=True, sort=False)
    if model_frames
    else pd.DataFrame()
  )
  correlations = add_fdr(
    pd.DataFrame(correlation_rows(model_level)),
    ["measurement", "event_stratum", "cancer"],
  )
  contrasts = add_fdr(
    pd.DataFrame(contrast_rows(model_level)),
    ["analysis", "cancer"],
  )
  status = pd.DataFrame(status_rows)

  output_root = resolve_path(args.output_root)
  evidence_root = output_root / "tables" / "supporting_evidence"
  paths = {
    "model_level": evidence_root / "model_level/model_level_expression_crispr_copy_number.tsv",
    "correlations": evidence_root / "expression/coexpression_by_event_group.tsv",
    "contrasts": evidence_root / "expression/compensation_and_dependency_contrasts.tsv",
    "status": evidence_root / "model_level/model_level_analysis_status.tsv",
  }
  for key, frame in (
    ("model_level", model_level),
    ("correlations", correlations),
    ("contrasts", contrasts),
    ("status", status),
  ):
    atomic_tsv(frame, paths[key])

  records = [
    output_record(
      "model_level_expression_crispr_copy_number",
      paths["model_level"],
      model_level,
      "Model-level values are measurements, not independent biological studies.",
    ),
    output_record(
      "coexpression_by_event_group",
      paths["correlations"],
      correlations,
      "Correlation is descriptive and is not automatic evidence of compensation.",
    ),
    output_record(
      "compensation_and_dependency_contrasts",
      paths["contrasts"],
      contrasts,
      "Compensation and conditional dependency are operational cancer-specific group comparisons.",
    ),
    output_record(
      "model_level_analysis_status",
      paths["status"],
      status,
      "Non-eligibility, mapping failure and insufficient samples are explicit.",
    ),
  ]
  manifest = pd.DataFrame([asdict(record) for record in records])
  manifest["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
  manifest_path = evidence_root / "model_level_supporting_evidence_manifest.tsv"
  atomic_tsv(manifest, manifest_path)
  for path in [*paths.values(), manifest_path]:
    if not path.exists() or path.stat().st_size == 0:
      raise RuntimeError(f"Model-level supporting output missing or empty: {path}")
  print(f"Model-level rows: {len(model_level):,}")
  print(f"Coexpression rows: {len(correlations):,}")
  print(f"Contrast rows: {len(contrasts):,}")
  print(f"Wrote {manifest_path}")


if __name__ == "__main__":
  main()
