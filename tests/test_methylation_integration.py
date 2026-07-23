from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from rses_onco.methylation import (
  pair_methylation_metrics,
  promoter_methylation_context_score,
)

ROOT = Path(__file__).resolve().parents[1]


def test_promoter_methylation_context_score_is_directional_and_bounded() -> None:
  assert promoter_methylation_context_score(0.8, 0.2) == pytest.approx(0.64)
  assert promoter_methylation_context_score(0.2, 0.8) == pytest.approx(0.04)
  assert promoter_methylation_context_score(1.5, -0.4) == pytest.approx(1.0)
  assert promoter_methylation_context_score(None, 0.2) is None


def test_pair_methylation_metrics_use_paired_cancer_samples() -> None:
  rows = []
  for index in range(12):
    sample_id = f"S{index:02d}"
    rows.extend([
      {
        "cancer": "colon",
        "sample_id": sample_id,
        "gene": "LOST",
        "promoter_beta": 0.72 + index * 0.005,
      },
      {
        "cancer": "colon",
        "sample_id": sample_id,
        "gene": "TARGET",
        "promoter_beta": 0.18 + index * 0.002,
      },
    ])
  frame = pd.DataFrame(rows)
  metrics = pair_methylation_metrics(
    frame,
    "LOST",
    "TARGET",
    "colon",
    min_samples=10,
  )
  assert metrics.evidence_status == "observed"
  assert metrics.n_samples == 12
  assert metrics.lost_median_beta is not None
  assert metrics.target_median_beta is not None
  assert metrics.lost_median_beta > metrics.target_median_beta
  assert metrics.promoter_methylation_context_score is not None
  assert 0 <= metrics.promoter_methylation_context_score <= 1
  assert metrics.paired_wilcoxon_p_value is not None


def ranking_row(*, tumor_event: float | None = 0.4) -> dict[str, object]:
  row: dict[str, object] = {
    "pair_id": "PAIR_A",
    "cancer": "colon",
    "score_comparability_group": "gene_pair",
    "functional_microniche_adjusted": 0.70,
    "component_tumor_event": tumor_event,
    "component_dependency": 0.60,
    "component_selectivity": 0.50,
    "component_expression_compensation": 0.80,
    "component_functional_relation": 0.90,
    "component_functional_microniche": 0.70,
    "component_validation_tractability": 0.30,
    "score_version": "RSES-Onco-expanded-v0.10.9",
    "expression_regulatory_semantics_version": (
      "eligibility-aware-wgcna-regulatory-v3"
    ),
  }
  for domain in (
    "tumor_event",
    "dependency",
    "selectivity",
    "expression_compensation",
    "functional_relation",
    "functional_microniche",
    "validation_tractability",
  ):
    row[f"eligible_component_{domain}"] = True
  return row


def methylation_row(score: float | None) -> dict[str, object]:
  return {
    "pair_id": "PAIR_A",
    "cancer": "colon",
    "evidence_status": "observed" if score is not None else "missing",
    "absence_reason": "" if score is not None else "no_eligible_probes",
    "n_samples": 25 if score is not None else 0,
    "lost_median_beta": 0.70 if score is not None else np.nan,
    "target_median_beta": 0.20 if score is not None else np.nan,
    "median_delta_beta": 0.50 if score is not None else np.nan,
    "methylation_spearman_rho": 0.10 if score is not None else np.nan,
    "promoter_methylation_context_score": score,
    "paired_wilcoxon_p_value": 0.001 if score is not None else np.nan,
    "paired_wilcoxon_q_value_bh": 0.01 if score is not None else np.nan,
    "paired_wilcoxon_q_value_bh_within_cancer": (
      0.005 if score is not None else np.nan
    ),
    "methylation_source": "GDC_SeSAMe_Methylation_Beta_Value",
    "methylation_workflow": "GDC_SeSAMe_Methylation_Array_Harmonization",
    "direct_gene_silencing_claim": False,
  }


def run_recompute(
  tmp_path: Path,
  ranking: pd.DataFrame,
  evidence: pd.DataFrame,
) -> pd.DataFrame:
  ranking_path = tmp_path / "ranking.tsv"
  evidence_path = tmp_path / "methylation.tsv"
  output_path = tmp_path / "output.tsv"
  ranking.to_csv(ranking_path, sep="\t", index=False)
  evidence.to_csv(evidence_path, sep="\t", index=False)
  subprocess.run([
    sys.executable,
    "scripts/recompute_rses_with_methylation.py",
    "--ranking",
    str(ranking_path),
    "--methylation-evidence",
    str(evidence_path),
    "--output",
    str(output_path),
  ], cwd=ROOT, check=True)
  return pd.read_csv(output_path, sep="\t", low_memory=False)


def test_methylation_is_internal_to_expression_compensation_weight(
  tmp_path: Path,
) -> None:
  result = run_recompute(
    tmp_path,
    pd.DataFrame([ranking_row()]),
    pd.DataFrame([methylation_row(0.50)]),
  )
  row = result.iloc[0]
  assert row["expression_compensation_expression_only"] == pytest.approx(0.80)
  assert row["component_expression_compensation"] == pytest.approx(0.71)
  assert row["expression_methylation_subcoverage"] == pytest.approx(1.0)
  assert row["coverage_adjusted_rses"] == pytest.approx(0.5428)
  assert row["score_version"] == "RSES-Onco-expanded-v0.11.1"
  assert row["methylation_score_role"] == "expression_compensation_sublayer"
  assert not bool(row["methylation_direct_silencing_claim"])


def test_missing_methylation_lowers_internal_coverage_without_zero_imputation(
  tmp_path: Path,
) -> None:
  result = run_recompute(
    tmp_path,
    pd.DataFrame([ranking_row()]),
    pd.DataFrame([methylation_row(None)]),
  )
  row = result.iloc[0]
  assert pd.isna(row["promoter_methylation_context_score"])
  assert row["expression_methylation_raw"] == pytest.approx(0.80)
  assert row["expression_methylation_subcoverage"] == pytest.approx(0.70)
  assert row["component_expression_compensation"] == pytest.approx(0.56)


def test_depmap_only_ranking_is_not_modified_by_tcga_methylation(
  tmp_path: Path,
) -> None:
  original = pd.DataFrame([ranking_row(tumor_event=None)])
  result = run_recompute(
    tmp_path,
    original,
    pd.DataFrame([methylation_row(0.50)]),
  )
  assert "promoter_methylation_context_score" not in result.columns
  assert result.iloc[0]["component_expression_compensation"] == pytest.approx(0.80)
  assert result.iloc[0]["score_version"] == "RSES-Onco-expanded-v0.10.9"


def test_methylation_acquisition_uses_official_gdc_annotation_ids() -> None:
  source = (
    ROOT / "scripts/download_gdc_methylation.py"
  ).read_text(encoding="utf-8")
  for uuid in (
    "5ce8ae8f-3386-4d12-9035-152742aa07e0",
    "e5182c42-bdc6-433e-9b4a-7b7c6696ce89",
    "021a2330-951d-474f-af24-1acd77e7664f",
  ):
    assert uuid in source
  assert "Methylation Beta Value" in source
  assert "Primary Tumor" in source


def test_pipeline_exposes_optional_methylation_modes_and_tables() -> None:
  pipeline = (
    ROOT / "scripts/resume_wgcna_regulatory_pipeline.sh"
  ).read_text(encoding="utf-8")
  materializer = (
    ROOT / "scripts/materialize_extended_supplementary_tables.py"
  ).read_text(encoding="utf-8")
  figures = (
    ROOT / "scripts/make_extended_supporting_figures.py"
  ).read_text(encoding="utf-8")
  assert "METHYLATION_MODE" in pipeline
  assert "methylation-only" in pipeline
  assert "recompute_rses_with_methylation.py" in pipeline
  for number in (45, 46, 47):
    assert f"Table_S{number}_" in materializer
  assert "pair_promoter_methylation_evidence.tsv" in figures
