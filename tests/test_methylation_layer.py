from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from rses_onco.methylation import (
  METHYLATION_SUBWEIGHTS,
  build_methylation_pair_metrics,
  read_promoter_methylation,
)
from scripts.integrate_methylation_regulatory_layer import (
  REGULATORY_METHYLATION_SUBWEIGHTS,
)

ROOT = Path(__file__).resolve().parents[1]


def synthetic_models() -> pd.DataFrame:
  return pd.DataFrame({
    "ModelID": [f"ACH-{index:06d}" for index in range(1, 13)],
    "OncotreeLineage": ["Bowel"] * 4 + ["Esophagus/Stomach"] * 4 + ["Lung"] * 4,
    "OncotreeSubtype": [
      "Colon Adenocarcinoma"
    ] * 4 + ["Stomach Adenocarcinoma"] * 4 + ["NSCLC"] * 4,
    "OncotreePrimaryDisease": [
      "Colorectal adenocarcinoma"
    ] * 4 + ["Gastric adenocarcinoma"] * 4 + ["Lung cancer"] * 4,
    "OncotreeCode": ["COAD"] * 4 + ["STAD"] * 4 + ["LUAD"] * 4,
  })


def test_methylation_weights_share_existing_regulatory_domain() -> None:
  assert np.isclose(sum(METHYLATION_SUBWEIGHTS.values()), 1.0)
  assert np.isclose(sum(REGULATORY_METHYLATION_SUBWEIGHTS.values()), 1.0)
  assert np.isclose(
    REGULATORY_METHYLATION_SUBWEIGHTS["promoter_methylation_context"],
    0.20,
  )


def test_promoter_methylation_reader_collapses_multiple_tss_features(
  tmp_path: Path,
) -> None:
  models = synthetic_models()
  frame = pd.DataFrame({
    "ModelID": models["ModelID"],
    "GENEA (1 1-1000)": [0.8, 0.7, 0.2, 0.1] * 3,
    "GENEB (2 1-1000)": [0.2, 0.3, 0.7, 0.8] * 3,
    "GENEB (2 2-1001)": [0.22, 0.32, 0.72, 0.82] * 3,
  })
  path = tmp_path / "methylation.csv"
  frame.to_csv(path, index=False)
  loaded = read_promoter_methylation(path, models)
  assert loaded.promoter_feature_counts["GENEB"] == 2
  assert {"ModelID", "GENEA", "GENEB"}.issubset(loaded.matrix.columns)
  assert loaded.diagnostics["missing_data_rule"] == "preserved_as_NA"
  assert loaded.matrix["GENEB"].between(0, 1).all()


def test_methylation_metrics_distinguish_missing_from_negative(
  tmp_path: Path,
) -> None:
  models = synthetic_models()
  frame = pd.DataFrame({
    "ModelID": models["ModelID"],
    "GENEA (1 1-1000)": [0.8, 0.7, 0.2, 0.1] * 3,
    "GENEB (2 1-1000)": [0.2, 0.3, 0.7, 0.8] * 3,
  })
  path = tmp_path / "methylation.csv"
  frame.to_csv(path, index=False)
  methylation = read_promoter_methylation(path, models)
  copy_number = pd.DataFrame({
    "ModelID": models["ModelID"],
    "GENEA": [0.1, 0.2, 1.0, 1.0] * 3,
  })
  pairs = pd.DataFrame([
    {"pair_id": "P1", "lost_gene": "GENEA", "target_gene": "GENEB"},
    {"pair_id": "P2", "lost_gene": "GENEA", "target_gene": "MISSING"},
  ])
  result = build_methylation_pair_metrics(
    methylation,
    copy_number,
    models,
    pairs,
    min_group_size=2,
  )
  observed = result.loc[result["pair_id"].eq("P1")]
  missing = result.loc[result["pair_id"].eq("P2")]
  assert observed["component_promoter_methylation_context"].notna().all()
  assert missing["component_promoter_methylation_context"].isna().all()
  assert missing["methylation_absence_reason"].str.contains("unmapped").all()


def test_integrator_recomputes_regulatory_coverage_without_double_counting(
  tmp_path: Path,
) -> None:
  models = synthetic_models()
  models_path = tmp_path / "Model.csv"
  models.to_csv(models_path, index=False)
  methylation = pd.DataFrame({
    "ModelID": models["ModelID"],
    "GENEA (1 1-1000)": [0.8, 0.7, 0.2, 0.1] * 3,
    "GENEB (2 1-1000)": [0.2, 0.3, 0.7, 0.8] * 3,
  })
  methylation_path = tmp_path / "methylation.csv"
  methylation.to_csv(methylation_path, index=False)
  copy_number = pd.DataFrame({
    "ModelID": models["ModelID"],
    "GENEA": [0.1, 0.2, 1.0, 1.0] * 3,
  })
  copy_path = tmp_path / "copy.csv"
  copy_number.to_csv(copy_path, index=False)
  candidates = pd.DataFrame([
    {"pair_id": "P1", "lost_gene": "GENEA", "target_gene": "GENEB"},
  ])
  candidate_path = tmp_path / "candidates.tsv"
  candidates.to_csv(candidate_path, sep="\t", index=False)
  evidence = pd.DataFrame([
    {
      "pair_id": "P1",
      "cancer": cancer,
      "regulatory_tf_association_divergence": 0.5,
      "regulatory_tf_expression_profile_divergence": 0.5,
      "regulatory_promoter_motif_divergence": 0.5,
      "promoter_evidence_type": "JASPAR_motif_prediction_not_direct_binding",
    }
    for cancer in ("colon", "stomach", "lung")
  ])
  input_path = tmp_path / "evidence.tsv"
  output_path = tmp_path / "output.tsv"
  metrics_path = tmp_path / "metrics.tsv"
  status_path = tmp_path / "status.json"
  evidence.to_csv(input_path, sep="\t", index=False)
  subprocess.run(
    [
      sys.executable,
      "scripts/integrate_methylation_regulatory_layer.py",
      "--methylation",
      str(methylation_path),
      "--copy-number",
      str(copy_path),
      "--models",
      str(models_path),
      "--candidates",
      str(candidate_path),
      "--input",
      str(input_path),
      "--output",
      str(output_path),
      "--metrics-output",
      str(metrics_path),
      "--status-output",
      str(status_path),
      "--min-group-size",
      "2",
    ],
    cwd=ROOT,
    check=True,
  )
  output = pd.read_csv(output_path, sep="\t")
  assert output["component_promoter_methylation_context"].notna().all()
  assert output["regulatory_network_coverage"].eq(1.0).all()
  assert output["regulatory_layer_version"].eq(
    "wgcna-promoter-methylation-regulatory-v3"
  ).all()


def test_resume_pipeline_integrates_methylation_as_optional_input() -> None:
  implementation = (
    ROOT / "scripts/resume_wgcna_regulatory_pipeline.sh"
  ).read_text(encoding="utf-8")
  assert "integrate_methylation_regulatory_layer.py" in implementation
  assert "METHYLATION=" in implementation
  assert "eligible-missing" in implementation
  assert "--methylation" in implementation
