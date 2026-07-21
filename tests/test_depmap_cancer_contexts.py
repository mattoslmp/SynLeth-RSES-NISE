from __future__ import annotations

import pandas as pd

from rses_onco.depmap import cancer_model_ids, cancer_model_mask


def _models() -> pd.DataFrame:
  return pd.DataFrame({
    "ModelID": [
      "ACH-COAD", "ACH-READ", "ACH-APPENDIX", "ACH-STAD",
      "ACH-ESCA", "ACH-LUNG", "ACH-BREAST",
    ],
    "OncotreeLineage": [
      "Bowel", "Bowel", "Bowel", "Esophagus/Stomach",
      "Esophagus/Stomach", "Lung", "Breast",
    ],
    "OncotreePrimaryDisease": [
      "Colorectal Adenocarcinoma", "Colorectal Adenocarcinoma",
      "Appendiceal Adenocarcinoma", "Esophagogastric Adenocarcinoma",
      "Esophagogastric Adenocarcinoma", "Non-Small Cell Lung Cancer",
      "Invasive Breast Carcinoma",
    ],
    "OncotreeSubtype": [
      "Colon Adenocarcinoma", "Rectal Adenocarcinoma",
      "Appendiceal Adenocarcinoma", "Stomach Adenocarcinoma",
      "Esophageal Adenocarcinoma", "Lung Adenocarcinoma",
      "Breast Invasive Ductal Carcinoma",
    ],
    "OncotreeCode": ["COAD", "READ", "APAD", "STAD", "ESCA", "LUAD", "IDC"],
    "DepmapModelType": ["COAD", "READ", "APAD", "STAD", "ESCA", "LUAD", "IDC"],
  })


def test_colorectal_context_uses_current_oncotree_labels() -> None:
  ids = set(cancer_model_ids(_models(), "colon"))
  assert ids == {"ACH-COAD", "ACH-READ"}


def test_gastric_context_excludes_esophageal_models() -> None:
  ids = set(cancer_model_ids(_models(), "stomach"))
  assert ids == {"ACH-STAD"}


def test_lung_context_is_lineage_based() -> None:
  mask = cancer_model_mask(_models(), "lung")
  assert _models().loc[mask, "ModelID"].tolist() == ["ACH-LUNG"]
