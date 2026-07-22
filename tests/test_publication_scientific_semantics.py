from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.publication_scientific_semantics import (
  add_display_pair_columns,
  compound_resolution,
  present_text,
  resolved_pair_label,
  resolved_pair_parts,
)


def test_present_text_rejects_nan_like_values() -> None:
  for value in (None, np.nan, pd.NA, "", "nan", "NA", "None", "null"):
    assert present_text(value) == ""
  assert present_text("CHEMBL123") == "CHEMBL123"


def test_composite_biomarker_is_not_rendered_as_nan_or_fake_gene() -> None:
  row = pd.Series({
    "analysis_lost_gene": np.nan,
    "lost_gene": np.nan,
    "lost_feature": "MSI/MMR deficiency",
    "analysis_target_gene": "WRN",
    "target_gene": "WRN",
  })
  lost, target, entity_type = resolved_pair_parts(row)
  assert lost == "MSI/MMR deficiency"
  assert target == "WRN"
  assert entity_type == "composite_feature"
  label = resolved_pair_label(row, width=80)
  assert label == "MSI/MMR deficiency ⇒ WRN"
  assert "nan" not in label.casefold()


def test_gene_loss_keeps_gene_arrow() -> None:
  row = pd.Series({
    "analysis_lost_gene": "MTAP",
    "lost_feature": "MTAP loss",
    "analysis_target_gene": "PRMT5",
  })
  assert resolved_pair_label(row, width=80) == "MTAP → PRMT5"


def test_drug_identifier_is_used_when_name_is_missing() -> None:
  row = pd.Series({
    "drug_name": np.nan,
    "drug_id": "CHEMBL12345",
    "drug_key": "CHEMBL12345",
  })
  display, resolution = compound_resolution(row)
  assert display == "CHEMBL12345"
  assert resolution == "compound_identifier"


def test_target_only_is_explicit() -> None:
  row = pd.Series({
    "drug_name": np.nan,
    "drug_id": np.nan,
    "drug_key": "TARGET_ONLY",
  })
  display, resolution = compound_resolution(row)
  assert display == "Target-level evidence only"
  assert resolution == "target_only"


def test_display_columns_preserve_provenance() -> None:
  frame = pd.DataFrame([
    {
      "analysis_lost_gene": np.nan,
      "lost_gene": np.nan,
      "lost_feature": "BRCA1/BRCA2 or HRD",
      "analysis_target_gene": "PARP1",
    }
  ])
  output = add_display_pair_columns(frame)
  assert output.loc[0, "display_lost_label"] == "BRCA1/BRCA2 or HRD"
  assert output.loc[0, "display_target_label"] == "PARP1"
  assert output.loc[0, "lost_entity_type"] == "composite_feature"
  assert output.loc[0, "display_pair_label"] == "BRCA1/BRCA2 or HRD ⇒ PARP1"
