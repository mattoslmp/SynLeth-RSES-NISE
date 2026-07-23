from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load(name: str, path: Path):
  spec = importlib.util.spec_from_file_location(name, path)
  module = importlib.util.module_from_spec(spec)
  assert spec.loader is not None
  spec.loader.exec_module(module)
  return module


INTEGRATE = load(
  "integrate_methylation_regulatory_layer",
  ROOT / "scripts/integrate_methylation_regulatory_layer.py",
)

from rses_onco import methylation as METHYLATION


def test_tcga_barcode_helpers():
  assert METHYLATION.tcga_sample_type(
    "TCGA-AA-1234-01A-01D"
  ) == "01"
  assert METHYLATION.tcga_sample_type(
    "TCGA-AA-1234-11A"
  ) == "11"
  assert METHYLATION.tcga_case_id(
    "TCGA-AA-1234-01A-01D"
  ) == "TCGA-AA-1234"


def test_position_rows_handles_columnar_xena_response():
  rows = METHYLATION.position_rows({
    "name": ["cg1", "cg2"],
    "chrom": ["chr1", "chr2"],
    "start": [100, 200],
  })
  assert rows[0]["name"] == "cg1"
  assert rows[1]["start"] == 200


def test_gene_probe_values_are_bounded():
  frame = METHYLATION.query_gene_frame(
    {"name": ["cg1"]},
    [[-0.2, 1.4]],
    ["TCGA-AA-1234-01A", "TCGA-AA-5678-11A"],
    "GENE1",
  )
  assert frame["beta_value"].between(0, 1).all()


def test_methylation_subscore_preserves_missingness():
  score = METHYLATION.coverage_subscore({
    "tumor_profile_divergence": 0.8,
    "tumor_normal_delta_divergence": None,
  })
  assert np.isclose(score["raw"], 0.8)
  assert np.isclose(score["coverage"], 0.70)
  assert np.isclose(score["adjusted"], 0.56)


def test_regulatory_methylation_weights_sum_to_one():
  assert np.isclose(
    sum(
      INTEGRATE.REGULATORY_METHYLATION_SUBWEIGHTS.values()
    ),
    1.0,
  )


def test_source_failure_preserves_original_three_component_score():
  score = INTEGRATE.eligibility_aware_subscore(
    {
      "tf_association_divergence": 0.4,
      "tf_expression_profile_divergence": 0.6,
      "promoter_motif_divergence": 0.2,
      "promoter_methylation_context": None,
    },
    INTEGRATE.REGULATORY_METHYLATION_SUBWEIGHTS,
    {
      "tf_association_divergence",
      "tf_expression_profile_divergence",
      "promoter_motif_divergence",
    },
  )
  old_score = 0.40 * 0.4 + 0.35 * 0.6 + 0.25 * 0.2
  assert np.isclose(score["raw"], old_score)
  assert np.isclose(score["coverage"], 1.0)
  assert np.isclose(score["adjusted"], old_score)


def test_available_but_missing_methylation_lowers_internal_coverage():
  score = INTEGRATE.eligibility_aware_subscore(
    {
      "tf_association_divergence": 0.4,
      "tf_expression_profile_divergence": 0.6,
      "promoter_motif_divergence": 0.2,
      "promoter_methylation_context": None,
    },
    INTEGRATE.REGULATORY_METHYLATION_SUBWEIGHTS,
    set(INTEGRATE.REGULATORY_METHYLATION_SUBWEIGHTS),
  )
  assert np.isclose(score["coverage"], 0.80)
  assert score["adjusted"] < score["raw"]


def test_repbase_is_explicitly_excluded():
  source = (
    ROOT / "scripts/acquire_tcga_nise_methylation.py"
  ).read_text(encoding="utf-8")
  assert '"repbase_used": False' in source
  assert "repetitive DNA sequence" in source
