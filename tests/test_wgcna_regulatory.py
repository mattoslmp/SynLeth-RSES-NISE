from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from scripts.build_wgcna_regulatory_layer import (
  REGULATORY_SUBWEIGHTS,
  WGCNA_SUBWEIGHTS,
  jaccard_divergence,
  weighted_subscore,
)
from scripts.download_ensembl_promoters import promoter_coordinates
from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  FUNCTIONAL_MICRONICHE_WEIGHTS,
)

ROOT = Path(__file__).resolve().parents[1]


def test_weighted_subscore_preserves_missingness_and_coverage() -> None:
  score = weighted_subscore(
    {
      "tom_divergence": 0.8,
      "module_divergence": None,
      "kme_divergence": 0.2,
    },
    WGCNA_SUBWEIGHTS,
  )
  expected_raw = (
    WGCNA_SUBWEIGHTS["tom_divergence"] * 0.8
    + WGCNA_SUBWEIGHTS["kme_divergence"] * 0.2
  ) / (
    WGCNA_SUBWEIGHTS["tom_divergence"]
    + WGCNA_SUBWEIGHTS["kme_divergence"]
  )
  expected_coverage = (
    WGCNA_SUBWEIGHTS["tom_divergence"]
    + WGCNA_SUBWEIGHTS["kme_divergence"]
  ) / sum(WGCNA_SUBWEIGHTS.values())
  assert np.isclose(score["raw"], expected_raw)
  assert np.isclose(score["coverage"], expected_coverage)
  assert np.isclose(score["adjusted"], expected_raw * expected_coverage)
  assert score["observed_subcomponents"] == 2


def test_regulator_jaccard_divergence() -> None:
  assert np.isclose(
    jaccard_divergence({"TF1", "TF2"}, {"TF2", "TF3"}) or 0,
    2 / 3,
  )
  assert jaccard_divergence(set(), set()) is None


def test_promoter_coordinates_respect_transcript_strand() -> None:
  forward = {
    "id": "ENSG1",
    "assembly_name": "GRCh38",
    "canonical_transcript": "ENST1.1",
    "Transcript": [{
      "id": "ENST1",
      "strand": 1,
      "start": 10000,
      "end": 12000,
      "seq_region_name": "1",
      "biotype": "protein_coding",
      "is_canonical": 1,
    }],
  }
  reverse = {
    "id": "ENSG2",
    "assembly_name": "GRCh38",
    "canonical_transcript": "ENST2.1",
    "Transcript": [{
      "id": "ENST2",
      "strand": -1,
      "start": 20000,
      "end": 24000,
      "seq_region_name": "2",
      "biotype": "protein_coding",
      "is_canonical": 1,
    }],
  }
  first = promoter_coordinates(forward, 2000, 500)
  second = promoter_coordinates(reverse, 2000, 500)
  assert first["tss"] == 10000
  assert first["promoter_start"] == 8000
  assert first["promoter_end"] == 10500
  assert first["region"] == "1:8000..10500:1"
  assert second["tss"] == 24000
  assert second["promoter_start"] == 23500
  assert second["promoter_end"] == 26000
  assert second["region"] == "2:23500..26000:-1"


def test_recompute_shares_existing_domain_weights(tmp_path: Path) -> None:
  ranking_row: dict[str, object] = {
    "pair_id": "PAIR1",
    "cancer": "lung",
    "score_comparability_group": "gene_pair",
    "microniche_expression_context": 0.20,
  }
  for domain in FUNCTIONAL_MICRONICHE_WEIGHTS:
    ranking_row.setdefault(f"microniche_{domain}", 0.40)
    ranking_row[f"eligible_microniche_{domain}"] = True
  for domain in EXPANDED_ONCO_WEIGHTS:
    ranking_row[f"component_{domain}"] = 0.50
    ranking_row[f"eligible_component_{domain}"] = True
  ranking = pd.DataFrame([ranking_row])
  evidence = pd.DataFrame([{
    "pair_id": "PAIR1",
    "cancer": "lung",
    "component_wgcna_expression_network": 0.80,
    "component_regulatory_network_composite": 0.60,
    "regulatory_tf_association_divergence": 0.40,
    "regulatory_tf_expression_profile_divergence": 0.70,
    "regulatory_promoter_motif_divergence": 0.50,
    "regulatory_network_raw": 0.55,
    "regulatory_network_coverage": 1.0,
    "promoter_evidence_type": "JASPAR_motif_prediction_not_direct_binding",
  }])
  ranking_path = tmp_path / "ranking.tsv"
  evidence_path = tmp_path / "evidence.tsv"
  output_path = tmp_path / "result.tsv"
  ranking.to_csv(ranking_path, sep="\t", index=False)
  evidence.to_csv(evidence_path, sep="\t", index=False)
  subprocess.run(
    [
      sys.executable,
      "-u",
      "scripts/recompute_rses_with_wgcna_regulatory.py",
      "--ranking",
      str(ranking_path),
      "--functional-evidence",
      str(evidence_path),
      "--output",
      str(output_path),
    ],
    cwd=ROOT,
    check=True,
  )
  result = pd.read_csv(output_path, sep="\t").iloc[0]
  assert np.isclose(result["pairwise_expression_context"], 0.20)
  assert np.isclose(result["wgcna_expression_network"], 0.80)
  assert np.isclose(result["expression_context_adjusted"], 0.50)
  assert np.isclose(result["microniche_expression_context"], 0.50)
  assert np.isclose(result["microniche_regulatory_network"], 0.60)
  assert result["scoring_semantics_version"] == "eligibility-aware-v1"
  assert (
    result["expression_regulatory_semantics_version"]
    == "eligibility-aware-wgcna-regulatory-v2"
  )
  assert not bool(result["direct_promoter_binding_claim"])
