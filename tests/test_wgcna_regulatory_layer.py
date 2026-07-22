from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from rses_onco.expanded import (
  EXPANDED_ONCO_WEIGHTS,
  FUNCTIONAL_MICRONICHE_WEIGHTS,
  coverage_aware_score,
)
from scripts.build_wgcna_regulatory_layer import (
  REGULATORY_SUBWEIGHTS,
  WGCNA_SUBWEIGHTS,
  cosine_divergence,
  jaccard_divergence,
  weighted_subscore,
)

ROOT = Path(__file__).resolve().parents[1]


def test_wgcna_subscore_preserves_missingness() -> None:
  score = weighted_subscore(
    {
      "tom_divergence": 0.8,
      "module_divergence": None,
      "kme_divergence": 0.4,
    },
    WGCNA_SUBWEIGHTS,
  )
  assert score["observed_subcomponents"] == 2
  assert 0 < score["coverage"] < 1
  assert np.isclose(
    score["adjusted"],
    score["raw"] * score["coverage"],
  )


def test_regulatory_sublayer_weights_sum_to_one() -> None:
  assert np.isclose(sum(REGULATORY_SUBWEIGHTS.values()), 1.0)
  assert np.isclose(sum(WGCNA_SUBWEIGHTS.values()), 1.0)


def test_set_and_profile_divergence_are_bounded() -> None:
  assert np.isclose(
    jaccard_divergence({"TF1", "TF2"}, {"TF2"}) or 0,
    0.5,
  )
  value = cosine_divergence(
    {"TF1": 1.0},
    {"TF1": -1.0},
  )
  assert value is not None
  assert np.isclose(value, 1.0)


def test_expression_pairwise_and_wgcna_share_one_domain_weight() -> None:
  score = coverage_aware_score(
    {
      "pairwise_expression_context": 0.8,
      "wgcna_expression_network": 0.2,
    },
    {
      "pairwise_expression_context": 0.5,
      "wgcna_expression_network": 0.5,
    },
  )
  assert np.isclose(score.adjusted_score, 0.5)
  assert np.isclose(score.coverage, 1.0)


def test_promoter_motifs_are_not_direct_binding_claims() -> None:
  scanner = (
    ROOT / "scripts/scan_promoter_motifs.py"
  ).read_text(encoding="utf-8")
  assert "direct_binding_claim" in scanner
  assert "= False" in scanner
  assert "not direct TF binding" in scanner
  validator = (
    ROOT / "scripts/validate_wgcna_regulatory_evidence.py"
  ).read_text(encoding="utf-8")
  assert "incorrectly claims direct promoter binding" in validator


def test_resume_pipeline_requires_actual_wgcna_and_fimo() -> None:
  pipeline = (
    ROOT / "scripts/resume_wgcna_regulatory_pipeline.sh"
  ).read_text(encoding="utf-8")
  builder = (
    ROOT / "scripts/build_wgcna_regulatory_layer.py"
  ).read_text(encoding="utf-8")
  implementation = pipeline + "\n" + builder
  required = {
    "run_wgcna_expression_network.R",
    "download_ensembl_promoters.py",
    "download_jaspar_core_vertebrates.py",
    "scan_promoter_motifs.py",
    "build_wgcna_regulatory_layer.py",
    "aggregate_wgcna_regulatory_layer.py",
    "recompute_rses_with_wgcna_regulatory.py",
    "Rscript",
    "fimo",
  }
  missing = sorted(
    value
    for value in required
    if value not in implementation
  )
  assert not missing, (
    f"WGCNA regulatory pipeline missing: {missing}"
  )


def test_consensus_aggregation_prevents_cancer_triplication(
  tmp_path: Path,
) -> None:
  base = pd.DataFrame({
    "pair_id": ["P1"],
    "component_localization": [0.4],
  })
  detailed = pd.DataFrame({
    "pair_id": ["P1", "P1", "P1"],
    "cancer": ["colon", "stomach", "lung"],
    "component_wgcna_expression_network": [0.2, 0.4, 0.6],
    "component_regulatory_network_composite": [0.3, 0.5, 0.7],
    "regulatory_network_coverage": [1.0, 1.0, 1.0],
    "wgcna_expression_network_coverage": [1.0, 1.0, 1.0],
  })
  base_path = tmp_path / "base.tsv"
  detailed_path = tmp_path / "detailed.tsv"
  output_path = tmp_path / "output.tsv"
  base.to_csv(base_path, sep="\t", index=False)
  detailed.to_csv(detailed_path, sep="\t", index=False)

  subprocess.run(
    [
      sys.executable,
      "scripts/aggregate_wgcna_regulatory_layer.py",
      "--base",
      str(base_path),
      "--cancer-specific",
      str(detailed_path),
      "--output",
      str(output_path),
    ],
    cwd=ROOT,
    check=True,
  )
  result = pd.read_csv(output_path, sep="\t")
  assert len(result) == 1
  assert np.isclose(
    result.loc[0, "component_wgcna_expression_network"],
    0.4,
  )
  assert result.loc[0, "consensus_cancers_observed"] == 3


def test_recompute_matches_cancer_specific_network_evidence(
  tmp_path: Path,
) -> None:
  rows = []
  for cancer in ("colon", "lung"):
    row: dict[str, object] = {
      "pair_id": "P1",
      "cancer": cancer,
      "score_comparability_group": "gene_pair",
      "microniche_expression_context": 0.8,
    }
    for domain in FUNCTIONAL_MICRONICHE_WEIGHTS:
      row.setdefault(f"microniche_{domain}", 0.5)
      row[f"eligible_microniche_{domain}"] = True
    for domain in EXPANDED_ONCO_WEIGHTS:
      row[f"component_{domain}"] = 0.5
      row[f"eligible_component_{domain}"] = True
    rows.append(row)
  ranking = pd.DataFrame(rows)
  evidence = pd.DataFrame([
    {
      "pair_id": "P1",
      "cancer": "colon",
      "component_wgcna_expression_network": 0.2,
      "component_regulatory_network_composite": 0.3,
      "regulatory_tf_association_divergence": 0.3,
      "regulatory_tf_expression_profile_divergence": 0.3,
      "regulatory_promoter_motif_divergence": 0.3,
      "regulatory_network_raw": 0.3,
      "regulatory_network_coverage": 1.0,
      "promoter_evidence_type": (
        "JASPAR_motif_prediction_not_direct_binding"
      ),
    },
    {
      "pair_id": "P1",
      "cancer": "lung",
      "component_wgcna_expression_network": 0.8,
      "component_regulatory_network_composite": 0.9,
      "regulatory_tf_association_divergence": 0.9,
      "regulatory_tf_expression_profile_divergence": 0.9,
      "regulatory_promoter_motif_divergence": 0.9,
      "regulatory_network_raw": 0.9,
      "regulatory_network_coverage": 1.0,
      "promoter_evidence_type": (
        "JASPAR_motif_prediction_not_direct_binding"
      ),
    },
  ])
  ranking_path = tmp_path / "ranking.tsv"
  evidence_path = tmp_path / "evidence.tsv"
  output_path = tmp_path / "output.tsv"
  ranking.to_csv(ranking_path, sep="\t", index=False)
  evidence.to_csv(evidence_path, sep="\t", index=False)
  subprocess.run(
    [
      sys.executable,
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
  result = pd.read_csv(output_path, sep="\t").set_index("cancer")
  assert np.isclose(result.loc["colon", "wgcna_expression_network"], 0.2)
  assert np.isclose(result.loc["lung", "wgcna_expression_network"], 0.8)
  assert (
    result.loc["lung", "coverage_adjusted_rses"]
    > result.loc["colon", "coverage_adjusted_rses"]
  )
  assert set(result["scoring_semantics_version"]) == {
    "eligibility-aware-v1"
  }
  assert set(result["expression_regulatory_semantics_version"]) == {
    "eligibility-aware-wgcna-regulatory-v2"
  }
