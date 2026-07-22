#!/usr/bin/env python3
"""Write supplementary methods for WGCNA and promoter-aware regulation."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  root = resolve_path(args.article_root)
  output_dir = root / "manuscript_assets/supplementary_methods"
  output_dir.mkdir(parents=True, exist_ok=True)

  text = """# WGCNA expression-network and promoter-aware regulatory layer

## Design principle

The layer does not add independent full-domain weights to evidence derived from the same expression matrix. Pairwise expression divergence and signed WGCNA share the existing expression-context domain. DoRothEA TF-target associations, TF-target expression consistency and JASPAR promoter motif predictions share the existing regulatory-network domain.

## WGCNA

For each cancer context, DepMap log2(TPM+1) expression is restricted to compatible models. All mapped candidate genes are retained and supplemented with cancer-specific highly variable genes. Genes require at least 80% observed values and non-zero variance. Remaining missing values are median-imputed only for network estimation, with the count recorded.

A signed WGCNA network is built with biweight midcorrelation. The soft-thresholding power is the first power with signed scale-free topology R2 at least 0.80 and negative slope; otherwise the finite power with the highest available R2 is selected. Signed adjacency, signed topological overlap, average-linkage clustering, dynamic tree cutting, module merging at eigengene dissimilarity 0.25, module membership and intramodular connectivity are calculated.

For each directed pair, the WGCNA subcomponent combines TOM divergence (0.40), module divergence (0.30) and absolute kME divergence (0.30). Cancer-specific values are retained in source tables and summarized by the median across cancers for the functional prior, preventing three cancer networks from being counted as independent evidence units.

## Expression-context integration

The expression-context domain uses equal internal weights for the pre-existing cancer-specific pairwise expression divergence and the consensus WGCNA subcomponent. Missing subcomponents remain missing and reduce internal coverage. The total expression-context weight in the functional microniche remains unchanged.

## Regulatory-network integration

The regulatory subcomponent combines DoRothEA regulator-set divergence (0.40), cancer-specific TF-target expression-profile divergence (0.35) and JASPAR promoter motif divergence (0.25). Ensembl canonical-transcript TSS windows are defined as 2,000 bp upstream and 500 bp downstream on the transcript strand. JASPAR 2026 CORE vertebrate non-redundant motifs are scanned with FIMO. Motif occurrence is a predicted cis-regulatory feature and is not direct TF binding, promoter occupancy, causal regulation or experimental validation.

No direct promoter-binding claim is generated unless a separate traceable ChIP-based source is added. Missing source evidence is not converted to zero. The total regulatory-network weight in the functional microniche remains unchanged.

## Interpretation boundary

WGCNA modules, coexpression, TF-target associations and promoter motifs are network-context evidence. They do not by themselves establish synthetic lethality, transcriptional compensation, direct TF binding, drug efficacy or clinical relevance.
"""
  (output_dir / "WGCNA_promoter_regulatory_layer.md").write_text(
    text,
    encoding="utf-8",
  )
  weights = pd.DataFrame([
    {"parent_domain": "expression_context", "subcomponent": "pairwise_expression_context", "internal_weight": 0.50},
    {"parent_domain": "expression_context", "subcomponent": "wgcna_expression_network", "internal_weight": 0.50},
    {"parent_domain": "regulatory_network", "subcomponent": "dorothea_regulator_divergence", "internal_weight": 0.40},
    {"parent_domain": "regulatory_network", "subcomponent": "tf_expression_profile_divergence", "internal_weight": 0.35},
    {"parent_domain": "regulatory_network", "subcomponent": "jaspar_promoter_motif_divergence", "internal_weight": 0.25},
  ])
  weights.to_csv(
    output_dir / "WGCNA_regulatory_internal_weights.tsv",
    sep="\t",
    index=False,
  )
  print(f"Wrote WGCNA/regulatory methods to {output_dir}")


if __name__ == "__main__":
  main()
