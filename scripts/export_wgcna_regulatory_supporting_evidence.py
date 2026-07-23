#!/usr/bin/env python3
"""Export WGCNA, promoter and TF regulatory source tables into article outputs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

ROOT = Path(__file__).resolve().parents[1]


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def sha256(path: Path) -> str:
  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
      digest.update(chunk)
  return digest.hexdigest()


def read_optional(path: Path) -> pd.DataFrame:
  if not path.exists() or path.stat().st_size == 0:
    return pd.DataFrame()
  try:
    return pd.read_csv(path, sep="\t", low_memory=False)
  except EmptyDataError:
    return pd.DataFrame()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temporary = path.with_suffix(path.suffix + ".tmp")
  frame.to_csv(temporary, sep="\t", index=False)
  temporary.replace(path)


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--article-root", default="article_outputs")
  args = parser.parse_args()
  article_root = resolve_path(args.article_root)
  output_dir = article_root / "tables/supporting_evidence/expression_regulatory"
  specifications = [
    (
      "cancer_specific_wgcna_pair_metrics",
      resolve_path(
        "data/processed/regulatory/wgcna/"
        "wgcna_pair_metrics_all_cancers.tsv"
      ),
      output_dir / "wgcna_pair_metrics_all_cancers.tsv",
      "Signed WGCNA TOM, modules and kME are expression-derived "
      "subcomponents and do not receive an independent full-domain weight.",
    ),
    (
      "wgcna_input_preparation",
      resolve_path(
        "data/processed/regulatory/wgcna/wgcna_input_preparation.tsv"
      ),
      output_dir / "wgcna_input_preparation.tsv",
      "All candidate genes plus cancer-specific highly variable genes "
      "were used; median imputation is reported only for network estimation.",
    ),
    (
      "wgcna_correlation_fallback_audit",
      resolve_path(
        "data/processed/regulatory/wgcna/"
        "wgcna_correlation_fallback_all_cancers.tsv"
      ),
      output_dir / "wgcna_correlation_fallback_all_cancers.tsv",
      "Biweight midcorrelation is primary; Pearson is used only for "
      "individual genes or module eigengenes with zero or non-finite MAD, "
      "and every affected entity is recorded.",
    ),
    (
      "wgcna_run_diagnostics",
      resolve_path(
        "data/processed/regulatory/wgcna/"
        "wgcna_run_diagnostics_all_cancers.tsv"
      ),
      output_dir / "wgcna_run_diagnostics_all_cancers.tsv",
      "Correlation settings, zero-MAD counts, fallback policy, selected "
      "power and module counts are reported per cancer.",
    ),
    (
      "promoter_tf_regulatory_pair_metrics",
      resolve_path(
        "data/processed/regulatory/promoter_tf_regulatory_pair_metrics.tsv"
      ),
      output_dir / "promoter_tf_regulatory_pair_metrics.tsv",
      "DoRothEA associations, TF-expression consistency and promoter "
      "motif predictions are distinct subcomponents.",
    ),
    (
      "tcga_gdc_methylation_pair_metrics",
      resolve_path(
        "data/processed/epigenetics/methylation/"
        "tcga_nise_methylation_pair_metrics.tsv"
      ),
      output_dir / "tcga_nise_methylation_pair_metrics.tsv",
      "GDC beta values accessed through UCSC Xena provide gene-associated "
      "methylation context; they are not direct proof of silencing or causality.",
    ),
    (
      "tcga_gdc_methylation_source_status",
      resolve_path(
        "data/processed/epigenetics/methylation/"
        "tcga_nise_methylation_source_status.tsv"
      ),
      output_dir / "tcga_nise_methylation_source_status.tsv",
      "Source, dataset, access date and technical availability are reported. "
      "Repbase is not used because it is a repeat-sequence library.",
    ),
    (
      "ensembl_canonical_promoters",
      resolve_path("data/raw/regulatory/ensembl_promoters.tsv"),
      output_dir / "ensembl_canonical_promoters.tsv",
      "Promoter coordinates are canonical-transcript TSS windows and are "
      "not proof of TF occupancy.",
    ),
    (
      "jaspar_promoter_motif_predictions",
      resolve_path(
        "data/processed/regulatory/jaspar_promoter_tf_summary.tsv"
      ),
      output_dir / "jaspar_promoter_tf_summary.tsv",
      "JASPAR/FIMO motif occurrence is predicted cis-regulatory support, "
      "not direct binding or causality.",
    ),
  ]
  rows = []
  for family, source, output, boundary in specifications:
    frame = read_optional(source)
    if frame.empty:
      raise FileNotFoundError(
        f"Missing or empty WGCNA/regulatory source: {source}"
      )
    atomic_tsv(frame, output)
    rows.append({
      "evidence_family": family,
      "source_path": str(source),
      "output_path": str(output),
      "rows": len(frame),
      "columns": len(frame.columns),
      "sha256": sha256(output),
      "status": "available",
      "interpretation_boundary": boundary,
    })
  manifest = pd.DataFrame(rows)
  manifest["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
  manifest_path = (
    output_dir / "wgcna_regulatory_supporting_evidence_manifest.tsv"
  )
  atomic_tsv(manifest, manifest_path)
  print(
    manifest[["evidence_family", "rows", "output_path"]]
      .to_string(index=False)
  )
  print(f"Wrote {manifest_path}")


if __name__ == "__main__":
  main()
