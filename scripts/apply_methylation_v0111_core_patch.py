#!/usr/bin/env python3
"""Apply core RSES-Onco v0.11.1 methylation integration patches."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
  return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
  (ROOT / path).write_text(text, encoding="utf-8")


def replace_once(path: str, old: str, new: str) -> None:
  text = read(path)
  if new in text:
    return
  if old not in text:
    raise RuntimeError(f"Patch anchor missing in {path}: {old[:100]!r}")
  write(path, text.replace(old, new, 1))


def patch_dependencies() -> None:
  replace_once("pyproject.toml", 'version = "0.11.0"', 'version = "0.11.1"')
  replace_once(
    "pyproject.toml",
    'description = "Coverage-aware all-NISE cancer vulnerability with WGCNA expression networks, promoter-aware TF regulation, evidence audit, robustness, pharmacology, structural atlas and scripted publication assets"',
    'description = "Coverage-aware all-NISE cancer vulnerability with WGCNA, TCGA/GDC methylation, promoter-aware TF regulation, evidence audit, robustness, pharmacology, structural atlas and scripted publication assets"',
  )
  replace_once(
    "pyproject.toml",
    '  "python-docx>=1.1",\n]',
    '  "python-docx>=1.1",\n  "xenaPython==1.0.14",\n]',
  )
  replace_once(
    "environment.yml",
    "  - pip:\n      - -e .",
    "  - pip:\n      - xenaPython==1.0.14\n      - -e .",
  )
  replace_once("CITATION.cff", "version: 0.11.0", "version: 0.11.1")
  replace_once(
    "CITATION.cff",
    "date-released: 2026-07-22",
    "date-released: 2026-07-23",
  )


def patch_resume_pipeline() -> None:
  path = "scripts/resume_wgcna_regulatory_pipeline.sh"
  replace_once(
    path,
    'PROMOTER_MOTIFS="${PROMOTER_MOTIFS:-$PROCESSED_DIR/regulatory/jaspar_promoter_tf_summary.tsv}"\nLOSS_THRESHOLD=',
    'PROMOTER_MOTIFS="${PROMOTER_MOTIFS:-$PROCESSED_DIR/regulatory/jaspar_promoter_tf_summary.tsv}"\nMETHYLATION_DIR="${METHYLATION_DIR:-$PROCESSED_DIR/epigenetics/methylation}"\nMETHYLATION_METRICS="${METHYLATION_METRICS:-$METHYLATION_DIR/tcga_nise_methylation_pair_metrics.tsv}"\nMETHYLATION_STRICT="${METHYLATION_STRICT:-0}"\nLOSS_THRESHOLD=',
  )
  replace_once(
    path,
    '  command -v fimo >/dev/null 2>&1 || {\n    echo "FIMO from the MEME suite is missing." >&2\n    return 1\n  }\n  bash -n scripts/resume_wgcna_regulatory_pipeline.sh',
    '  command -v fimo >/dev/null 2>&1 || {\n    echo "FIMO from the MEME suite is missing." >&2\n    return 1\n  }\n  python -c "import xenaPython; print(\'xenaPython methylation client: OK\')"\n  bash -n scripts/resume_wgcna_regulatory_pipeline.sh',
  )
  replace_once(
    path,
    '  log_stage "Build cancer-specific signed WGCNA and TF regulatory evidence"',
    '''  log_stage "Acquire TCGA/GDC candidate-gene methylation through UCSC Xena"
  methylation_command=(
    python -u scripts/acquire_tcga_nise_methylation.py
      --candidates "$CANDIDATES"
      --output-dir "$METHYLATION_DIR"
      --min-samples 20
  )
  if [[ "$METHYLATION_STRICT" == "1" ]]; then
    methylation_command+=(--strict)
  fi
  run_logged "$LOG_DIR/03b_tcga_methylation.log" "${methylation_command[@]}"
  require_file "$METHYLATION_METRICS"

  log_stage "Build cancer-specific signed WGCNA and TF regulatory evidence"''',
  )
  replace_once(
    path,
    '      --promoter-motifs "$PROMOTER_MOTIFS" \\\n      --output "$CANCER_SPECIFIC_EVIDENCE"\n\n  log_stage "Build pair-level consensus for source compatibility"',
    '      --promoter-motifs "$PROMOTER_MOTIFS" \\\n      --output "$CANCER_SPECIFIC_EVIDENCE"\n\n  log_stage "Integrate methylation into the regulatory microniche domain"\n  run_logged "$LOG_DIR/04b_methylation_regulatory_integration.log" \\\n    python -u scripts/integrate_methylation_regulatory_layer.py \\\n      --regulatory-evidence "$CANCER_SPECIFIC_EVIDENCE" \\\n      --methylation-metrics "$METHYLATION_METRICS" \\\n      --output "$CANCER_SPECIFIC_EVIDENCE"\n\n  log_stage "Build pair-level consensus for source compatibility"',
  )


def patch_aggregate() -> None:
  path = "scripts/aggregate_wgcna_regulatory_layer.py"
  replace_once(
    path,
    '    "regulatory_promoter_motif_divergence",\n    "regulatory_lost_regulator_count",',
    '    "regulatory_promoter_motif_divergence",\n    "regulatory_promoter_methylation_context",\n    "component_promoter_methylation_context",\n    "methylation_context_raw",\n    "methylation_context_coverage",\n    "methylation_primary_tumor_overlap_n",\n    "methylation_pair_spearman_rho",\n    "methylation_pair_median_absolute_beta_difference",\n    "methylation_tumor_normal_delta_divergence",\n    "regulatory_lost_regulator_count",',
  )
  replace_once(
    path,
    '    "JASPAR_promoter_motif_prediction"\n  )',
    '    "JASPAR_promoter_motif_prediction_plus_TCGA_GDC_methylation"\n  )',
  )
  replace_once(
    path,
    '"wgcna-promoter-regulatory-v2"',
    '"wgcna-promoter-methylation-regulatory-v3"',
  )


def patch_recompute() -> None:
  path = "scripts/recompute_rses_with_wgcna_regulatory.py"
  replace_once(
    path,
    '"eligibility-aware-wgcna-regulatory-v3"',
    '"eligibility-aware-wgcna-regulatory-methylation-v4"',
  )
  replace_once(
    path,
    '      "regulatory_promoter_motif_divergence": pair_evidence.get(\n        "regulatory_promoter_motif_divergence"\n      ),',
    '      "regulatory_promoter_motif_divergence": pair_evidence.get(\n        "regulatory_promoter_motif_divergence"\n      ),\n      "regulatory_promoter_methylation_context": pair_evidence.get(\n        "regulatory_promoter_methylation_context"\n      ),\n      "methylation_context_raw": pair_evidence.get(\n        "methylation_context_raw"\n      ),\n      "methylation_context_subcoverage": pair_evidence.get(\n        "methylation_context_coverage"\n      ),\n      "methylation_primary_tumor_overlap_n": pair_evidence.get(\n        "methylation_primary_tumor_overlap_n"\n      ),\n      "methylation_source_status": pair_evidence.get(\n        "methylation_source_status"\n      ),',
  )
  replace_once(
    path,
    '"score_version": "RSES-Onco-expanded-v0.10.9"',
    '"score_version": "RSES-Onco-expanded-v0.11.1"',
  )
  replace_once(
    path,
    '        "0.40*DoRothEA_regulator_divergence + "\n        "0.35*TF_expression_profile_divergence + "\n        "0.25*JASPAR_promoter_motif_divergence, coverage-adjusted "',
    '        "0.32*DoRothEA_regulator_divergence + "\n        "0.28*TF_expression_profile_divergence + "\n        "0.20*JASPAR_promoter_motif_divergence + "\n        "0.20*TCGA_GDC_methylation_context, eligibility- and coverage-adjusted "',
  )
  replace_once(
    path,
    '      "direct_promoter_binding_claim": False,',
    '      "direct_promoter_binding_claim": False,\n      "direct_methylation_silencing_claim": False,',
  )


def patch_support_export() -> None:
  path = "scripts/export_wgcna_regulatory_supporting_evidence.py"
  replace_once(
    path,
    '    (\n      "ensembl_canonical_promoters",',
    '''    (
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
      "ensembl_canonical_promoters",''',
  )


def patch_validator() -> None:
  path = "scripts/validate_wgcna_regulatory_evidence.py"
  replace_once(
    path,
    '"eligibility-aware-wgcna-regulatory-v3"',
    '"eligibility-aware-wgcna-regulatory-methylation-v4"',
  )
  replace_once(
    path,
    '"regulatory_promoter_motif_divergence",\n      "regulatory_network_subcoverage",',
    '"regulatory_promoter_motif_divergence",\n      "regulatory_promoter_methylation_context",\n      "methylation_context_subcoverage",\n      "regulatory_network_subcoverage",',
  )
  replace_once(
    path,
    'if score_versions != {"RSES-Onco-expanded-v0.10.9"}:',
    'if score_versions != {"RSES-Onco-expanded-v0.11.1"}:',
  )
  replace_once(
    path,
    '"expected RSES-Onco-expanded-v0.10.9"',
    '"expected RSES-Onco-expanded-v0.11.1"',
  )
  replace_once(
    path,
    '      "regulatory_promoter_motif_divergence",\n      "regulatory_network_subcoverage",',
    '      "regulatory_promoter_motif_divergence",\n      "regulatory_promoter_methylation_context",\n      "methylation_context_subcoverage",\n      "regulatory_network_subcoverage",',
  )
  replace_once(
    path,
    '    "jaspar_promoter_motif_predictions",\n  }',
    '    "jaspar_promoter_motif_predictions",\n    "tcga_gdc_methylation_pair_metrics",\n    "tcga_gdc_methylation_source_status",\n  }',
  )


def main() -> None:
  patch_dependencies()
  patch_resume_pipeline()
  patch_aggregate()
  patch_recompute()
  patch_support_export()
  patch_validator()
  print("Applied core methylation integration patch.")


if __name__ == "__main__":
  main()
