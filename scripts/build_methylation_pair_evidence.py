#!/usr/bin/env python3
"""Build cancer-specific promoter-methylation evidence for directed candidate pairs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
  sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from rses_onco.methylation import pair_methylation_metrics
from rses_onco.utils import bh_adjust, canonical_gene_name

CANCERS = ("colon", "stomach", "lung")
SIMPLE_GENE = re.compile(r"^[A-Za-z0-9-]+$")


def resolve_path(value: str) -> Path:
  path = Path(value)
  return path if path.is_absolute() else ROOT / path


def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--candidates",
    default="data/processed/expanded_candidate_universe.tsv",
  )
  parser.add_argument(
    "--gene-sample",
    default="data/processed/methylation/gdc_promoter_methylation_gene_sample.tsv",
  )
  parser.add_argument(
    "--output",
    default="data/processed/methylation/pair_promoter_methylation_evidence.tsv",
  )
  parser.add_argument("--min-samples", type=int, default=10)
  args = parser.parse_args()

  candidates = pd.read_csv(resolve_path(args.candidates), sep="\t", low_memory=False)
  gene_sample = pd.read_csv(resolve_path(args.gene_sample), sep="\t", low_memory=False)
  rows: list[dict[str, object]] = []
  generated_at = datetime.now(timezone.utc).isoformat()

  for candidate in candidates.to_dict("records"):
    pair_id = str(candidate.get("pair_id"))
    lost = canonical_gene_name(candidate.get("lost_gene"))
    target = canonical_gene_name(candidate.get("target_gene"))
    simple_pair = bool(
      lost
      and target
      and SIMPLE_GENE.fullmatch(lost)
      and SIMPLE_GENE.fullmatch(target)
    )
    for cancer in CANCERS:
      flag = candidate.get(cancer, 1)
      try:
        included = int(float(flag)) == 1
      except (TypeError, ValueError):
        included = True
      if not included:
        continue
      if not simple_pair:
        rows.append({
          "pair_id": pair_id,
          "cancer": cancer,
          "lost_gene": lost,
          "target_gene": target,
          "n_samples": 0,
          "evidence_status": "not_eligible",
          "absence_reason": "composite_event_not_executable_as_gene_pair_methylation_analysis",
          "methylation_source": "GDC_SeSAMe_Methylation_Beta_Value",
          "methylation_score_role": "expression_compensation_sublayer",
          "generated_at_utc": generated_at,
          "direct_gene_silencing_claim": False,
        })
        continue
      metrics = pair_methylation_metrics(
        gene_sample,
        lost,
        target,
        cancer,
        min_samples=args.min_samples,
      )
      rows.append({
        "pair_id": pair_id,
        **metrics.__dict__,
        "methylation_source": "GDC_SeSAMe_Methylation_Beta_Value",
        "methylation_workflow": "GDC_SeSAMe_Methylation_Array_Harmonization",
        "methylation_score_role": "expression_compensation_sublayer",
        "methylation_formula": "lost_gene_median_promoter_beta * (1 - target_gene_median_promoter_beta)",
        "promoter_definition": "GDC GENCODE-v36 promoter/TSS probes; -2000 to +500 bp or annotated promoter/TSS group",
        "generated_at_utc": generated_at,
        "direct_gene_silencing_claim": False,
      })

  result = pd.DataFrame(rows)
  if result.empty:
    raise RuntimeError("No candidate/cancer rows were generated for methylation evidence")
  result["paired_wilcoxon_q_value_bh"] = np.nan
  result["paired_wilcoxon_q_value_bh_within_cancer"] = np.nan
  observed = result["paired_wilcoxon_p_value"].notna() if "paired_wilcoxon_p_value" in result else pd.Series(False, index=result.index)
  if observed.any():
    result.loc[observed, "paired_wilcoxon_q_value_bh"] = bh_adjust(
      result.loc[observed, "paired_wilcoxon_p_value"]
    )
    result.loc[observed, "paired_wilcoxon_q_value_bh_within_cancer"] = (
      result.loc[observed]
      .groupby("cancer")["paired_wilcoxon_p_value"]
      .transform(lambda values: bh_adjust(values))
    )
  result = result.sort_values(
    ["cancer", "evidence_status", "promoter_methylation_context_score", "pair_id"],
    ascending=[True, True, False, True],
    na_position="last",
  )
  output = resolve_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  temporary = output.with_suffix(output.suffix + ".tmp")
  result.to_csv(temporary, sep="\t", index=False)
  temporary.replace(output)
  print(f"Wrote cancer-specific promoter methylation evidence: {output} ({len(result):,} rows)")
  print(result["evidence_status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
  main()
