from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_expression_completion_uses_na_sentinels_not_zero(
  tmp_path: Path,
) -> None:
  coordinates = tmp_path / "coordinates.tsv"
  values = tmp_path / "values.tsv"
  summary = tmp_path / "summary.tsv"
  pd.DataFrame({"gene": ["GENEA", "GENEB"]}).to_csv(
    coordinates,
    sep="\t",
    index=False,
  )
  pd.DataFrame([
    {
      "ModelID": "ACH-1",
      "gene": "GENEA",
      "expression_log2_tpm_plus_1": 2.5,
      "cancer": "colon",
      "source_file": "expression.csv",
    },
  ]).to_csv(values, sep="\t", index=False)
  pd.DataFrame([
    {
      "cancer": "colon",
      "gene": "GENEA",
      "n_models": 1,
      "observed_values": 1,
      "median_expression": 2.5,
      "mean_expression": 2.5,
      "q25_expression": 2.5,
      "q75_expression": 2.5,
      "unit": "log2(TPM+1)",
      "source_file": "expression.csv",
    },
  ]).to_csv(summary, sep="\t", index=False)

  subprocess.run([
    sys.executable,
    "scripts/complete_genomic_circos_expression_summary.py",
    "--coordinates",
    str(coordinates),
    "--model-values",
    str(values),
    "--summary",
    str(summary),
  ], cwd=ROOT, check=True)

  completed_summary = pd.read_csv(
    summary,
    sep="\t",
    low_memory=False,
  )
  assert len(completed_summary) == 6
  missing_summary = completed_summary.loc[
    completed_summary["observed_values"].eq(0)
  ]
  assert not missing_summary.empty
  assert missing_summary["median_expression"].isna().all()
  assert missing_summary["evidence_status"].eq(
    "gene_or_context_unavailable_in_expression_matrix"
  ).all()

  completed_values = pd.read_csv(values, sep="\t", low_memory=False)
  sentinel = completed_values.loc[
    completed_values["is_measurement"].astype(str).str.casefold().eq("false")
  ]
  assert not sentinel.empty
  assert sentinel["expression_log2_tpm_plus_1"].isna().all()
  assert set(completed_values["gene"]) == {"GENEA", "GENEB"}
