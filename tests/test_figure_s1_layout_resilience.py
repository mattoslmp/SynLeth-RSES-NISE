from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rses_onco.publication import audit_figure_layout
from scripts.make_main_figures_resilient import read_optional
from scripts.publication_figure_s1_resilience import build_coverage_figure


def coverage_frame(
  labels: list[str],
  column: str,
) -> pd.DataFrame:
  count = len(labels)
  return pd.DataFrame({
    column: labels,
    "coverage_fraction": np.linspace(0.08, 0.96, count),
    "coverage_label": [
      f"{index + 1}/{count}"
      for index in range(count)
    ],
  })


def test_figure_s1_long_mechanistic_labels_pass_strict_layout() -> None:
  domain = coverage_frame(
    [
      f"Evidence domain {index} with complete eligibility description"
      for index in range(9)
    ],
    "domain_label",
  )
  source = coverage_frame(
    [
      f"Source {index} with a deliberately descriptive provenance label"
      for index in range(11)
    ],
    "evidence_source",
  )

  cancer_rows = []
  for domain_family in (
    "RSES-Onco",
    "Functional microniche",
    "Expression and regulatory evidence",
    "Validation and tractability",
    "Biochemical and structural evidence",
    "Conditional genetic phenotype evidence",
  ):
    for cancer in ("colon", "stomach", "lung"):
      cancer_rows.append({
        "domain_family": domain_family,
        "cancer": cancer,
        "coverage_fraction": 0.62,
        "coverage_label": "62/100",
      })
  cancer = pd.DataFrame(cancer_rows)

  classes = pd.DataFrame({
    "domain_family": [
      "Functional microniche and regulatory-network evidence"
    ] * 24,
    "mechanistic_class": [
      (
        "non_homologous_isofunctional_enzyme_or_"
        f"homologous_paralog_mechanistic_class_{index:02d}"
      )
      for index in range(24)
    ],
    "coverage_fraction": np.linspace(0.15, 0.98, 24),
    "coverage_label": [f"{index + 1}/24" for index in range(24)],
  })

  reasons = pd.DataFrame({
    "evidence_state": [
      "missing",
      "insufficient_sample",
      "technical_failure",
      "not_eligible",
    ] * 6,
    "absence_reason": [
      (
        "candidate_gene_or_context_is_unavailable_for_the_"
        f"declared_evidence_source_reason_{index:02d}"
      )
      for index in range(24)
    ],
    "rows": np.arange(1, 25),
  })

  fig, source_data = build_coverage_figure(
    domain,
    source,
    cancer,
    classes,
    reasons,
  )
  audit = audit_figure_layout(fig, "Figure_S1")
  try:
    assert audit.warnings == ()
    assert audit.status == "pass"
    assert set(source_data["panel"].dropna()) == {
      "A_domain",
      "B_source",
      "C_cancer",
      "D_class",
      "E_missingness_reason",
    }
    assert fig.get_size_inches()[1] > 17.0
  finally:
    plt.close(fig)


def test_main_figure_reader_handles_mixed_identifier_types_without_dtype_warning(
  tmp_path: Path,
) -> None:
  path = tmp_path / "mixed.tsv"
  pd.DataFrame({
    "activity": [1, "enzyme activity", None, 4],
    "ec_number": ["1.1.1.1", 2, None, "3.6.1.1"],
    "drug_name": [None, "compound A", 17, "compound B"],
    "coverage_adjusted_rses": [0.5, 0.4, 0.3, 0.2],
  }).to_csv(path, sep="\t", index=False)

  with warnings.catch_warnings():
    warnings.simplefilter("error", pd.errors.DtypeWarning)
    frame = read_optional(path)

  assert len(frame) == 4
  assert set(frame.columns) == {
    "activity",
    "ec_number",
    "drug_name",
    "coverage_adjusted_rses",
  }
