from __future__ import annotations

from pathlib import Path

from rses_onco.expanded import EXPANDED_ONCO_WEIGHTS

ROOT = Path(__file__).resolve().parents[1]


def test_expanded_scoring_records_eligibility_and_diagnostics() -> None:
  script = (
    ROOT / "scripts/run_expanded_rses_onco.py"
  ).read_text(encoding="utf-8")
  required = {
    'SCORING_SEMANTICS_VERSION = "eligibility-aware-v1"',
    "eligible_component_",
    "eligible_microniche_",
    "score_comparability_group",
    "dependency_absence_reason",
    "expression_compensation_absence_reason",
    "eligible_domain_weight",
  }
  missing = sorted(value for value in required if value not in script)
  assert not missing, f"expanded scoring script missing semantics: {missing}"


def test_audit_rejects_pre_eligibility_ranking() -> None:
  script = (
    ROOT / "scripts/build_publication_evidence_audit_complete.py"
  ).read_text(encoding="utf-8")
  assert "eligibility-aware-v1" in script
  for domain in EXPANDED_ONCO_WEIGHTS:
    assert f"eligible_component_{{domain}}" in script or (
      "eligible_component_{domain}" in script
    )
  assert "Rerun the expanded scoring stage" in script


def test_missing_and_noneligible_rules_are_documented() -> None:
  methods = (
    ROOT / "scripts/build_publication_methods_documentation.py"
  ).read_text(encoding="utf-8")
  assert "A missing component is not converted to zero" in methods
  assert "A non-eligible component does not enter" in methods
  assert "Technical source failure is not biological negative evidence" in methods
