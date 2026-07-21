from rses_onco.score import rses_onco_score


def test_missing_domains_reduce_coverage_not_observed_score():
  complete = rses_onco_score({
    "tumor_event": 1.0,
    "dependency": 1.0,
    "selectivity": 1.0,
    "expression_compensation": 1.0,
    "functional_relation": 1.0,
    "validation_tractability": 1.0,
  })
  partial = rses_onco_score({
    "tumor_event": 1.0,
    "dependency": 1.0,
    "selectivity": None,
    "expression_compensation": None,
    "functional_relation": 1.0,
    "validation_tractability": 1.0,
  })
  assert complete.observed_score == 1.0
  assert partial.observed_score == 1.0
  assert partial.coverage < complete.coverage
  assert partial.adjusted_score < complete.adjusted_score
