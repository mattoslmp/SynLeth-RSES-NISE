import numpy as np

from rses_onco.utils import bh_adjust, canonical_gene_name


def test_canonical_gene_name():
  assert canonical_gene_name("PRMT5 (10419)") == "PRMT5"


def test_bh_adjust_monotonic_and_bounded():
  adjusted = bh_adjust([0.01, 0.04, 0.03, np.nan])
  assert np.nanmax(adjusted) <= 1
  assert np.isfinite(adjusted[:3]).all()
