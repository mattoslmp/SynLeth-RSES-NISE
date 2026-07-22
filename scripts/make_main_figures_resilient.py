#!/usr/bin/env python3
"""Run main-figure generation with axis-off audit resilience."""
from __future__ import annotations

import scripts.make_main_figures as target
from scripts.publication_layout_resilience import save_figure_triplet


target.save_figure_triplet = save_figure_triplet


if __name__ == "__main__":
  target.main()
