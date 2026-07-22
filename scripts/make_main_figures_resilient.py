#!/usr/bin/env python3
"""Run main-figure generation with axis-off audit resilience."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

import scripts.make_main_figures as target
from scripts.publication_layout_resilience import save_figure_triplet


target.save_figure_triplet = save_figure_triplet


if __name__ == "__main__":
  target.main()
