# Publication layout recovery (v0.10.5)

## Problem

A complete real-data run reached figure generation and strict layout validation.
Two issues were observed:

1. axis-off schematic panels retained latent Matplotlib tick-label artists and
   triggered false `tick_outside_figure` warnings;
2. Figure S2 contained long wrapped activity labels whose rendered line count
   exceeded the height estimated only from the number of rows, producing a real
   `tick_overlap:0:y` failure.

## Resolution

The publication workflow now uses resilient figure entry points:

- `scripts/make_main_figures_resilient.py`
- `scripts/make_supplementary_figures_resilient.py`
- `scripts/publication_layout_resilience.py`

Axis-off panels have only their latent ticks cleared before auditing. Visible
axes remain under the original strict clipping and overlap checks.

Horizontal supplementary bar charts now calculate height from the total number
of wrapped text lines, use explicit numeric y positions and preserve readable
9-point tick labels with additional padding.

## Resume from cached analyses

The pharmacology, structural and article tables do not need to be recomputed.
Run:

```bash
STRICT_LAYOUT=1 bash scripts/run_publication_pipeline.sh assets-only
```

This regenerates priorities, tables, all figures, workbook, manifests and final
validation using cached API and structural products.

## Expected result

- 8 main figures;
- 32 supplementary figures;
- 120 PNG/PDF/SVG image files;
- 4 main tables;
- 18 supplementary tables;
- all layout audits pass.
