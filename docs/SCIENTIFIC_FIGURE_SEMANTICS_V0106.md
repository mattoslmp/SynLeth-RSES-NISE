# Scientific figure semantics recovery (v0.10.6)

## Corrections

### Figure 3

Composite biomarkers are not genes. When `analysis_lost_gene` and `lost_gene` are
missing but `lost_feature` is present, the figure now displays the real feature
using `⇒` instead of printing `nan` or inventing a gene symbol.

Examples:

- `MSI/MMR deficiency ⇒ WRN`
- `BRCA1/BRCA2 or HRD ⇒ PARP1`
- `High alkylation-damage state ⇒ NTHL1`

Gene-to-gene hypotheses retain `→`, such as `MTAP → MAT2A`.

### Figure 6

The all-target discovery panel now reports the exact FDR-supported rows, including
`delta_effect`, nominal P value, within-loss/cancer FDR q value, `n_loss` and
`n_intact`. A single discovery remains a valid single discovery; the figure does
not manufacture additional targets. The class panel displays both median and
maximum scores plus the number of unique directions.

### Figure 7

Missing `drug_name` no longer masks a valid `drug_id`. Compound display precedence
is:

1. `drug_name` plus identifier when available;
2. `drug_id`;
3. non-`TARGET_ONLY` `drug_key`;
4. explicit `Target-level evidence only`.

The main table includes only compound-resolved hypotheses, with cancer, pair or
composite context, compound identifier, source databases and score. Target-only
evidence remains in source and supplementary tables but is not mislabeled as a
compound hypothesis. The full pharmacology universe is represented as a density
plot to avoid unreadable overplotting.

## Recovery

After updating to v0.10.6, cached pharmacology and structural data can be reused:

```bash
STRICT_LAYOUT=1 bash scripts/run_publication_pipeline.sh assets-only
```

Do not rerun ChEMBL, AlphaFold or the complete functional-evidence acquisition for
these display corrections.
