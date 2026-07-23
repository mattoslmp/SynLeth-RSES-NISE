#!/usr/bin/env python3
"""Idempotently synchronize canonical documentation with PR #25 methylation semantics."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BEGIN = "<!-- BEGIN PROMOTER METHYLATION V0.11.1 -->"
END = "<!-- END PROMOTER METHYLATION V0.11.1 -->"


def read(path: str) -> str:
  return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
  (ROOT / path).write_text(text.rstrip() + "\n", encoding="utf-8")


def replace_block(text: str, heading: str, body: str) -> str:
  block = f"{BEGIN}\n\n## {heading}\n\n{body.strip()}\n\n{END}"
  if BEGIN in text and END in text:
    prefix, remainder = text.split(BEGIN, 1)
    _, suffix = remainder.split(END, 1)
    return prefix.rstrip() + "\n\n" + block + suffix
  return text.rstrip() + "\n\n" + block + "\n"


def patch_protocol() -> None:
  path = "docs/END_TO_END_ARTICLE_PROTOCOL.md"
  text = read(path)
  replacements = {
    "# RSES-Onco v0.11.0:": "# RSES-Onco v0.11.1:",
    "Repository and publication framework: 0.11.0": (
      "Repository and publication framework: 0.11.1"
    ),
    "Scientific score: RSES-Onco-expanded-v0.10.9": (
      "Scientific score: RSES-Onco-expanded-v0.10.10"
    ),
    "WGCNA/regulatory semantics: eligibility-aware-wgcna-regulatory-v3": (
      "Expression/regulatory semantics: eligibility-aware-wgcna-regulatory-v3"
    ),
    'SESSION="rses_v0110_': 'SESSION="rses_v0111_',
    "run_rses_v0110_": "run_rses_v0111_",
    "last_rses_v0110_": "last_rses_v0111_",
    "pre_v0110_": "pre_v0111_",
  }
  for old, new in replacements.items():
    text = text.replace(old, new)
  methylation_link = (
    "- [`METHYLATION_DATA_AND_SCORING_V0111.md`]"
    "(METHYLATION_DATA_AND_SCORING_V0111.md)"
  )
  if methylation_link not in text:
    anchor = (
      "- [`DATA_ACQUISITION_AND_REPRODUCTION_V0110.md`]"
      "(DATA_ACQUISITION_AND_REPRODUCTION_V0110.md)"
    )
    text = text.replace(anchor, anchor + "\n" + methylation_link)
  body = r"""
Promoter methylation is an epigenetic subcomponent of the existing functional-microniche regulatory-network domain. It is not a new top-level RSES-Onco domain and therefore cannot receive an additional independent global weight.

The supported sources are the DepMap custom-download dataset `Methylation (1kb upstream TSS)` and the traceable historical `CCLE_RRBS_TSS1kb_20181022` matrix. The resume workflow searches under `DEPMAP_DIR`, or the user may provide an explicit file:

```bash
export DEPMAP_DIR="$NEW/data/raw/depmap"
export METHYLATION="$DEPMAP_DIR/Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv"
```

The regulatory-network internal weights are:

```text
DoRothEA TF-association divergence       0.32
TF-expression-profile divergence        0.28
JASPAR/FIMO promoter-motif divergence    0.20
promoter-methylation context             0.20
```

The methylation context itself combines pair promoter-profile divergence (0.50) and conditional target-promoter hypomethylation in lost-gene-loss versus intact models (0.50). Missing methylation remains NA and lowers regulatory subcoverage; it is not converted to zero.

Run the complete recalculation after setting `METHYLATION`:

```bash
MPLBACKEND=Agg \
PYTHONUNBUFFERED=1 \
bash scripts/resume_wgcna_regulatory_pipeline.sh resume-regulatory
```

The final ranking must report:

```text
score_version=RSES-Onco-expanded-v0.10.10
methylation_semantics_version=promoter-methylation-context-v1
regulatory_layer_version=wgcna-promoter-methylation-regulatory-v3
```

Supplementary Figure S52 and the exact source-data table report the methylation layer. Beta-values are association evidence and are not direct proof of silencing, derepression, compensation or synthetic lethality.
"""
  write(path, replace_block(text, "Promoter methylation layer", body))


def patch_acquisition() -> None:
  path = "docs/DATA_ACQUISITION_AND_REPRODUCTION_V0110.md"
  text = read(path)
  body = r"""
The promoter methylation input is acquired separately from the four required DepMap matrices. Use the DepMap custom-download resource named `Methylation (1kb upstream TSS)` or the historical CCLE RRBS file `CCLE_RRBS_TSS1kb_20181022.txt.gz`.

Accepted default names under `DEPMAP_DIR` are:

```text
Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv
Methylation_1kb_upstream_TSS.csv
CCLE_RRBS_TSS1kb_20181022.txt.gz
CCLE_RRBS_TSS1kb_20181022.txt
```

A non-default path is supplied with:

```bash
export METHYLATION="/absolute/path/to/promoter_methylation.csv"
```

Validate that values are beta-like ratios in [0,1], that model or cell-line identifiers map through `Model.csv`, and that candidate genes have promoter/TSS features when expected. The parser accepts ModelID-row, long ModelID/gene/value and historical feature-row layouts. Multiple features assigned to one gene are collapsed by the median and feature counts are retained.

Generate the methylation-aware regulatory layer with:

```bash
python -u scripts/integrate_methylation_regulatory_layer.py \
  --methylation "$METHYLATION" \
  --copy-number "$DEPMAP_DIR/OmicsCNGeneWGS.csv" \
  --models "$DEPMAP_DIR/Model.csv" \
  --candidates data/processed/expanded_candidate_universe.tsv \
  --input data/processed/regulatory/expanded_pair_functional_evidence_by_cancer_pre_methylation.tsv \
  --output data/processed/regulatory/expanded_pair_functional_evidence_by_cancer.tsv
```

The canonical resume workflow performs this stage automatically. Source absence, unmapped genes, unmapped models and insufficient groups are recorded explicitly rather than converted to biological zero. Detailed formulas and accepted layouts are documented in [`METHYLATION_DATA_AND_SCORING_V0111.md`](METHYLATION_DATA_AND_SCORING_V0111.md).
"""
  write(path, replace_block(text, "Promoter methylation data acquisition", body))


def patch_supplement() -> None:
  path = "supplementary/Supplementary_Methods_RSES_Onco_v0110.md"
  text = read(path)
  body = r"""
Promoter methylation was evaluated from the DepMap custom-download dataset `Methylation (1kb upstream TSS)` or the historical `CCLE_RRBS_TSS1kb_20181022` matrix. Model or cell-line identifiers were mapped through `Model.csv`. Multiple promoter/TSS measurements assigned to the same gene were collapsed by the median within each model, while missing values remained missing.

For a directional pair with lost gene $L$ and target gene $T$, methylation profile divergence combined Spearman correlation divergence and the median absolute beta-value difference:

$$
D_{profile}=\operatorname{mean}\left(\frac{1-\rho_{L,T}}{2},\operatorname{clip}\left(\frac{\widetilde{|\beta_L-\beta_T|}}{0.25},0,1\right)\right)
$$

Models were also stratified by the same copy-number loss threshold used in dependency and expression analyses. Conditional target hypomethylation support was:

$$
H_T=\operatorname{clip}\left(\frac{-(\widetilde{\beta}_{T,loss}-\widetilde{\beta}_{T,intact})}{0.25},0,1\right)
$$

The methylation context combined the two available terms with equal internal weights:

$$
M=0.50D_{profile}+0.50H_T
$$

with eligibility-aware subcoverage when either term was missing. The regulatory-network internal composition was 0.32 DoRothEA TF-association divergence, 0.28 TF-expression-profile divergence, 0.20 JASPAR/FIMO promoter-motif divergence and 0.20 promoter-methylation context. The global regulatory-network weight and all seven top-level RSES-Onco weights were unchanged.

Methylation values were treated as association evidence. They were not interpreted automatically as direct causal silencing, promoter occupancy, compensatory transcription, synthetic lethality, therapeutic validation or clinical efficacy.

### Promoter methylation references

1. DepMap Project. Custom downloads and cell-line molecular data resources.
2. Ghandi M, Huang FW, Jané-Valbuena J, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature*. 2019;569:503-508.
3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300.
"""
  write(path, replace_block(text, "Promoter methylation methods", body))


def patch_manuscript() -> None:
  path = "manuscript/RSES_Onco_intro_methods_draft_v0110.md"
  text = read(path)
  body = r"""
Promoter methylation was incorporated as an epigenetic subcomponent of the regulatory-network microniche domain. The analysis supported the DepMap `Methylation (1kb upstream TSS)` custom-download matrix and the historical CCLE RRBS TSS1kb matrix. Multiple promoter features assigned to the same gene were collapsed by the median within each model.

For each candidate pair and cancer lineage, the analysis quantified divergence between the lost- and target-gene promoter-methylation profiles and tested whether target-promoter methylation was lower in models with loss of the origin gene than in intact models. The latter comparison used a Mann-Whitney U test, with Benjamini-Hochberg adjustment globally and within cancer context. The methylation context combined pair-profile divergence and conditional target hypomethylation with equal internal weights and explicit missing-data coverage.

Methylation shared the existing regulatory-network domain with DoRothEA TF associations, TF-expression consistency and JASPAR/FIMO promoter-motif evidence. Its inclusion did not alter the seven top-level RSES-Onco weights. Promoter methylation was considered association evidence and not direct proof of transcriptional silencing or causal compensation.
"""
  write(path, replace_block(text, "Promoter methylation analysis", body))


def main() -> None:
  patch_protocol()
  patch_acquisition()
  patch_supplement()
  patch_manuscript()
  print("Methylation-aware canonical documentation synchronized.")


if __name__ == "__main__":
  main()
