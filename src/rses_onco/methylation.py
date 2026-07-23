from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from .utils import canonical_gene_name


@dataclass(frozen=True)
class MethylationPairMetrics:
  lost_gene: str
  target_gene: str
  cancer: str
  n_samples: int
  lost_median_beta: float | None
  target_median_beta: float | None
  median_delta_beta: float | None
  methylation_spearman_rho: float | None
  promoter_methylation_context_score: float | None
  paired_wilcoxon_p_value: float | None
  evidence_status: str
  absence_reason: str


def promoter_methylation_context_score(
  lost_beta: float | None,
  target_beta: float | None,
) -> float | None:
  """Return directional promoter-methylation support in [0, 1].

  The score is high only when the hypothesized lost gene is promoter-methylated
  and the putative backup/target remains promoter-hypomethylated. Beta values are
  not interpreted as proof of silencing; they provide an epigenetic-context
  subcomponent that must be interpreted together with expression.
  """
  if lost_beta is None or target_beta is None:
    return None
  if not np.isfinite(lost_beta) or not np.isfinite(target_beta):
    return None
  lost = float(np.clip(lost_beta, 0.0, 1.0))
  target = float(np.clip(target_beta, 0.0, 1.0))
  return float(lost * (1.0 - target))


def pair_methylation_metrics(
  gene_sample: pd.DataFrame,
  lost_gene: str,
  target_gene: str,
  cancer: str,
  min_samples: int = 10,
) -> MethylationPairMetrics:
  required = {"cancer", "sample_id", "gene", "promoter_beta"}
  missing = sorted(required - set(gene_sample.columns))
  if missing:
    raise ValueError(
      "Promoter methylation gene-sample table lacks columns: "
      + ", ".join(missing)
    )

  lost = canonical_gene_name(lost_gene)
  target = canonical_gene_name(target_gene)
  work = gene_sample.loc[
    gene_sample["cancer"].astype(str).eq(str(cancer))
    & gene_sample["gene"].astype(str).isin({lost, target}),
    ["sample_id", "gene", "promoter_beta"],
  ].copy()
  work["promoter_beta"] = pd.to_numeric(
    work["promoter_beta"],
    errors="coerce",
  )
  work = work.dropna(subset=["promoter_beta"])
  if work.empty:
    return MethylationPairMetrics(
      lost,
      target,
      cancer,
      0,
      None,
      None,
      None,
      None,
      None,
      None,
      "missing",
      "genes_or_cancer_context_absent_from_methylation_matrix",
    )

  pivot = work.pivot_table(
    index="sample_id",
    columns="gene",
    values="promoter_beta",
    aggfunc="median",
  )
  if lost not in pivot.columns or target not in pivot.columns:
    return MethylationPairMetrics(
      lost,
      target,
      cancer,
      0,
      None,
      None,
      None,
      None,
      None,
      None,
      "missing",
      "one_gene_missing_after_probe_to_gene_aggregation",
    )
  paired = pivot[[lost, target]].dropna()
  if len(paired) < min_samples:
    return MethylationPairMetrics(
      lost,
      target,
      cancer,
      len(paired),
      None,
      None,
      None,
      None,
      None,
      None,
      "insufficient_sample",
      f"paired_methylation_samples_below_{min_samples}",
    )

  lost_median = float(paired[lost].median())
  target_median = float(paired[target].median())
  delta = float((paired[lost] - paired[target]).median())
  rho = float(
    spearmanr(
      paired[lost],
      paired[target],
      nan_policy="omit",
    ).statistic
  )
  if not np.isfinite(rho):
    rho = None

  differences = paired[lost] - paired[target]
  p_value: float | None = None
  if np.any(np.abs(differences.to_numpy(dtype=float)) > 0):
    result = wilcoxon(
      differences,
      alternative="greater",
      zero_method="wilcox",
      method="auto",
    )
    if np.isfinite(result.pvalue):
      p_value = float(result.pvalue)

  return MethylationPairMetrics(
    lost,
    target,
    cancer,
    len(paired),
    lost_median,
    target_median,
    delta,
    rho,
    promoter_methylation_context_score(lost_median, target_median),
    p_value,
    "observed",
    "",
  )
