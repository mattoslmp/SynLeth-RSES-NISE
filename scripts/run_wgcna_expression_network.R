#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: run_wgcna_expression_network.R ",
    "EXPRESSION_TSV PAIRS_TSV OUTPUT_DIR CANCER"
  )
}

expression_path <- args[[1]]
pairs_path <- args[[2]]
output_dir <- args[[3]]
cancer <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression <- read.delim(
  expression_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!"ModelID" %in% names(expression)) {
  stop("WGCNA input lacks ModelID")
}
if (nrow(expression) < 20) {
  stop(sprintf("WGCNA requires at least 20 samples; observed %d", nrow(expression)))
}

model_ids <- as.character(expression$ModelID)
datExpr <- as.data.frame(expression[, setdiff(names(expression), "ModelID"), drop = FALSE])
for (column in names(datExpr)) {
  datExpr[[column]] <- as.numeric(datExpr[[column]])
}
rownames(datExpr) <- model_ids

quality <- goodSamplesGenes(datExpr, verbose = 0)
if (!quality$allOK) {
  datExpr <- datExpr[quality$goodSamples, quality$goodGenes, drop = FALSE]
  model_ids <- rownames(datExpr)
}
if (nrow(datExpr) < 20 || ncol(datExpr) < 50) {
  stop(
    sprintf(
      "Insufficient WGCNA matrix after QC: %d samples x %d genes",
      nrow(datExpr),
      ncol(datExpr)
    )
  )
}

correlation_method <- "bicor"
correlation_max_p_outliers <- 0.10
correlation_pearson_fallback <- "individual"
correlation_options <- list(
  use = "p",
  maxPOutliers = correlation_max_p_outliers,
  pearsonFallback = correlation_pearson_fallback
)
kme_correlation_options <- paste0(
  "use = 'p', maxPOutliers = ",
  correlation_max_p_outliers,
  ", pearsonFallback = '",
  correlation_pearson_fallback,
  "'"
)

powers <- c(1:12, seq(14, 20, by = 2))
soft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = correlation_method,
  corOptions = correlation_options,
  verbose = 0
)
fit <- soft$fitIndices
eligible_power <- fit$Power[
  is.finite(fit$SFT.R.sq)
  & fit$SFT.R.sq >= 0.80
  & fit$slope < 0
]
if (length(eligible_power)) {
  power <- min(eligible_power)
  power_rule <- "first_signed_scale_free_R2_at_least_0.80_with_negative_slope"
} else {
  valid <- which(is.finite(fit$SFT.R.sq) & fit$Power >= 4)
  if (!length(valid)) {
    power <- 6
    power_rule <- "fallback_power_6_no_finite_fit"
  } else {
    best <- valid[which.max(fit$SFT.R.sq[valid])]
    power <- fit$Power[[best]]
    power_rule <- "maximum_available_signed_scale_free_R2"
  }
}

adjacency_matrix <- adjacency(
  datExpr,
  power = power,
  type = "signed",
  corFnc = correlation_method,
  corOptions = correlation_options
)
TOM <- TOMsimilarity(adjacency_matrix, TOMType = "signed", verbose = 0)
dissTOM <- 1 - TOM
gene_tree <- hclust(as.dist(dissTOM), method = "average")
minimum_module_size <- max(20, min(40, floor(ncol(datExpr) / 10)))
dynamic_labels <- cutreeDynamic(
  dendro = gene_tree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = minimum_module_size
)
dynamic_colors <- labels2colors(dynamic_labels)
merged <- mergeCloseModules(
  datExpr,
  dynamic_colors,
  cutHeight = 0.25,
  verbose = 0
)
module_colors <- merged$colors
MEs <- orderMEs(merged$newMEs)

kme <- signedKME(
  datExpr,
  MEs,
  outputColumnName = "kME",
  corFnc = correlation_method,
  corOptions = kme_correlation_options
)
connectivity <- intramodularConnectivity(adjacency_matrix, module_colors)
module_names <- sub("^ME", "", names(MEs))
module_column <- match(module_colors, module_names)
kme_own <- rep(NA_real_, length(module_colors))
for (index in seq_along(module_colors)) {
  if (!is.na(module_column[[index]])) {
    kme_own[[index]] <- kme[index, module_column[[index]]]
  }
}

genes <- colnames(datExpr)
gene_mad <- vapply(
  datExpr,
  function(values) suppressWarnings(mad(values, na.rm = TRUE)),
  numeric(1)
)
gene_zero_mad <- !is.finite(gene_mad) | gene_mad == 0
module_eigengene_mad <- vapply(
  MEs,
  function(values) suppressWarnings(mad(values, na.rm = TRUE)),
  numeric(1)
)
module_eigengene_zero_mad <- (
  !is.finite(module_eigengene_mad)
  | module_eigengene_mad == 0
)

gene_table <- data.frame(
  cancer = cancer,
  gene = genes,
  module = module_colors,
  kME_own = kme_own,
  kTotal = connectivity$kTotal,
  kWithin = connectivity$kWithin,
  expression_mad = unname(gene_mad[genes]),
  pearson_fallback_expected = unname(gene_zero_mad[genes]),
  stringsAsFactors = FALSE
)
write.table(
  gene_table,
  file.path(output_dir, "wgcna_gene_modules.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

eigengenes <- data.frame(ModelID = rownames(MEs), MEs, check.names = FALSE)
write.table(
  eigengenes,
  file.path(output_dir, "wgcna_module_eigengenes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pairs <- read.delim(pairs_path, stringsAsFactors = FALSE)
required_pair_columns <- c("pair_id", "lost_gene", "target_gene")
if (!all(required_pair_columns %in% names(pairs))) {
  stop("Pair input lacks pair_id, lost_gene or target_gene")
}
gene_lookup <- split(gene_table, gene_table$gene)
rows <- lapply(seq_len(nrow(pairs)), function(i) {
  record <- pairs[i, , drop = FALSE]
  lost <- as.character(record$lost_gene)
  target <- as.character(record$target_gene)
  lost_index <- match(lost, genes)
  target_index <- match(target, genes)
  if (is.na(lost_index) || is.na(target_index)) {
    return(data.frame(
      cancer = cancer,
      pair_id = record$pair_id,
      lost_gene = lost,
      target_gene = target,
      wgcna_status = "gene_missing_after_wgcna_qc",
      wgcna_tom_similarity = NA_real_,
      wgcna_tom_divergence = NA_real_,
      wgcna_same_module = NA,
      wgcna_module_divergence = NA_real_,
      wgcna_lost_module = NA_character_,
      wgcna_target_module = NA_character_,
      wgcna_lost_kME = NA_real_,
      wgcna_target_kME = NA_real_,
      wgcna_kME_divergence = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  lost_row <- gene_lookup[[lost]][1, , drop = FALSE]
  target_row <- gene_lookup[[target]][1, , drop = FALSE]
  same_module <- identical(
    as.character(lost_row$module),
    as.character(target_row$module)
  ) && as.character(lost_row$module) != "grey"
  tom_similarity <- TOM[lost_index, target_index]
  data.frame(
    cancer = cancer,
    pair_id = record$pair_id,
    lost_gene = lost,
    target_gene = target,
    wgcna_status = "available",
    wgcna_tom_similarity = tom_similarity,
    wgcna_tom_divergence = 1 - tom_similarity,
    wgcna_same_module = same_module,
    wgcna_module_divergence = ifelse(same_module, 0, 1),
    wgcna_lost_module = lost_row$module,
    wgcna_target_module = target_row$module,
    wgcna_lost_kME = lost_row$kME_own,
    wgcna_target_kME = target_row$kME_own,
    wgcna_kME_divergence = abs(abs(lost_row$kME_own) - abs(target_row$kME_own)),
    stringsAsFactors = FALSE
  )
})
pair_table <- do.call(rbind, rows)
write.table(
  pair_table,
  file.path(output_dir, "wgcna_pair_metrics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  fit,
  file.path(output_dir, "wgcna_soft_threshold_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

fallback_audit <- rbind(
  data.frame(
    cancer = cancer,
    entity_type = "gene",
    entity = genes,
    mad = unname(gene_mad[genes]),
    pearson_fallback_expected = unname(gene_zero_mad[genes]),
    fallback_reason = ifelse(
      unname(gene_zero_mad[genes]),
      "zero_or_nonfinite_MAD",
      ""
    ),
    primary_correlation = correlation_method,
    fallback_correlation = "pearson",
    pearson_fallback_policy = correlation_pearson_fallback,
    max_p_outliers = correlation_max_p_outliers,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cancer = cancer,
    entity_type = "module_eigengene",
    entity = names(MEs),
    mad = unname(module_eigengene_mad[names(MEs)]),
    pearson_fallback_expected = unname(
      module_eigengene_zero_mad[names(MEs)]
    ),
    fallback_reason = ifelse(
      unname(module_eigengene_zero_mad[names(MEs)]),
      "zero_or_nonfinite_MAD",
      ""
    ),
    primary_correlation = correlation_method,
    fallback_correlation = "pearson",
    pearson_fallback_policy = correlation_pearson_fallback,
    max_p_outliers = correlation_max_p_outliers,
    stringsAsFactors = FALSE
  )
)
write.table(
  fallback_audit,
  file.path(output_dir, "wgcna_correlation_fallback.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

diagnostics <- data.frame(
  cancer = cancer,
  samples = nrow(datExpr),
  genes = ncol(datExpr),
  selected_power = power,
  power_selection_rule = power_rule,
  network_type = "signed",
  correlation = correlation_method,
  correlation_policy = (
    "bicor_primary_with_individual_Pearson_fallback_for_zero_MAD_only"
  ),
  max_p_outliers = correlation_max_p_outliers,
  pearson_fallback = correlation_pearson_fallback,
  signed_kme_correlation = correlation_method,
  signed_kme_max_p_outliers = correlation_max_p_outliers,
  signed_kme_pearson_fallback = correlation_pearson_fallback,
  zero_mad_gene_count = sum(gene_zero_mad),
  zero_mad_module_eigengene_count = sum(module_eigengene_zero_mad),
  pearson_fallback_entity_count = (
    sum(gene_zero_mad) + sum(module_eigengene_zero_mad)
  ),
  tom_type = "signed",
  min_module_size = minimum_module_size,
  merge_cut_height = 0.25,
  module_count_excluding_grey = length(setdiff(unique(module_colors), "grey")),
  wgcna_version = as.character(packageVersion("WGCNA")),
  stringsAsFactors = FALSE
)
write.table(
  diagnostics,
  file.path(output_dir, "wgcna_run_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
