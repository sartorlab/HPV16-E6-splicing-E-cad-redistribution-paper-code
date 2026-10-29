#!/usr/bin/env Rscript
# Differential expression for S7A-F and the DEG inputs to S8A-E.  Inputs: gene_counts.tsv,
# sample_metadata.tsv, gene_annotation.tsv and the six comparison directions in
# comparison_order.tsv.  Outputs: results/de/res_*.tsv, sample size factors and a DE summary.
# Model: Bailey F. Garb's DESeq2 ~cell_line analysis, with three independent biological replicates
# per condition. Replicate numbers do not imply pairing.

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[[1]], mustWork = TRUE)
set.seed(20261016)

# Preserve missing statistics as NA, including independently filtered genes.
write_tsv <- function(table, path) {
  fwrite(table, path, sep = "\t", na = "NA", quote = "auto")
}

# 1. Align the 12 count columns with the experimental conditions.
samples <- as.data.frame(fread(file.path(project, "config/sample_metadata.tsv")))
raw <- fread(file.path(project, "inputs/gene_counts.tsv"))
stopifnot(identical(names(raw)[-1], samples$count_column), !anyDuplicated(raw$gene_id))
counts <- as.matrix(raw[, -1])
rownames(counts) <- raw$gene_id
colnames(counts) <- samples$sample_id
stopifnot(nrow(counts) == 60675, ncol(counts) == 12, all(is.finite(counts)), all(counts >= 0), all(counts ==
  floor(counts)))
storage.mode(counts) <- "integer"
rownames(samples) <- samples$sample_id
samples$cell_line <- factor(samples$cell_line, levels = c("A", "B", "C", "D"))
annotation <- as.data.frame(fread(file.path(project, "resources/gene_annotation.tsv"), na.strings = ""))
annotation <- annotation[match(rownames(counts), annotation$ensembl_gene_id), ]
stopifnot(identical(annotation$ensembl_gene_id, rownames(counts)))

# 2. Fit one model to genes with >=10 counts in >=3 samples. There is no batch term, pairing term,
# fold-change shrinkage or contrast-specific count filtering.
keep <- rowSums(counts >= 10) >= 3
stopifnot(sum(keep) == 17725L)
dds <- DESeqDataSetFromMatrix(counts[keep, ], samples, design = ~cell_line)
dds <- DESeq(dds, parallel = FALSE)

# 3. Compute the six Wald contrasts. For comparisons to PCD, preserve the PCD-first fitted
# contrast, then change only LFC and Wald-statistic signs.  Positive displayed LFC always means
# greater expression in the numerator.
comparisons <- as.data.frame(fread(file.path(project, "config/comparison_order.tsv")))
condition_codes <- c(PCD = "A", E6I = "B", FL = "C", FLm = "D")
summaries <- list()
for (i in seq_len(nrow(comparisons))) {
  comparison <- comparisons[i, ]
  numerator <- unname(condition_codes[comparison$numerator])
  denominator <- unname(condition_codes[comparison$denominator])
  if (comparison$orientation == -1L) {
    temporary <- numerator
    numerator <- denominator
    denominator <- temporary
  }
  # alpha=0.1 is retained for DESeq2 independent filtering. Reported DEGs use the separate
  # thresholds padj<0.05 and absolute log2 fold change>1.
  result <- results(dds, contrast = c("cell_line", numerator, denominator))
  table <- cbind(annotation[keep, c("ensembl_gene_id", "gene_symbol")], as.data.frame(result))
  rownames(table) <- NULL
  table$log2FoldChange <- table$log2FoldChange * comparison$orientation
  table$stat <- table$stat * comparison$orientation
  table$comparison <- comparison$comparison
  table$comparison_order <- comparison$comparison_order
  table$numerator <- comparison$numerator
  table$denominator <- comparison$denominator
  filename <- sprintf("res_%02d_%s.tsv", comparison$comparison_order, comparison$comparison)
  write_tsv(table, file.path(project, "results/de", filename))
  summaries[[i]] <- data.frame(comparison = comparison$comparison, modeled_genes = nrow(table), DEGs = sum(!is.na(table$padj) &
    table$padj < 0.05 & abs(table$log2FoldChange) > 1), missing_pvalue = sum(is.na(table$pvalue)),
    missing_padj = sum(is.na(table$padj)), results_alpha = metadata(result)$alpha, independent_filter_threshold = unname(metadata(result)$filterThreshold))
}

# 4. Export normalization and filtering summaries needed to interpret the tables.
write_tsv(do.call(rbind, summaries), file.path(project, "results/de/DE_summary.tsv"))
write_tsv(data.frame(sample_id = samples$sample_id, size_factor = sizeFactors(dds)), file.path(project,
  "results/de/sample_size_factors.tsv"))
message("Completed six DESeq2 contrasts from 12 samples; 17725 modeled genes.")
