# Normalize batch-corrected counts using Bailey Garb's edgeR procedure.
# The full gene background is needed for the library sizes and normalization
# factors, even though only selected gene signatures are scored downstream.
module_root <- Sys.getenv("E6_PATIENT_MODULE_ROOT")
if (!nzchar(module_root)) {
  stop("Run this module through patient_scores/run.py.")
}
source(file.path(module_root, "scripts", "common.R"))
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("Install the R package edgeR.")
}

counts <- read.csv(file.path(output_dir, "batch_corrected_counts_n230.csv"),
                   check.names = FALSE)

# Sort by symbol before filtering. Retain genes with CPM >=1 in at least three
# patients, then recompute library sizes on the retained counts.
counts <- counts[order(counts[[1]]), , drop = FALSE]
dge <- edgeR::DGEList(as.matrix(counts[, -1, drop = FALSE]), genes = counts[[1]])
expressed <- rowSums(edgeR::cpm(dge) >= 1) >= 3
dge <- dge[expressed, ]
dge$samples$lib.size <- colSums(dge$counts)

# Use upper-quartile factors and edgeR's original default prior count for
# logCPM. Retain the first duplicate symbol only after normalization.
dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
log_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
rownames(log_cpm) <- dge$genes$genes
colnames(log_cpm) <- names(counts)[-1]
log_cpm <- log_cpm[!duplicated(rownames(log_cpm)), , drop = FALSE]
stopifnot(ncol(log_cpm) == 230, all(is.finite(log_cpm)))
write.csv(log_cpm, file.path(output_dir, "log_cpm_n230.csv"))
cat("Calculated upperquartile logCPM for", nrow(log_cpm), "genes.\n")
