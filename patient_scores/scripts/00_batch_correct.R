# Apply Bailey Garb's zero-count filter and human-gene ComBat_seq correction.
# Input: the complete four-source count join from 00_prepare_counts.py.
# Output: batch_corrected_counts_n230.csv, with gene rows and patient columns.
module_root <- Sys.getenv("E6_PATIENT_MODULE_ROOT")
if (!nzchar(module_root)) {
  stop("Run this module through patient_scores/run.py.")
}
source(file.path(module_root, "scripts", "common.R"))
if (!requireNamespace("sva", quietly = TRUE)) {
  stop("Install the R package sva.")
}

count_file <- file.path(output_dir, "uncorrected_joined_counts_n230.csv")
counts <- read.csv(count_file, check.names = FALSE)
cohort <- read.csv(file.path(module_root, "data", "cohort_n230.csv"),
                   check.names = FALSE)
sample_ids <- names(counts)[-1]
stopifnot(length(sample_ids) == 230,
          !anyDuplicated(sample_ids),
          !anyDuplicated(cohort$samples),
          setequal(sample_ids, cohort$samples))

# Scale the original <50-zero-patients threshold from 210 to 230 patients.
# The resulting strict inequality is <55, not <=55. Keep duplicate gene rows
# at this stage; Bailey removed duplicate symbols after logCPM normalization.
count_matrix <- as.matrix(counts[, -1, drop = FALSE])
stopifnot(is.numeric(count_matrix), all(is.finite(count_matrix)),
          all(count_matrix >= 0))
zero_threshold <- round(50 * 230 / 210)
keep <- rowSums(count_matrix == 0) < zero_threshold
count_matrix <- count_matrix[keep, , drop = FALSE]

# Match the eight sequencing batches to count-column order. No biological
# covariates are supplied to ComBat_seq, matching the original analysis.
batch <- cohort$batch[match(sample_ids, cohort$samples)]
stopifnot(!anyNA(batch), length(unique(batch)) == 8)
adjusted <- sva::ComBat_seq(count_matrix, batch = batch)
rownames(adjusted) <- counts[[1]][keep]
write.csv(adjusted, file.path(output_dir, "batch_corrected_counts_n230.csv"))
cat("Completed ComBat_seq:", nrow(adjusted), "gene rows and",
    ncol(adjusted), "patients.\n")
