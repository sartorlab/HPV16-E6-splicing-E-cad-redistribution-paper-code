# Calculate the eleven patient scores needed for Figure 7B-F.
# Scientific methods and fixed gene lists follow Bailey Garb's analysis.
# For every gene, expression is ranked ACROSS THE 230 PATIENTS. Scores then
# combine those gene-specific patient ranks. Patients are never ranked by
# comparing different genes within their own expression profile.
module_root <- Sys.getenv("E6_PATIENT_MODULE_ROOT")
if (!nzchar(module_root)) {
  stop("Run this module through patient_scores/run.py.")
}
source(file.path(module_root, "scripts", "common.R"))

expr <- as.matrix(read.csv(file.path(output_dir, "log_cpm_n230.csv"),
                            row.names = 1, check.names = FALSE))
refdir <- file.path(module_root, "references", "gene_lists")
stopifnot(ncol(expr) == 230, !anyDuplicated(rownames(expr)),
          !anyDuplicated(colnames(expr)), all(is.finite(expr)))

rank_score <- function(genes, signs = rep(1, length(genes))) {
  # Select gene rows in the supplied signature order, reversing expression
  # where a negative direction is specified. apply(..., 1, rank) ranks each
  # gene across patients, producing a patient-by-gene rank matrix.
  selected <- expr[match(genes, rownames(expr)), , drop = FALSE]
  selected <- selected * signs
  ranks <- apply(selected, 1, rank)

  # Preserve base::rank defaults: average tied ranks; missing values receive
  # final sequential ranks. Consequently missing genes retain the original
  # patient-order-dependent behavior rather than being silently discarded.
  # Sum gene ranks per patient, then z-standardize those sums across patients.
  as.numeric(scale(rowSums(ranks)))
}

influence_score <- function(filename) {
  # The E6 percent-FL signature contains positively and negatively associated
  # genes. Rank each gene across patients with minimum ranks for ties. For
  # negative genes, rank negated expression; leave missing values as NA.
  genes <- read.csv(file.path(refdir, filename), check.names = FALSE)
  up <- as.character(genes$Gene.Symbol[genes$Direction == "up"])
  down <- as.character(genes$Gene.Symbol[genes$Direction == "down"])
  positive_ranks <- apply(
    expr[match(up, rownames(expr)), , drop = FALSE], 1,
    rank, ties.method = "min", na.last = "keep"
  )
  negative_ranks <- apply(
    -expr[match(down, rownames(expr)), , drop = FALSE], 1,
    rank, ties.method = "min", na.last = "keep"
  )
  ranks <- cbind(positive_ranks, negative_ranks)

  # Bailey's influence formula combines gene ranks within each patient as
  # (sum - mean) / SD. This is not the across-patient z-score used above.
  (rowSums(ranks, na.rm = TRUE) - rowMeans(ranks, na.rm = TRUE)) /
    apply(ranks, 1, sd, na.rm = TRUE)
}

read_two_column_symbols <- function(filename) {
  # These original reference files store gene symbols in their second column.
  read.table(file.path(refdir, filename), header = FALSE)[, 2]
}

# Calculate in the expression matrix's original patient order. In the EMT
# list, the first 18 genes are positive and the final four are negative.
emt <- scan(file.path(refdir, "emt_marker.txt"), what = "", quiet = TRUE)
out <- data.frame(samples = colnames(expr), check.names = FALSE)
out[["e6%FL_influence_scores"]] <- influence_score("e6%FL_influence.txt")
out$mes_score <- rank_score(read.table(
  file.path(refdir, "mesenchymal_cell_differentiation_ensembl_genes.txt")
)[, 1])
out$tcell_score <- rank_score(read_two_column_symbols("T_cell_activation_genes.txt"))
out$cell_adhesion <- rank_score(scan(
  file.path(refdir, "go_cell_adhesion.txt"), what = "", quiet = TRUE
))
out$go_resp_electron_transport <- rank_score(read.delim(
  file.path(refdir, "GO_0022904_ox_phos.txt")
)$x)
out$p_emt_zt <- rank_score(scan(
  file.path(refdir, "pemt_genes.txt"), what = "", quiet = TRUE
))
out$kegg_ox_phos_score <- rank_score(read.csv(
  file.path(refdir, "oxidative_phos_kegg.txt")
)$x)
out$kerat_score <- rank_score(read_two_column_symbols("keratinization_genes.txt"))
out$emt_score <- rank_score(emt, c(rep(1, 18), rep(-1, 4)))
out$bcell_score <- rank_score(read_two_column_symbols("B_cell_activation_genes.txt"))
out$kegg_focal_adhesion <- rank_score(read.csv(
  file.path(refdir, "kegg_focal_adhesion.csv"), header = FALSE
)[, 1])

# Sort only the finished score table for readable, stable patient ordering.
stopifnot(all(is.finite(as.matrix(out[, -1, drop = FALSE]))))
out <- out[order(out$samples), , drop = FALSE]
write.csv(out, file.path(output_dir, "patient_scores_n230.csv"), row.names = FALSE)
cat("Calculated eleven patient scores.\n")
