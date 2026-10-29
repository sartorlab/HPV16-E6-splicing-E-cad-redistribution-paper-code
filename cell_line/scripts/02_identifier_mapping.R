#!/usr/bin/env Rscript
# Select one Ensembl representative per Entrez gene for overrepresentation tests.  Inputs: first DE
# table (all contrasts have identical all-sample baseMean) and resources/ensembl_entrez_mapping.tsv
# from the Ensembl109 annotation.  Output: results/mapping/identifier_mapping.tsv, consumed by
# 03_ora.R.  Assumption: representatives depend on expression, never on significance.

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[[1]], mustWork = TRUE)

# Read identifiers as strings so their formatting cannot be altered by inference.
read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = c("NA",
    ""))
}
write_tsv <- function(table, path) {
  write.table(table, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
reference <- read_tsv(file.path(project, "results/de/res_01_E6IvsFL.tsv"))
reference$baseMean <- as.numeric(reference$baseMean)
dir.create(file.path(project, "results/mapping"), showWarnings = FALSE)

# 1. Reject missing, nonnumeric or ambiguous Entrez mappings. A gene with multiple mappings is
# excluded rather than assigning a favorable pathway ID.
mapping <- read_tsv(file.path(project, "resources/ensembl_entrez_mapping.tsv"))
mapping <- unique(mapping[, c("ensembl_gene_id", "entrez_id")])
mapping <- mapping[!is.na(mapping$entrez_id) & grepl("^[0-9]+$", mapping$entrez_id), ]
by_ensembl <- split(mapping$entrez_id, mapping$ensembl_gene_id)
candidate_ids <- lapply(reference$ensembl_gene_id, function(id) sort(unique(by_ensembl[[id]])))
selected <- data.frame(ensembl_gene_id = reference$ensembl_gene_id, gene_symbol = reference$gene_symbol,
  baseMean = reference$baseMean, n_entrez_mappings = lengths(candidate_ids), entrez_candidates = vapply(candidate_ids,
    paste, character(1), collapse = ";"), stringsAsFactors = FALSE)
selected$entrez_id <- vapply(candidate_ids, function(ids) {
  if (length(ids) == 1)
    ids else NA_character_
}, character(1))
selected$mapping_status <- ifelse(selected$n_entrez_mappings == 0, "unmapped", ifelse(selected$n_entrez_mappings >
  1, "ambiguous_ensembl_to_entrez", "candidate"))

# 2. For multiple Ensembl genes sharing one Entrez ID, select the greatest all-sample baseMean;
# exact ties are broken by lexicographic Ensembl identifier.
eligible <- which(selected$mapping_status == "candidate")
eligible <- eligible[order(-selected$baseMean[eligible], selected$ensembl_gene_id[eligible], method = "radix")]
representatives <- eligible[!duplicated(selected$entrez_id[eligible])]
selected$mapping_status[eligible] <- "duplicate_entrez_nonrepresentative"
selected$mapping_status[representatives] <- "representative"
representative_ids <- setNames(selected$ensembl_gene_id[representatives], selected$entrez_id[representatives])
selected$representative_ensembl <- unname(representative_ids[selected$entrez_id])
write_tsv(selected, file.path(project, "results/mapping/identifier_mapping.tsv"))
