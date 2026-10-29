#!/usr/bin/env Rscript
# Overrepresentation analysis for Figure S8A-E.  Inputs: six DE tables, fixed
# identifier_mapping.tsv, and exact GO-BP, Hallmark and KEGG membership/metadata tables. Only five
# contrasts are tested.  Outputs: all_ORA_results.tsv (all tested pathways, including zero
# overlaps) and ORA_run_summary.tsv (background/query sizes and counts by test family).
# Assumptions: independent up/down DEG lists, finite-padj background genes, collection-specific
# universes, and BH correction within each test family.
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
out <- file.path(project, "results", "ora")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character",
  na.strings = "NA")
write_tsv <- function(x, f) utils::write.table(x, f, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
# 1. Read the five manuscript contrasts and versioned reference definitions.
cfg <- read_tsv(file.path(project, "config", "comparison_order.tsv"))
cfg$comparison_order <- as.integer(cfg$comparison_order)
cfg <- cfg[order(cfg$comparison_order), ]
cfg$file_stem <- sprintf("%02d_%s", cfg$comparison_order, cfg$comparison)
stopifnot(identical(cfg$comparison_order, 1:6))
cfg <- cfg[cfg$comparison_order != 4L, ]
collection_order <- c("GO_BP", "HALLMARK", "KEGG_PATHWAY")
membership <- read_tsv(file.path(project, "resources", "gene_sets", "gene_set_membership.tsv"))
set_meta <- read_tsv(file.path(project, "resources", "gene_sets", "gene_set_metadata.tsv"))
mapping <- read_tsv(file.path(project, "results", "mapping", "identifier_mapping.tsv"))
representatives <- mapping[mapping$mapping_status == "representative", ]
stopifnot(!anyDuplicated(representatives$entrez_id), !anyDuplicated(representatives$ensembl_gene_id))
symbols <- setNames(representatives$gene_symbol, representatives$entrez_id)
all_results <- list()
summaries <- list()
minimum <- 15L
maximum <- 500L
family_index <- 0L
# 2. For each comparison, form significant up/down lists from DESeq2 adjusted P.  Independently
# filtered padj=NA genes are not eligible for either query or universe.
for (i in seq_len(nrow(cfg))) {
  c <- cfg[i, ]
  de <- read_tsv(file.path(project, "results", "de", paste0("res_", c$file_stem, ".tsv")))
  stopifnot(all(de$comparison == c$comparison))
  de$padj <- as.numeric(de$padj)
  de$log2FoldChange <- as.numeric(de$log2FoldChange)
  map <- mapping[match(de$ensembl_gene_id, mapping$ensembl_gene_id), ]
  valid_padj <- is.finite(de$padj) & de$padj >= 0 & de$padj <= 1
  up <- valid_padj & is.finite(de$log2FoldChange) & de$padj < 0.05 & de$log2FoldChange > 1
  down <- valid_padj & is.finite(de$log2FoldChange) & de$padj < 0.05 & de$log2FoldChange < -1
  chosen <- map$mapping_status == "representative"
  base_universe <- sort(unique(map$entrez_id[valid_padj & chosen]))
  queries <- list(up = sort(unique(map$entrez_id[up & chosen])), down = sort(unique(map$entrez_id[down &
    chosen])))
  for (coll in collection_order) {
    gm <- membership[membership$collection == coll, ]
    cm <- set_meta[set_meta$collection == coll, ]
    # Each collection defines its own annotated background. Restrict each pathway to that
    # background before applying the 15-500 gene size limits.
    source_sets <- lapply(split(gm$entrez_id, gm$set_id), unique)
    universe <- intersect(base_universe, unique(gm$entrez_id))
    universe <- sort(universe)
    N <- length(universe)
    stopifnot(N > 0L)
    sets <- lapply(source_sets, intersect, universe)
    sizes <- lengths(sets)
    sizes_for_meta <- unname(sizes[cm$set_id])
    sizes_for_meta[is.na(sizes_for_meta)] <- 0L
    sa <- data.frame(comparison_order = c$comparison_order, contrast = c$comparison, collection = coll,
      set_id = cm$set_id, set_name = cm$set_name, n_source_genes = unname(lengths(source_sets)[cm$set_id]),
      universe_size = N, set_size = sizes_for_meta, size_status = ifelse(sizes_for_meta < minimum,
        "below_min_size", ifelse(sizes_for_meta > maximum, "above_max_size", "tested")), stringsAsFactors = FALSE)
    sa$n_source_genes[is.na(sa$n_source_genes)] <- 0L
    tested_ids <- sa$set_id[sa$size_status == "tested"]
    stopifnot(length(tested_ids) > 0L)
    for (direction in c("up", "down")) {
      family_index <- family_index + 1L
      query <- intersect(queries[[direction]], universe)
      n <- length(query)
      message("ORA ", family_index, "/30: ", c$comparison, " / ", coll, " / ", direction, "; background=",
        N, ", query=", n, ", tested sets=", length(tested_ids))
      # Hypergeometric notation: N=background, M=pathway, n=query, k=overlap.  Test P[X >= k] and
      # include zero-overlap sets (P=1) in BH correction.
      hits <- lapply(sets[tested_ids], function(gs) sort(intersect(gs, query)))
      k <- lengths(hits)
      M <- unname(sizes[tested_ids])
      stopifnot(all(k <= n), all(k <= M), all(N - M - n + k >= 0))
      p <- if (n == 0L)
        rep(1, length(k)) else stats::phyper(k - 1L, M, N - M, n, lower.tail = FALSE)
      padj <- stats::p.adjust(p, method = "BH")
      md <- cm[match(tested_ids, cm$set_id), ]
      # Descriptive enrichment/odds statistics accompany the hypothesis test.  The odds ratio is
      # the sample 2x2 ratio, not Fisher's conditional MLE.
      odds_denominator <- (n - k) * (M - k)
      odds_numerator <- k * (N - M - n + k)
      odds <- odds_numerator/odds_denominator
      odds[is.nan(odds) | n == 0L] <- NA_real_
      rr <- data.frame(comparison_order = c$comparison_order, contrast = c$comparison, comparison_label = c$label,
        numerator = c$numerator, denominator = c$denominator, collection = coll, set_id = tested_ids,
        set_name = md$set_name, source_id = md$source_id, original_source = md$original_source,
        version = md$version, direction = direction, gene_count = as.integer(k), query_size = n,
        set_size = M, universe_size = N, gene_ratio = if (n > 0L)
          k/n else rep(NA_real_, length(k)), background_ratio = M/N, fold_enrichment = if (n > 0L)
          (k/n)/(M/N) else rep(NA_real_, length(k)), odds_ratio = odds, pvalue = p, padj = padj, hit_entrez = vapply(hits,
          paste, character(1), collapse = ";"), hit_symbols = vapply(hits, function(ids) {
          sy <- unname(symbols[ids])
          sy[is.na(sy) | !nzchar(sy)] <- ids[is.na(sy) | !nzchar(sy)]
          paste(sy, collapse = ";")
        }, character(1)), status = if (n == 0L)
          rep("no_query_genes", length(k)) else ifelse(k == 0L, "zero_hits", "tested"), stringsAsFactors = FALSE)
      rr <- rr[order(rr$padj, rr$pvalue, rr$set_id), ]
      all_results[[paste(c$comparison, coll, direction)]] <- rr
      summaries[[paste(c$comparison, coll, direction)]] <- data.frame(comparison_order = c$comparison_order,
        contrast = c$comparison, collection = coll, direction = direction, n_modeled_ensembl_genes = nrow(de),
        n_finite_padj_genes = sum(valid_padj), n_full_DEG_direction = sum(if (direction == "up") up else down),
        n_frozen_representatives_finite_padj = length(base_universe), n_mapped_DEG_before_collection = length(queries[[direction]]),
        universe_size = N, query_size = n, n_source_sets = nrow(cm), n_size_eligible_tests = nrow(rr),
        n_zero_hit_tests = sum(rr$gene_count == 0L), n_FDR05 = sum(rr$padj < 0.05), n_zero_pvalues = sum(rr$pvalue ==
          0), query_status = if (n == 0L)
          "no_query_genes" else "available", stringsAsFactors = FALSE)

    }
  }
}
# 3. Export complete results, before plotting selects the top ten per direction.
combined <- do.call(rbind, all_results)
rownames(combined) <- NULL
summary <- do.call(rbind, summaries)
rownames(summary) <- NULL
stopifnot(identical(unique(combined$contrast), cfg$comparison), nrow(summary) == 30L)
write_tsv(combined, file.path(out, "all_ORA_results.tsv"))
write_tsv(summary, file.path(out, "ORA_run_summary.tsv"))

message("Completed five ORA comparisons (30 multiple-testing families).")
