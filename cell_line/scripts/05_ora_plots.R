# Draw Figure S8A-E from complete overrepresentation results.
# Inputs: results/ora/all_ORA_results.tsv and comparison_order.tsv.
# Outputs: a combined GO-BP/Hallmark/KEGG PDF and its displayed-pathway TSV.
# FIGURE_PANEL selects S8A-S8E; empty means all five comparisons.
# Up/down bars show signed DEG counts; color shows signed -log2(adjusted P).
# Each collection contributes up to ten significant pathways per direction.

# Preserve physical row spacing while allowing enough width for full pathway names.
render_reference_enrichment <- function(project) {
  method <- "ORA"
  main_subdir <- ""
  project <- normalizePath(project, mustWork = TRUE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gtable)
  })
  source(file.path(project, "scripts/plot_export_helpers.R"), local = TRUE)
  initialize_manuscript_export(project)
  read_tsv <- function(p) read.delim(p, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
  write_tsv <- function(x, p) write.table(x, p, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  figdir <- file.path(project, "figures")
  data_dir <- file.path(figdir, "source_data")
  specs <- read_tsv(file.path(project, "config/comparison_order.tsv"))
  specs <- specs[specs$comparison_order != 4, ]
  all <- read_tsv(file.path(project, "results/ora/all_ORA_results.tsv"))
  collections <- c("GO_BP", "HALLMARK", "KEGG_PATHWAY")
  titles <- c(GO_BP = "GO_BP", HALLMARK = "HALLMARK", KEGG_PATHWAY = "KEGG")
  stopifnot(setequal(unique(all$contrast), specs$comparison), all(is.finite(all$padj)))
  stopifnot(all(all$direction %in% c("up", "down")))
  all$plot_direction <- ifelse(all$direction == "up", 1, -1)
  all$bar_gene_count <- all$gene_count
  all$signed_gene_count <- all$plot_direction * all$bar_gene_count
  all$signed_minus_log2_FDR <- all$plot_direction * (-log2(pmax(all$padj, .Machine$double.xmin)))
  # Convert MSigDB identifiers to readable names; preserve official KEGG case.
  label_text <- function(x) {
    original <- x
    from_msigdb <- grepl("^(GOBP_|HALLMARK_)", x)
    x <- tolower(gsub("_", " ", sub("^(GOBP_|HALLMARK_)", "", x)))
    acronyms <- c(dna = "DNA", rna = "RNA", mrna = "mRNA", rrna = "rRNA", trna = "tRNA", myc = "MYC",
      e2f = "E2F", g2m = "G2M", mtorc1 = "MTORC1", tnfa = "TNFA", nfkb = "NFKB", kras = "KRAS",
      il2 = "IL2", il6 = "IL6", stat3 = "STAT3", stat5 = "STAT5", jak = "JAK", uv = "UV", pi3k = "PI3K",
      akt = "AKT", mtor = "MTOR", ecm = "ECM", mapk = "MAPK", atp = "ATP", nad = "NAD", p53 = "p53",
      isr = "ISR")
    for (w in names(acronyms)) x <- gsub(paste0("\\b", w, "\\b"), acronyms[[w]], x, perl = TRUE)
    # KEGG already supplies display-ready names: preserve IgSF CAM, Type II,
    # PI3K-Akt and other official capitalization instead of lowercasing them.
    x[!from_msigdb] <- original[!from_msigdb]
    x <- gsub("[[:space:]]+", " ", trimws(x))
    # Sentence case must not change standard mixed-case RNA abbreviations.
    capitalise <- !grepl("^(mRNA|rRNA|tRNA)\\b", x, perl = TRUE)
    if (any(capitalise))
      substr(x[capitalise], 1, 1) <- toupper(substr(x[capitalise], 1, 1))
    x
  }
  label_pt <- 8
  row_pitch_mm <- 3.5
  bar_fraction <- 0.88
  panel_width_mm <- 25
  font <- .manuscript_export_state$font_family
  # Shared typography and physical margins for the three database components.
  theme_ref <- function() theme_classic(base_size = label_pt, base_family = font) + theme(text = element_text(family = font,
    colour = "black"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.22),
    axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.22), axis.ticks.length = unit(0.9,
      "mm"), axis.text = element_text(size = label_pt, colour = "black"), axis.text.y = element_text(margin = margin(r = 1.5),
      lineheight = 1), axis.title = element_text(size = label_pt + 0.5), axis.title.x = element_text(margin = margin(t = 3)),
    axis.title.y = element_text(margin = margin(r = 4)), plot.title = element_text(size = 10, face = "bold",
      hjust = 0.5, margin = margin(b = 4)), plot.margin = margin(2, 2, 2, 2, "mm"), legend.position = "right",
    legend.title = element_text(size = 7, lineheight = 1), legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(3, "mm"), legend.box.margin = margin(0,
      0, 0, 0))
  # Use the same natural SVG font metrics as the production export when measuring
  # each complete label. Fix the plotting rectangle in physical units, then let
  # its labels and legend determine the total width without wrapping/shrinking.
  metric_svg <- tempfile(fileext = ".svg")
  svglite::svglite(metric_svg, width = 16, height = 5)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    unlink(metric_svg)
  }, add = TRUE)
  # Reserve a fixed plotting rectangle; let full labels determine total width.
  measure_grob <- function(p, height_mm) {
    g <- ggplotGrob(p)
    pane <- which(g$layout$name == "panel")
    stopifnot(length(pane) == 1)
    g$widths[g$layout$l[pane]] <- unit(panel_width_mm, "mm")
    g$heights[g$layout$t[pane]] <- unit(height_mm, "mm")
    list(g = g, width = convertWidth(sum(g$widths), "mm", valueOnly = TRUE), height = convertHeight(sum(g$heights),
      "mm", valueOnly = TRUE))
  }
  selected_panel <- Sys.getenv("FIGURE_PANEL", "")
  panel_ids <- paste0("S8", LETTERS[1:5])
  indices <- if (nzchar(selected_panel))
    match(selected_panel, panel_ids) else seq_len(nrow(specs))
  stopifnot(!anyNA(indices))
  for (i in indices) {
    spec <- specs[i, ]
    cmp <- spec$comparison
    x <- all[all$contrast == cmp, ]
    # 1. Select the ten lowest-FDR pathways in each direction/database.
    picks <- list()
    for (coll in collections) {
      z <- x[x$collection == coll & x$padj < 0.05, ]
      up <- z[z$plot_direction > 0, ]
      down <- z[z$plot_direction < 0, ]
      up <- head(up[order(up$padj, up$pvalue, up$set_id), ], 10)
      down <- head(down[order(down$padj, down$pvalue, down$set_id), ], 10)
      z <- rbind(up, down[rev(seq_len(nrow(down))), ])
      z$display_row <- seq_len(nrow(z))
      z$display_label <- label_text(z$set_name)
      z$plot_variant <- rep("compact", nrow(z))
      picks[[coll]] <- z

    }
    selected <- do.call(rbind, picks)
    rownames(selected) <- NULL
    nslots <- max(c(6L, vapply(picks, nrow, integer(1))))
    panel_height <- row_pitch_mm * nslots
    selected$plot_y <- nslots + 1 - selected$display_row
    selected$bar_height_rows <- rep(bar_fraction, nrow(selected))
    selected$row_pitch_mm <- rep(row_pitch_mm, nrow(selected))
    # 2. Share axis/color scales across the three components of this comparison.
    bar_range <- range(c(0, selected$signed_gene_count))
    if (diff(bar_range) == 0)
      bar_range <- c(-1, 1)
    xlimits <- bar_range + c(-1, 1) * max(diff(bar_range) * 0.035, 0.5)
    # Small, fixed-width bar areas need a hard tick limit: pretty(n) can return
    # more than n values and crowd adjacent three-digit labels together.
    xbreaks <- (scales::breaks_pretty(n = 3))(xlimits)
    xbreaks <- sort(unique(c(0, xbreaks[xbreaks >= xlimits[1] & xbreaks <= xlimits[2]])))
    if (length(xbreaks) > 3)
      xbreaks <- unique(xbreaks[c(1, which.min(abs(xbreaks)), length(xbreaks))])
    colour_limit <- max(5, ceiling(max(c(0, abs(selected$signed_minus_log2_FDR)))/5) * 5)
    panel_grobs <- list()
    for (coll in collections) {
      z <- selected[selected$collection == coll, ]
      p <- ggplot(z, aes(signed_gene_count, plot_y, fill = signed_minus_log2_FDR)) + geom_col(orientation = "y",
        width = bar_fraction, colour = "black", linewidth = 0.22) + scale_y_continuous(breaks = z$plot_y,
        labels = z$display_label, limits = c(0.5, nslots + 0.5), expand = c(0, 0)) + scale_x_continuous(limits = xlimits,
        breaks = if (nrow(z))
          xbreaks else NULL, expand = c(0, 0)) + scale_fill_gradientn(colours = c("#2166AC", "#4393C3",
        "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), limits = c(-colour_limit,
        colour_limit), breaks = c(-colour_limit, 0, colour_limit), name = "Signed\n-log2(FDR)") +
        guides(fill = guide_colourbar(direction = "vertical", barwidth = unit(2.7, "mm"), barheight = unit(17,
          "mm"), title.position = "top", ticks.colour = "white", frame.colour = NA)) + labs(x = "Gene Count",
        y = "Function", title = titles[[coll]]) + theme_ref()
      p <- if (nrow(z))
        p + geom_vline(xintercept = 0, linewidth = 0.22) else p + annotate("text", x = mean(xlimits), y = (nslots + 1)/2, label = "No significant\npathways",
        family = font, size = 7/ggplot2::.pt, lineheight = 1.05)
      g <- measure_grob(p, panel_height)
      panel_grobs[[coll]] <- g
    }
    # An empty database has no x tick labels. Align every component row so its
    # title and frame remain level with nonempty neighboring database panels.
    aligned_heights <- do.call(unit.pmax, lapply(panel_grobs, function(g) g$g$heights))
    for (coll in collections) {
      z <- selected[selected$collection == coll, ]
      g <- panel_grobs[[coll]]
      g$g$heights <- aligned_heights
      g$height <- convertHeight(sum(aligned_heights), "mm", valueOnly = TRUE)
      panel_grobs[[coll]] <- g

    }
    # 3. Combine the components without truncating labels or creating extra figures.
    widths <- vapply(panel_grobs, function(g) g$width, numeric(1))
    heights <- vapply(panel_grobs, function(g) g$height, numeric(1))
    title_height <- 6
    combined <- gtable(widths = unit(widths, "mm"), heights = unit(c(title_height, max(heights)),
      "mm"))
    combined <- gtable_add_grob(combined, textGrob(spec$label, gp = gpar(fontfamily = font, fontface = "bold",
      fontsize = 10)), t = 1, l = 1, r = 3, clip = "off")
    for (j in seq_along(panel_grobs)) combined <- gtable_add_grob(combined, panel_grobs[[j]]$g, t = 2,
      l = j, clip = "off")
    report_cmp <- paste(spec$numerator, "vs", spec$denominator, sep = "_")
    stem <- paste0("FigureS8", LETTERS[[i]])
    manuscript_export(combined, stem, sum(widths), sum(c(title_height, max(heights))), figdir = figdir,
      subdir = main_subdir, description = paste(method, spec$label, "GO_BP / HALLMARK / KEGG; top 10 significant pathways per direction; single-line reference layout."))
    write_tsv(selected, file.path(data_dir, paste0(stem, ".tsv")))
    message(sprintf("Rendered %s %s: %.1f x %.1f mm; %d rows maximum", method, cmp, sum(widths),
      title_height + max(heights), nslots))
  }

}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  render_reference_enrichment(args[[1]])
}
