#!/usr/bin/env Rscript
# Draw Figure S7A-F from complete DESeq2 result tables; no model is refitted.  Inputs:
# results/de/res_*.tsv, comparison_order.tsv, and an installed font.  Outputs: one PDF plus
# plotted-gene and label-coordinate TSVs per selected panel.  FIGURE_PANEL can select S7A-S7F;
# empty means all six panels.  Only finite LFC/padj genes are plotted. DEGs require padj<0.05 and
# |LFC|>1.  Label the top 20 DEGs by padj, P, decreasing |LFC|, then Ensembl identifier.

# Render fixed-size volcanoes while retaining the input point coordinates.
render_volcanoes <- function(project) {
  project <- normalizePath(project, mustWork = TRUE)
  if (!nzchar(Sys.getenv("TZ")))
    Sys.setenv(TZ = "UTC")
  suppressPackageStartupMessages(library(ggplot2))
  source(file.path(project, "scripts/plot_export_helpers.R"), local = TRUE)
  initialize_manuscript_export(project)
  figdir <- file.path(project, "figures")
  datadir <- file.path(figdir, "source_data")
  read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
  write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE,
    na = "NA")
  specs <- read_tsv(file.path(project, "config/comparison_order.tsv"))
  stopifnot(identical(as.integer(specs$comparison_order), 1:6))
  selected_panel <- Sys.getenv("FIGURE_PANEL", "")
  if (nzchar(selected_panel))
    stopifnot(selected_panel %in% paste0("S7", LETTERS[1:6]))
  palette <- c(Up = "#B2182B", Down = "#2166AC", N.S. = "#BDBDBD")
  width_mm <- 100
  height_mm <- 100
  gene_label_pt <- 8.5
  leader_width_mm <- 0.12
  for (i in seq_len(nrow(specs))) {
    comp <- specs$comparison[[i]]
    panel_id <- paste0("S7", c("A", "B", "C", "F", "D", "E")[[i]])
    if (nzchar(selected_panel) && panel_id != selected_panel)
      next
    stem <- paste0("Figure", panel_id)
    tab <- read_tsv(file.path(project, "results/de", sprintf("res_%02d_%s.tsv", i, comp)))
    stopifnot(all(tab$comparison == comp), all(tab$comparison_order == i))
    # 1. Select finite points and define up/down categories without changing NA values.
    valid <- is.finite(tab$padj) & is.finite(tab$log2FoldChange)
    tab <- tab[valid, , drop = FALSE]
    stopifnot(!anyDuplicated(tab$ensembl_gene_id), all(tab$padj >= 0 & tab$padj <= 1))
    tab$minus_log10_padj <- -log10(pmax(tab$padj, .Machine$double.xmin))
    tab$category <- "N.S."
    tab$category[tab$padj < 0.05 & tab$log2FoldChange > 1] <- "Up"
    tab$category[tab$padj < 0.05 & tab$log2FoldChange < -1] <- "Down"
    tab$direction <- unname(c(Up = "Higher in numerator", Down = "Higher in denominator", N.S. = "Below threshold")[tab$category])
    # 2. Rank the significant genes for deterministic top-20 labeling.
    sig <- tab[tab$category != "N.S.", , drop = FALSE]
    sig <- sig[order(sig$padj, sig$pvalue, -abs(sig$log2FoldChange), sig$ensembl_gene_id), ]
    label_rows <- head(sig, 20L)
    label_rows$significance_rank <- seq_len(nrow(label_rows))
    label_rows$display_label <- ifelse(is.na(label_rows$gene_symbol) | !nzchar(label_rows$gene_symbol) |
      label_rows$gene_symbol == ".", label_rows$ensembl_gene_id, label_rows$gene_symbol)
    duplicate_symbol <- duplicated(label_rows$display_label) | duplicated(label_rows$display_label,
      fromLast = TRUE)
    label_rows$display_label[duplicate_symbol] <- paste0(label_rows$display_label[duplicate_symbol],
      " (", label_rows$ensembl_gene_id[duplicate_symbol], ")")
    tab$is_labeled <- tab$ensembl_gene_id %in% label_rows$ensembl_gene_id
    tab$significance_rank <- match(tab$ensembl_gene_id, label_rows$ensembl_gene_id)
    tab$point_colour <- unname(palette[tab$category])
    # Draw significant genes last, keeping this exact drawing order in the TSV.
    tab <- tab[order(match(tab$category, c("N.S.", "Down", "Up"))), ]
    counts <- table(factor(tab$category, levels = names(palette)))
    legend_labels <- paste0(names(palette), " (", formatC(as.integer(counts), format = "d", big.mark = ","),
      ")")
    names(legend_labels) <- names(palette)
    # 3. Plot all eligible genes; the dashed lines mark the DEG thresholds.
    v <- ggplot(tab, aes(log2FoldChange, minus_log10_padj)) + geom_point(aes(colour = category),
      size = 0.55, alpha = 0.7) + geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25,
      colour = "#888888") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25,
      colour = "#888888") + scale_colour_manual(values = palette, breaks = names(palette), labels = legend_labels,
      drop = FALSE) + scale_x_continuous(limits = range(tab$log2FoldChange), oob = scales::oob_keep,
      expand = expansion(mult = c(0.16, 0.16)), n.breaks = 5) + scale_y_continuous(limits = range(tab$minus_log10_padj),
      oob = scales::oob_keep, expand = expansion(mult = c(0.025, 0.15)), n.breaks = 5) + labs(title = specs$label[[i]],
      x = "Log2 fold change", y = "−Log10(FDR)") + guides(colour = guide_legend(nrow = 1, byrow = TRUE,
      override.aes = list(size = 2.2, alpha = 1))) + theme_classic(base_family = .manuscript_export_state$font_family,
      base_size = 10.5) + theme(text = element_text(family = .manuscript_export_state$font_family,
      colour = "black"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1, "mm"), axis.text = element_text(size = 10.5, colour = "black"),
      axis.title = element_text(size = 11, colour = "black"), axis.title.x = element_text(margin = margin(t = 3)),
      axis.title.y = element_text(margin = margin(r = 3)), plot.title = element_text(size = 13,
        face = "bold", hjust = 0.5, margin = margin(b = 3)), plot.margin = margin(4, 5, 4, 4),
      legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 10),
      legend.key = element_blank(), legend.key.size = grid::unit(2.8, "mm"), legend.key.spacing.x = grid::unit(1.8,
        "mm"), legend.margin = margin(0, 0, 1, 0), legend.box.spacing = grid::unit(1, "mm"))
    panel <- ggplot_build(v)$layout$panel_params[[1]]
    # Measure the actual plot rectangle, then place labels using conservative system-font glyph
    # bounds. This avoids unresolved overlaps in dense one-sided top-20 lists without changing any
    # gene's scatter coordinates.
    staging <- tempfile(pattern = ".volcano_layout_", tmpdir = figdir)
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    probe <- file.path(staging, "panel.svg")
    ggsave(probe, v, device = svglite::svglite, width = width_mm, height = height_mm, units = "mm",
      bg = "white")
    xml <- xml2::read_xml(probe)
    first_point <- xml2::xml_find_first(xml, "//d1:circle")
    group <- xml2::xml_find_first(first_point, "ancestor::d1:g[@clip-path][1]")
    clip_id <- sub("^url\\(#(.*)\\)$", "\\1", xml2::xml_attr(group, "clip-path"))
    clip <- xml2::xml_find_first(xml, paste0("//d1:clipPath[@id='", clip_id, "']/d1:rect"))
    panel_width <- as.numeric(xml2::xml_attr(clip, "width"))
    panel_height <- as.numeric(xml2::xml_attr(clip, "height"))
    stopifnot(is.finite(panel_width), is.finite(panel_height), panel_width > 0, panel_height > 0)
    layout_input <- list(panel_width_pt = panel_width, panel_height_pt = panel_height, x_range = panel$x.range,
      y_range = panel$y.range, font_path = .manuscript_export_state$font_path, font_size_pt = gene_label_pt,
      labels = label_rows[, c("ensembl_gene_id", "display_label", "log2FoldChange", "minus_log10_padj",
        "category")], points = data.frame(x = tab$log2FoldChange, y = tab$minus_log10_padj))
    input_path <- file.path(staging, "input.json")
    output_path <- file.path(staging, "output.json")
    jsonlite::write_json(layout_input, input_path, auto_unbox = TRUE, digits = 16)
    .manuscript_run(Sys.getenv("PYTHON", "python3"), c(file.path(project, "scripts/volcano_label_layout.py"),
      input_path, output_path), "Deterministic volcano label placement")
    placed <- jsonlite::read_json(output_path, simplifyVector = TRUE)
    unlink(staging, recursive = TRUE)
    stopifnot(identical(placed$labels$ensembl_gene_id, label_rows$ensembl_gene_id))
    stopifnot(all(c("label_x", "label_y", "leader_x", "leader_y", "leader_start_x", "leader_start_y",
      "draw_leader", "leader_length_pt", "leader_ambiguous") %in% names(placed$labels)), !anyNA(placed$labels$draw_leader))
    label_rows$label_x <- placed$labels$label_x
    label_rows$label_y <- placed$labels$label_y
    label_rows$leader_x <- placed$labels$leader_x
    label_rows$leader_y <- placed$labels$leader_y
    label_rows$leader_start_x <- placed$labels$leader_start_x
    label_rows$leader_start_y <- placed$labels$leader_start_y
    label_rows$draw_leader <- placed$labels$draw_leader
    label_rows$leader_length_pt <- placed$labels$leader_length_pt
    label_rows$leader_ambiguous <- placed$labels$leader_ambiguous
    label_rows$label_colour <- unname(palette[label_rows$category])
    # Freeze the original coordinate range so displaced labels cannot rescale the underlying
    # volcano. Leaders stop at the text box, not its center.  Borderless white backing keeps
    # crossing leaders/cutoffs out of glyphs.
    v <- v + coord_cartesian(xlim = panel$x.range, ylim = panel$y.range, expand = FALSE) + geom_segment(data = label_rows[label_rows$draw_leader,
      , drop = FALSE], aes(x = leader_start_x, y = leader_start_y, xend = leader_x, yend = leader_y),
      colour = "#999999", linewidth = leader_width_mm, inherit.aes = TRUE) + geom_label(data = label_rows,
      aes(x = label_x, y = label_y, label = display_label, colour = category), family = .manuscript_export_state$font_family,
      size = gene_label_pt/ggplot2::.pt, show.legend = FALSE, fill = "white", border.colour = NA,
      linewidth = 0, label.padding = grid::unit(0.1, "mm"), label.r = grid::unit(0, "mm"), hjust = 0.5,
      vjust = 0.5)
    # 4. Export the panel and the exact scientific values represented by its marks.
    manuscript_export(v, stem, width_mm, height_mm, figdir = figdir, description = paste0(specs$label[[i]],
      "; padj < 0.05 and |log2FC| > 1; top 20 total DEGs ranked by padj, pvalue, decreasing |log2FC|, then Ensembl ID; Up/Down/N.S. counts partition finite plotted genes."))
    write_tsv(tab, file.path(datadir, paste0(stem, ".tsv")))
    write_tsv(label_rows, file.path(datadir, paste0(stem, "_labels.tsv")))
    message("Rendered ", stem, ": ", paste(legend_labels, collapse = ", "), "; labels = ", nrow(label_rows))
  }

}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  render_volcanoes(if (length(args))
    args[[1]] else ".")
}
