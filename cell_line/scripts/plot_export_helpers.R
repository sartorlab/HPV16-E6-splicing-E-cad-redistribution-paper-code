# Shared vector export for cell-line manuscript panels.
# Input: a ggplot/gtable object plus physical dimensions in millimeters.
# Output: one PDF with text converted to paths, preserving the selected font.
# Requires installed DejaVu Sans (or PLOT_FONT), rsvg-convert and Ghostscript.
# SVG text is shaped naturally before outlining; temporary source files are removed.

.manuscript_export_state <- new.env(parent = emptyenv())

# Resolve the font and required executables once for this rendering process.
initialize_manuscript_export <- function(project) {
  project <- normalizePath(project, mustWork = TRUE)
  for (pkg in c("ggplot2", "systemfonts", "svglite")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Required export package missing: ", pkg)
  }
  family <- Sys.getenv("PLOT_FONT", "DejaVu Sans")
  resolved <- systemfonts::match_fonts(family)
  available <- systemfonts::system_fonts()
  if (!family %in% available$family)
    stop("Install requested PLOT_FONT: ", family)
  font_path <- resolved$path
  gs <- unname(Sys.which("gs"))
  rsvg_candidates <- unname(Sys.which("rsvg-convert"))
  if (!nzchar(gs) || !nzchar(rsvg_candidates))
    stop("Figure export needs Ghostscript gs and librsvg rsvg-convert on PATH.")
  .manuscript_export_state$project <- project
  .manuscript_export_state$gs <- gs
  .manuscript_export_state$rsvg <- rsvg_candidates[[1L]]
  # No font program is bundled or embedded into SVG source.
  .manuscript_export_state$font_face_css <- ""
  .manuscript_export_state$font_family <- family
  .manuscript_export_state$font_path <- font_path
  invisible(list(project = project, font_family = family, font_path = font_path))
}

# Guard against calling an exporter before the run/font state is initialized.
.manuscript_require_initialized <- function() {
  if (is.null(.manuscript_export_state$project)) {
    stop("Call initialize_manuscript_export(project) before exporting figures.")
  }
}

# Execute a graphics command with quoted arguments and informative failures.
.manuscript_run <- function(command, args, step) {
  # system2 invokes a shell: quote every argument, including paths with spaces.
  result <- suppressWarnings(system2(command, args = vapply(args, shQuote, character(1)), stdout = TRUE,
    stderr = TRUE))
  status <- attr(result, "status")
  if (is.null(status))
    status <- 0L
  if (status != 0L) {
    stop(step, " failed (exit ", status, "):\n", paste(tail(result, 50L), collapse = "\n"))
  }
  invisible(result)
}

# Keep vector geometry while replacing glyph programs with their outlines.
.manuscript_gs_pdf_args <- function(output, input) {
  c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7", "-dAutoRotatePages=/None",
    "-dNoOutputFonts", paste0("-sOutputFile=", output), input)
}

# Finalize a complete staged output only after the graphics command succeeds.
.manuscript_promote <- function(source, destination) {
  if (!file.exists(source) || file.info(source)$size <= 0) {
    stop("Export did not produce a nonempty file: ", source)
  }
  if (!file.rename(source, destination))
    stop("Cannot finalize export: ", destination)
  invisible(destination)
}

# Export one panel at its specified physical size; no alternate formats are emitted.
manuscript_export <- function(p, stem, width_mm, height_mm, figdir = NULL, subdir = "", description = "") {
  .manuscript_require_initialized()
  if (length(stem) != 1L || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", stem)) {
    stop("Figure stem must be a single filename without directories.")
  }
  if (length(subdir) != 1L || grepl("(^/|(^|/)\\.\\.($|/)|\\\\)", subdir)) {
    stop("Figure subdir must be a relative path within figdir.")
  }
  if (length(width_mm) != 1L || length(height_mm) != 1L || !is.finite(width_mm) || !is.finite(height_mm) ||
    width_mm <= 0 || height_mm <= 0) {
    stop("Figure dimensions must be positive millimeters.")
  }
  if (is.null(figdir))
    figdir <- file.path(.manuscript_export_state$project, "figures")
  destination <- file.path(figdir, subdir)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(pattern = paste0(".", stem, "_export_"), tmpdir = destination)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  source_pdf <- file.path(staging, "editable.pdf")
  source_svg <- file.path(staging, "editable.svg")
  outlined_pdf <- file.path(staging, "outlined.pdf")
  ggplot2::ggsave(source_svg, plot = p, device = svglite::svglite, width = width_mm, height = height_mm,
    units = "mm", bg = "white")
  svg_text <- paste(readLines(source_svg, warn = FALSE), collapse = "\n")
  # Let the shaping engine use natural font advances, without stretch-to-fit.
  svg_text <- gsub("\\s+(textLength|lengthAdjust)=[\"'][^\"']*[\"']", "", svg_text, perl = TRUE)
  embedded_style <- paste0("\n<style type=\"text/css\"><![CDATA[\n", .manuscript_export_state$font_face_css,
    "]]></style>\n")
  svg_text <- sub("(<svg\\b[^>]*>)", paste0("\\1", embedded_style), svg_text, perl = TRUE)
  writeLines(svg_text, source_svg, useBytes = TRUE)
  .manuscript_run(.manuscript_export_state$rsvg, c("--format=pdf", "--output", source_pdf, source_svg),
    "Natural SVG text shaping and editable PDF export")
  .manuscript_run(.manuscript_export_state$gs, .manuscript_gs_pdf_args(outlined_pdf, source_pdf), "PDF text outlining")
  .manuscript_promote(outlined_pdf, file.path(destination, paste0(stem, ".pdf")))
  invisible(file.path(destination, paste0(stem, ".pdf")))
}
