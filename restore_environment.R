#!/usr/bin/env Rscript
# Restore the R packages needed by the analysis and plotting scripts.
# Usage: Rscript --vanilla restore_environment.R /path/to/this/repository
# The resulting .R_library_path lets every Python entry point find that library.

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
if (as.character(getRversion()) != "4.4.3") {
    stop("Use R 4.4.3 for the recorded environment.")
}

# renv.lock fixes both package versions and their CRAN/Bioconductor sources.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}
renv::restore(project = root, lockfile = file.path(root, "renv.lock"), prompt = FALSE)
writeLines(renv::paths$library(project = root), file.path(root, ".R_library_path"))
message("R packages restored; wrote the local library path for the figure scripts.")
