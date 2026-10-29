# Shared path initialization for the patient analysis scripts.
# run.py passes the module path in an environment variable so that paths with
# spaces work without parsing R's encoded --file command-line argument.
module_root <- Sys.getenv("E6_PATIENT_MODULE_ROOT")
if (!nzchar(module_root)) {
  stop("Run this module through patient_scores/run.py.")
}
module_root <- normalizePath(module_root, mustWork = TRUE)

# Each stage receives the same output directory. Intermediate CSV tables are
# read from and written to this directory, leaving the bundled inputs unchanged.
args <- commandArgs(TRUE)
if (length(args) != 1) {
  stop("The runner must supply exactly one output directory.")
}
output_dir <- args[1]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
