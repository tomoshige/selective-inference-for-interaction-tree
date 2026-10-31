# Run from any directory; outputs are written beside this script.
script <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script)) setwd(dirname(normalizePath(sub("^--file=", "", script[1L]))))
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript reproduce.R quick|simulation|actg|all [--workers=1]"
if (!length(args) || !args[1L] %in% c("quick", "simulation", "actg", "all"))
  stop(usage, call. = FALSE)
stage <- args[1L]
workers <- 1L
for (a in args[-1L]) {
  if (!grepl("^--workers=[0-9]+$", a)) stop(usage, call. = FALSE)
  workers <- as.integer(sub("^--workers=", "", a))
}
if (is.na(workers) || workers < 1L) stop("Use at least one worker.", call. = FALSE)
if (.Platform$OS.type == "windows" && workers > 1L)
  stop("Use --workers=1 on Windows.", call. = FALSE)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
if (stage %in% c("quick", "simulation", "all")) {
  required <- c("data.table", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the R packages listed in README.md: ",
                            paste(missing, collapse = ", "), call. = FALSE)
}
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run <- function(file, ...) {
  cat("Running", file, "\n")
  status <- system2(rscript, shQuote(c("--vanilla", normalizePath(file), ...)))
  if (status != 0L) stop("Execution failed: ", file, call. = FALSE)
}
simulation <- function(quick = FALSE) {
  mode <- if (quick) "quick" else "main"
  extra <- if (quick) c("--reps=2", "--width=2") else character()
  run("simulation/run.R", paste0("--stage=", mode),
      paste0("--workers=", if (quick) 1L else workers), extra)
  run("simulation/summarize.R", "--stage", mode)
  run("simulation/figures.R", "--stage", mode)
}
actg <- function() {
  run("actg/run.R")
  run("actg/figures.R")
}
if (stage == "quick") {
  simulation(TRUE)
  actg()
} else if (stage == "simulation") {
  simulation()
} else if (stage == "actg") {
  actg()
} else {
  simulation()
  actg()
}
cat("Complete. Results are in results/ and figures/ next to reproduce.R.\n")
