## ============================================================
## scripts/make_outputs.R -- tables/figures from saved results
##
## Generates every summary, table and figure from the stored
## replication-level file without rerunning the Monte Carlo
##
## Usage
##   source("scripts/make_outputs.R")            # production results
##   MODE <- "smoke"; source("scripts/make_outputs.R")
## ============================================================

local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  cand <- c(file.path("scripts", "_common.R"), "_common.R",
            if (length(fa)) file.path(dirname(fa[1]), "_common.R"))
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("cannot find scripts/_common.R; ",
                       "run from the project root")
  source(hit)
})
ROOT <- locate_root()
SRC  <- source_project(ROOT)

if (!exists("MODE"))         MODE         <- "production"
if (!exists("RUN_TABLES"))   RUN_TABLES   <- TRUE
if (!exists("RUN_FIGURE"))   RUN_FIGURE   <- TRUE
if (!exists("RUN_CHECKS"))   RUN_CHECKS   <- TRUE
if (!exists("PLRD_WINDOWS")) PLRD_WINDOWS <- list("PLRD" = NULL)

cfg <- make_config(mode = MODE, plrd_windows = PLRD_WINDOWS,
                   out_root = file.path(ROOT, "output"))
if (file.exists(cfg$files$metadata)) {
  meta <- readRDS(cfg$files$metadata)
  cfg$reps    <- meta$config$reps
  cfg$methods <- meta$config$methods
  cfg$designs <- meta$config$designs
}

res <- load_results(cfg$files$results)
produce_outputs(cfg, res,
                run_tables = RUN_TABLES,
                run_figure = RUN_FIGURE,
                run_checks = RUN_CHECKS)
