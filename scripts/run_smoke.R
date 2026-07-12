## ============================================================
## scripts/run_smoke.R
##
## Exercises every stage end to end with tiny reps All output are
## written under output/smoke/ with "_smoke"-tagged filenames and
## a "% SMOKE OUTPUT" banner inside every table, so it can never
## be mistaken for (or overwrite) production results.
##
## Pre-assignable parameters, e.g.  SMOKE_REPS <- 5; source(...).
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

if (!exists("SMOKE_REPS"))    SMOKE_REPS    <- 20
if (!exists("SMOKE_DESIGNS")) SMOKE_DESIGNS <- DESIGNS   # all seven
if (!exists("RUN_AK"))        RUN_AK        <- TRUE
if (!exists("RUN_UNIFORM"))   RUN_UNIFORM   <- TRUE
if (!exists("MAX_CORES"))     MAX_CORES     <- 30L
if (!exists("PLRD_WINDOWS"))  PLRD_WINDOWS  <- list("PLRD" = NULL)

cat(sprintf("** SMOKE mode: %d reps, n = 500 only; outputs are pipeline tests, not results **\n",
            SMOKE_REPS))

cfg <- make_config(
  mode         = "smoke",
  reps         = SMOKE_REPS,
  ns_default   = c(500),
  unif_ns      = c(500),
  run_uniform  = RUN_UNIFORM,
  run_ak       = RUN_AK,
  designs      = SMOKE_DESIGNS,
  plrd_windows = PLRD_WINDOWS,
  max_cores    = MAX_CORES,
  out_root     = file.path(ROOT, "output")
)

res <- run_mc(cfg, SRC)
produce_outputs(cfg, res)

cat("\nSmoke run complete. Inspect output/smoke/logs/validation_report_smoke.txt\n",
    "and the checklist in README.md before launching production.\n", sep = "")
