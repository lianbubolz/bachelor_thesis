## ============================================================
## scripts/run_mc.R -- PRODUCTION Monte Carlo
##
## Usage
##   Rscript scripts/run_mc.R          from the project root, or
##   source("scripts/run_mc.R")        from an R session.
##
## e.g.  REPS <- 2000; source("scripts/run_mc.R").
## Output goes to output/production/ 
## ============================================================

## --- find and source scripts/_common.R ------------------
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

if (!exists("REPS"))          REPS          <- 10000
if (!exists("NS"))            NS            <- c(500)  # Beta(2,4) grid
if (!exists("UNIF_NS"))       UNIF_NS       <- c(500)
if (!exists("SEED"))          SEED          <- 42
if (!exists("RUN_AK"))        RUN_AK        <- TRUE
if (!exists("RUN_UNIFORM"))   RUN_UNIFORM   <- TRUE
if (!exists("MAX_CORES"))     MAX_CORES     <- 30L
if (!exists("NS_OVERRIDE"))   NS_OVERRIDE   <- list()
if (!exists("PLRD_WINDOWS"))  PLRD_WINDOWS  <- list("PLRD" = NULL)
if (!exists("RUN_TABLES"))    RUN_TABLES    <- TRUE
if (!exists("RUN_FIGURE"))    RUN_FIGURE    <- TRUE
if (!exists("RUN_CHECKS"))    RUN_CHECKS    <- TRUE

## --- sample-size extension ------------------------------------------
## The engine supports any grid. For the n-extension of the three
## focal designs, pre-assign (or uncomment):
##   NS_OVERRIDE <- list(pl_smooth = c(250, 500, 1000, 2000),
##                       pl_local  = c(250, 500, 1000, 2000),
##                       gap       = c(250, 500, 1000, 2000))
##
## --- PLRD window alternatives ----------------------------------------
## The package's full-range default is the main specification.
## Shrinking-window, if wanted, e.g.:
##   PLRD_WINDOWS <- list("PLRD"         = NULL,
##                        "PLRD (w=0.5)" = 0.5,
##                        "PLRD (w=2h)"  = function(dat, ctx) 2 * ctx[["h_l"]])
## ----------------------------------------------------------------------

cfg <- make_config(
  mode         = "production",
  reps         = REPS,
  ns_default   = NS,
  ns_override  = NS_OVERRIDE,
  unif_ns      = UNIF_NS,
  run_uniform  = RUN_UNIFORM,
  run_ak       = RUN_AK,
  plrd_windows = PLRD_WINDOWS,
  master_seed  = SEED,
  max_cores    = MAX_CORES,
  out_root     = file.path(ROOT, "output")
)

res <- run_mc(cfg, SRC)
produce_outputs(cfg, res,
                run_tables = RUN_TABLES,
                run_figure = RUN_FIGURE,
                run_checks = RUN_CHECKS)

cat("\nDone. Key outputs under ", cfg$paths$base, ":\n",
    "  raw/       replication-level results (+ checkpoint during runs)\n",
    "  processed/ mc_summary, plrd_bound_diag, conditional coverage,\n",
    "             pop_projections.txt\n",
    "  tables/    tab_main, tab_alln, tab_het, tab_bound(_het),\n",
    "             tab_safeguard, tab_unif, tab_fits\n",
    "  figures/   fig_weights.pdf (+ info)\n",
    "  logs/      run_metadata, failures.csv, warnings_summary.csv,\n",
    "             validation_report.txt\n", sep = "")
