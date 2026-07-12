## ============================================================
## config.R -- run configuration, output paths, schema version
##
## Everything the simulation depends on is collected in a single
## config object produced by make_config().
## ============================================================

## Bump this whenever the raw-results schema changes. Loaders refuse
## files whose stored schema does not match (see summarize_results.R),
## so stale result files cannot silently feed new tables.
SCHEMA_VERSION <- "v6"

## Canonical design and method orderings 
DESIGNS      <- c("lee", "lm", "pl_smooth", "pl_local", "gap",
                  "var_jump", "var_flat")
MAIN_DESIGNS <- c("lee", "lm", "pl_smooth", "pl_local", "gap")
HET_DESIGNS  <- c("var_jump", "var_flat")

DLAB <- c(lee = "Lee", lm = "Ludwig--Miller", pl_smooth = "PL, smooth",
          pl_local = "PL, local feature", gap = "Curvature gap",
          var_jump = "Variance jump", var_flat = "Variance flat")
DLAB_SHORT <- c(lee = "Lee", lm = "Ludwig--Miller", pl_smooth = "PL, smooth",
                pl_local = "PL, local", gap = "Curvature gap",
                var_jump = "Variance jump", var_flat = "Variance flat")

## LaTeX labels; unknown methods (e.g. extra PLRD windows) fall back
## to their plain name via mlab().
MLAB_BASE <- c("LL (h_MSE)"  = "LL ($h_{\\mathrm{MSE}}$)",
               "RBC (h_MSE)" = "\\RBC{} ($h_{\\mathrm{MSE}}$)",
               "RBC (h_CER)" = "\\RBC{} ($h_{\\mathrm{CER}}$)",
               "PLRD"        = "\\PLRD{}",
               "AK honest"   = "AK honest")
mlab <- function(m) ifelse(m %in% names(MLAB_BASE), MLAB_BASE[m], m)

## ------------------------------------------------------------
## make_config(): 
## ------------------------------------------------------------
## mode          "production" or "smoke". Smoke forces small REPS by
##               default and redirects ALL output below out_root/smoke.
## reps          Monte Carlo replications per (dist, design, n) cell.
## ns_default    sample-size grid for the Beta(2,4) block, e.g. c(500)
##               or c(250, 500, 1000, 2000)
## ns_override   named list, per-design n-grid override for the Beta
##               block only (mirrors NS_OVERRIDE in v5), e.g.
##               list(pl_smooth = c(250, 500, 1000, 2000),
##                    pl_local  = c(250, 500, 1000, 2000),
##                    gap       = c(250, 500, 1000, 2000)).
## unif_ns       n-grid for the uniform-X appendix block.
## run_uniform   include the uniform-X block?
## run_ak        include the RDHonest arm (data-driven ROT bound)?
## designs       which designs to run (default: all seven).
## plrd_windows  named list of PLRD window specifications. Names are
##               the method labels; values are
##                 NULL               package default (full data range),
##                 a positive number  passed as max.window, or
##                 function(dat, ctx) returning a positive number,
##                                    where ctx holds the realized
##                                    rdrobust bandwidths (h_l, h_r,
##                                    b_l, b_r) of the mserd fit.
##               The first entry is the main specification. Example
##               shrinking-window alternatives:
##                 list("PLRD"          = NULL,
##                      "PLRD (w=0.5)"  = 0.5,
##                      "PLRD (w=2h)"   = function(dat, ctx) 2 * ctx["h_l"])
## master_seed   master seed; per-task seeds are derived from it and
##               from the task identifiers (see simulation_engine.R),
##               so results are independent of the worker count.
## n_cores       NULL = auto (min(max_cores, physical cores - 1)).
## tol           numerical tolerances used by the per-replication
##               validation:
##                 rec       weight-reconstruction / decomposition
##                           identities (relative to max(1, |est|)),
##                 balance   exact balance identities (absolute),
##                 plrd_xmom soft threshold for PLRD's pooled x-moment
##                           (only reported, discretization-limited).
## out_root      root of the output tree.
make_config <- function(mode        = c("production", "smoke"),
                        reps        = NULL,
                        ns_default  = c(500),
                        ns_override = list(),
                        unif_ns     = c(500),
                        run_uniform = TRUE,
                        run_ak      = TRUE,
                        designs     = DESIGNS,
                        plrd_windows = list("PLRD" = NULL),
                        master_seed = 42,
                        max_cores   = 30L,
                        n_cores     = NULL,
                        tol = list(rec = 1e-8, balance = 1e-8,
                                   plrd_xmom = 5e-3),
                        out_root    = "output",
                        fig = list(design = "pl_smooth", n = 500,
                                   seed = 42, xlim = c(-1, 1))) {
  mode <- match.arg(mode)
  if (is.null(reps)) reps <- if (mode == "smoke") 20L else 10000L

  stopifnot(all(designs %in% DESIGNS),
            length(plrd_windows) >= 1,
            !is.null(names(plrd_windows)),
            all(nzchar(names(plrd_windows))),
            reps >= 1)

  methods <- c("LL (h_MSE)", "RBC (h_MSE)", "RBC (h_CER)",
               names(plrd_windows),
               if (run_ak) "AK honest")

  base <- file.path(out_root, mode)
  tag  <- if (mode == "smoke") "_smoke" else ""
  paths <- list(
    base      = base,
    raw       = file.path(base, "raw"),
    processed = file.path(base, "processed"),
    tables    = file.path(base, "tables"),
    figures   = file.path(base, "figures"),
    logs      = file.path(base, "logs")
  )
  files <- list(
    results    = file.path(paths$raw, sprintf("mc_results_%s%s.rds",
                                              SCHEMA_VERSION, tag)),
    checkpoint = file.path(paths$raw, sprintf("mc_checkpoint_%s%s.rds",
                                              SCHEMA_VERSION, tag)),
    metadata   = file.path(paths$logs, sprintf("run_metadata%s.rds", tag)),
    metadata_txt = file.path(paths$logs, sprintf("run_metadata%s.txt", tag)),
    failures   = file.path(paths$logs, sprintf("failures%s.csv", tag)),
    warnings   = file.path(paths$logs, sprintf("warnings_summary%s.csv", tag)),
    validation = file.path(paths$logs, sprintf("validation_report%s.txt", tag)),
    summary    = file.path(paths$processed, sprintf("mc_summary_%s%s.csv",
                                                    SCHEMA_VERSION, tag)),
    bound_diag = file.path(paths$processed,
                           sprintf("plrd_bound_diag_%s%s.csv",
                                   SCHEMA_VERSION, tag)),
    cond_cov   = file.path(paths$processed,
                           sprintf("plrd_conditional_coverage_%s%s.csv",
                                   SCHEMA_VERSION, tag)),
    pop_proj   = file.path(paths$processed, sprintf("pop_projections%s.txt", tag))
  )

  cfg <- list(mode = mode, schema = SCHEMA_VERSION,
              reps = as.integer(reps),
              ns_default = as.integer(ns_default),
              ns_override = ns_override,
              unif_ns = as.integer(unif_ns),
              run_uniform = isTRUE(run_uniform),
              run_ak = isTRUE(run_ak),
              designs = designs,
              methods = methods,
              plrd_windows = plrd_windows,
              master_seed = as.integer(master_seed),
              max_cores = as.integer(max_cores),
              n_cores = n_cores,
              tol = tol,
              out_root = out_root,
              paths = paths, files = files, fig = fig)
  class(cfg) <- "mc_config"
  cfg
}

## Create the output tree for this config.
ensure_output_dirs <- function(cfg) {
  for (p in cfg$paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(cfg)
}

## Marker inserted at the top of every smoke-generated table.
smoke_banner <- function(cfg) {
  if (cfg$mode != "smoke") return(character(0))
  c("% ------------------------------------------------------------",
    sprintf("%% SMOKE OUTPUT (%d reps) -- pipeline test only, NOT results.",
            cfg$reps),
    "% Regenerate from a production run before using in the thesis.",
    "% ------------------------------------------------------------")
}

## Number of workers actually used.
resolve_cores <- function(cfg) {
  if (!is.null(cfg$n_cores)) return(max(1L, as.integer(cfg$n_cores)))
  nc <- parallel::detectCores(logical = FALSE)
  if (is.na(nc)) nc <- parallel::detectCores(logical = TRUE)
  if (is.na(nc)) nc <- 1L
  min(cfg$max_cores, max(1L, nc - 1L))
}
