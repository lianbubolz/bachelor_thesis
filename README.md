# Thesis Monte Carlo: RBC vs PLRD inference

## Layout

```
R/                      library code (no side effects on source)
  config.R              run configuration, paths, schema version
  dgp_definitions.R     the seven designs, both running-variable laws
  estimator_rdrobust.R  strict wrapper: LL + RBC via rdrobust
  estimator_plrd.R      strict wrapper: PLRD (+ window resolution)
  estimator_rdhonest.R  strict wrapper: AK honest (data-driven ROT)
  weight_reconstruction.R  exact LL/RBC weights at realized bandwidths
  diagnostics.R         conditional bias, decomposition, balance, bounds
  simulation_engine.R   task seeds, one_rep(), run_mc(), metadata
  summarize_results.R   per-cell summaries, PLRD extras
  validate_results.R    post-run checks -> validation_report.txt
  latex_tables.R        all .tex writers + console diagnostics
  figures.R             fig_weights.pdf + info file
scripts/
  run_smoke.R            end-to-end test with small REPS -> output/smoke/
  run_mc.R              production run             -> output/production/
  make_outputs.R        generate tables/figures from saved results
output/
  smoke/…               "_smoke"-tagged files, SMOKE banner in tables
  production/…          raw/ processed/ tables/ figures/ logs/
```

## Quick start

```r

## 1. smoke run: all seven designs, 20 reps, n = 500, ~seconds
Rscript scripts/run_smoke.R
##   -> work through the checklist below

## 2. production (Windows/PSOCK, potentially long runtime)
Rscript scripts/run_mc.R

## 3. tables/figures again later
Rscript scripts/make_outputs.R
```

All scripts also work via `source()` from an R session but pre-assign parameters first e.g.

```r
REPS <- 2000
NS          <- c(500, 1000)
RUN_UNIFORM <- FALSE
PLRD_WINDOWS <- list("PLRD" = NULL, "PLRD (w=0.5)" = 0.5,
                     "PLRD (w=2h)" = function(dat, ctx) 2 * ctx[["h_l"]])
source("scripts/run_mc.R")
```

When changing parameters, restart R or rm(list = ls()) the working environment before running the scripts again.

### n-grid extension

```r
NS_OVERRIDE <- list(pl_smooth = c(250, 500, 1000, 2000),
                    pl_local  = c(250, 500, 1000, 2000),
                    gap       = c(250, 500, 1000, 2000))
```

## Smoke-run verification checklist

Run `Rscript scripts/run_smoke.R`, then check, in order:

1. Console: `estimator probes ok` appeared before the loop (all three
   wrappers exercised once on synthetic data.
2. `output/smoke/logs/validation_report_smoke.txt` ends in
   `ALL CHECKS PASSED`. It verifies completeness (rows = reps for
   every cell x method), per-cell seed uniqueness, success-flag and
   cover/length consistency, weight identities (expect max relative
   gaps ~1e-15), exact balance identities (~1e-15), PLRD pooled
   moments under the 5e-3 screen.
3. `logs/failures_smoke.csv` is empty (header only) or every row is
   explained; failures also print per method in the validation report.
4. `tables/` contains tab_main, tab_alln, tab_bound, tab_safeguard,
   tab_fits (+ tab_unif, tab_het, tab_bound_het when the uniform
   block / variance pair ran), each starting with the SMOKE banner.
5. `figures/fig_weights_info.txt`: reconstruction gaps ~1e-15,
   `sum(gamma*Y) - tau.hat` ~1e-16, balance deviation ~1e-15.


Only then launch `scripts/run_mc.R`.

## Notes

- Coverage and interval-length summaries use successful fits only,
  `tab_fits.tex` reports success counts and failure rates per cell.
