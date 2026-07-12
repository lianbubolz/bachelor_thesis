## ============================================================
## latex_tables.R -- LaTeX table writers 
## ============================================================

fmt_num <- function(x, d = 3) {
  ifelse(!is.finite(x), "--", sprintf(paste0("%.", d, "f"), x))
}
fmt_pct <- function(x) ifelse(!is.finite(x), "--", sprintf("%.1f", 100 * x))
fmt_cnt <- function(x) ifelse(!is.finite(x), "--", sprintf("%d", as.integer(x)))

## provenance + smoke marker for every table
table_header <- function(cfg) {
  c(smoke_banner(cfg),
    sprintf("%% generated %s from %s (schema %s)",
            format(Sys.time(), tz = "UTC", usetz = TRUE),
            basename(cfg$files$results), SCHEMA_VERSION))
}

write_table <- function(lines, file, cfg) {
  writeLines(c(table_header(cfg), lines), file)
  cat("wrote ", file, "\n", sep = "")
}

## summary-cell accessor
make_cell <- function(summ) {
  function(dist_use, dg, n_use, m, col) {
    r <- summ[summ$dist == dist_use & summ$design == dg &
                summ$n == n_use & summ$method == m, col]
    if (length(r) == 0) NA_real_ else r[1]
  }
}

## ---------------- EC/IL table at one n (tab_main, tab_unif) --------
write_ecil <- function(summ, designs, methods, dist_use, n_use,
                       file, cfg) {
  if (!any(summ$dist == dist_use & summ$n == n_use)) {
    message("skip ", file, " (no rows for dist = ", dist_use,
            ", n = ", n_use, ")")
    return(invisible(NULL))
  }
  cell <- make_cell(summ)
  k <- length(methods)
  heads <- paste(sprintf("\\multicolumn{2}{c}{%s}", mlab(methods)),
                 collapse = " & ")
  cmid  <- paste(sprintf("\\cmidrule(lr){%d-%d}",
                         seq(2, by = 2, length.out = k),
                         seq(3, by = 2, length.out = k)), collapse = "")
  sub   <- paste(c("", rep(c("EC", "IL"), k)), collapse = " & ")
  body <- character(0)
  for (dg in designs) {
    cells <- character(0)
    for (m in methods) {
      cells <- c(cells,
                 fmt_pct(cell(dist_use, dg, n_use, m, "cover")),
                 fmt_num(cell(dist_use, dg, n_use, m, "med_len"), 3))
    }
    body <- c(body, sprintf("%s \\\\",
                            paste(c(DLAB[dg], cells), collapse = " & ")))
  }
  out <- c(sprintf("\\begin{tabular}{l %s}",
                   paste(rep("cc", k), collapse = " ")),
           "\\toprule",
           sprintf(" & %s \\\\", heads),
           cmid,
           sprintf("%s \\\\", sub),
           "\\midrule",
           body,
           "\\bottomrule", "\\end{tabular}")
  write_table(out, file, cfg)
}

## ---------------- all designs x all n x all arms --------------------
write_alln <- function(summ, methods, file, cfg) {
  cell <- make_cell(summ)
  ns_beta <- sort(unique(summ$n[summ$dist == "beta"]))
  k <- length(methods)
  heads <- paste(sprintf("\\multicolumn{2}{c}{%s}", mlab(methods)),
                 collapse = " & ")
  cmid  <- paste(sprintf("\\cmidrule(lr){%d-%d}",
                         seq(3, by = 2, length.out = k),
                         seq(4, by = 2, length.out = k)), collapse = "")
  sub   <- paste(c("", "$n$", rep(c("EC", "IL"), k)), collapse = " & ")
  body <- character(0)
  first <- TRUE
  for (dg in intersect(DESIGNS, unique(summ$design))) {
    if (!first) body <- c(body, "\\addlinespace")
    first <- FALSE
    for (n_use in ns_beta) {
      if (!any(summ$dist == "beta" & summ$design == dg &
                 summ$n == n_use)) next
      cells <- character(0)
      for (m in methods) {
        cells <- c(cells,
                   fmt_pct(cell("beta", dg, n_use, m, "cover")),
                   fmt_num(cell("beta", dg, n_use, m, "med_len"), 3))
      }
      body <- c(body, sprintf("%s \\\\",
                              paste(c(DLAB[dg], as.integer(n_use), cells),
                                    collapse = " & ")))
    }
  }
  out <- c(sprintf("\\begin{tabular}{l r %s}",
                   paste(rep("cc", k), collapse = " ")),
           "\\toprule",
           sprintf(" & & %s \\\\", heads),
           cmid,
           sprintf("%s \\\\", sub),
           "\\midrule",
           body,
           "\\bottomrule", "\\end{tabular}")
  write_table(out, file, cfg)
}

## ---------------- heteroskedasticity pair with IL ratio -------------
write_het <- function(summ, methods, file, cfg) {
  cell <- make_cell(summ)
  ns_beta <- sort(unique(summ$n[summ$dist == "beta" &
                                  summ$design %in% HET_DESIGNS]))
  if (!length(ns_beta)) {
    message("skip ", file, " (no variance-pair rows)")
    return(invisible(NULL))
  }
  k <- length(methods)
  heads <- paste(sprintf("\\multicolumn{2}{c}{%s}", mlab(methods)),
                 collapse = " & ")
  cmid  <- paste(sprintf("\\cmidrule(lr){%d-%d}",
                         seq(3, by = 2, length.out = k),
                         seq(4, by = 2, length.out = k)), collapse = "")
  sub   <- paste(c("", "$n$", rep(c("EC", "IL"), k)), collapse = " & ")
  body <- character(0)
  first <- TRUE
  for (n_use in ns_beta) {
    if (!first) body <- c(body, "\\addlinespace")
    first <- FALSE
    for (dg in HET_DESIGNS) {
      cells <- character(0)
      for (m in methods) {
        cells <- c(cells,
                   fmt_pct(cell("beta", dg, n_use, m, "cover")),
                   fmt_num(cell("beta", dg, n_use, m, "med_len"), 3))
      }
      body <- c(body, sprintf("%s \\\\",
                              paste(c(DLAB[dg], as.integer(n_use), cells),
                                    collapse = " & ")))
    }
    ratio_cells <- sapply(methods, function(m) {
      rj <- cell("beta", "var_jump", n_use, m, "med_len")
      rf <- cell("beta", "var_flat", n_use, m, "med_len")
      r  <- if (is.finite(rj) && is.finite(rf) && rf > 0) rj / rf else NA_real_
      sprintf("\\multicolumn{2}{c}{%s}", fmt_num(r, 2))
    })
    body <- c(body, sprintf("%s \\\\",
                            paste(c("IL ratio (jump/flat)",
                                    as.integer(n_use), ratio_cells),
                                  collapse = " & ")))
  }
  out <- c(sprintf("\\begin{tabular}{l r %s}",
                   paste(rep("cc", k), collapse = " ")),
           "\\toprule",
           sprintf(" & & %s \\\\", heads),
           cmid,
           sprintf("%s \\\\", sub),
           "\\midrule",
           body,
           "\\bottomrule", "\\end{tabular}")
  write_table(out, file, cfg)
}

## ---------------- PLRD bound vs realized conditional bias -----------
## median |exact conditional bias|, quantiles of their per-replication
## ratio, the bound-exceedance rate, and coverage.
write_bound <- function(summ, designs, file, cfg, method = "PLRD") {
  cell <- make_cell(summ)
  ns_beta <- sort(unique(summ$n[summ$dist == "beta"]))
  body <- character(0)
  for (dg in designs) {
    for (n_use in ns_beta) {
      if (!any(summ$dist == "beta" & summ$design == dg &
                 summ$n == n_use & summ$method == method)) next
      body <- c(body, sprintf(
        "%s & %d & %s & %s & %s & %s & %s & %s \\\\",
        DLAB_SHORT[dg], as.integer(n_use),
        fmt_num(cell("beta", dg, n_use, method, "med_bhat"), 3),
        fmt_num(cell("beta", dg, n_use, method, "mean_abs_cond_bias"), 3),
        fmt_num(cell("beta", dg, n_use, method, "ratio_q50"), 2),
        fmt_num(cell("beta", dg, n_use, method, "ratio_q90"), 2),
        fmt_pct(cell("beta", dg, n_use, method, "exceed_rate")),
        fmt_pct(cell("beta", dg, n_use, method, "cover"))))
    }
  }
  out <- c("\\begin{tabular}{lrcccccc}", "\\toprule",
           paste("design & $n$ & med.\\ $\\hat{b}$ &",
                 "mean $\\abs{\\widehat{\\mathrm{bias}}_n}$ &",
                 "$q_{50}$ ratio & $q_{90}$ ratio &",
                 "exceed (\\%) & EC (\\%) \\\\"),
           "\\midrule", body, "\\bottomrule", "\\end{tabular}",
           "% ratio = |conditional bias| / reported bound, per replication;",
           "% exceed = share of replications with |conditional bias| > bound.")
  write_table(out, file, cfg)
}

## ---------------- safeguard diagnostics with counts ------------------
write_safeguard <- function(summ, ccov, file, cfg, method = "PLRD") {
  cell <- make_cell(summ)
  ns_beta <- sort(unique(summ$n[summ$dist == "beta"]))
  get_cc <- function(dg, n_use, flag) {
    r <- ccov[ccov$dist == "beta" & ccov$design == dg &
                ccov$n == n_use & ccov$method == method &
                ccov$pretest == flag, , drop = FALSE]
    if (nrow(r) == 0) return(c(NA_real_, NA_real_))
    c(r$cov[1], r$cnt[1])
  }
  fmt_cc <- function(v) {
    if (!is.finite(v[1])) return("--")
    sprintf("%s {\\scriptsize(%d)}", fmt_pct(v[1]), as.integer(v[2]))
  }
  body <- character(0)
  first <- TRUE
  for (dg in intersect(DESIGNS, unique(summ$design))) {
    if (!first) body <- c(body, "\\addlinespace")
    first <- FALSE
    for (n_use in ns_beta) {
      if (!any(summ$dist == "beta" & summ$design == dg &
                 summ$n == n_use & summ$method == method)) next
      body <- c(body, sprintf(
        "%s & %d & %s & %s & %s & %s & %s \\\\",
        DLAB[dg], as.integer(n_use),
        fmt_pct(cell("beta", dg, n_use, method, "pretest_rate")),
        fmt_cc(get_cc(dg, n_use, 1)),
        fmt_cc(get_cc(dg, n_use, 0)),
        fmt_num(cell("beta", dg, n_use, "RBC (h_MSE)", "mean_h"), 3),
        fmt_num(cell("beta", dg, n_use, "RBC (h_MSE)", "mean_rho"), 2)))
    }
  }
  out <- c("\\begin{tabular}{lrccccc}", "\\toprule",
           paste("design & $n$ & reject (\\%) & EC, refit & EC, no refit &",
                 "mean $\\hat{h}_{\\mathrm{MSE}}$ & mean $\\hat{\\rho}$ \\\\"),
           "\\midrule", body, "\\bottomrule", "\\end{tabular}",
           "% cell counts of the conditional-coverage estimates in parentheses.")
  write_table(out, file, cfg)
}

## ---------------- successful fits and failure rates ------------------
write_fits <- function(summ, methods, file, cfg) {
  cell <- make_cell(summ)
  ns_beta <- sort(unique(summ$n[summ$dist == "beta"]))
  k <- length(methods)
  heads <- paste(sprintf("\\multicolumn{2}{c}{%s}", mlab(methods)),
                 collapse = " & ")
  cmid  <- paste(sprintf("\\cmidrule(lr){%d-%d}",
                         seq(3, by = 2, length.out = k),
                         seq(4, by = 2, length.out = k)), collapse = "")
  sub   <- paste(c("", "$n$", rep(c("$R_{\\mathrm{ok}}$", "fail \\%"), k)),
                 collapse = " & ")
  body <- character(0)
  first <- TRUE
  for (dg in intersect(DESIGNS, unique(summ$design))) {
    if (!first) body <- c(body, "\\addlinespace")
    first <- FALSE
    for (n_use in ns_beta) {
      if (!any(summ$dist == "beta" & summ$design == dg &
                 summ$n == n_use)) next
      cells <- character(0)
      for (m in methods) {
        cells <- c(cells,
                   fmt_cnt(cell("beta", dg, n_use, m, "n_success")),
                   fmt_num(100 * cell("beta", dg, n_use, m, "fail_rate"), 2))
      }
      body <- c(body, sprintf("%s \\\\",
                              paste(c(DLAB[dg], as.integer(n_use), cells),
                                    collapse = " & ")))
    }
  }
  out <- c(sprintf("\\begin{tabular}{l r %s}",
                   paste(rep("cc", k), collapse = " ")),
           "\\toprule",
           sprintf(" & & %s \\\\", heads),
           cmid,
           sprintf("%s \\\\", sub),
           "\\midrule", body,
           "\\bottomrule", "\\end{tabular}",
           "% coverage and interval-length summaries use the successful",
           "% fits; failed fits are logged in logs/failures*.csv.")
  write_table(out, file, cfg)
}

## ---------------- write everything ----------------------------------
write_all_tables <- function(res, cfg) {
  summ <- summarize_results(res)
  ccov <- plrd_conditional_coverage(res)
  tdir <- cfg$paths$tables
  main_methods <- cfg$methods

  write_ecil(summ, intersect(MAIN_DESIGNS, cfg$designs), main_methods,
             "beta", 500, file.path(tdir, "tab_main.tex"), cfg)
  if (any(summ$dist == "unif")) {
    write_ecil(summ, intersect(DESIGNS, cfg$designs), main_methods,
               "unif", 500, file.path(tdir, "tab_unif.tex"), cfg)
  } else message("skip tab_unif.tex (no uniform block in results)")

  write_alln(summ, main_methods, file.path(tdir, "tab_alln.tex"), cfg)
  if (all(HET_DESIGNS %in% cfg$designs)) {
    write_het(summ, main_methods, file.path(tdir, "tab_het.tex"), cfg)
    write_bound(summ, HET_DESIGNS,
                file.path(tdir, "tab_bound_het.tex"), cfg)
  }
  write_bound(summ, intersect(MAIN_DESIGNS, cfg$designs),
              file.path(tdir, "tab_bound.tex"), cfg)
  write_safeguard(summ, ccov, file.path(tdir, "tab_safeguard.tex"), cfg)
  write_fits(summ, main_methods, file.path(tdir, "tab_fits.tex"), cfg)

  invisible(summ)
}

## ---------------- console diagnostics  ---------------------
print_console_diagnostics <- function(res, summ) {
  cat("\n--- fits per cell (successes / total) ---\n")
  agg <- summ[, c("dist", "design", "n", "method", "n_success",
                  "n_total", "fail_rate")]
  print(utils::head(agg[order(-agg$fail_rate), ], 15), row.names = FALSE)

  cat("\n--- safeguard rejection and conditional coverage (PLRD, Beta) ---\n")
  ccov <- plrd_conditional_coverage(res)
  cb <- ccov[ccov$dist == "beta", , drop = FALSE]
  if (nrow(cb)) print(cb, digits = 3, row.names = FALSE)

  cat("\n--- PLRD bound diagnostics (Beta) ---\n")
  pb <- summ[summ$dist == "beta" & grepl("^PLRD", summ$method),
             c("design", "n", "method", "med_bhat", "mean_abs_cond_bias",
               "ratio_q50", "ratio_q90", "exceed_rate", "pretest_rate")]
  if (nrow(pb)) print(pb[order(match(pb$design, DESIGNS), pb$n), ],
                      digits = 3, row.names = FALSE)

  cat("\n--- mean realized h_MSE and rho (rdrobust, Beta) ---\n")
  hb <- summ[summ$dist == "beta" & summ$method == "RBC (h_MSE)",
             c("design", "n", "mean_h", "mean_rho")]
  if (nrow(hb)) print(hb[order(match(hb$design, DESIGNS), hb$n), ],
                      digits = 3, row.names = FALSE)
  invisible(NULL)
}
