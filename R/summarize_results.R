## ============================================================
## summarize_results.R
## ============================================================

REQUIRED_COLS <- c("dist", "design", "n", "rep", "task_seed", "method",
                   "status", "error", "n_warn", "warn_text", "n_msg",
                   "success", "est", "lo", "hi", "se", "tau", "cover",
                   "len", "h_l", "h_r", "b_l", "b_r", "rho", "pretest",
                   "bhat", "Bhat", "plrd_window", "cond_bias", "noise",
                   "wt_gap", "decomp_gap", "bal_dev", "bal_gx",
                   "bal_gx2", "ok_weights", "ok_balance", "bound_ratio",
                   "bound_exceed", "bound_exceed_z")

load_results <- function(path) {
  if (!file.exists(path)) {
    stop("results file not found: ", path,
         "\nRun scripts/run_mc.R (or run_smoke.R) first.")
  }
  res <- readRDS(path)
  sv <- attr(res, "schema_version")
  if (is.null(sv) || !identical(sv, SCHEMA_VERSION)) {
    stop("results file '", path, "' has schema '",
         if (is.null(sv)) "<none>" else sv, "' but this code expects '",
         SCHEMA_VERSION, "'; rerun the Monte Carlo.")
  }
  miss <- setdiff(REQUIRED_COLS, names(res))
  if (length(miss)) {
    stop("results file is missing columns: ", paste(miss, collapse = ", "))
  }
  res
}

## Quantile helper that tolerates all-NA input.
qna <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, p, names = FALSE, type = 7))
}
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
med_na  <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)

## ------------------------------------------------------------
## summarize_results(): one row per dist x design x n x method
## ------------------------------------------------------------
summarize_results <- function(res) {
  key <- interaction(res$dist, res$design, res$n, res$method,
                     drop = TRUE, sep = "\r")
  parts <- split(res, key)

  out <- lapply(parts, function(d) {
    s <- d[d$success, , drop = FALSE]
    n_total <- nrow(d)
    n_succ  <- nrow(s)
    err     <- s$est - s$tau
    p_cov   <- mean_na(s$cover)
    data.frame(
      dist = d$dist[1], design = d$design[1], n = d$n[1],
      method = d$method[1],
      n_total = n_total, n_success = n_succ,
      n_fail = n_total - n_succ,
      fail_rate = (n_total - n_succ) / n_total,
      cover = p_cov,
      cover_mcse = if (n_succ > 0 && is.finite(p_cov))
        sqrt(p_cov * (1 - p_cov) / n_succ) else NA_real_,
      mc_bias = mean_na(err),
      mc_sd   = if (n_succ > 1) stats::sd(err) else NA_real_,
      rmse    = if (n_succ > 0) sqrt(mean(err^2)) else NA_real_,
      med_len = med_na(s$len),
      mean_len = mean_na(s$len),
      mean_se = mean_na(s$se),
      med_se  = med_na(s$se),
      mean_h  = mean_na(s$h_l),
      mean_rho = mean_na(s$rho),
      ## exact weight-based diagnostics (LL / RBC / PLRD)
      mean_cond_bias = mean_na(s$cond_bias),
      mean_abs_cond_bias = mean_na(abs(s$cond_bias)),
      sd_cond_bias = if (sum(is.finite(s$cond_bias)) > 1)
        stats::sd(s$cond_bias[is.finite(s$cond_bias)]) else NA_real_,
      ok_weights_rate = mean_na(as.numeric(s$ok_weights)),
      ok_balance_rate = mean_na(as.numeric(s$ok_balance)),
      warn_rate = mean(d$n_warn > 0),
      ## PLRD-specific (NA elsewhere)
      pretest_rate = mean_na(s$pretest),
      med_bhat = med_na(s$bhat),
      mean_Bhat = mean_na(s$Bhat),
      ratio_q50 = qna(s$bound_ratio, 0.50),
      ratio_q90 = qna(s$bound_ratio, 0.90),
      ratio_q99 = qna(s$bound_ratio, 0.99),
      ratio_max = if (any(is.finite(s$bound_ratio)))
        max(s$bound_ratio, na.rm = TRUE) else NA_real_,
      exceed_rate = mean_na(s$bound_exceed),
      exceed_z_q99 = qna(s$bound_exceed_z, 0.99),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$dist, match(out$design, DESIGNS), out$n,
            match(out$method, unique(res$method))), ]
}

## ------------------------------------------------------------
## PLRD conditional coverage by safeguard decision, with counts.
## plrd_methods defaults to every arm that carries a pretest flag.
## ------------------------------------------------------------
plrd_conditional_coverage <- function(res, min_reps = 10) {
  p <- res[res$success & !is.na(res$pretest), , drop = FALSE]
  if (nrow(p) == 0L) {
    return(data.frame(dist = character(0), design = character(0),
                      n = integer(0), method = character(0),
                      pretest = numeric(0), cov = numeric(0),
                      cnt = integer(0), reliable = logical(0)))
  }
  agg <- merge(
    stats::aggregate(list(cov = p$cover),
                     by = list(dist = p$dist, design = p$design,
                               n = p$n, method = p$method,
                               pretest = p$pretest), FUN = mean),
    stats::aggregate(list(cnt = p$cover),
                     by = list(dist = p$dist, design = p$design,
                               n = p$n, method = p$method,
                               pretest = p$pretest), FUN = length),
    by = c("dist", "design", "n", "method", "pretest"))
  agg$reliable <- agg$cnt >= min_reps
  agg[order(agg$dist, match(agg$design, DESIGNS), agg$n,
            agg$method, agg$pretest), ]
}

## ------------------------------------------------------------
## Replication-level PLRD bound diagnostics table (long form),
## one row per dist x design x n x PLRD arm.
## ------------------------------------------------------------
plrd_bound_summary <- function(res) {
  p <- res[grepl("^PLRD", res$method), , drop = FALSE]
  if (nrow(p) == 0L) return(NULL)
  summ <- summarize_results(p)
  keep <- c("dist", "design", "n", "method", "n_success", "n_fail",
            "med_bhat", "mean_abs_cond_bias", "ratio_q50", "ratio_q90",
            "ratio_q99", "ratio_max", "exceed_rate", "exceed_z_q99",
            "pretest_rate", "cover")
  summ[, keep]
}

## Write all processed summaries for a config.
write_summaries <- function(res, cfg) {
  summ <- summarize_results(res)
  utils::write.csv(summ, cfg$files$summary, row.names = FALSE)
  cat("wrote ", cfg$files$summary, "\n", sep = "")

  bd <- plrd_bound_summary(res)
  if (!is.null(bd)) {
    utils::write.csv(bd, cfg$files$bound_diag, row.names = FALSE)
    cat("wrote ", cfg$files$bound_diag, "\n", sep = "")
  }

  ccov <- plrd_conditional_coverage(res)
  utils::write.csv(ccov, cfg$files$cond_cov, row.names = FALSE)
  cat("wrote ", cfg$files$cond_cov, "\n", sep = "")

  invisible(list(summary = summ, bound = bd, cond_cov = ccov))
}
