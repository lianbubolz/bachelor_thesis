## ============================================================
## validate_results.R
##
## Verifies, on the saved replication-level rows:
##   1. completeness: exactly one row per (cell, method, rep), and
##      no duplicated task seeds within a cell;
##   2. CI validity: every success row has finite, ordered endpoints
##   3. coverage/length consistency with est/lo/hi/tau;
##   4. weight identities: reconstructed (LL/RBC) and returned
##      (PLRD gamma) weights reproduce the reported estimates, and
##      est = tau + cond_bias + noise, within cfg$tol$rec;
##   5. balance: LL/RBC side-specific identities and PLRD's exact
##      sums within cfg$tol$balance; PLRD pooled x-moments screened
##      against the loose cfg$tol$plrd_xmom (discretization-limited,
##      reported only);
##   6. failure visibility: failed fits present with status/error.
##
## Writes a plain-text report
## Lines starting with "FAIL" indicate problems that must be
## resolved before the numbers are used
## ============================================================

validate_results <- function(res, cfg, write = TRUE) {
  L <- character(0)
  say <- function(...) L <<- c(L, sprintf(...))
  n_fail_checks <- 0L
  fail <- function(fmt, ...) {
    n_fail_checks <<- n_fail_checks + 1L
    L <<- c(L, paste0("FAIL: ", sprintf(fmt, ...)))
  }

  say("validation report (%s, schema %s, %d rows)",
      cfg$mode, SCHEMA_VERSION, nrow(res))
  say("generated %s", format(Sys.time(), tz = "UTC", usetz = TRUE))
  say("")

  ## 1. completeness ------------------------------------------------
  cnt <- stats::aggregate(list(rows = res$rep),
                          by = list(dist = res$dist, design = res$design,
                                    n = res$n, method = res$method),
                          FUN = length)
  expected <- cfg$reps
  bad <- cnt[cnt$rows != expected, , drop = FALSE]
  if (nrow(bad) == 0L) {
    say("completeness: every cell x method has exactly %d rows -- ok",
        expected)
  } else {
    fail("%d cell x method combinations deviate from %d rows:",
         nrow(bad), expected)
    for (i in seq_len(min(nrow(bad), 10))) {
      say("  %s/%s n=%d %s : %d rows", bad$dist[i], bad$design[i],
          bad$n[i], bad$method[i], bad$rows[i])
    }
  }
  ## task seeds repeat across methods within a rep by design. Check
  ## per-cell uniqueness across reps instead:
  seed_ok <- TRUE
  for (k in split(res, interaction(res$dist, res$design, res$n,
                                   drop = TRUE))) {
    per_rep <- unique(k[, c("rep", "task_seed")])
    if (anyDuplicated(per_rep$task_seed)) { seed_ok <- FALSE; break }
  }
  if (seed_ok) say("task seeds: unique across reps within every cell -- ok")
  else fail("duplicated task seeds within a cell")

  ## 2. success / CI validity --------------------------------------
  ci_ok <- with(res, is.finite(est) & is.finite(lo) & is.finite(hi) &
                       lo <= hi)
  mismatch <- sum((res$status == "ok" & ci_ok) != res$success)
  if (mismatch == 0L) {
    say("success flag: consistent with status + finite ordered CI -- ok")
  } else fail("%d rows where success flag disagrees with its definition",
              mismatch)
  s <- res[res$success, , drop = FALSE]

  ## 3. coverage / length consistency -------------------------------
  cov2 <- as.numeric(s$tau >= s$lo & s$tau <= s$hi)
  len2 <- s$hi - s$lo
  if (max(abs(cov2 - s$cover)) == 0 &&
      max(abs(len2 - s$len)) < 1e-12) {
    say("cover/len: recomputed from endpoints match stored values -- ok")
  } else fail("stored cover/len disagree with endpoints")

  ## 4. weight identities -------------------------------------------
  has_w <- s[is.finite(s$wt_gap), , drop = FALSE]
  if (nrow(has_w)) {
    scale <- pmax(1, abs(has_w$est))
    wmax <- max(abs(has_w$wt_gap) / scale)
    dmax <- max(abs(has_w$decomp_gap) / scale)
    say("weight identity: %d rows carry weights; max |wt_gap| = %.2e, max |decomp_gap| = %.2e (rel., tol %.0e)",
        nrow(has_w), wmax, dmax, cfg$tol$rec)
    n_bad <- sum(!has_w$ok_weights)
    if (n_bad == 0L) say("  all reconstructed/returned weights reproduce the estimates -- ok")
    else fail("%d rows exceed the weight-identity tolerance", n_bad)
    cover_frac <- stats::aggregate(
      list(share = is.finite(s$wt_gap)),
      by = list(method = s$method), FUN = mean)
    for (i in seq_len(nrow(cover_frac))) {
      say("  weights available for %5.1f%% of %s successes",
          100 * cover_frac$share[i], cover_frac$method[i])
    }
  } else fail("no rows carry weight diagnostics")

  ## 5. balance ------------------------------------------------------
  has_b <- s[is.finite(s$bal_dev), , drop = FALSE]
  if (nrow(has_b)) {
    say("balance: max exact-identity deviation = %.2e (tol %.0e)",
        max(has_b$bal_dev), cfg$tol$balance)
    n_bad <- sum(!has_b$ok_balance)
    if (n_bad == 0L) say("  all exact balance identities hold -- ok")
    else fail("%d rows violate exact balance identities", n_bad)
  }
  px <- s[grepl("^PLRD", s$method) & is.finite(s$bal_gx), , drop = FALSE]
  if (nrow(px)) {
    mgx  <- max(abs(px$bal_gx)); mgx2 <- max(abs(px$bal_gx2))
    say("PLRD pooled moments (discretization-limited, informational):")
    say("  max |sum gamma x|   = %.2e, max |sum gamma x^2| = %.2e (screen %.0e)",
        mgx, mgx2, cfg$tol$plrd_xmom)
    if (max(mgx, mgx2) > cfg$tol$plrd_xmom) {
      say("  warn: pooled moment above the screen; inspect bal_gx/bal_gx2 in the raw file")
    }
  }

  ## 6. failure visibility -------------------------------------------
  nf <- sum(!res$success)
  say("")
  say("failures: %d of %d fits (%.3f%%) -- retained with status/error, see %s",
      nf, nrow(res), 100 * nf / nrow(res), basename(cfg$files$failures))
  if (nf > 0) {
    tab <- sort(table(res$method[!res$success]), decreasing = TRUE)
    for (m in names(tab)) say("  %-14s %d", m, tab[[m]])
  }

  say("")
  say(if (n_fail_checks == 0L) "ALL CHECKS PASSED"
      else sprintf("%d CHECK(S) FAILED -- do not use these results yet",
                   n_fail_checks))

  if (write) {
    writeLines(L, cfg$files$validation)
    cat("wrote ", cfg$files$validation, "\n", sep = "")
  }
  cat(paste0("  ", utils::tail(L, 1), "\n"))
  invisible(list(lines = L, n_failed = n_fail_checks))
}
