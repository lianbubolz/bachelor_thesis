## ============================================================
## simulation_engine.R -- task seeding, one_rep(), run_mc()
## ============================================================

## ------------------------------------------------------------
## Reproducible task seeds, independent of the worker count
##
## Every (dist, design, n, rep) task receives a deterministic seed
## derived from those identifiers and the master seed via a
## polynomial string hash. Within a cell, seeds are consecutive
## offsets of the cell hash, hence distinct. Because each task
## seeds its own RNG, results are identical under serial and
## parallel execution and under any core count, and single tasks
## can be reproduced in isolation.
## ------------------------------------------------------------
hash_string <- function(s) {
  h <- 0
  for (cc in utf8ToInt(s)) h <- (h * 31 + cc) %% 2147483647
  h
}

task_seed <- function(master_seed, dist, design, n, rep) {
  cell <- hash_string(paste("plrd-mc", SCHEMA_VERSION, master_seed,
                            dist, design, n, sep = "|"))
  as.integer((cell + rep - 1) %% 2147483646) + 1L
}

set_task_seed <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
}

## ------------------------------------------------------------
## Result schema: one row per requested (task, method), ALWAYS.
## ------------------------------------------------------------
result_row <- function(dist, design, n, rep, seed, method,
                       status = "error", error = "",
                       n_warn = 0L, warn_text = "", n_msg = 0L,
                       est = NA_real_, lo = NA_real_, hi = NA_real_,
                       se = NA_real_, tau = NA_real_,
                       h_l = NA_real_, h_r = NA_real_,
                       b_l = NA_real_, b_r = NA_real_,
                       pretest = NA_real_, bhat = NA_real_,
                       Bhat = NA_real_, plrd_window = NA_real_,
                       cond_bias = NA_real_, noise = NA_real_,
                       wt_gap = NA_real_, decomp_gap = NA_real_,
                       bal_dev = NA_real_, bal_gx = NA_real_,
                       bal_gx2 = NA_real_,
                       ok_weights = NA, ok_balance = NA,
                       bound_ratio = NA_real_, bound_exceed = NA_real_,
                       bound_exceed_z = NA_real_) {
  ci_ok  <- is.finite(est) && is.finite(lo) && is.finite(hi) && lo <= hi
  success <- (status == "ok") && ci_ok
  if (status == "ok" && !ci_ok && !nzchar(error)) {
    error <- "invalid CI (non-finite or lo > hi)"
  }
  data.frame(
    dist = dist, design = design, n = as.integer(n),
    rep = as.integer(rep), task_seed = as.integer(seed),
    method = method, status = status, error = error,
    n_warn = as.integer(n_warn), warn_text = warn_text,
    n_msg = as.integer(n_msg), success = success,
    est = est, lo = lo, hi = hi, se = se, tau = tau,
    cover = if (success) as.numeric(tau >= lo && tau <= hi) else NA_real_,
    len   = if (success) hi - lo else NA_real_,
    h_l = h_l, h_r = h_r, b_l = b_l, b_r = b_r,
    rho = if (is.finite(h_l) && is.finite(b_l) && b_l > 0) h_l / b_l
          else NA_real_,
    pretest = pretest, bhat = bhat, Bhat = Bhat,
    plrd_window = plrd_window,
    cond_bias = cond_bias, noise = noise,
    wt_gap = wt_gap, decomp_gap = decomp_gap,
    bal_dev = bal_dev, bal_gx = bal_gx, bal_gx2 = bal_gx2,
    ok_weights = ok_weights, ok_balance = ok_balance,
    bound_ratio = bound_ratio, bound_exceed = bound_exceed,
    bound_exceed_z = bound_exceed_z,
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------
## one_rep(): a single Monte Carlo task
##
## Returns exactly length(cfg$methods) rows -- estimator failures
## produce a status = "error" row instead of vanishing, and captured
## warning/message texts travel with the row.
## ------------------------------------------------------------
one_rep <- function(rep, dist, design, n, cfg) {
  seed <- task_seed(cfg$master_seed, dist, design, n, rep)
  set_task_seed(seed)

  dgp <- get_dgp(dist, design)
  dat <- dgp$gen(n)

  tol_rec <- cfg$tol$rec
  tol_bal <- cfg$tol$balance
  rows <- vector("list", length(cfg$methods)); names(rows) <- cfg$methods
  base_args <- list(dist = dist, design = design, n = n, rep = rep,
                    seed = seed, tau = dat$tau)

  row_from <- function(method, cap, fill = list()) {
    args <- c(base_args, list(method = method,
                              status = cap$status, error = cap$error,
                              n_warn = length(cap$warns),
                              warn_text = trunc_msg(cap$warns),
                              n_msg = length(cap$msgs)),
              fill)
    do.call(result_row, args)
  }
  ## attach a reconstruction/diagnostic note to an existing row
  add_note <- function(row, note) {
    row$n_warn <- row$n_warn + 1L
    row$warn_text <- trunc_msg(c(row$warn_text[nzchar(row$warn_text)], note))
    row
  }

  ## ---- rdrobust, MSE-optimal bandwidths -> LL and RBC (h_MSE) ----
  cm  <- capture_conditions(fit_rdrobust(dat$y, dat$x, "mserd"))
  ctx <- c(h_l = NA_real_, h_r = NA_real_, b_l = NA_real_, b_r = NA_real_)

  if (cm$status == "ok") {
    v <- cm$value
    ctx[] <- c(v$h_l, v$h_r, v$b_l, v$b_r)
    wt <- capture_conditions(
      build_rd_weights(dat$x, v$h_l, v$h_r, v$b_l, v$b_r))

    fill_ll <- list(est = v$conv$est, lo = v$conv$lo, hi = v$conv$hi,
                    se = v$conv$se,
                    h_l = v$h_l, h_r = v$h_r, b_l = v$b_l, b_r = v$b_r)
    fill_rb <- list(est = v$rob$est, lo = v$rob$lo, hi = v$rob$hi,
                    se = v$rob$se,
                    h_l = v$h_l, h_r = v$h_r, b_l = v$b_l, b_r = v$b_r)
    if (wt$status == "ok") {
      d_ll <- weight_diagnostics(wt$value$ll,  dat, v$conv$est, tol_rec)
      d_rb <- weight_diagnostics(wt$value$rbc, dat, v$rob$est,  tol_rec)
      b_ll <- balance_llrbc(wt$value$max_dev_ll,  tol_bal)
      b_rb <- balance_llrbc(wt$value$max_dev_rbc, tol_bal)
      fill_ll <- c(fill_ll, d_ll, b_ll)
      fill_rb <- c(fill_rb, d_rb, b_rb)
    }
    rows[["LL (h_MSE)"]]  <- row_from("LL (h_MSE)",  cm, fill_ll)
    rows[["RBC (h_MSE)"]] <- row_from("RBC (h_MSE)", cm, fill_rb)
    if (wt$status != "ok") {
      note <- paste0("[weight reconstruction failed] ", wt$error)
      rows[["LL (h_MSE)"]]  <- add_note(rows[["LL (h_MSE)"]],  note)
      rows[["RBC (h_MSE)"]] <- add_note(rows[["RBC (h_MSE)"]], note)
    }
  } else {
    rows[["LL (h_MSE)"]]  <- row_from("LL (h_MSE)",  cm)
    rows[["RBC (h_MSE)"]] <- row_from("RBC (h_MSE)", cm)
  }

  ## ---- rdrobust, CER-optimal bandwidths -> RBC (h_CER) ----
  cc <- capture_conditions(fit_rdrobust(dat$y, dat$x, "cerrd"))
  if (cc$status == "ok") {
    v <- cc$value
    wt <- capture_conditions(
      build_rd_weights(dat$x, v$h_l, v$h_r, v$b_l, v$b_r))
    fill <- list(est = v$rob$est, lo = v$rob$lo, hi = v$rob$hi,
                 se = v$rob$se,
                 h_l = v$h_l, h_r = v$h_r, b_l = v$b_l, b_r = v$b_r)
    if (wt$status == "ok") {
      fill <- c(fill,
                weight_diagnostics(wt$value$rbc, dat, v$rob$est, tol_rec),
                balance_llrbc(wt$value$max_dev_rbc, tol_bal))
    }
    rows[["RBC (h_CER)"]] <- row_from("RBC (h_CER)", cc, fill)
    if (wt$status != "ok") {
      rows[["RBC (h_CER)"]] <- add_note(
        rows[["RBC (h_CER)"]],
        paste0("[weight reconstruction failed] ", wt$error))
    }
  } else {
    rows[["RBC (h_CER)"]] <- row_from("RBC (h_CER)", cc)
  }

  ## ---- PLRD, one arm per configured window ----
  for (wname in names(cfg$plrd_windows)) {
    spec <- cfg$plrd_windows[[wname]]
    res  <- capture_conditions({
      mw <- resolve_plrd_window(spec, dat, ctx)
      fit_plrd(dat$y, dat$x, max_window = mw)
    })
    if (res$status == "ok") {
      v  <- res$value
      dg <- weight_diagnostics(v$gamma, dat, v$est, tol_rec)
      bl <- balance_plrd(v$gamma, dat$x, dat$w, tol_bal)
      bd <- plrd_bound_diag(dg$cond_bias, v$bhat, v$se)
      rows[[wname]] <- row_from(
        wname, res,
        c(list(est = v$est, lo = v$lo, hi = v$hi, se = v$se,
               pretest = v$pretest, bhat = v$bhat, Bhat = v$Bhat,
               plrd_window = v$window),
          dg, bl, bd))
    } else {
      rows[[wname]] <- row_from(wname, res)
    }
  }

  ## ---- RDHonest, data-driven rule-of-thumb bound ----
  if (cfg$run_ak) {
    ca <- capture_conditions(fit_rdhonest(dat$y, dat$x))
    if (ca$status == "ok") {
      v <- ca$value
      rows[["AK honest"]] <- row_from(
        "AK honest", ca,
        list(est = v$est, lo = v$lo, hi = v$hi, se = v$se,
             bhat = v$bhat, Bhat = v$M, h_l = v$h, h_r = v$h))
    } else {
      rows[["AK honest"]] <- row_from("AK honest", ca)
    }
  }

  out <- do.call(rbind, rows[cfg$methods])
  rownames(out) <- NULL
  out
}

## ------------------------------------------------------------
## Cell/block enumeration
## ------------------------------------------------------------
build_cells <- function(cfg) {
  cells <- list()
  add <- function(dist, design, ns) {
    for (n in ns) {
      cells[[length(cells) + 1L]] <<-
        list(dist = dist, design = design, n = as.integer(n))
    }
  }
  for (dg in cfg$designs) {
    ns <- if (dg %in% names(cfg$ns_override)) cfg$ns_override[[dg]]
          else cfg$ns_default
    add("beta", dg, ns)
  }
  if (cfg$run_uniform) for (dg in cfg$designs) add("unif", dg, cfg$unif_ns)
  cells
}

## ------------------------------------------------------------
## run_mc(): the full Monte Carlo
## ------------------------------------------------------------
run_mc <- function(cfg, src_dir) {
  ensure_output_dirs(cfg)
  t_start <- Sys.time()

  for (pkg in c("rdrobust", "plrd", if (cfg$run_ak) "RDHonest")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("package '", pkg, "' is required but not installed")
    }
  }

  ## fail-fast probe: exercise every wrapper once before the long run
  {
    set_task_seed(1L)
    xx <- draw_x_beta(400)
    yy <- 1 + xx + 0.25 * (xx >= 0) + 0.2 * rnorm(400)
    for (probe in list(
      quote(fit_rdrobust(yy, xx, "mserd")),
      quote(fit_plrd(yy, xx)),
      if (cfg$run_ak) quote(fit_rdhonest(yy, xx)))) {
      if (is.null(probe)) next
      pr <- capture_conditions(eval(probe))
      if (pr$status != "ok") {
        stop("estimator probe failed before the long run: ",
             deparse(probe), " -> ", pr$error)
      }
    }
    cat("estimator probes ok\n")
  }

  n_cores <- resolve_cores(cfg)
  use_par <- n_cores > 1L
  cl <- NULL
  if (use_par) {
    cat(sprintf("Starting PSOCK cluster with %d workers...\n", n_cores))
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(sd, run_ak) {
      suppressPackageStartupMessages({
        library(rdrobust); library(plrd)
        if (run_ak) library(RDHonest)
      })
      for (f in sort(list.files(sd, pattern = "\\.R$", full.names = TRUE))) {
        source(f)
      }
      NULL
    }, sd = src_dir, run_ak = cfg$run_ak)
  } else {
    cat("Running serially with 1 worker.\n")
  }

  cells <- build_cells(cfg)
  all_res <- vector("list", length(cells))
  cell_times <- numeric(length(cells))

  for (ci in seq_along(cells)) {
    cell <- cells[[ci]]
    cat(sprintf("[%s] dist %-4s design %-10s n = %-5d (%d/%d) ... ",
                format(Sys.time(), "%H:%M"), cell$dist, cell$design,
                cell$n, ci, length(cells)))
    t0 <- Sys.time()

    seeds <- vapply(seq_len(cfg$reps), function(r) {
      task_seed(cfg$master_seed, cell$dist, cell$design, cell$n, r)
    }, integer(1))
    if (anyDuplicated(seeds)) {
      stop("duplicated task seeds within cell ", cell$dist, "/",
           cell$design, "/", cell$n)
    }

    res <- if (use_par) {
      parallel::parLapply(cl, seq_len(cfg$reps), one_rep,
                          dist = cell$dist, design = cell$design,
                          n = cell$n, cfg = cfg)
    } else {
      lapply(seq_len(cfg$reps), one_rep,
             dist = cell$dist, design = cell$design,
             n = cell$n, cfg = cfg)
    }
    res <- do.call(rbind, res)
    stopifnot(nrow(res) == cfg$reps * length(cfg$methods))
    all_res[[ci]] <- res

    ## checkpoint so a crash late in a long run loses at most one cell
    saveRDS(do.call(rbind, all_res[seq_len(ci)]), cfg$files$checkpoint)

    cell_times[ci] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("done (%.1f min, %d/%d fits ok)\n", cell_times[ci],
                sum(res$success), nrow(res)))
  }

  all_res <- do.call(rbind, all_res)
  rownames(all_res) <- NULL
  attr(all_res, "schema_version") <- SCHEMA_VERSION
  attr(all_res, "mode")    <- cfg$mode
  attr(all_res, "created") <- format(Sys.time(), tz = "UTC",
                                     usetz = TRUE)
  saveRDS(all_res, cfg$files$results)
  if (file.exists(cfg$files$checkpoint)) file.remove(cfg$files$checkpoint)
  cat("saved ", cfg$files$results, " with ", nrow(all_res), " rows\n",
      sep = "")

  write_run_metadata(cfg, src_dir, t_start, Sys.time(), n_cores,
                     cells, cell_times)
  write_failure_log(all_res, cfg)
  write_warning_summary(all_res, cfg)

  invisible(all_res)
}

## ------------------------------------------------------------
## Reproducibility information
## ------------------------------------------------------------
describe_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return("not installed")
  d <- utils::packageDescription(pkg)
  bits <- c(Version = d$Version %||% NA,
            Packaged = d$Packaged %||% NA,
            RemoteSha = d$RemoteSha %||% d$GithubSHA1 %||% NA,
            Built = d$Built %||% NA)
  paste(names(bits)[!is.na(bits)], bits[!is.na(bits)],
        sep = "=", collapse = "; ")
}

## config with window functions made printable
serializable_cfg <- function(cfg) {
  cfg$plrd_windows <- lapply(cfg$plrd_windows, function(w) {
    if (is.null(w)) "NULL (package default: full range)"
    else if (is.function(w)) paste(deparse(w), collapse = " ")
    else w
  })
  cfg
}

write_run_metadata <- function(cfg, src_dir, t_start, t_end, n_cores,
                               cells, cell_times) {
  src_files <- sort(list.files(src_dir, pattern = "\\.R$",
                               full.names = TRUE))
  meta <- list(
    schema_version = SCHEMA_VERSION,
    mode = cfg$mode,
    started = format(t_start, tz = "UTC", usetz = TRUE),
    finished = format(t_end, tz = "UTC", usetz = TRUE),
    minutes_total = as.numeric(difftime(t_end, t_start, units = "mins")),
    r_version = R.version.string,
    platform = R.version$platform,
    n_cores = n_cores,
    master_seed = cfg$master_seed,
    seed_rule = paste("task_seed = 1 + (hash31('plrd-mc|schema|master|",
                      "dist|design|n') + rep - 1) mod (2^31 - 2);",
                      "Mersenne-Twister / Inversion / Rejection"),
    packages = list(rdrobust = describe_pkg("rdrobust"),
                    plrd = describe_pkg("plrd"),
                    RDHonest = describe_pkg("RDHonest")),
    config = serializable_cfg(cfg),
    cells = data.frame(
      dist = vapply(cells, `[[`, "", "dist"),
      design = vapply(cells, `[[`, "", "design"),
      n = vapply(cells, `[[`, 0L, "n"),
      minutes = round(cell_times, 2)),
    source_md5 = tools::md5sum(src_files)
  )
  saveRDS(meta, cfg$files$metadata)

  txt <- c(
    sprintf("run metadata (%s, schema %s)", cfg$mode, SCHEMA_VERSION),
    sprintf("started  : %s", meta$started),
    sprintf("finished : %s (%.1f min)", meta$finished, meta$minutes_total),
    sprintf("R        : %s on %s", meta$r_version, meta$platform),
    sprintf("cores    : %d", n_cores),
    sprintf("reps     : %d", cfg$reps),
    sprintf("master seed : %d", cfg$master_seed),
    sprintf("seed rule   : %s", meta$seed_rule),
    "packages :",
    sprintf("  rdrobust : %s", meta$packages$rdrobust),
    sprintf("  plrd     : %s", meta$packages$plrd),
    sprintf("  RDHonest : %s", meta$packages$RDHonest),
    sprintf("methods  : %s", paste(cfg$methods, collapse = ", ")),
    sprintf("designs  : %s", paste(cfg$designs, collapse = ", ")),
    sprintf("beta n grid : %s%s", paste(cfg$ns_default, collapse = ", "),
            if (length(cfg$ns_override))
              paste0("  (overrides: ",
                     paste(names(cfg$ns_override),
                           vapply(cfg$ns_override, paste, "",
                                  collapse = "/"),
                           sep = "=", collapse = "; "), ")") else ""),
    sprintf("uniform block : %s (n = %s)", cfg$run_uniform,
            paste(cfg$unif_ns, collapse = ", ")),
    "plrd windows :",
    paste0("  ", names(cfg$plrd_windows), " = ",
           unlist(serializable_cfg(cfg)$plrd_windows)),
    "source md5 :",
    sprintf("  %s  %s", meta$source_md5, basename(names(meta$source_md5)))
  )
  writeLines(txt, cfg$files$metadata_txt)
  cat("wrote ", cfg$files$metadata_txt, "\n", sep = "")
}

write_failure_log <- function(res, cfg) {
  bad <- res[!res$success,
             c("dist", "design", "n", "rep", "task_seed", "method",
               "status", "error", "n_warn", "warn_text"), drop = FALSE]
  utils::write.csv(bad, cfg$files$failures, row.names = FALSE)
  cat(sprintf("wrote %s (%d failed fits of %d)\n",
              cfg$files$failures, nrow(bad), nrow(res)))
}

write_warning_summary <- function(res, cfg) {
  w <- res[res$n_warn > 0 | res$n_msg > 0, , drop = FALSE]
  if (nrow(w) == 0L) {
    utils::write.csv(
      data.frame(dist = character(0), design = character(0),
                 method = character(0), warn_text = character(0),
                 count = integer(0)),
      cfg$files$warnings, row.names = FALSE)
    cat("wrote ", cfg$files$warnings, " (no warnings/messages)\n", sep = "")
    return(invisible(NULL))
  }
  w$warn_text[!nzchar(w$warn_text)] <- "(message only, no warning text)"
  agg <- stats::aggregate(list(count = w$rep),
                          by = list(dist = w$dist, design = w$design,
                                    method = w$method,
                                    warn_text = w$warn_text),
                          FUN = length)
  agg <- agg[order(-agg$count), ]
  utils::write.csv(agg, cfg$files$warnings, row.names = FALSE)
  cat("wrote ", cfg$files$warnings, "\n", sep = "")
}
