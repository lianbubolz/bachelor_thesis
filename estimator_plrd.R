## ============================================================
## estimator_plrd.R -- wrapper around plrd
##
## github.com/ghoshadi/plrd (v0.0.1):
##   * plrd(Y, X, threshold, max.window = NULL, alpha = 0.95,
##          seed = 42, ...): max.window defaults to the full data
##     range; alpha is the COVERAGE level; seed governs only the
##     internal cross-fitting fold split (drawn under
##     withr::with_seed, so it does not use the RNG
##     stream and is deterministic given the data).
##   * Returned fields used here: tau.hat, ci.lower, ci.upper,
##     max.bias (= b-hat), sampling.se (= s-hat),
##     Lipschitz.constant (= B-hat, original units),
##     diff.curvatures (LOGICAL safeguard flag: TRUE = the ANOVA
##     F-test of partial linearity at signif.curvature = 0.001
##     rejected and the side-specific-curvature program was used),
##     gamma (length(Y), input order, folds already combined with
##     the 1/2 factors; satisfies sum(gamma * Y) = tau.hat and
##     sums to +1 / -1 on the treated / control side; zero outside
##     max.window).
##
## The window is configurable: the package's full-range default is
## the main specification; shrinking-window alternatives are run
## by passing a positive max_window (see config$plrd_windows).
## ============================================================

fit_plrd <- function(y, x, max_window = NULL, alpha = 0.95, seed = 42) {
  if (!is.null(max_window)) {
    max_window <- strict_scalar(max_window, "plrd max_window")
    if (max_window <= 0) stop("plrd max_window must be positive")
  }

  f <- plrd::plrd(Y = y, X = x, threshold = 0,
                  max.window = max_window,
                  alpha = alpha, seed = seed)

  need <- c("tau.hat", "ci.lower", "ci.upper", "max.bias",
            "sampling.se", "Lipschitz.constant", "diff.curvatures",
            "gamma", "max.window")
  if (!all(need %in% names(f))) {
    stop("plrd return fields changed; missing: ",
         paste(setdiff(need, names(f)), collapse = ", "),
         "; found: ", paste(names(f), collapse = ", "))
  }

  gamma <- as.numeric(f$gamma)
  if (length(gamma) != length(y) || anyNA(gamma)) {
    stop("plrd gamma is not aligned with the sample (length ",
         length(gamma), " vs n = ", length(y), ")")
  }
  pretest <- f$diff.curvatures
  if (!is.logical(pretest) || length(pretest) != 1L || is.na(pretest)) {
    stop("plrd diff.curvatures is not a single logical")
  }

  list(est  = strict_scalar(f$tau.hat, "plrd tau.hat"),
       lo   = strict_scalar(f$ci.lower, "plrd ci.lower"),
       hi   = strict_scalar(f$ci.upper, "plrd ci.upper"),
       se   = strict_scalar(f$sampling.se, "plrd sampling.se"),
       bhat = strict_scalar(f$max.bias, "plrd max.bias"),
       Bhat = strict_scalar(f$Lipschitz.constant, "plrd Lipschitz.constant"),
       pretest = as.numeric(pretest),
       gamma = gamma,
       window = strict_scalar(f$max.window, "plrd max.window"))
}

resolve_plrd_window <- function(spec, dat, ctx) {
  if (is.null(spec)) return(NULL)
  if (is.function(spec)) {
    w <- spec(dat, ctx)
    return(strict_scalar(w, "plrd window function result"))
  }
  strict_scalar(spec, "plrd window spec")
}
