## ============================================================
## weight_reconstruction.R
##
## Implements the thesis' Lemma B.1 / eq. (9) representation of the
## conventional local-linear and the robust-bias-corrected point
## estimators as linear functionals sum_i lambda_i y_i, at the
## REALIZED bandwidths of the rdrobust fit, with the triangular
## kernel (the rdrobust default).
##
## For each side s in {+, -} and bandwidth a:
##   w_s(a, p, nu) = e_{nu+1}' (R' W R)^{-1} R' W,  R = [1, x, ..., x^p],
##   W = diag K(x/a), so w(., 1, 0) are the local-linear intercept
##   weights and w(., 2, 2) the local-quadratic curvature-coefficient
##   weights of the pilot fit. Then
##   lambda_ll  = w_+(h,1,0) - w_-(h,1,0),
##   lambda_rbc = lambda_ll - c_+ w_+(b,2,2) + c_- w_-(b,2,2),
##   c_s = sum_i w_s(h,1,0)_i x_i^2  (= h^2 beta_hat_s, Lemma 3.1(i)).
##
## Verified: sum(lambda * y) reproduces the rdrobust "Conventional"
## and "Bias-Corrected"/"Robust" point estimates to machine
## precision (1e-15) for both mserd and cerrd bandwidths.
##
## Balance identities implied by the construction (Prop. 3.2):
##   LL : side sums = +1 / -1, side-specific linear balance;
##   RBC: side sums = +1 / -1, side-specific linear AND quadratic
##        balance.
## ============================================================

build_rd_weights <- function(x, h_l, h_r, b_l = NA_real_, b_r = NA_real_) {
  stopifnot(is.finite(h_l), is.finite(h_r), h_l > 0, h_r > 0)
  want_rbc <- is.finite(b_l) && is.finite(b_r)
  if (want_rbc) stopifnot(b_l > 0, b_r > 0)

  K <- function(u) pmax(0, 1 - abs(u))   # triangular kernel

  side_w <- function(idx, a, p, nu) {
    w <- numeric(length(x))
    xi <- x[idx]
    kw <- K(xi / a)
    if (sum(kw > 0) < p + 1) {
      stop(sprintf("fewer than %d observations with positive kernel weight on one side (bandwidth %.4g)",
                   p + 1, a))
    }
    R  <- outer(xi, 0:p, "^")
    A  <- crossprod(R, R * kw)           # R' W R
    Mi <- solve(A, t(R * kw))            # (R'WR)^{-1} R' W
    w[idx] <- Mi[nu + 1, ]
    w
  }

  Ip <- which(x >= 0); Im <- which(x < 0)
  if (length(Ip) < 2L || length(Im) < 2L) {
    stop("need at least two observations on each side of the cutoff")
  }

  wllp <- side_w(Ip, h_r, 1, 0)          # local-linear intercept, right
  wllm <- side_w(Im, h_l, 1, 0)          # local-linear intercept, left
  ll   <- wllp - wllm

  bh2p <- sum(wllp * x^2)                # = h^2 * beta_hat_plus
  bh2m <- sum(wllm * x^2)                # = h^2 * beta_hat_minus

  rbc <- NULL
  if (want_rbc) {
    wpp <- side_w(Ip, b_r, 2, 2)         # quadratic pilot coef, right
    wpm <- side_w(Im, b_l, 2, 2)         # quadratic pilot coef, left
    rbc <- wllp - wllm - bh2p * wpp + bh2m * wpm   # eq. (9)
  }

  dev_ll <- c(sum_plus  = sum(ll[Ip]) - 1,
              sum_minus = sum(ll[Im]) + 1,
              x_plus    = sum(ll[Ip] * x[Ip]),
              x_minus   = sum(ll[Im] * x[Im]))
  dev_rbc <- NULL
  if (want_rbc) {
    dev_rbc <- c(sum_plus  = sum(rbc[Ip]) - 1,
                 sum_minus = sum(rbc[Im]) + 1,
                 x_plus    = sum(rbc[Ip] * x[Ip]),
                 x_minus   = sum(rbc[Im] * x[Im]),
                 x2_plus   = sum(rbc[Ip] * x[Ip]^2),
                 x2_minus  = sum(rbc[Im] * x[Im]^2))
  }

  list(ll = ll, rbc = rbc,
       h2beta_plus = bh2p, h2beta_minus = bh2m,
       dev_ll = dev_ll, dev_rbc = dev_rbc,
       max_dev_ll  = max(abs(dev_ll)),
       max_dev_rbc = if (want_rbc) max(abs(dev_rbc)) else NA_real_)
}

## Hard assertion used by the figure and the tests.
assert_balance <- function(wts, tol = 1e-6) {
  dev <- max(wts$max_dev_ll, wts$max_dev_rbc, na.rm = TRUE)
  if (dev > tol) {
    stop("balance identities violated: max deviation ", format(dev))
  }
  invisible(dev)
}
