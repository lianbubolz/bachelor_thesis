## ============================================================
## diagnostics.R
##
## For any arm whose point estimate is linear in y with known
## weights lambda, and a sample carrying the true conditional
## quantities (mu, sigma, tau), the estimate decomposes EXACTLY as
##
##   est = tau + cond_bias + noise,
##   cond_bias = sum_i lambda_i mu_i - tau,
##   noise     = sum_i lambda_i (y_i - mu_i).
##
## wt_gap    = sum(lambda * y) - est  checks that the weights
##             reproduce the reported estimate (reconstruction for
##             LL/RBC and the returned gamma for PLRD).
## decomp_gap = est - (tau + cond_bias + noise)
## Both are compared against tol_rec relative to max(1, |est|).
## ============================================================

weight_diagnostics <- function(lambda, dat, est, tol_rec = 1e-8) {
  eps       <- dat$y - dat$mu
  cond_bias <- sum(lambda * dat$mu) - dat$tau
  noise     <- sum(lambda * eps)
  wt_gap    <- sum(lambda * dat$y) - est
  decomp_gap <- est - (dat$tau + cond_bias + noise)
  scale <- max(1, abs(est))
  list(cond_bias = cond_bias, noise = noise,
       wt_gap = wt_gap, decomp_gap = decomp_gap,
       ok_weights = is.finite(wt_gap) &&
         abs(wt_gap) <= tol_rec * scale &&
         abs(decomp_gap) <= tol_rec * scale)
}

## Balance for LL / RBC: exact side-specific identities from the
## reconstruction (see weight_reconstruction.R).
balance_llrbc <- function(max_dev, tol_bal = 1e-8) {
  list(bal_dev = max_dev,
       bal_gx = NA_real_, bal_gx2 = NA_real_,
       ok_balance = is.finite(max_dev) && max_dev <= tol_bal)
}

## Balance for PLRD: sum(gamma) = 0 and side sums +1/-1 are exact
## by construction and are checked against tol_bal.
## The pooled linear and quadratic moments
## sum(gamma * x), sum(gamma * x^2) hold only up to the package's
## grid discretization (400 bins) and post-hoc renormalization
balance_plrd <- function(gamma, x, w, tol_bal = 1e-8) {
  dev <- c(total     = sum(gamma),
           sum_plus  = sum(gamma[w == 1]) - 1,
           sum_minus = sum(gamma[w == 0]) + 1)
  list(bal_dev = max(abs(dev)),
       bal_gx  = sum(gamma * x),
       bal_gx2 = sum(gamma * x^2),
       ok_balance = max(abs(dev)) <= tol_bal)
}

## PLRD bias-bound diagnostics: absolute conditional bias against
## the reported worst-case bound b-hat, per replication.
##   bound_ratio    |cond_bias| / b-hat
##   bound_exceed   1{|cond_bias| > b-hat}
##   bound_exceed_z (|cond_bias| - b-hat) / s-hat  (negative while
##                  inside the bound; the exceedance in SE units)
plrd_bound_diag <- function(cond_bias, bhat, shat) {
  ab <- abs(cond_bias)
  list(bound_ratio    = if (is.finite(bhat) && bhat > 0) ab / bhat else NA_real_,
       bound_exceed   = as.numeric(is.finite(bhat) && ab > bhat),
       bound_exceed_z = if (is.finite(shat) && shat > 0) (ab - bhat) / shat
                        else NA_real_)
}
