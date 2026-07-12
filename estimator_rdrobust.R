## ============================================================
## estimator_rdrobust.R -- wrapper around rdrobust
##
##   "LL (h_MSE)"   conventional local-linear estimate and CI
##                  ("Conventional" rows),
##   "RBC (...)"    bias-corrected estimate with the robust CI
##                  ("Robust" rows; the Bias-Corrected and Robust
##                  point estimates coincide by construction).
##
## package defaults (rdrobust 4.1.0):
##   c = 0, p = 1, q = 2, deriv = 0, kernel = "triangular",
##   vce = "nn", nnmatch = 3, level = 95, masspoints = "adjust",
##   scalepar = 1, scaleregul = 1.
## masspoints = "adjust" is the package default 
## with a continuous running variable it is equivalent to "off".
## ============================================================

fit_rdrobust <- function(y, x, bwselect = c("mserd", "cerrd")) {
  bwselect <- match.arg(bwselect)
  f <- rdrobust::rdrobust(y, x, c = 0, p = 1, q = 2, deriv = 0,
                          kernel = "triangular", bwselect = bwselect,
                          vce = "nn", nnmatch = 3, level = 95,
                          masspoints = "adjust",
                          scalepar = 1, scaleregul = 1)

  need_rows <- c("Conventional", "Bias-Corrected", "Robust")
  for (fld in c("coef", "se", "ci", "bws")) {
    if (is.null(f[[fld]])) stop("rdrobust object missing $", fld)
  }
  if (!all(need_rows %in% rownames(f$coef)) ||
      !all(need_rows %in% rownames(f$ci)) ||
      !all(need_rows %in% rownames(f$se))) {
    stop("rdrobust coef/se/ci rows changed; found: ",
         paste(rownames(f$coef), collapse = ", "))
  }
  if (!all(c("h", "b") %in% rownames(f$bws)) ||
      !all(c("left", "right") %in% colnames(f$bws))) {
    stop("rdrobust bws layout changed; found rows [",
         paste(rownames(f$bws), collapse = ","), "] cols [",
         paste(colnames(f$bws), collapse = ","), "]")
  }

  h_l <- strict_scalar(f$bws["h", "left"],  "h (left)")
  h_r <- strict_scalar(f$bws["h", "right"], "h (right)")
  b_l <- strict_scalar(f$bws["b", "left"],  "b (left)")
  b_r <- strict_scalar(f$bws["b", "right"], "b (right)")
  if (min(h_l, h_r, b_l, b_r) <= 0) {
    stop("rdrobust returned a non-positive bandwidth")
  }

  arm <- function(row) {
    list(est = strict_scalar(f$coef[row, 1], paste0(row, " coef")),
         lo  = strict_scalar(f$ci[row, 1],   paste0(row, " ci lower")),
         hi  = strict_scalar(f$ci[row, 2],   paste0(row, " ci upper")),
         se  = strict_scalar(f$se[row, 1],   paste0(row, " se")))
  }

  list(conv = arm("Conventional"),
       rob  = arm("Robust"),
       h_l = h_l, h_r = h_r, b_l = b_l, b_r = b_r)
}
