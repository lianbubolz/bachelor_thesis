## ============================================================
## estimator_rdhonest.R -- strict wrapper around RDHonest
##
## Uses the data-driven rule-of-thumb smoothness bound: M
## is left unspecified, so RDHonest applies the Armstrong & Kolesar
## (2020) ROT
##
## Arguments (RDHonest defaults):
##   cutoff = 0, kern = "triangular", opt.criterion = "MSE",
##   se.method = "nn", J = 3, alpha = 0.05, sclass = "H".
## ============================================================

fit_rdhonest <- function(y, x) {
  f <- RDHonest::RDHonest(y ~ x, data = data.frame(y = y, x = x),
                          cutoff = 0, kern = "triangular",
                          opt.criterion = "MSE", se.method = "nn",
                          J = 3, alpha = 0.05, sclass = "H")

  cf <- f$coefficients
  if (is.null(cf) || !is.data.frame(cf) || nrow(cf) < 1L) {
    stop("RDHonest object has no coefficients data frame")
  }
  need <- c("estimate", "std.error", "maximum.bias",
            "conf.low", "conf.high", "bandwidth", "M")
  if (!all(need %in% colnames(cf))) {
    stop("RDHonest coefficient columns changed; missing: ",
         paste(setdiff(need, colnames(cf)), collapse = ", "),
         "; found: ", paste(colnames(cf), collapse = ", "))
  }

  h <- strict_scalar(cf[1, "bandwidth"], "RDHonest bandwidth")
  list(est  = strict_scalar(cf[1, "estimate"],  "RDHonest estimate"),
       lo   = strict_scalar(cf[1, "conf.low"],  "RDHonest conf.low"),
       hi   = strict_scalar(cf[1, "conf.high"], "RDHonest conf.high"),
       se   = strict_scalar(cf[1, "std.error"], "RDHonest std.error"),
       bhat = strict_scalar(cf[1, "maximum.bias"], "RDHonest maximum.bias"),
       M    = strict_scalar(cf[1, "M"], "RDHonest M"),
       h    = h)
}
