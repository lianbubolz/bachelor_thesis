## ============================================================
## dgp_definitions.R -- data-generating processes
##
## Every design is defined through explicit side-specific mean
## functions mu0(x) (control) and mu1(x) (treated) plus a
## conditional-sd function sd_fun(x). A generated sample carries
## the true conditional quantities needed for exact replication-
## level diagnostics:
##   x      running variable
##   w      treatment indicator, w = 1{x >= 0}
##   y      outcome
##   mu     realized conditional mean  mu_{w_i}(x_i)
##   sigma  realized conditional sd    sd_fun(x_i)
##   tau    true cutoff effect         mu1(0) - mu0(0)
## ============================================================

## ---------------- design parameters ----------------
DGP_PARAMS <- list(
  A_LOC     = 0.045,  # local-feature amplitude (pl_local)
  S_LOC     = 0.06,   # local-feature scale, ~0.5 * mean(h_MSE)
  DELTA_GAP = 16,     # curvature gap mu1'' - mu0'' (design gap)
  SD_BASE   = 0.1295, # Lee / Ludwig-Miller noise sd (CCT 2014)
  SD_PL     = 0.25,   # noise sd for the PL family
  SD_JUMP_L = 0.15,   # var_jump: sd left of the cutoff
  SD_JUMP_R = 0.45    # var_jump: sd right of the cutoff
)

## ---------------- running-variable draws and densities ----------------
draw_x_beta <- function(n) 2 * stats::rbeta(n, 2, 4) - 1
draw_x_unif <- function(n) stats::runif(n, -1, 1)

## Densities on [-1, 1] (normalized; used by the population checks
## and by the var_flat cross-check in the tests).
dens_beta <- function(x) 0.625 * (x + 1) * (1 - x)^3
dens_unif <- function(x) rep(0.5, length(x))

## P(X < 0) under each design density. Beta(2,4): pbeta(1/2, 2, 4)
## = 13/16 exactly; uniform: 1/2.
P_LEFT <- c(beta = 13 / 16, unif = 1 / 2)

## ------------------------------------------------------------
##  E[sigma^2(X)] = P(X<0) * sd_L^2 + P(X>=0) * sd_R^2,
## under the relevant running-variable distribution. 
## sqrt((sd_L^2 + sd_R^2)/2) = 0.335, for uniform X.
## Under Beta(2,4) (P(X<0) = 13/16) the matched value is
## sqrt(13/16 * 0.15^2 + 3/16 * 0.45^2) = sqrt(0.05625) = 0.23717.
## ------------------------------------------------------------
flat_sd <- function(dist) {
  stopifnot(dist %in% names(P_LEFT))
  p <- P_LEFT[[dist]]
  sqrt(p * DGP_PARAMS$SD_JUMP_L^2 + (1 - p) * DGP_PARAMS$SD_JUMP_R^2)
}

## ---------------- CCT (2014, S.3.1) quintics ----------------
## Model 1: Lee (2008) calibration.
lee_mu0 <- function(x) 0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5
lee_mu1 <- function(x) 0.52 + 0.84*x - 3.00*x^2 +  7.99*x^3 -  9.01*x^4 + 3.56*x^5
## Model 2: Ludwig-Miller (2007) calibration.
lm_mu0  <- function(x) 3.71 +  2.30*x +  3.28*x^2 +  1.45*x^3 +  0.23*x^4 + 0.03*x^5
lm_mu1  <- function(x) 0.26 + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5

## local feature used by pl_local
pl_feature <- function(x) {
  A <- DGP_PARAMS$A_LOC; S <- DGP_PARAMS$S_LOC
  A * sin(pi * x / S) * exp(-x^2 / (2 * (2 * S)^2))
}

## ---------------- DGP constructor ----------------
make_dgp <- function(mu0, mu1, sd_fun, draw, label) {
  force(mu0); force(mu1); force(sd_fun); force(draw); force(label)
  tau <- mu1(0) - mu0(0)
  gen <- function(n) {
    x     <- draw(n)
    w     <- as.numeric(x >= 0)
    mu    <- ifelse(w == 1, mu1(x), mu0(x))
    sigma <- rep_len(sd_fun(x), n)          # constants recycle to length n
    y     <- mu + sigma * stats::rnorm(n)
    list(x = x, w = w, y = y, mu = mu, sigma = sigma, tau = tau)
  }
  list(gen = gen, mu0 = mu0, mu1 = mu1, sd_fun = sd_fun,
       tau = tau, label = label)
}

## ---------------- design registry ----------------
## dist is "beta" or "unif";
make_dgp_list <- function(dist = c("beta", "unif")) {
  dist <- match.arg(dist)
  draw <- switch(dist, beta = draw_x_beta, unif = draw_x_unif)
  P    <- DGP_PARAMS
  sd_flat <- flat_sd(dist)

  pl_mu0  <- function(x) 0.5 + 1.27*x + 4*x^2 + 2*x^3
  pl_mu1  <- function(x) pl_mu0(x) + 0.25
  loc_mu0 <- function(x) 1 + x + pl_feature(x)
  loc_mu1 <- function(x) loc_mu0(x) + 0.25
  gap_mu0 <- function(x) 0.5 + 1.27*x + 2*x^2
  gap_mu1 <- function(x) gap_mu0(x) + 0.25 + 0.5*x + 0.5*P$DELTA_GAP*x^2

  list(
    ## D1: canonical applied baseline (tau = 0.04)
    lee = make_dgp(lee_mu0, lee_mu1,
                   sd_fun = function(x) P$SD_BASE, draw = draw,
                   label = "lee"),
    ## D2: strong-curvature stress test (tau = -3.45)
    lm = make_dgp(lm_mu0, lm_mu1,
                  sd_fun = function(x) P$SD_BASE, draw = draw,
                  label = "lm"),
    ## D3: PL true, modest smoothness (tau = 0.25)
    pl_smooth = make_dgp(pl_mu0, pl_mu1,
                         sd_fun = function(x) P$SD_PL, draw = draw,
                         label = "pl_smooth"),
    ## D4': PL true, odd local feature at the selectors' blind scale
    pl_local = make_dgp(loc_mu0, loc_mu1,
                        sd_fun = function(x) P$SD_PL, draw = draw,
                        label = "pl_local"),
    ## D5': pure quadratic PL violation, no cubic base content
    gap = make_dgp(gap_mu0, gap_mu1,
                   sd_fun = function(x) P$SD_PL, draw = draw,
                   label = "gap"),
    ## D6: PL true, variance jump at the cutoff (appendix)
    var_jump = make_dgp(pl_mu0, pl_mu1,
                        sd_fun = function(x) ifelse(x < 0, P$SD_JUMP_L,
                                                    P$SD_JUMP_R),
                        draw = draw, label = "var_jump"),
    ## D6': homoskedastic, sd matched to E[sigma^2(X)] of
    ## var_jump under this running-variable distribution
    var_flat = make_dgp(pl_mu0, pl_mu1,
                        sd_fun = function(x) sd_flat,
                        draw = draw, label = "var_flat")
  )
}

get_dgp <- function(dist, design) {
  reg <- make_dgp_list(dist)
  if (!design %in% names(reg)) stop("unknown design: ", design)
  reg[[design]]
}
