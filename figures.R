## ============================================================
## figures.R -- weight-vector figure (single draw)
##
## LL / RBC weights are reconstructed at the realized rdrobust
## bandwidths, PLRD weights are the packages gamma.
## ============================================================

make_weight_figure <- function(cfg) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("figure skipped: ggplot2 not installed")
    return(invisible(NULL))
  }
  fig <- cfg$fig
  pdf_file  <- file.path(cfg$paths$figures, "fig_weights.pdf")
  info_file <- file.path(cfg$paths$figures, "fig_weights_info.txt")

  set_task_seed(fig$seed)
  dgp <- get_dgp("beta", fig$design)
  dat <- dgp$gen(fig$n)

  rr <- fit_rdrobust(dat$y, dat$x, "mserd")
  wts <- build_rd_weights(dat$x, rr$h_l, rr$h_r, rr$b_l, rr$b_r)
  assert_balance(wts, tol = 1e-6)
  gap_ll  <- sum(wts$ll  * dat$y) - rr$conv$est
  gap_rbc <- sum(wts$rbc * dat$y) - rr$rob$est

  pl <- capture_conditions(fit_plrd(dat$y, dat$x))
  if (pl$status != "ok") {
    message("figure skipped: plrd failed on the figure draw: ", pl$error)
    return(invisible(NULL))
  }
  pv <- pl$value
  gam <- pv$gamma
  tau_dev <- sum(gam * dat$y) - pv$est
  out_share <- sum(abs(gam[abs(dat$x) > max(abs(fig$xlim))])) / sum(abs(gam))

  pd <- rbind(
    data.frame(x = dat$x, w = wts$ll,  method = "LL (h_MSE)"),
    data.frame(x = dat$x, w = wts$rbc, method = "RBC (h_MSE)"),
    data.frame(x = dat$x, w = gam,     method = "PLRD")
  )
  pd$method <- factor(pd$method,
                      levels = c("LL (h_MSE)", "RBC (h_MSE)", "PLRD"))

  h_fig <- rr$h_l; b_fig <- rr$b_l
  p <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = w, colour = method,
                                        shape = method)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                        linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-h_fig, h_fig), linetype = "dashed",
                        colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-b_fig, b_fig), linetype = "dotted",
                        colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.1, alpha = 0.75) +
    ggplot2::coord_cartesian(xlim = fig$xlim) +
    ggplot2::labs(x = "running variable x", y = "weight",
                  colour = NULL, shape = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(pdf_file, p, width = 6.3, height = 3.2)
  cat("wrote ", pdf_file, "\n", sep = "")

  info <- c(
    if (cfg$mode == "smoke") "** SMOKE-mode figure: pipeline test only **",
    sprintf("design = %s, n = %d, seed = %d (single draw)",
            fig$design, fig$n, fig$seed),
    sprintf("realized h_MSE = %.4f, b = %.4f, rho = h/b = %.3f",
            h_fig, b_fig, h_fig / b_fig),
    sprintf("h^2 beta_hat: plus = %.5f, minus = %.5f",
            wts$h2beta_plus, wts$h2beta_minus),
    sprintf("reconstruction: sum(ll*y) - conv est = %.2e, sum(rbc*y) - bc est = %.2e",
            gap_ll, gap_rbc),
    sprintf("rbc/ll balance identities: max abs deviation = %.2e",
            max(wts$max_dev_ll, wts$max_dev_rbc)),
    sprintf("plrd: sum(gamma*Y) - tau.hat = %.2e (should be ~0)", tau_dev),
    sprintf("plrd: diff.curvatures (safeguard fired) = %s",
            as.logical(pv$pretest)),
    sprintf("plrd: B-hat = %.4f, b-hat = %.5f, s-hat = %.5f",
            pv$Bhat, pv$bhat, pv$se),
    sprintf("plrd: exact conditional bias = %+.5f, |bias|/b-hat = %.3f",
            sum(gam * dat$mu) - dat$tau,
            abs(sum(gam * dat$mu) - dat$tau) / pv$bhat),
    sprintf("share of |gamma| outside plotted xlim [%.2f, %.2f]: %.4f",
            fig$xlim[1], fig$xlim[2], out_share),
    "vertical lines: dashed = +/- h_MSE, dotted = +/- b (pilot);",
    "PLRD window is the full data range (package default)."
  )
  writeLines(info, info_file)
  cat("wrote ", info_file, "\n", sep = "")
  cat(paste0("  ", info, "\n"), sep = "")
  invisible(list(pdf = pdf_file, info = info_file))
}
