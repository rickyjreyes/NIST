#!/usr/bin/env Rscript
# model_comparison.R
# ---------------------------------------------------------------------------
# Compare two nested Poisson models on the Fe II binned spectrum:
#   M0: smooth null            log(mu) = log(B) + poly(ell, degree)
#   M1: smooth + log-periodic  ... + a*cos(k*ell) + b*sin(k*ell)   at selected k
#
# Reports log-likelihood, parameter count, AIC, AICc, BIC, deltaAIC, deltaBIC,
# deviance difference, and a held-out predictive log-likelihood from a blocked
# split with k LOCKED from the training block.
#
# IMPORTANT caveat (written into the output and report): because k is chosen by
# a frequency scan, naive AIC/BIC at the selected k do NOT fully account for the
# look-elsewhere search. The scan-global bootstrap p-value remains the principal
# multiplicity correction. M0 is a smooth statistical baseline, NOT a complete
# physical atomic model.
#
# Writes:
#   tables_r/statistical_audit/model_comparison.csv
#   figures_r/statistical_audit/model_comparison.png (+ fig11)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

poisson_loglik <- function(y, mu) sum(stats::dpois(y, pmax(mu, EPS), log = TRUE))

aic_bic <- function(loglik, npar, n) {
  aic <- -2 * loglik + 2 * npar
  aicc <- if (n - npar - 1 > 0) aic + (2 * npar * (npar + 1)) / (n - npar - 1) else NA_real_
  bic <- -2 * loglik + npar * log(n)
  c(aic = aic, aicc = aicc, bic = bic)
}

# fit_models(): fit M0 and M1 at a given k on (ell,y,baseline). Returns metrics.
fit_models <- function(ell, y, baseline, degree, k) {
  X0 <- design_poly(ell, degree)
  f0 <- fit_poisson_loglinear(y, baseline, X0)
  X1 <- cbind(X0, cos(k * ell), sin(k * ell))
  f1 <- fit_poisson_loglinear(y, baseline, X1)
  n <- length(y)
  p0 <- ncol(X0); p1 <- ncol(X1)
  ll0 <- poisson_loglik(y, f0$mu); ll1 <- poisson_loglik(y, f1$mu)
  m0 <- aic_bic(ll0, p0, n); m1 <- aic_bic(ll1, p1, n)
  list(
    ll0 = ll0, ll1 = ll1, p0 = p0, p1 = p1, n = n,
    deviance0 = f0$deviance, deviance1 = f1$deviance,
    m0 = m0, m1 = m1, mu0 = f0$mu, mu1 = f1$mu)
}

# heldout_predictive(): blocked split (lower half train / upper half test).
# k is estimated on the training block and LOCKED, then both models are refit
# on training and scored by predictive Poisson log-likelihood on the test block.
heldout_predictive <- function(ell, y, baseline, degree, k_grid) {
  n <- length(y)
  mid <- floor(n / 2)
  tr <- seq_len(mid); te <- seq.int(mid + 1L, n)
  # lock k on training
  sk <- scan_k(ell[tr], y[tr], baseline[tr], k_grid, degree)
  k_lock <- sk$best$k_best
  X0tr <- design_poly(ell[tr], degree)
  X1tr <- cbind(X0tr, cos(k_lock * ell[tr]), sin(k_lock * ell[tr]))
  f0 <- fit_poisson_loglinear(y[tr], baseline[tr], X0tr)
  f1 <- fit_poisson_loglinear(y[tr], baseline[tr], X1tr)
  # predict on test using training betas (design built on test ell, centered by
  # the training transform is approximated by recomputing design on test ell)
  predict_mu <- function(beta, ell_te, base_te, degree, k = NULL) {
    X0 <- design_poly(ell_te, degree)
    X <- if (is.null(k)) X0 else cbind(X0, cos(k * ell_te), sin(k * ell_te))
    eta <- pmin(pmax(as.numeric(X %*% beta), -10), 10)
    pmax(base_te * exp(eta), EPS)
  }
  mu0_te <- predict_mu(f0$beta, ell[te], baseline[te], degree)
  mu1_te <- predict_mu(f1$beta, ell[te], baseline[te], degree, k_lock)
  list(k_lock = k_lock,
       ll0_test = poisson_loglik(y[te], mu0_te),
       ll1_test = poisson_loglik(y[te], mu1_te))
}

main <- function(argv = commandArgs(TRUE)) {
  cfg <- default_audit_config()
  root <- audit_repo_root()
  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              "Fe", 2L, "wavenumber")$lines
  k_grid <- audit_k_grid(cfg)
  res <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)
  k <- res$best$k_best
  fm <- fit_models(res$ell, res$y, res$baseline, cfg$degree, k)
  ho <- heldout_predictive(res$ell, res$y, res$baseline, cfg$degree, k_grid)

  tab <- data.frame(
    model = c("M0_smooth_null", "M1_smooth_plus_logperiodic"),
    loglik = c(fm$ll0, fm$ll1),
    n_params = c(fm$p0, fm$p1),
    n_bins = fm$n,
    deviance = c(fm$deviance0, fm$deviance1),
    AIC = c(fm$m0["aic"], fm$m1["aic"]),
    AICc = c(fm$m0["aicc"], fm$m1["aicc"]),
    BIC = c(fm$m0["bic"], fm$m1["bic"]),
    heldout_loglik_test = c(ho$ll0_test, ho$ll1_test),
    stringsAsFactors = FALSE)
  tab$deltaAIC <- tab$AIC - min(tab$AIC)
  tab$deltaBIC <- tab$BIC - min(tab$BIC)
  tab$selected_k <- k
  tab$heldout_k_lock <- ho$k_lock
  tab$deviance_difference <- fm$deviance0 - fm$deviance1
  tab$note <- "k selected by scan; AIC/BIC do not fully account for look-elsewhere"
  write_table(tab, file.path(root, "tables_r/statistical_audit/model_comparison.csv"))

  plotdf <- data.frame(
    criterion = rep(c("AIC", "BIC"), each = 2),
    model = rep(c("M0", "M1"), 2),
    value = c(tab$AIC, tab$BIC))
  p <- ggplot2::ggplot(plotdf, ggplot2::aes(model, value, fill = model)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::facet_wrap(~criterion, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(M0 = "#999999", M1 = species_colour("Fe")),
                               guide = "none") +
    ggplot2::labs(
      title = "Model comparison: smooth null (M0) vs smooth + log-periodic (M1)",
      subtitle = sprintf("deltaAIC(M1)=%.1f  deltaBIC(M1)=%.1f  held-out loglik gain=%.2f (k locked from training)",
                         tab$deltaAIC[2], tab$deltaBIC[2],
                         ho$ll1_test - ho$ll0_test),
      x = NULL, y = "information criterion (lower is better)",
      caption = "M0 is a smooth statistical baseline, not a complete physical model.") +
    theme_audit()
  save_fig(p, file.path(root, "figures_r/statistical_audit/model_comparison.png"))
  save_fig(p, file.path(root, "figures_r/statistical_audit/fig11_model_comparison.png"))

  cat(sprintf("[model_comparison] deltaAIC=%.2f deltaBIC=%.2f heldout_gain=%.3f\n",
              tab$deltaAIC[2], tab$deltaBIC[2], ho$ll1_test - ho$ll0_test))
  invisible(tab)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("model_comparison\\.R$", .invoked_file)) {
  main()
}
