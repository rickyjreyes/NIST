#!/usr/bin/env Rscript
# run_holdout_replication.R
# ---------------------------------------------------------------------------
# Blocked holdout validation in log-wavenumber space.  k is learned on TRAINING
# bins and LOCKED before the test block is evaluated.  Confirmatory test blocks
# are never rescanned to choose k; exploratory rescans remain labelled.
#
# Two fixed-k null calibrations are reported:
#   fixed_k_test_p          : historical calibration with test baseline fixed
#   fixed_k_test_p_refit    : adversarial calibration that re-estimates the
#                             smoothing baseline inside every null replicate
# The conservative p-value is max(the two) and drives final holdout claims.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

poisson_loglik <- function(y, mu) sum(stats::dpois(y, pmax(mu, EPS), log = TRUE))

fit_block <- function(ell, y, baseline, degree, k) {
  X0 <- design_poly(ell, degree)
  f0 <- fit_poisson_loglinear(y, baseline, X0)
  X1 <- cbind(X0, cos(k * ell), sin(k * ell))
  f1 <- fit_poisson_loglinear(y, baseline, X1)
  b <- f1$beta; a_coef <- b[length(b) - 1L]; b_coef <- b[length(b)]
  list(deltaD = f0$deviance - f1$deviance,
       amplitude = sqrt(a_coef^2 + b_coef^2),
       phase = atan2(-b_coef, a_coef),
       ll0 = poisson_loglik(y, f0$mu), ll1 = poisson_loglik(y, f1$mu))
}

# Parametric fixed-k test.  The observed statistic always uses the observed
# test-block baseline.  Under refit_baseline=TRUE, every synthetic null dataset
# gets its own freshly estimated smoothing baseline before M0/M1 are fitted.
fixed_k_test_p <- function(ell, y, baseline, degree, k, B,
                           sigma = NULL, refit_baseline = FALSE) {
  if (isTRUE(refit_baseline) && (is.null(sigma) || !is.finite(sigma)))
    stop("sigma is required when refit_baseline=TRUE")
  obs <- fit_block(ell, y, baseline, degree, k)$deltaD
  X0 <- design_poly(ell, degree)
  mu0 <- fit_poisson_loglinear(y, baseline, X0)$mu
  seeds <- sample.int(.Machine$integer.max, B)
  vals <- vapply(seq_len(B), function(i) {
    set.seed(seeds[i]); y0 <- rpois(length(mu0), mu0)
    base0 <- if (isTRUE(refit_baseline))
      pmax(gaussian_filter_nearest(y0, sigma), EPS) else baseline
    fit_block(ell, y0, base0, degree, k)$deltaD
  }, numeric(1))
  list(p = emp_p(sum(vals >= obs), B), obs = obs, B = B,
       refit_baseline = isTRUE(refit_baseline))
}

phase_diff <- function(a, b) abs(((a - b + pi) %% (2 * pi)) - pi)

one_holdout <- function(design, ell, y, sigma, degree, k_grid, tr, te, B) {
  base_tr <- pmax(gaussian_filter_nearest(y[tr], sigma), EPS)
  base_te <- pmax(gaussian_filter_nearest(y[te], sigma), EPS)

  # Training-only frequency selection.
  sk <- scan_k(ell[tr], y[tr], base_tr, k_grid, degree)
  k_lock <- sk$best$k_best
  ftr <- fit_block(ell[tr], y[tr], base_tr, degree, k_lock)

  # Confirmatory test at locked k only.
  fte <- fit_block(ell[te], y[te], base_te, degree, k_lock)
  p_fixed <- fixed_k_test_p(ell[te], y[te], base_te, degree, k_lock, B,
                            sigma = sigma, refit_baseline = FALSE)
  p_refit <- fixed_k_test_p(ell[te], y[te], base_te, degree, k_lock, B,
                            sigma = sigma, refit_baseline = TRUE)
  p_cons <- max(p_fixed$p, p_refit$p)

  # Exploratory only: test-block rescan cannot rescue confirmatory failure.
  sk_te <- scan_k(ell[te], y[te], base_te, k_grid, degree)

  data.frame(
    design = design, n_train = length(tr), n_test = length(te),
    k_lock = k_lock,
    test_deltaD_lockedk = fte$deltaD,
    test_amplitude_lockedk = fte$amplitude,
    test_phase_lockedk = fte$phase,
    train_phase = ftr$phase,
    phase_diff = phase_diff(fte$phase, ftr$phase),
    predictive_loglik_gain = fte$ll1 - fte$ll0,
    direction_consistent = (fte$deltaD > 0),
    fixed_k_test_p = p_fixed$p,
    fixed_k_test_p_refit = p_refit$p,
    conservative_fixed_k_test_p = p_cons,
    test_B = B,
    exploratory_rescan_k = sk_te$best$k_best,
    exploratory_rescan_deltaD = sk_te$best$deltaD,
    stringsAsFactors = FALSE)
}

run_holdout <- function(cfg, bins = cfg$bins_primary, B = 500L, kfold = 5L) {
  root <- audit_repo_root()
  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              "Fe", 2L, "wavenumber")$lines
  k_grid <- audit_k_grid(cfg)
  binned <- build_binned(lines, bins)
  ell <- binned$ell; y <- binned$count
  n <- length(y); half <- floor(n / 2)
  sigma <- cfg$baseline_sigma; degree <- cfg$degree

  rows <- list()
  rows[[1]] <- one_holdout("lower_train_upper_test", ell, y, sigma, degree, k_grid,
                           seq_len(half), seq.int(half + 1L, n), B)
  rows[[2]] <- one_holdout("upper_train_lower_test", ell, y, sigma, degree, k_grid,
                           seq.int(half + 1L, n), seq_len(half), B)

  blk <- cut(seq_len(n), breaks = 4, labels = FALSE)
  tr3 <- which(blk %in% c(1, 3)); te3 <- which(blk %in% c(2, 4))
  rows[[3]] <- one_holdout("alternating_blocks", ell, y, sigma, degree, k_grid, tr3, te3, B)

  fold <- cut(seq_len(n), breaks = kfold, labels = FALSE)
  for (f in seq_len(kfold)) {
    te <- which(fold == f); tr <- which(fold != f)
    rows[[length(rows) + 1L]] <- one_holdout(sprintf("kfold_%d_of_%d", f, kfold),
                                             ell, y, sigma, degree, k_grid, tr, te, B)
  }
  do.call(rbind, rows)
}

plot_holdout <- function(ho) {
  ho$design <- factor(ho$design, levels = rev(ho$design))
  ggplot2::ggplot(ho, ggplot2::aes(test_deltaD_lockedk, design,
                                   fill = direction_consistent)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("p_cons=%.3f", conservative_fixed_k_test_p)),
                       hjust = -0.1, size = 3) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = species_colour("Fe"), `FALSE` = "#999999"),
                               name = "direction consistent") +
    ggplot2::expand_limits(x = max(ho$test_deltaD_lockedk, na.rm = TRUE) * 1.28) +
    ggplot2::labs(
      title = "Blocked holdout validation (k locked from training)",
      subtitle = "Conservative fixed-k p = max(fixed-baseline, baseline-refit null calibration).",
      x = "test-block deltaD at locked k", y = NULL,
      caption = "No confirmatory test-block rescan. Exploratory rescans are table diagnostics only.") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  B <- 500L
  i <- which(argv == "--null-n"); if (length(i) == 1L) B <- as.integer(argv[i + 1L])
  cfg <- default_audit_config()
  setup_rng(cfg$seed)
  root <- audit_repo_root()
  ho <- run_holdout(cfg, B = B)
  write_table(ho, file.path(root, "tables_r/statistical_audit/holdout_results.csv"))
  save_fig(plot_holdout(ho), file.path(root, "figures_r/statistical_audit/holdout_validation.png"))
  save_fig(plot_holdout(ho), file.path(root, "figures_r/statistical_audit/fig09_holdout_validation.png"))
  cat(sprintf("[holdout] %d designs; %d direction-consistent; median conservative p=%.3f\n",
              nrow(ho), sum(ho$direction_consistent),
              stats::median(ho$conservative_fixed_k_test_p)))
  invisible(ho)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_holdout_replication\\.R$", .invoked_file)) main()
