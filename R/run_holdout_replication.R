#!/usr/bin/env Rscript
# run_holdout_replication.R
# ---------------------------------------------------------------------------
# Blocked holdout validation in log-wavenumber space (NOT random per-line
# splits alone). Designs:
#   1. lower-range discovery / upper-range validation
#   2. upper-range discovery / lower-range validation
#   3. alternating contiguous blocks
#   4. repeated blocked K-fold
#
# CRITICAL confirmatory rule: k is estimated on TRAINING bins and LOCKED before
# the test block is touched. The test block is NOT rescanned in confirmatory
# mode. An exploratory rescan of the test block is reported separately and
# clearly labelled.
#
# Writes:
#   tables_r/statistical_audit/holdout_results.csv
#   figures_r/statistical_audit/holdout_validation.png (+ fig09)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

poisson_loglik <- function(y, mu) sum(stats::dpois(y, pmax(mu, EPS), log = TRUE))

# fit_block(): fit M0 and M1 (at a given k) on a set of bins; return metrics.
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

# fixed_k_test_p(): parametric null p-value for the test block deltaD AT the
# locked k (no rescan). y0 ~ Poisson(mu0_test); deltaD at fixed k.
fixed_k_test_p <- function(ell, y, baseline, degree, k, B) {
  obs <- fit_block(ell, y, baseline, degree, k)$deltaD
  X0 <- design_poly(ell, degree)
  mu0 <- fit_poisson_loglinear(y, baseline, X0)$mu
  seeds <- sample.int(.Machine$integer.max, B)
  vals <- vapply(seq_len(B), function(i) {
    set.seed(seeds[i]); y0 <- rpois(length(mu0), mu0)
    fit_block(ell, y0, baseline, degree, k)$deltaD
  }, numeric(1))
  list(p = emp_p(sum(vals >= obs), B), obs = obs, B = B)
}

# circular phase difference in [0, pi]
phase_diff <- function(a, b) {
  d <- abs(((a - b + pi) %% (2 * pi)) - pi); d
}

# one_holdout(): train on `tr` bins, lock k, validate on `te` bins.
one_holdout <- function(design, ell, y, sigma, degree, k_grid, tr, te, B) {
  base_tr <- pmax(gaussian_filter_nearest(y[tr], sigma), EPS)
  base_te <- pmax(gaussian_filter_nearest(y[te], sigma), EPS)
  # lock k on training
  sk <- scan_k(ell[tr], y[tr], base_tr, k_grid, degree)
  k_lock <- sk$best$k_best
  ftr <- fit_block(ell[tr], y[tr], base_tr, degree, k_lock)
  # confirmatory: locked-k on test
  fte <- fit_block(ell[te], y[te], base_te, degree, k_lock)
  ptest <- fixed_k_test_p(ell[te], y[te], base_te, degree, k_lock, B)
  # exploratory: rescan test
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
    fixed_k_test_p = ptest$p, test_B = B,
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
  # 1. lower train / upper test
  rows[[1]] <- one_holdout("lower_train_upper_test", ell, y, sigma, degree, k_grid,
                           seq_len(half), seq.int(half + 1L, n), B)
  # 2. upper train / lower test
  rows[[2]] <- one_holdout("upper_train_lower_test", ell, y, sigma, degree, k_grid,
                           seq.int(half + 1L, n), seq_len(half), B)
  # 3. alternating contiguous blocks (4 blocks: train on 1,3 ; test on 2,4)
  blk <- cut(seq_len(n), breaks = 4, labels = FALSE)
  tr3 <- which(blk %in% c(1, 3)); te3 <- which(blk %in% c(2, 4))
  rows[[3]] <- one_holdout("alternating_blocks", ell, y, sigma, degree, k_grid, tr3, te3, B)
  # 4. repeated blocked K-fold
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
    ggplot2::geom_text(ggplot2::aes(label = sprintf("p=%.3f", fixed_k_test_p)),
                       hjust = -0.1, size = 3) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = species_colour("Fe"), `FALSE` = "#999999"),
                               name = "direction consistent") +
    ggplot2::expand_limits(x = max(ho$test_deltaD_lockedk, na.rm = TRUE) * 1.25) +
    ggplot2::labs(
      title = "Blocked holdout validation (k locked from training)",
      subtitle = "Test-block deviance improvement at the locked k; fixed-k empirical p shown.",
      x = "test-block deltaD at locked k", y = NULL,
      caption = "Confirmatory: no test-block rescan. Exploratory rescans reported in the table only.") +
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
  cat(sprintf("[holdout] %d designs; %d with consistent direction; median fixed-k p=%.3f\n",
              nrow(ho), sum(ho$direction_consistent), stats::median(ho$fixed_k_test_p)))
  invisible(ho)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_holdout_replication\\.R$", .invoked_file)) {
  main()
}
