#!/usr/bin/env Rscript
# calibrate_false_positive_rate.R
# ---------------------------------------------------------------------------
# Synthetic-null calibration: does the full scan-global discovery procedure
# produce nominal false-positive rates?
#
# Design (tractable nested Monte Carlo): all synthetic datasets are generated
# under the SAME fitted smooth null mu0, so the null distribution of the
# scan-global max statistic is shared. We therefore:
#   1. fit mu0 on the Fe II primary binned spectrum;
#   2. build ONE reference null-max distribution of size null_n (the inner null);
#   3. draw calibration_n INDEPENDENT synthetic datasets y* ~ Poisson(mu0),
#      compute each scan-global max statistic, and its empirical p-value against
#      the reference null using p = (1 + r)/(1 + null_n);
#   4. estimate the false-positive rate = mean(p <= alpha) at alpha in
#      {0.10, 0.05, 0.01}, with an exact binomial CI.
#
# Same frequency grid and selection rule are used throughout. We do NOT claim
# exact calibration when the Monte Carlo CI is wide.
#
# Writes:
#   tables_r/statistical_audit/null_calibration.csv
#   figures_r/statistical_audit/null_calibration.png (+ fig13)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

calibrate <- function(cfg, calibration_n, null_n, parallel = FALSE) {
  root <- audit_repo_root()
  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              "Fe", 2L, "wavenumber")$lines
  k_grid <- audit_k_grid(cfg)
  ref <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)
  mu0 <- ref$mu0; ell <- ref$ell; baseline <- ref$baseline; degree <- cfg$degree

  # inner reference null
  nd <- null_distribution(ell, ref$y, baseline, mu0, k_grid, degree, null_n, parallel = parallel)
  ref_null <- nd$max_vals

  # outer synthetic datasets (independent of inner draws)
  seeds <- sample.int(.Machine$integer.max, calibration_n)
  one <- function(i) {
    set.seed(seeds[i]); y <- rpois(length(mu0), mu0)
    stat <- scan_k(ell, y, baseline, k_grid, degree)$best$deltaD
    emp_p(sum(ref_null >= stat), null_n)
  }
  pvals <- unlist(audit_lapply(seq_len(calibration_n), one, parallel = parallel))

  alphas <- c(0.10, 0.05, 0.01)
  rows <- lapply(alphas, function(a) {
    obs <- sum(pvals <= a); n <- length(pvals)
    ci <- binom_ci(obs, n)
    data.frame(
      nominal_alpha = a, n_outer_sims = n, inner_null_n = null_n,
      observed_count = obs, expected_count = a * n,
      observed_fpr = obs / n,
      ci_lo = ci["lower"], ci_hi = ci["upper"],
      calibration_ratio = (obs / n) / a,
      mc_uncertainty = sqrt(a * (1 - a) / n),
      compatible = (a >= ci["lower"] && a <= ci["upper"]),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

plot_calibration <- function(cal) {
  ggplot2::ggplot(cal, ggplot2::aes(nominal_alpha, observed_fpr)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), width = 0.004,
                           colour = species_colour("Fe")) +
    ggplot2::geom_point(size = 3, colour = species_colour("Fe")) +
    ggplot2::scale_x_continuous(limits = c(0, 0.12)) +
    ggplot2::scale_y_continuous(limits = c(0, max(0.12, max(cal$ci_hi, na.rm = TRUE)))) +
    ggplot2::labs(
      title = "Synthetic-null calibration of the scan-global procedure",
      subtitle = "Observed false-positive rate vs nominal alpha (exact binomial CI). Dashed = ideal.",
      x = "nominal alpha", y = "observed false-positive rate",
      caption = "Compatibility judged by whether nominal alpha lies within the Monte Carlo CI.") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  calibration_n <- 2000L; null_n <- 2000L
  i <- which(argv == "--calibration-n"); if (length(i) == 1L) calibration_n <- as.integer(argv[i + 1L])
  j <- which(argv == "--null-n"); if (length(j) == 1L) null_n <- as.integer(argv[j + 1L])
  parallel <- any(argv == "--parallel")
  cfg <- default_audit_config()
  if (parallel) configure_parallel(TRUE, "auto")
  setup_rng(cfg$seed)
  root <- audit_repo_root()
  cal <- calibrate(cfg, calibration_n, null_n, parallel)
  write_table(cal, file.path(root, "tables_r/statistical_audit/null_calibration.csv"))
  save_fig(plot_calibration(cal), file.path(root, "figures_r/statistical_audit/null_calibration.png"))
  save_fig(plot_calibration(cal), file.path(root, "figures_r/statistical_audit/fig13_null_calibration.png"))
  cat("[calibration] observed FPR vs nominal:\n")
  print(cal[, c("nominal_alpha", "observed_fpr", "ci_lo", "ci_hi", "compatible")])
  invisible(cal)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("calibrate_false_positive_rate\\.R$", .invoked_file)) {
  main()
}
