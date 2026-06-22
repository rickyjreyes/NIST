#!/usr/bin/env Rscript
# run_injection_recovery.R
# ---------------------------------------------------------------------------
# Injection-recovery power analysis. A known log-periodic mode is injected into
# synthetic data built from the fitted Fe II smooth baseline:
#
#   mu_true = mu0 * exp(A * cos(k_inj * ell - phi0))
#   y* ~ Poisson(mu_true)
#
# Amplitude grid: 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10.
# Injected frequencies: Fe reference peak, a lower control, a higher control,
# and an off-grid value between target windings.
#
# For each (amplitude, frequency) cell we report, with binomial CIs:
#   - detection-anywhere probability (global p <= 0.05 vs the A=0 reference null)
#   - correct-region recovery probability (k_best within tol of k_inj)
#   - globally-significant correct-region probability (both)
#   - median recovered k and recovered-amplitude bias
#   - false localisation to another peak.
#
# Writes:
#   tables_r/statistical_audit/injection_recovery.csv
#   figures_r/statistical_audit/injection_recovery_power.png (+ fig14)
#   figures_r/statistical_audit/injection_frequency_bias.png
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

run_injection <- function(cfg, inject_n, null_n, alpha = 0.05, parallel = FALSE) {
  root <- audit_repo_root()
  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              "Fe", 2L, "wavenumber")$lines
  k_grid <- audit_k_grid(cfg)
  ref <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)
  mu0 <- ref$mu0; ell <- ref$ell; baseline <- ref$baseline; degree <- cfg$degree
  k_ref <- ref$best$k_best

  # reference null under A = 0 (shared across cells)
  nd <- null_distribution(ell, ref$y, baseline, mu0, k_grid, degree, null_n, parallel = parallel)
  ref_null <- nd$max_vals
  thresh_stat <- stats::quantile(ref_null, 1 - alpha)   # detection threshold

  freqs <- c(fe_reference = k_ref, lower_control = 12.0,
             higher_control = 55.0, off_grid = k_ref * 1.013)
  amps <- cfg$injection_amplitudes
  phi0 <- 0.6
  tol <- cfg$tol_primary

  grid <- expand.grid(freq_name = names(freqs), amplitude = amps,
                      stringsAsFactors = FALSE)
  cell <- function(i) {
    fn <- grid$freq_name[i]; A <- grid$amplitude[i]; k_inj <- freqs[[fn]]
    seeds <- sample.int(.Machine$integer.max, inject_n)
    rec <- lapply(seq_len(inject_n), function(r) {
      set.seed(seeds[r])
      mu_true <- pmax(mu0 * exp(A * cos(k_inj * ell - phi0)), EPS)
      y <- rpois(length(mu_true), mu_true)
      sk <- scan_k(ell, y, baseline, k_grid, degree)
      kb <- sk$best$k_best
      c(stat = sk$best$deltaD, kb = kb, amp = sk$best$amplitude,
        detect = as.numeric(sk$best$deltaD >= thresh_stat),
        correct = as.numeric(abs(kb - k_inj) / k_inj <= tol))
    })
    M <- do.call(rbind, rec)
    n <- nrow(M)
    det <- sum(M[, "detect"]); cor <- sum(M[, "correct"])
    sig_cor <- sum(M[, "detect"] == 1 & M[, "correct"] == 1)
    false_loc <- sum(M[, "detect"] == 1 & M[, "correct"] == 0)
    ci_det <- binom_ci(det, n); ci_cor <- binom_ci(sig_cor, n)
    data.frame(
      freq_name = fn, k_injected = k_inj, amplitude = A, n_sims = n,
      detection_prob = det / n, det_ci_lo = ci_det["lower"], det_ci_hi = ci_det["upper"],
      correct_region_prob = cor / n,
      sig_correct_prob = sig_cor / n, sigcor_ci_lo = ci_cor["lower"], sigcor_ci_hi = ci_cor["upper"],
      false_localisation_prob = false_loc / n,
      median_recovered_k = stats::median(M[, "kb"]),
      freq_bias = stats::median(M[, "kb"]) - k_inj,
      amplitude_bias = stats::median(M[, "amp"]) - A,
      stringsAsFactors = FALSE)
  }
  out <- audit_lapply(seq_len(nrow(grid)), cell, parallel = parallel)
  res <- do.call(rbind, out)
  res$alpha <- alpha; res$null_n <- null_n
  res
}

plot_power <- function(inj) {
  ggplot2::ggplot(inj, ggplot2::aes(amplitude, sig_correct_prob, colour = freq_name)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = sigcor_ci_lo, ymax = sigcor_ci_hi, fill = freq_name),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = 3, colour = "grey40") +
    viridis::scale_colour_viridis(discrete = TRUE, end = 0.85, name = "injected freq") +
    viridis::scale_fill_viridis(discrete = TRUE, end = 0.85, guide = "none") +
    ggplot2::labs(
      title = "Injection-recovery power: globally-significant correct-region detection",
      subtitle = "Probability of detecting AND localising the injected mode vs injected amplitude.",
      x = "injected amplitude", y = "power", caption = "Dotted line = 0.8 power reference.") +
    theme_audit()
}

plot_freq_bias <- function(inj) {
  d <- inj[inj$amplitude > 0, ]
  ggplot2::ggplot(d, ggplot2::aes(amplitude, freq_bias, colour = freq_name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
    viridis::scale_colour_viridis(discrete = TRUE, end = 0.85, name = "injected freq") +
    ggplot2::labs(title = "Injection-recovery frequency bias",
                  subtitle = "Median recovered k minus injected k vs amplitude.",
                  x = "injected amplitude", y = "frequency bias (recovered - injected)") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  inject_n <- 500L; null_n <- 1000L
  i <- which(argv == "--injection-n"); if (length(i) == 1L) inject_n <- as.integer(argv[i + 1L])
  j <- which(argv == "--null-n"); if (length(j) == 1L) null_n <- as.integer(argv[j + 1L])
  parallel <- any(argv == "--parallel")
  cfg <- default_audit_config(fast = is_fast(argv))
  if (parallel) configure_parallel(TRUE, "auto")
  setup_rng(cfg$seed)
  root <- audit_repo_root()
  inj <- run_injection(cfg, inject_n, null_n, parallel = parallel)
  write_table(inj, file.path(root, "tables_r/statistical_audit/injection_recovery.csv"))
  save_fig(plot_power(inj), file.path(root, "figures_r/statistical_audit/injection_recovery_power.png"))
  save_fig(plot_power(inj), file.path(root, "figures_r/statistical_audit/fig14_injection_recovery_power.png"))
  save_fig(plot_freq_bias(inj), file.path(root, "figures_r/statistical_audit/injection_frequency_bias.png"))
  cat("[injection] power at top amplitude (fe_reference):\n")
  fe <- inj[inj$freq_name == "fe_reference", ]
  print(fe[, c("amplitude", "detection_prob", "sig_correct_prob")])
  invisible(inj)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_injection_recovery\\.R$", .invoked_file)) {
  main()
}
