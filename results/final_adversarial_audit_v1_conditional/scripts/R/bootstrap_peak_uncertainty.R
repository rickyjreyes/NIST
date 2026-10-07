#!/usr/bin/env Rscript
# bootstrap_peak_uncertainty.R
# ---------------------------------------------------------------------------
# Bootstrap uncertainty and peak estimates for the canonical Fe II scan.
#
# Two RESAMPLING SCHEMES are kept strictly distinct:
#
#   A. Parametric NULL bootstrap  (null_distribution() in audit_utils.R)
#      - y0 ~ Poisson(mu0) under the fitted smooth null
#      - used ONLY for significance / calibration, NOT for the CI of the
#        alternative-model peak.
#
#   B. DATA bootstrap (this module, data_bootstrap_peaks())
#      - resample the retained unique lines WITH replacement, rebin, rebuild
#        the baseline, and rerun the full scan
#      - preserves the line-position data-generating structure
#      - used for uncertainty and stability of the detected peak.
#
# Writes:
#   tables_r/statistical_audit/peak_estimates.csv
#   tables_r/statistical_audit/bootstrap_peak_draws.csv
#   tables_r/statistical_audit/peak_confidence_intervals.csv
#   figures_r/statistical_audit/bootstrap_peak_distributions.png
#   figures_r/statistical_audit/fig05_bootstrap_peak_distribution.png
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

# peak_estimate_row(): peak-derived quantities for one reference analysis.
peak_estimate_row <- function(res, analysis_id) {
  best <- res$best
  scan <- res$scan
  k <- best$k_best
  # separation from the second-best LOCAL peak
  dist <- max(1L, length(scan$deltaD) %/% 80L)
  pk_idx <- find_peaks_distance(scan$deltaD, dist) + 1L
  if (length(pk_idx) == 0L) pk_idx <- which.max(scan$deltaD)
  pk <- scan[pk_idx, , drop = FALSE]
  pk <- pk[order(-pk$deltaD), , drop = FALSE]
  k_best_idx <- which.min(abs(pk$k - k))
  second_k <- if (nrow(pk) >= 2L) pk$k[setdiff(seq_len(nrow(pk)), k_best_idx)][1] else NA_real_
  data.frame(
    analysis_id = analysis_id,
    k_best = k,
    omega_best = k,                       # angular-frequency convention
    delta_log_x = delta_log_x(k),
    scale_ratio = scale_ratio(k),
    amplitude = best$amplitude,
    phase = best$phase,
    deltaD = best$deltaD,
    D_base = best$D_base,
    D_harmonic = best$D_harmonic,
    n_obs = res$n_obs,
    n_lines = res$n_lines,
    bins = res$bins,
    peak_rank = 1L,
    second_best_k = second_k,
    separation_from_second = if (is.na(second_k)) NA_real_ else abs(k - second_k),
    stringsAsFactors = FALSE
  )
}

# data_bootstrap_peaks(): scheme B. Returns one row per resample.
data_bootstrap_peaks <- function(lines, cfg, k_ref, bins = cfg$bins_primary,
                                 B = cfg$bootstrap_n, parallel = FALSE) {
  k_grid <- audit_k_grid(cfg)
  n <- nrow(lines)
  k_lo <- cfg$k_min; k_hi <- cfg$k_max
  seeds <- sample.int(.Machine$integer.max, B)  # per-draw seeds (worker-independent)
  one <- function(b) {
    set.seed(seeds[b])
    idx <- sample.int(n, n, replace = TRUE)
    bl <- lines[idx, , drop = FALSE]
    # rebin on resampled positions; sort needed for histogram range
    bl <- bl[order(bl$ell), , drop = FALSE]
    res <- tryCatch(run_scan_analysis(bl, bins, k_grid, cfg$degree, cfg$baseline_sigma),
                    error = function(e) NULL)
    if (is.null(res)) {
      return(data.frame(draw = b, k_best = NA_real_, deltaD = NA_real_,
                        amplitude = NA_real_, phase = NA_real_,
                        delta_log_x = NA_real_, scale_ratio = NA_real_,
                        converged = FALSE, boundary = NA, in_region = NA,
                        stringsAsFactors = FALSE))
    }
    kb <- res$best$k_best
    boundary <- (abs(kb - k_lo) < 1e-9) || (abs(kb - k_hi) < 1e-9)
    data.frame(
      draw = b, k_best = kb, deltaD = res$best$deltaD,
      amplitude = res$best$amplitude, phase = res$best$phase,
      delta_log_x = delta_log_x(kb), scale_ratio = scale_ratio(kb),
      converged = TRUE, boundary = boundary,
      in_region = in_peak_region(kb, k_ref, "relative", cfg$tol_primary),
      stringsAsFactors = FALSE)
  }
  draws <- audit_lapply(seq_len(B), one, parallel = parallel)
  out <- do.call(rbind, draws)
  rownames(out) <- NULL
  out
}

# summarise_ci(): percentile summary of a bootstrap vector, with multimodality
# diagnostics. Returns a one-row data.frame.
summarise_ci <- function(x, observed, quantity) {
  v <- x[is.finite(x)]
  n_valid <- length(v)
  if (n_valid < 2L) {
    return(data.frame(quantity = quantity, observed = observed,
                      median = NA, mean = NA, sd = NA,
                      ci_lo = NA, ci_hi = NA, iqr = NA, bias = NA,
                      n_valid = n_valid, multimodal = NA, stringsAsFactors = FALSE))
  }
  # crude multimodality flag: > 1 well-separated density mode
  dens <- stats::density(v)
  pk <- find_peaks_distance(dens$y, max(1L, length(dens$y) %/% 20L))
  big <- sum(dens$y[pk + 1L] > 0.15 * max(dens$y))
  data.frame(
    quantity = quantity, observed = observed,
    median = stats::median(v), mean = mean(v), sd = stats::sd(v),
    ci_lo = unname(stats::quantile(v, 0.025)),
    ci_hi = unname(stats::quantile(v, 0.975)),
    iqr = stats::IQR(v),
    bias = mean(v) - observed,
    n_valid = n_valid,
    multimodal = big > 1L,
    stringsAsFactors = FALSE)
}

# plot_bootstrap_distributions(): faceted histograms of the bootstrap draws.
plot_bootstrap_distributions <- function(draws, ref) {
  d <- draws[draws$converged & is.finite(draws$k_best), , drop = FALSE]
  long <- rbind(
    data.frame(quantity = "k_best (rad/ell)", value = d$k_best),
    data.frame(quantity = "log-period 2pi/k", value = d$delta_log_x),
    data.frame(quantity = "scale ratio exp(2pi/k)", value = d$scale_ratio),
    data.frame(quantity = "amplitude", value = d$amplitude)
  )
  vlines <- data.frame(
    quantity = c("k_best (rad/ell)", "log-period 2pi/k",
                 "scale ratio exp(2pi/k)", "amplitude"),
    value = c(ref$k_best, ref$delta_log_x, ref$scale_ratio, ref$amplitude))
  ggplot2::ggplot(long, ggplot2::aes(value)) +
    ggplot2::geom_histogram(bins = 50, fill = species_colour("Fe"), colour = "white",
                            linewidth = 0.1) +
    ggplot2::geom_vline(data = vlines, ggplot2::aes(xintercept = value),
                        linetype = 2, colour = "grey20") +
    ggplot2::facet_wrap(~quantity, scales = "free") +
    ggplot2::labs(
      title = "Fe II data-bootstrap peak distributions",
      subtitle = sprintf("%d resamples; dashed line = reference estimate. Scheme B (data resampling).",
                         nrow(d)),
      x = NULL, y = "count") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  B <- NULL
  if (length(argv) > 0L) {
    i <- which(argv == "--bootstrap-n")
    if (length(i) == 1L) B <- as.integer(argv[i + 1L])
  }
  cfg <- default_audit_config(bootstrap_n = if (is.null(B)) 2000L else B)
  setup_rng(cfg$seed)
  root <- audit_repo_root()

  raw <- read_nist_csv(file.path(root, "data/Fe_lines.csv"))
  cl <- clean_lines_source(raw, "Fe", 2L, "wavenumber")
  k_grid <- audit_k_grid(cfg)
  ref <- run_scan_analysis(cl$lines, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)
  k_ref <- ref$best$k_best

  # --- peak estimates -------------------------------------------------------
  pe <- peak_estimate_row(ref, sprintf("fe_ion2_wn_bin%d", cfg$bins_primary))
  write_table(pe, file.path(root, "tables_r/statistical_audit/peak_estimates.csv"))

  # --- data bootstrap (scheme B) -------------------------------------------
  cat(sprintf("[bootstrap] data resampling B=%d (scheme B)\n", cfg$bootstrap_n))
  draws <- data_bootstrap_peaks(cl$lines, cfg, k_ref, parallel = FALSE)
  write_table(draws, file.path(root, "tables_r/statistical_audit/bootstrap_peak_draws.csv"))

  # --- confidence intervals -------------------------------------------------
  ci <- rbind(
    summarise_ci(draws$k_best,      pe$k_best,      "k_best"),
    summarise_ci(draws$delta_log_x, pe$delta_log_x, "delta_log_x"),
    summarise_ci(draws$scale_ratio, pe$scale_ratio, "scale_ratio"),
    summarise_ci(draws$amplitude,   pe$amplitude,   "amplitude"),
    summarise_ci(draws$phase,       pe$phase,       "phase"),
    summarise_ci(draws$deltaD,      pe$deltaD,      "deltaD")
  )
  ci$n_failed <- sum(!draws$converged)
  ci$failure_rate <- mean(!draws$converged)
  write_table(ci, file.path(root, "tables_r/statistical_audit/peak_confidence_intervals.csv"))

  p <- plot_bootstrap_distributions(draws, pe)
  save_fig(p, file.path(root, "figures_r/statistical_audit/bootstrap_peak_distributions.png"))
  save_fig(p, file.path(root, "figures_r/statistical_audit/fig05_bootstrap_peak_distribution.png"))

  # fig02: canonical Fe II scan curve with the selected peak and bootstrap k CI
  kci <- ci[ci$quantity == "k_best", ]
  scan_df <- ref$scan
  p2 <- ggplot2::ggplot(scan_df, ggplot2::aes(k, deltaD)) +
    ggplot2::annotate("rect", xmin = kci$ci_lo, xmax = kci$ci_hi, ymin = -Inf, ymax = Inf,
                      fill = species_colour("Fe"), alpha = 0.15) +
    ggplot2::geom_line(colour = "grey25", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = pe$k_best, linetype = 2, colour = species_colour("Fe")) +
    ggplot2::labs(
      title = "Fe II log-cosine scan curve and peak stability",
      subtitle = sprintf("Selected k=%.3f (log-period 2pi/k=%.4f, scale ratio %.4f); shaded = bootstrap 95%% CI for k.",
                         pe$k_best, pe$delta_log_x, pe$scale_ratio),
      x = "k (rad / ln cm^-1)", y = "deltaD = D_base - D_harmonic",
      caption = "Reference analysis at the primary bin count.") +
    theme_audit()
  save_fig(p2, file.path(root, "figures_r/statistical_audit/fig02_feii_scan_stability.png"))

  cat(sprintf("[bootstrap] valid=%d/%d  in-region(%.0f%% tol)=%.1f%%  multimodal-k=%s\n",
              sum(draws$converged), nrow(draws), 100 * cfg$tol_primary,
              100 * mean(draws$in_region, na.rm = TRUE),
              ci$multimodal[ci$quantity == "k_best"]))
  invisible(list(peak = pe, ci = ci, draws = draws))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("bootstrap_peak_uncertainty\\.R$", .invoked_file)) {
  main()
}
