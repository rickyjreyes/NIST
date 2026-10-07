#!/usr/bin/env Rscript
# run_resolution_mode_diagnostics.R
# ---------------------------------------------------------------------------
# FOLLOW-UP diagnostic for the Fe II resolution-dependent mode transition.
#
# This module is deliberately separate from run_final_adversarial_audit.R.
# It MUST NOT overwrite, rescue, or reinterpret the frozen adversarial verdicts.
# In particular, the original 60:240 bin-grid stability result remains the
# primary robustness test and its pass/fail criterion is unchanged.
#
# Questions addressed here:
#   1. Does the observed k~9.6 -> k~31.3 transition persist when Gaussian
#      smoothing is held approximately fixed in physical ell=log(wavenumber)
#      width rather than fixed in bin units?
#   2. Are k=9.602325620315224 and the canonical Fe primary k simultaneously
#      present as competing local modes even when only one is the global winner?
#   3. Does the previously frozen GWTC k=9.602325620315224 survive a direct
#      fixed-frequency null test in Fe and neighbouring ion-II datasets?
#
# The GWTC frequency existed before this NIST follow-up, but noticing its match
# to the coarse-bin NIST branch is post-hoc. Results are therefore exploratory
# cross-domain follow-up evidence, not original preregistration.
#
# Writes:
#   tables_r/statistical_audit/resolution_mode_dense.csv
#   tables_r/statistical_audit/resolution_mode_top_peaks.csv
#   tables_r/statistical_audit/external_gwtc_k_fixed_tests.csv
#   figures_r/statistical_audit/resolution_mode_transition.png
#   figures_r/statistical_audit/resolution_mode_strength.png
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) {
  source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this),
                   "audit_utils.R"))
}

GWTC_K_FROZEN <- 9.602325620315224

matched_sigma <- function(bins, reference_bins = 160L, reference_sigma = 6.0) {
  as.numeric(reference_sigma) * as.numeric(bins) / as.numeric(reference_bins)
}

local_peak_indices <- function(deltaD) {
  x <- as.numeric(deltaD)
  n <- length(x)
  if (n == 0L) return(integer())
  if (n == 1L) return(1L)
  x[!is.finite(x)] <- -Inf
  idx <- integer()
  if (x[1] > x[2]) idx <- c(idx, 1L)
  if (n > 2L) {
    mid <- which(x[2:(n - 1L)] >= x[1:(n - 2L)] &
                   x[2:(n - 1L)] > x[3:n]) + 1L
    idx <- c(idx, mid)
  }
  if (x[n] >= x[n - 1L]) idx <- c(idx, n)
  unique(idx[is.finite(x[idx])])
}

extract_top_peaks <- function(scan, n_top = 5L) {
  idx <- local_peak_indices(scan$deltaD)
  if (length(idx) == 0L) idx <- which.max(scan$deltaD)
  out <- data.frame(
    k = scan$k[idx],
    deltaD = scan$deltaD[idx],
    stringsAsFactors = FALSE
  )
  out <- out[order(out$deltaD, decreasing = TRUE), , drop = FALSE]
  out <- head(out, max(1L, as.integer(n_top)))
  out$peak_rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

classify_peak_branch <- function(k, k_external = GWTC_K_FROZEN, k_primary,
                                 tol_relative = 0.02) {
  if (!is.finite(k)) return("other")
  if (abs(k - k_external) / k_external <= tol_relative) return("external_gwtc_9p602")
  if (abs(k - k_primary) / k_primary <= tol_relative) return("primary_fe")
  "other"
}

scan_value_at_k <- function(res, k_target) {
  idx <- nearest_k_index(res$scan$k, k_target)
  c(k_grid = res$scan$k[idx], deltaD = res$scan$deltaD[idx])
}

resolution_scan_one <- function(lines, bins, sigma_bins, smoothing_mode, cfg,
                                k_primary, k_external = GWTC_K_FROZEN,
                                n_top = 5L) {
  k_grid <- audit_k_grid(cfg)
  res <- run_scan_analysis(lines, bins, k_grid, cfg$degree, sigma_bins)
  pext <- scan_value_at_k(res, k_external)
  ppri <- scan_value_at_k(res, k_primary)
  bin_width <- stats::median(res$binned$edge_hi - res$binned$edge_lo)
  best <- res$best$deltaD
  summary <- data.frame(
    smoothing_mode = smoothing_mode,
    bins = as.integer(bins),
    sigma_bins = as.numeric(sigma_bins),
    sigma_ell_approx = as.numeric(sigma_bins) * bin_width,
    n_lines = res$n_lines,
    k_best = res$best$k_best,
    deltaD_best = best,
    best_branch = classify_peak_branch(res$best$k_best, k_external, k_primary,
                                       cfg$tol_primary),
    k_external_target = k_external,
    k_external_grid = unname(pext["k_grid"]),
    deltaD_external = unname(pext["deltaD"]),
    external_to_best = unname(pext["deltaD"]) / best,
    k_primary_target = k_primary,
    k_primary_grid = unname(ppri["k_grid"]),
    deltaD_primary = unname(ppri["deltaD"]),
    primary_to_best = unname(ppri["deltaD"]) / best,
    stringsAsFactors = FALSE
  )
  peaks <- extract_top_peaks(res$scan, n_top)
  peaks$smoothing_mode <- smoothing_mode
  peaks$bins <- as.integer(bins)
  peaks$sigma_bins <- as.numeric(sigma_bins)
  peaks$sigma_ell_approx <- as.numeric(sigma_bins) * bin_width
  peaks$branch <- vapply(peaks$k, classify_peak_branch, character(1),
                         k_external = k_external, k_primary = k_primary,
                         tol_relative = cfg$tol_primary)
  peaks <- peaks[, c("smoothing_mode", "bins", "sigma_bins", "sigma_ell_approx",
                     "peak_rank", "k", "deltaD", "branch")]
  list(summary = summary, peaks = peaks)
}

run_dense_resolution_map <- function(lines, cfg, k_primary,
                                     bins_dense = seq(40L, 300L, by = 5L),
                                     n_top = 5L, parallel = FALSE) {
  jobs <- do.call(rbind, lapply(as.integer(bins_dense), function(b) {
    data.frame(
      bins = c(b, b),
      smoothing_mode = c("fixed_sigma_bins", "matched_sigma_ell"),
      sigma_bins = c(cfg$baseline_sigma,
                     matched_sigma(b, cfg$bins_primary, cfg$baseline_sigma)),
      stringsAsFactors = FALSE
    )
  }))
  one <- function(i) {
    g <- jobs[i, ]
    resolution_scan_one(lines, g$bins, g$sigma_bins, g$smoothing_mode,
                        cfg, k_primary, GWTC_K_FROZEN, n_top)
  }
  ans <- audit_lapply(seq_len(nrow(jobs)), one, parallel = parallel, seed = cfg$seed)
  list(
    summary = do.call(rbind, lapply(ans, `[[`, "summary")),
    peaks = do.call(rbind, lapply(ans, `[[`, "peaks"))
  )
}

fixed_frequency_null <- function(res, k_target, B, seed,
                                 refit_baseline = FALSE, parallel = FALSE) {
  if (B <= 0L) return(list(observed = NA_real_, p = NA_real_, tail = NA_integer_))
  observed <- scan_k(res$ell, res$y, res$baseline, c(k_target), res$degree)$best$deltaD
  setup_rng(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(i) {
    set.seed(seeds[i])
    y0 <- stats::rpois(length(res$mu0), res$mu0)
    baseline0 <- if (isTRUE(refit_baseline)) {
      pmax(gaussian_filter_nearest(y0, res$baseline_sigma), EPS)
    } else {
      res$baseline
    }
    scan_k(res$ell, y0, baseline0, c(k_target), res$degree)$best$deltaD
  }
  sims <- unlist(audit_lapply(seq_len(B), one, parallel = parallel, seed = seed),
                 use.names = FALSE)
  tail <- sum(sims >= observed)
  list(observed = observed, p = emp_p(tail, B), tail = as.integer(tail))
}

fixed_external_target_one <- function(raw, species, source, cfg, B,
                                      parallel = FALSE) {
  lines <- clean_lines_source(raw, species, 2L, source)$lines
  k_grid <- audit_k_grid(cfg)
  res <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree,
                           cfg$baseline_sigma)
  fixed <- fixed_frequency_null(res, GWTC_K_FROZEN, B, cfg$seed,
                                refit_baseline = FALSE, parallel = parallel)
  refit <- fixed_frequency_null(res, GWTC_K_FROZEN, B, cfg$seed + 101L,
                                refit_baseline = TRUE, parallel = parallel)
  data.frame(
    species = species,
    ion = 2L,
    source = source,
    bins = cfg$bins_primary,
    sigma = cfg$baseline_sigma,
    degree = cfg$degree,
    n_lines = nrow(lines),
    external_target = "GWTC_frozen_k",
    k_fixed = GWTC_K_FROZEN,
    k_selected_full_scan = res$best$k_best,
    selected_deltaD = res$best$deltaD,
    fixed_deltaD = fixed$observed,
    fixed_to_selected_ratio = fixed$observed / res$best$deltaD,
    fixed_baseline_p = fixed$p,
    fixed_baseline_tail = fixed$tail,
    refit_baseline_p = refit$p,
    refit_baseline_tail = refit$tail,
    conservative_fixed_p = max(fixed$p, refit$p, na.rm = TRUE),
    null_n = as.integer(B),
    posthoc_cross_domain_followup = TRUE,
    stringsAsFactors = FALSE
  )
}

run_external_target_tests <- function(root, cfg, B = 5000L, parallel = FALSE) {
  species <- c("Fe", "Cr", "Mn", "Co", "Ni", "Ti")
  rows <- list()
  j <- 1L
  for (sp in species) {
    path <- file.path(root, "data", sprintf("%s_lines.csv", sp))
    raw <- read_nist_csv(path)
    sources <- if (sp == "Fe") c("wavenumber", "observed", "ritz") else "wavenumber"
    for (src in sources) {
      rows[[j]] <- fixed_external_target_one(raw, sp, src, cfg, B, parallel)
      cat(sprintf("[external-k] %s/%s k=%.10f p_fixed=%.4g p_refit=%.4g p_cons=%.4g\n",
                  sp, src, GWTC_K_FROZEN,
                  rows[[j]]$fixed_baseline_p,
                  rows[[j]]$refit_baseline_p,
                  rows[[j]]$conservative_fixed_p))
      j <- j + 1L
    }
  }
  do.call(rbind, rows)
}

plot_resolution_transition <- function(dense, k_primary) {
  ggplot2::ggplot(dense, ggplot2::aes(bins, k_best, linetype = smoothing_mode)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(yintercept = GWTC_K_FROZEN, linetype = 3) +
    ggplot2::geom_hline(yintercept = k_primary, linetype = 2) +
    ggplot2::labs(
      title = "Fe II resolution-dependent selected mode",
      subtitle = "Original fixed-sigma-bin pipeline versus approximately fixed ell-space smoothing width.",
      x = "bin count", y = "selected k (rad / ln cm^-1)",
      linetype = "smoothing") +
    theme_audit()
}

plot_mode_strength <- function(dense) {
  long <- rbind(
    data.frame(smoothing_mode = dense$smoothing_mode, bins = dense$bins,
               target = "GWTC k=9.6023256", ratio = dense$external_to_best),
    data.frame(smoothing_mode = dense$smoothing_mode, bins = dense$bins,
               target = "Fe primary k", ratio = dense$primary_to_best)
  )
  ggplot2::ggplot(long, ggplot2::aes(bins, ratio, linetype = target)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~smoothing_mode) +
    ggplot2::geom_hline(yintercept = 1, linetype = 3) +
    ggplot2::labs(
      title = "Strength of the two tracked modes relative to the winning peak",
      subtitle = "A ratio near 1 means the tracked mode is competitive with the global winner.",
      x = "bin count", y = "deltaD(target) / deltaD(best)", linetype = "target") +
    theme_audit()
}

arg_int <- function(argv, flag, default) {
  i <- which(argv == flag)
  if (length(i) == 0L) return(as.integer(default))
  as.integer(argv[i[1] + 1L])
}

main <- function(argv = commandArgs(TRUE)) {
  parallel <- any(argv == "--parallel")
  fast <- is_fast(argv)
  if (parallel) configure_parallel(TRUE, "auto")
  cfg <- default_audit_config(fast = fast)
  root <- audit_repo_root()

  min_bin <- arg_int(argv, "--min-bin", if (fast) 60L else 40L)
  max_bin <- arg_int(argv, "--max-bin", if (fast) 220L else 300L)
  step_bin <- arg_int(argv, "--step-bin", if (fast) 20L else 5L)
  fixed_null_n <- arg_int(argv, "--fixed-null-n", if (fast) 25L else 5000L)
  bins_dense <- seq(min_bin, max_bin, by = step_bin)

  raw_fe <- read_nist_csv(file.path(root, "data/Fe_lines.csv"))
  lines_fe <- clean_lines_source(raw_fe, "Fe", 2L, "wavenumber")$lines
  k_primary <- run_scan_analysis(lines_fe, cfg$bins_primary, audit_k_grid(cfg),
                                 cfg$degree, cfg$baseline_sigma)$best$k_best

  cat(sprintf("[resolution] Fe primary k=%.10f; frozen external GWTC k=%.10f\n",
              k_primary, GWTC_K_FROZEN))
  dense <- run_dense_resolution_map(lines_fe, cfg, k_primary, bins_dense,
                                    n_top = 5L, parallel = parallel)
  write_table(dense$summary,
              file.path(root, "tables_r/statistical_audit/resolution_mode_dense.csv"))
  write_table(dense$peaks,
              file.path(root, "tables_r/statistical_audit/resolution_mode_top_peaks.csv"))
  save_fig(plot_resolution_transition(dense$summary, k_primary),
           file.path(root, "figures_r/statistical_audit/resolution_mode_transition.png"))
  save_fig(plot_mode_strength(dense$summary),
           file.path(root, "figures_r/statistical_audit/resolution_mode_strength.png"))

  ext <- run_external_target_tests(root, cfg, fixed_null_n, parallel)
  write_table(ext,
              file.path(root, "tables_r/statistical_audit/external_gwtc_k_fixed_tests.csv"))

  cat(sprintf("[resolution] wrote %d dense resolution rows and %d top-peak rows\n",
              nrow(dense$summary), nrow(dense$peaks)))
  cat(sprintf("[resolution] external fixed-k tests=%d, B=%d; results are exploratory follow-up\n",
              nrow(ext), fixed_null_n))
  invisible(list(dense = dense, external = ext, k_primary = k_primary))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_resolution_mode_diagnostics\\.R$", .invoked_file)) {
  main()
}
