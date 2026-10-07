#!/usr/bin/env Rscript
# run_resolution_invariant_protocol_v1.R
# ---------------------------------------------------------------------------
# Prospective runner for config/resolution_invariant_protocol_v1.json.
#
# This analysis is intentionally separate from the canonical final audit. It
# cannot overwrite the historical fixed-sigma result or the 120/160/200 gate.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
source(file.path(R_DIR, "audit_utils.R"))

.diag <- new.env(parent = environment())
sys.source(file.path(R_DIR, "run_resolution_mode_diagnostics.R"), envir = .diag)

load_protocol <- function(root) {
  p <- file.path(root, "config", "resolution_invariant_protocol_v1.json")
  if (!file.exists(p)) stop("missing frozen protocol: ", p)
  jsonlite::fromJSON(p, simplifyVector = TRUE)
}

protocol_sigma <- function(bins, protocol) {
  as.numeric(protocol$smoothing$reference_sigma_bins) * as.numeric(bins) /
    as.numeric(protocol$smoothing$reference_bins)
}

one_resolution_observed <- function(lines, bins, protocol) {
  cfg <- default_audit_config(fast = FALSE)
  cfg$k_min <- protocol$scan$k_min
  cfg$k_max <- protocol$scan$k_max
  cfg$n_k <- protocol$scan$n_k
  cfg$degree <- protocol$scan$degree
  k_grid <- audit_k_grid(cfg)
  sigma <- protocol_sigma(bins, protocol)
  res <- run_scan_analysis(lines, bins, k_grid, cfg$degree, sigma)
  target <- as.numeric(protocol$target_k)
  target_fit <- scan_k(res$ell, res$y, res$baseline, c(target), cfg$degree)$best
  peaks <- .diag$extract_top_peaks(res$scan, protocol$scan$top_local_maxima_to_track)
  tol <- as.numeric(protocol$target_relative_tolerance)
  winner_ok <- abs(res$best$k_best - target) / target <= tol
  top5_ok <- any(abs(peaks$k - target) / target <= tol)
  target_rank <- {
    z <- which(abs(peaks$k - target) / target <= tol)
    if (length(z)) min(peaks$peak_rank[z]) else NA_integer_
  }
  list(
    res = res,
    target_fit = target_fit,
    peaks = peaks,
    row = data.frame(
      bins = as.integer(bins),
      sigma_bins = sigma,
      n_lines = res$n_lines,
      k_best = res$best$k_best,
      deltaD_best = res$best$deltaD,
      best_amplitude = res$best$amplitude,
      best_phase = res$best$phase,
      target_k = target,
      target_deltaD = target_fit$deltaD,
      target_amplitude = target_fit$amplitude,
      target_phase = target_fit$phase,
      target_to_best = target_fit$deltaD / res$best$deltaD,
      winner_in_target_region = winner_ok,
      target_in_top5 = top5_ok,
      target_rank = target_rank,
      stringsAsFactors = FALSE))
}

null_one_resolution <- function(obs, protocol, parallel = FALSE) {
  B <- as.integer(protocol$null$replicates_per_bin)
  seed <- as.integer(protocol$null$seed) + as.integer(obs$row$bins[1])
  target <- as.numeric(protocol$target_k)
  tol <- as.numeric(protocol$target_relative_tolerance)
  res <- obs$res
  sigma <- obs$row$sigma_bins[1]
  k_grid <- res$scan$k
  degree <- as.integer(protocol$scan$degree)

  setup_rng(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(i) {
    set.seed(seeds[i])
    y0 <- stats::rpois(length(res$mu0), res$mu0)
    baseline0 <- pmax(gaussian_filter_nearest(y0, sigma), EPS)
    sk <- scan_k(res$ell, y0, baseline0, k_grid, degree)
    c(max_deltaD = sk$best$deltaD,
      winner_k = sk$best$k_best,
      winner_target = as.numeric(abs(sk$best$k_best - target) / target <= tol))
  }
  M <- do.call(rbind, audit_lapply(seq_len(B), one, parallel = parallel, seed = seed))
  threshold <- unname(stats::quantile(M[, "max_deltaD"], 1 - protocol$decision_rules$alpha,
                                      names = FALSE, type = 8))
  tail <- sum(M[, "max_deltaD"] >= obs$row$deltaD_best[1])
  significant_target_lock <- M[, "max_deltaD"] >= threshold & M[, "winner_target"] == 1
  summary <- data.frame(
    bins = obs$row$bins[1],
    sigma_bins = sigma,
    null_n = B,
    observed_max_deltaD = obs$row$deltaD_best[1],
    global_tail = as.integer(tail),
    global_p = emp_p(tail, B),
    alpha_threshold_deltaD = threshold,
    unconditional_target_winner_rate = mean(M[, "winner_target"] == 1),
    significant_target_lock_rate = mean(significant_target_lock),
    stringsAsFactors = FALSE)
  list(summary = summary, draws = M, threshold = threshold)
}

inject_one_resolution <- function(obs, null, protocol, parallel = FALSE) {
  n <- as.integer(protocol$injection$replicates_per_bin_per_amplitude)
  ratios <- as.numeric(protocol$injection$amplitude_ratios)
  target <- as.numeric(protocol$target_k)
  tol <- as.numeric(protocol$target_relative_tolerance)
  degree <- as.integer(protocol$scan$degree)
  res <- obs$res
  sigma <- obs$row$sigma_bins[1]
  k_grid <- res$scan$k
  A_obs <- obs$row$target_amplitude[1]
  phi <- obs$row$target_phase[1]
  base_seed <- as.integer(protocol$injection$seed) + as.integer(obs$row$bins[1]) * 100L

  one_ratio <- function(j) {
    ratio <- ratios[j]
    A <- ratio * A_obs
    setup_rng(base_seed + j)
    seeds <- sample.int(.Machine$integer.max, n)
    one <- function(i) {
      set.seed(seeds[i])
      mu_true <- pmax(res$mu0 * exp(A * cos(target * res$ell - phi)), EPS)
      y <- stats::rpois(length(mu_true), mu_true)
      baseline <- pmax(gaussian_filter_nearest(y, sigma), EPS)
      sk <- scan_k(res$ell, y, baseline, k_grid, degree)
      correct <- abs(sk$best$k_best - target) / target <= tol
      c(stat = sk$best$deltaD,
        k_best = sk$best$k_best,
        amplitude = sk$best$amplitude,
        detect = as.numeric(sk$best$deltaD >= null$threshold),
        correct = as.numeric(correct))
    }
    M <- do.call(rbind, audit_lapply(seq_len(n), one, parallel = parallel,
                                    seed = base_seed + j))
    data.frame(
      bins = obs$row$bins[1],
      sigma_bins = sigma,
      amplitude_ratio = ratio,
      injected_amplitude = A,
      n_sims = n,
      detection_prob = mean(M[, "detect"] == 1),
      correct_region_prob = mean(M[, "correct"] == 1),
      sig_correct_prob = mean(M[, "detect"] == 1 & M[, "correct"] == 1),
      false_location_prob = mean(M[, "detect"] == 1 & M[, "correct"] == 0),
      median_recovered_k = stats::median(M[, "k_best"]),
      k_bias = stats::median(M[, "k_best"]) - target,
      amplitude_bias = stats::median(M[, "amplitude"]) - A,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, lapply(seq_along(ratios), one_ratio))
}

classify_protocol <- function(observed, nulls, injections, protocol) {
  dr <- protocol$decision_rules
  winner_retention <- mean(observed$winner_in_target_region)
  top5_retention <- mean(observed$target_in_top5)
  zero <- injections[abs(injections$amplitude_ratio) < 1e-12, , drop = FALSE]
  one <- injections[abs(injections$amplitude_ratio - 1) < 1e-12, , drop = FALSE]

  zero_ok <- nrow(zero) == nrow(observed) &&
    max(zero$sig_correct_prob, na.rm = TRUE) <= dr$zero_signal_false_lock_rate_max
  one_ok <- nrow(one) == nrow(observed) &&
    min(one$sig_correct_prob, na.rm = TRUE) >= dr$observed_amplitude_correct_recovery_min
  null_lock_ok <- max(nulls$significant_target_lock_rate, na.rm = TRUE) <=
    dr$zero_signal_false_lock_rate_max

  pass <- winner_retention >= dr$winner_target_retention_min &&
    top5_retention >= dr$target_top5_retention_min && zero_ok && one_ok && null_lock_ok

  data.frame(
    winner_target_retention = winner_retention,
    target_top5_retention = top5_retention,
    max_null_significant_target_lock_rate = max(nulls$significant_target_lock_rate),
    max_zero_signal_sig_correct_rate = max(zero$sig_correct_prob),
    min_one_x_sig_correct_rate = min(one$sig_correct_prob),
    classification = if (pass) "PASS" else "FAIL",
    stringsAsFactors = FALSE)
}

main <- function(argv = commandArgs(TRUE)) {
  root <- audit_repo_root()
  protocol <- load_protocol(root)
  if (!identical(protocol$status, "FROZEN_BEFORE_FINAL_RUN"))
    stop("protocol is not marked frozen")
  parallel <- any(argv == "--parallel")
  if (parallel) configure_parallel(TRUE, "auto")

  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              protocol$species, protocol$ion, protocol$source)$lines
  bins <- as.integer(protocol$bin_grid)

  observed_list <- lapply(bins, function(b) one_resolution_observed(lines, b, protocol))
  observed <- do.call(rbind, lapply(observed_list, `[[`, "row"))
  peaks <- do.call(rbind, lapply(seq_along(observed_list), function(i) {
    p <- observed_list[[i]]$peaks
    p$bins <- observed_list[[i]]$row$bins[1]
    p$sigma_bins <- observed_list[[i]]$row$sigma_bins[1]
    p$in_target_region <- abs(p$k - protocol$target_k) / protocol$target_k <=
      protocol$target_relative_tolerance
    p[, c("bins", "sigma_bins", "peak_rank", "k", "deltaD", "in_target_region")]
  }))

  null_list <- lapply(observed_list, null_one_resolution,
                      protocol = protocol, parallel = parallel)
  nulls <- do.call(rbind, lapply(null_list, `[[`, "summary"))
  injections <- do.call(rbind, lapply(seq_along(observed_list), function(i)
    inject_one_resolution(observed_list[[i]], null_list[[i]], protocol, parallel)))

  verdict <- classify_protocol(observed, nulls, injections, protocol)

  td <- file.path(root, "tables_r/statistical_audit")
  write_table(observed, file.path(td, "resolution_invariant_v1_observed.csv"))
  write_table(peaks, file.path(td, "resolution_invariant_v1_top_peaks.csv"))
  write_table(nulls, file.path(td, "resolution_invariant_v1_null_summary.csv"))
  write_table(injections, file.path(td, "resolution_invariant_v1_injection_summary.csv"))
  write_table(verdict, file.path(td, "resolution_invariant_v1_verdict.csv"))

  machine <- list(
    schema_version = "resolution_invariant_protocol_v1_result",
    protocol = "config/resolution_invariant_protocol_v1.json",
    git_commit = tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"),
                                  stdout = TRUE)[1], error = function(e) NA_character_),
    target_k = protocol$target_k,
    bin_grid = bins,
    historical_fixed_sigma_reference_retention =
      protocol$reporting$original_fixed_sigma_reference_retention,
    verdict = as.list(verdict[1, , drop = FALSE]))
  jsonlite::write_json(machine,
                       file.path(root, "resolution_invariant_protocol_v1_result.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = NA)

  cat("\nRESOLUTION-INVARIANT PROTOCOL V1\n")
  print(verdict, row.names = FALSE)
  invisible(list(observed = observed, peaks = peaks, nulls = nulls,
                 injections = injections, verdict = verdict))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_resolution_invariant_protocol_v1\\.R$", .invoked_file)) main()
