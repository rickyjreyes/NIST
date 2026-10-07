#!/usr/bin/env Rscript
# run_final_contract_checks.R
# ---------------------------------------------------------------------------
# Contract-specific final checks required by the final NIST adversarial audit.
#
# This module does not search for a better statistic or retune the detector.
# It reuses the frozen scanner/configuration and adds only the explicitly
# requested final checks that were not represented as first-class artifacts:
#   * exact Fe II / Co II 120/160/200 reproduction table
#   * data-bootstrap stability for BOTH Fe II and Co II
#   * canonical-null candidate-location distributions for BOTH species
#   * null p-value calibration at 10%, 5%, 1%, and 0.1%
#   * injection/recovery at 0, 0.25x, 0.5x, 1.0x, and 1.5x observed amplitude
#   * deterministic sorting and leave-decile-out influence diagnostics
#
# Co II is never pooled with Fe II. Its target is its own frozen primary scan
# result and its fixed-frequency p-value is explicitly labelled descriptive
# because the target was scan-selected from the same species data.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
if (!exists("emp_p", mode = "function")) source(file.path(R_DIR, "audit_utils.R"))

parse_contract_args <- function(argv) {
  d <- list(`bootstrap-n` = 5000L, `null-n` = 5000L,
            `calibration-n` = 10000L, `injection-n` = 2000L,
            seed = 20260517L, parallel = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!key %in% names(d)) { i <- i + 1L; next }
    if (key == "parallel") {
      nxt <- if (i + 1L <= length(argv)) tolower(argv[i + 1L]) else NA_character_
      if (!is.na(nxt) && nxt %in% c("true", "false")) {
        d[[key]] <- nxt == "true"; i <- i + 2L
      } else {
        d[[key]] <- TRUE; i <- i + 1L
      }
    } else {
      if (i + 1L > length(argv)) stop("missing value for --", key)
      d[[key]] <- as.integer(argv[i + 1L]); i <- i + 2L
    }
  }
  for (k in c("bootstrap-n", "null-n", "calibration-n", "injection-n")) {
    if (!is.finite(d[[k]]) || d[[k]] <= 0L) stop("invalid --", k)
  }
  d
}

species_data <- function(root, species) {
  p <- file.path(root, sprintf("data/%s_lines.csv", species))
  if (!file.exists(p)) stop("missing frozen input: ", p)
  raw <- read_nist_csv(p)
  cl <- clean_lines_source(raw, species, 2L, "wavenumber")
  if (!isTRUE(cl$flow$reconciles)) stop(species, " cleaning flow does not reconcile")
  list(path = p, raw = raw, lines = cl$lines, flow = cl$flow)
}

observed_species <- function(lines, species, cfg) {
  kg <- audit_k_grid(cfg)
  rows <- lapply(c(120L, 160L, 200L), function(b) {
    r <- run_scan_analysis(lines, b, kg, cfg$degree, cfg$baseline_sigma)
    data.frame(
      species = species, ion = 2L, source = "wavenumber", bins = b,
      n_lines = r$n_lines, k_best = r$best$k_best,
      delta_log_x = delta_log_x(r$best$k_best),
      scale_ratio = scale_ratio(r$best$k_best),
      n_obs = r$n_obs, amplitude = r$best$amplitude, phase = r$best$phase,
      deltaD = r$best$deltaD, D_base = r$best$D_base,
      D_harmonic = r$best$D_harmonic,
      k_min = cfg$k_min, k_max = cfg$k_max,
      n_scanned_frequencies = cfg$n_k, baseline_sigma = cfg$baseline_sigma,
      degree = cfg$degree, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  kref <- tab$k_best[tab$bins == cfg$bins_primary][1]
  tab$primary_k <- kref
  tab$relative_k_error <- abs(tab$k_best - kref) / kref
  tab$compatible_2pct <- tab$relative_k_error <= cfg$tol_primary
  tab
}

bootstrap_species <- function(lines, species, cfg, B, seed, parallel = FALSE) {
  kg <- audit_k_grid(cfg)
  ref <- run_scan_analysis(lines, cfg$bins_primary, kg, cfg$degree, cfg$baseline_sigma)
  kref <- ref$best$k_best
  setup_rng(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(i) {
    set.seed(seeds[i])
    idx <- sample.int(nrow(lines), nrow(lines), replace = TRUE)
    rr <- tryCatch(run_scan_analysis(lines[idx, , drop = FALSE],
                                    cfg$bins_primary, kg, cfg$degree,
                                    cfg$baseline_sigma),
                   error = function(e) NULL)
    if (is.null(rr)) {
      return(data.frame(species = species, draw = i, converged = FALSE,
                        k_best = NA_real_, deltaD = NA_real_, amplitude = NA_real_,
                        phase = NA_real_, in_reference_region = NA,
                        exceeds_observed_stat = NA, stringsAsFactors = FALSE))
    }
    data.frame(
      species = species, draw = i, converged = TRUE,
      k_best = rr$best$k_best, deltaD = rr$best$deltaD,
      amplitude = rr$best$amplitude, phase = rr$best$phase,
      in_reference_region = abs(rr$best$k_best - kref) / kref <= cfg$tol_primary,
      exceeds_observed_stat = rr$best$deltaD >= ref$best$deltaD,
      stringsAsFactors = FALSE)
  }
  draws <- do.call(rbind, audit_lapply(seq_len(B), one, parallel = parallel))
  ok <- draws$converged %in% TRUE
  kval <- draws$k_best[ok & is.finite(draws$k_best)]
  dval <- draws$deltaD[ok & is.finite(draws$deltaD)]
  aval <- draws$amplitude[ok & is.finite(draws$amplitude)]
  summary <- data.frame(
    species = species, B = B, observed_k = kref,
    observed_deltaD = ref$best$deltaD, observed_amplitude = ref$best$amplitude,
    n_converged = sum(ok), failure_rate = mean(!ok),
    candidate_retention_2pct = mean(draws$in_reference_region, na.rm = TRUE),
    statistic_exceedance_count = sum(draws$exceeds_observed_stat %in% TRUE, na.rm = TRUE),
    statistic_exceedance_fraction = mean(draws$exceeds_observed_stat, na.rm = TRUE),
    k_median = stats::median(kval), k_ci_lo = unname(stats::quantile(kval, 0.025)),
    k_ci_hi = unname(stats::quantile(kval, 0.975)),
    deltaD_median = stats::median(dval),
    deltaD_ci_lo = unname(stats::quantile(dval, 0.025)),
    deltaD_ci_hi = unname(stats::quantile(dval, 0.975)),
    amplitude_median = stats::median(aval),
    amplitude_ci_lo = unname(stats::quantile(aval, 0.025)),
    amplitude_ci_hi = unname(stats::quantile(aval, 0.975)),
    note = "bootstrap exceedance fraction is a stability diagnostic, not a null p-value",
    stringsAsFactors = FALSE)
  list(draws = draws, summary = summary, reference = ref)
}

null_species <- function(lines, species, cfg, B, seed, parallel = FALSE) {
  kg <- audit_k_grid(cfg)
  primary <- run_scan_analysis(lines, cfg$bins_primary, kg, cfg$degree, cfg$baseline_sigma)
  kref <- primary$best$k_best
  summaries <- list(); draws_all <- list(); ref_by_bin <- list()
  bins_vec <- c(120L, 160L, 200L)
  for (j in seq_along(bins_vec)) {
    b <- bins_vec[j]
    ref <- run_scan_analysis(lines, b, kg, cfg$degree, cfg$baseline_sigma)
    ref_by_bin[[as.character(b)]] <- ref
    kidx <- nearest_k_index(kg, kref)
    obs_fixed <- ref$scan$deltaD[kidx]
    setup_rng(seed + 1000L * j)
    seeds <- sample.int(.Machine$integer.max, B)
    one <- function(i) {
      set.seed(seeds[i])
      y0 <- stats::rpois(length(ref$mu0), ref$mu0)
      sk <- scan_k(ref$ell, y0, ref$baseline, kg, cfg$degree)
      data.frame(
        species = species, bins = b, draw = i,
        max_deltaD = sk$best$deltaD, k_at_max = sk$best$k_best,
        fixed_k_deltaD = sk$scan$deltaD[kidx],
        stringsAsFactors = FALSE)
    }
    dr <- do.call(rbind, audit_lapply(seq_len(B), one, parallel = parallel))
    gt <- sum(dr$max_deltaD >= ref$best$deltaD)
    ft <- sum(dr$fixed_k_deltaD >= obs_fixed)
    summaries[[as.character(b)]] <- data.frame(
      species = species, bins = b, B = B,
      observed_global_stat = ref$best$deltaD,
      observed_k_best = ref$best$k_best,
      fixed_k = kref, observed_fixed_k_stat = obs_fixed,
      global_exceedances = gt, global_p = emp_p(gt, B),
      fixed_exceedances = ft, fixed_frequency_p = emp_p(ft, B),
      mc_resolution = resolution_floor(B),
      null_max = max(dr$max_deltaD),
      null_q99 = unname(stats::quantile(dr$max_deltaD, 0.99)),
      null_candidate_k_median = stats::median(dr$k_at_max),
      candidate_origin = paste0(
        "species primary 160-bin scan-selected k; fixed-frequency p is descriptive ",
        "and is not a preregistered independent-frequency significance"),
      stringsAsFactors = FALSE)
    draws_all[[as.character(b)]] <- dr
  }
  list(summary = do.call(rbind, summaries),
       draws = do.call(rbind, draws_all),
       primary = primary, refs = ref_by_bin)
}

calibrate_species <- function(lines, species, cfg, calibration_n, null_n,
                              seed, parallel = FALSE) {
  kg <- audit_k_grid(cfg)
  ref <- run_scan_analysis(lines, cfg$bins_primary, kg, cfg$degree, cfg$baseline_sigma)

  setup_rng(seed + 20001L)
  inner_seeds <- sample.int(.Machine$integer.max, null_n)
  inner_one <- function(i) {
    set.seed(inner_seeds[i])
    y0 <- stats::rpois(length(ref$mu0), ref$mu0)
    scan_k(ref$ell, y0, ref$baseline, kg, cfg$degree)$best$deltaD
  }
  ref_null <- unlist(audit_lapply(seq_len(null_n), inner_one, parallel = parallel),
                     use.names = FALSE)

  setup_rng(seed + 30001L)
  outer_seeds <- sample.int(.Machine$integer.max, calibration_n)
  outer_one <- function(i) {
    set.seed(outer_seeds[i])
    y0 <- stats::rpois(length(ref$mu0), ref$mu0)
    stat <- scan_k(ref$ell, y0, ref$baseline, kg, cfg$degree)$best$deltaD
    data.frame(species = species, draw = i, statistic = stat,
               empirical_p = emp_p(sum(ref_null >= stat), null_n),
               stringsAsFactors = FALSE)
  }
  pdraws <- do.call(rbind, audit_lapply(seq_len(calibration_n), outer_one,
                                        parallel = parallel))
  alphas <- c(0.10, 0.05, 0.01, 0.001)
  rows <- lapply(alphas, function(a) {
    n <- nrow(pdraws); obs <- sum(pdraws$empirical_p <= a)
    ci <- binom_ci(obs, n)
    data.frame(
      species = species, nominal_alpha = a, n_outer_sims = n,
      inner_null_n = null_n, observed_count = obs,
      expected_count = a * n, observed_fpr = obs / n,
      ci_lo = unname(ci["lower"]), ci_hi = unname(ci["upper"]),
      compatible = a >= ci["lower"] && a <= ci["upper"],
      stringsAsFactors = FALSE)
  })
  list(summary = do.call(rbind, rows), pvalues = pdraws, reference_null = ref_null)
}

inject_species <- function(lines, species, cfg, inject_n, null_n, seed,
                           parallel = FALSE) {
  kg <- audit_k_grid(cfg)
  ref <- run_scan_analysis(lines, cfg$bins_primary, kg, cfg$degree, cfg$baseline_sigma)
  kref <- ref$best$k_best
  observed_amp <- abs(ref$best$amplitude)
  ratios <- c(0, 0.25, 0.50, 1.00, 1.50)
  phi0 <- 0.6

  setup_rng(seed + 40001L)
  inner_seeds <- sample.int(.Machine$integer.max, null_n)
  inner_one <- function(i) {
    set.seed(inner_seeds[i])
    y0 <- stats::rpois(length(ref$mu0), ref$mu0)
    scan_k(ref$ell, y0, ref$baseline, kg, cfg$degree)$best$deltaD
  }
  ref_null <- unlist(audit_lapply(seq_len(null_n), inner_one, parallel = parallel),
                     use.names = FALSE)
  threshold <- unname(stats::quantile(ref_null, 0.95))

  rows <- list()
  for (j in seq_along(ratios)) {
    ratio <- ratios[j]
    A <- ratio * observed_amp
    setup_rng(seed + 50000L + 1000L * j)
    seeds <- sample.int(.Machine$integer.max, inject_n)
    one <- function(i) {
      set.seed(seeds[i])
      mu <- pmax(ref$mu0 * exp(A * cos(kref * ref$ell - phi0)), EPS)
      y <- stats::rpois(length(mu), mu)
      sk <- scan_k(ref$ell, y, ref$baseline, kg, cfg$degree)
      kb <- sk$best$k_best
      c(stat = sk$best$deltaD, k = kb, amp = sk$best$amplitude,
        detect = as.numeric(sk$best$deltaD >= threshold),
        correct = as.numeric(abs(kb - kref) / kref <= cfg$tol_primary))
    }
    M <- do.call(rbind, audit_lapply(seq_len(inject_n), one, parallel = parallel))
    det <- sum(M[, "detect"] == 1)
    cor <- sum(M[, "correct"] == 1)
    sigcor <- sum(M[, "detect"] == 1 & M[, "correct"] == 1)
    false_loc <- sum(M[, "detect"] == 1 & M[, "correct"] == 0)
    dci <- binom_ci(det, inject_n); sci <- binom_ci(sigcor, inject_n)
    rows[[j]] <- data.frame(
      species = species, amplitude_ratio = ratio, injected_amplitude = A,
      observed_amplitude = observed_amp, k_injected = kref,
      n_sims = inject_n, null_n = null_n, threshold_alpha = 0.05,
      detection_prob = det / inject_n,
      det_ci_lo = unname(dci["lower"]), det_ci_hi = unname(dci["upper"]),
      correct_region_prob = cor / inject_n,
      sig_correct_prob = sigcor / inject_n,
      sigcor_ci_lo = unname(sci["lower"]), sigcor_ci_hi = unname(sci["upper"]),
      false_location_prob = false_loc / inject_n,
      median_recovered_k = stats::median(M[, "k"]),
      k_bias = stats::median(M[, "k"]) - kref,
      amplitude_bias = stats::median(M[, "amp"]) - A,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

artifact_species <- function(lines, flow, species, cfg) {
  kg <- audit_k_grid(cfg)
  ref <- run_scan_analysis(lines, cfg$bins_primary, kg, cfg$degree, cfg$baseline_sigma)
  rev_ref <- run_scan_analysis(lines[nrow(lines):1L, , drop = FALSE],
                               cfg$bins_primary, kg, cfg$degree, cfg$baseline_sigma)
  sort_k_diff <- abs(ref$best$k_best - rev_ref$best$k_best)
  sort_stat_diff <- abs(ref$best$deltaD - rev_ref$best$deltaD)

  influence <- lapply(0:9, function(fold) {
    keep <- (seq_len(nrow(lines)) %% 10L) != fold
    rr <- run_scan_analysis(lines[keep, , drop = FALSE], cfg$bins_primary, kg,
                            cfg$degree, cfg$baseline_sigma)
    data.frame(
      species = species, left_out_decile = fold,
      retained_lines = sum(keep), k_best = rr$best$k_best,
      deltaD = rr$best$deltaD, amplitude = rr$best$amplitude,
      in_reference_region = abs(rr$best$k_best - ref$best$k_best) /
        ref$best$k_best <= cfg$tol_primary,
      stringsAsFactors = FALSE)
  })
  influence <- do.call(rbind, influence)
  checks <- data.frame(
    species = species,
    raw_rows = flow$n_raw,
    retained_lines = flow$retained,
    duplicates_removed = flow$n_duplicates,
    missing_wavenumber = flow$missing_wavenumber,
    nonpositive_removed = flow$excl_nonpositive,
    cleaning_reconciles = isTRUE(flow$reconciles),
    sorting_k_abs_diff = sort_k_diff,
    sorting_stat_abs_diff = sort_stat_diff,
    sorting_invariant = sort_k_diff <= 1e-12 && sort_stat_diff <= 1e-8,
    influence_reference_fraction = mean(influence$in_reference_region),
    fft_applicable = FALSE,
    fft_note = "canonical scanner is direct harmonic regression over an explicit k grid; no FFT normalization is used",
    stringsAsFactors = FALSE)
  list(checks = checks, influence = influence)
}

plot_contract <- function(root, observed, bootstrap_draws, null_draws, calibration, injection) {
  fd <- file.path(root, "figures_r/statistical_audit")
  dir.create(fd, showWarnings = FALSE, recursive = TRUE)

  prim <- observed[observed$bins == 160, c("species", "deltaD")]
  p1 <- ggplot2::ggplot(null_draws[null_draws$bins == 160, ],
                        ggplot2::aes(max_deltaD, fill = species)) +
    ggplot2::geom_histogram(bins = 60, alpha = 0.55, position = "identity") +
    ggplot2::geom_vline(data = prim, ggplot2::aes(xintercept = deltaD, colour = species),
                        linewidth = 0.8, linetype = 2) +
    ggplot2::facet_wrap(~species, scales = "free_y") +
    ggplot2::labs(title = "Observed statistic versus canonical null",
                  x = "scan-global max deltaD", y = "count") + theme_audit()
  save_fig(p1, file.path(fd, "final_observed_vs_null.png"))

  p2 <- ggplot2::ggplot(bootstrap_draws[bootstrap_draws$converged %in% TRUE, ],
                        ggplot2::aes(k_best, fill = species)) +
    ggplot2::geom_histogram(bins = 60, alpha = 0.55, position = "identity") +
    ggplot2::facet_wrap(~species, scales = "free_y") +
    ggplot2::labs(title = "Bootstrap candidate-location stability",
                  x = "selected k", y = "count") + theme_audit()
  save_fig(p2, file.path(fd, "final_bootstrap_candidate_stability.png"))

  p3 <- ggplot2::ggplot(calibration,
                        ggplot2::aes(nominal_alpha, observed_fpr, colour = species)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), width = 0.002) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Null p-value calibration", x = "nominal alpha",
                  y = "empirical false-positive rate") + theme_audit()
  save_fig(p3, file.path(fd, "final_calibration_curve.png"))

  p4 <- ggplot2::ggplot(injection,
                        ggplot2::aes(amplitude_ratio, sig_correct_prob, colour = species)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = 3) +
    ggplot2::labs(title = "Relative-amplitude injection/recovery",
                  x = "injected amplitude / observed amplitude",
                  y = "globally significant correct-location recovery") + theme_audit()
  save_fig(p4, file.path(fd, "final_injection_recovery_power.png"))

  p5 <- ggplot2::ggplot(observed, ggplot2::aes(bins, k_best, colour = species)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
    ggplot2::labs(title = "Predeclared 120/160/200-bin robustness",
                  x = "bins", y = "selected k") + theme_audit()
  save_fig(p5, file.path(fd, "final_binning_robustness.png"))

  p6 <- ggplot2::ggplot(observed[observed$bins == 160, ],
                        ggplot2::aes(species, k_best, colour = species)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::labs(title = "Fe II and Co II are separate species-level results",
                  subtitle = "No pooling or post-hoc significance combination",
                  x = NULL, y = "primary selected k") + theme_audit()
  save_fig(p6, file.path(fd, "final_fe_co_comparison.png"))
}

main <- function(argv = commandArgs(TRUE)) {
  args <- parse_contract_args(argv)
  root <- audit_repo_root()
  td <- file.path(root, "tables_r/statistical_audit")
  dir.create(td, showWarnings = FALSE, recursive = TRUE)

  if (isTRUE(args$parallel)) configure_parallel(TRUE, "auto")
  cfg <- default_audit_config(seed = args$seed,
                              bootstrap_n = args[["bootstrap-n"]],
                              null_n = args[["null-n"]],
                              calibration_n = args[["calibration-n"]],
                              injection_n = args[["injection-n"]])

  sp <- c("Fe", "Co")
  dat <- setNames(lapply(sp, function(s) species_data(root, s)), sp)

  obs <- do.call(rbind, lapply(sp, function(s) observed_species(dat[[s]]$lines, s, cfg)))
  rownames(obs) <- NULL
  write_table(obs, file.path(td, "final_species_observed_results.csv"))

  bs <- lapply(seq_along(sp), function(i)
    bootstrap_species(dat[[sp[i]]]$lines, sp[i], cfg, args[["bootstrap-n"]],
                      args$seed + 100000L * i, args$parallel))
  bdraw <- do.call(rbind, lapply(bs, `[[`, "draws"))
  bsum <- do.call(rbind, lapply(bs, `[[`, "summary"))
  write_table(bdraw, file.path(td, "final_species_bootstrap_draws.csv"))
  write_table(bsum, file.path(td, "final_species_bootstrap_summary.csv"))

  ns <- lapply(seq_along(sp), function(i)
    null_species(dat[[sp[i]]]$lines, sp[i], cfg, args[["null-n"]],
                 args$seed + 200000L * i, args$parallel))
  ndraw <- do.call(rbind, lapply(ns, `[[`, "draws"))
  nsum <- do.call(rbind, lapply(ns, `[[`, "summary"))
  write_table(ndraw, file.path(td, "final_species_null_draws.csv"))
  write_table(nsum, file.path(td, "final_species_null_summary.csv"))
  null_design <- data.frame(
    null_model = "canonical_fitted_smooth_poisson_fixed_baseline",
    preserves = paste("species-specific line count expectation after canonical binning;",
                      "observed wavelength/wavenumber coverage; bin grid; smooth density baseline;",
                      "scan range and frequency grid; normalization and selection rule"),
    destroys = paste("coherent log-periodic phase/structure beyond the fitted smooth baseline;",
                     "individual catalog line positions are not retained in the synthetic count draw"),
    caveat = paste("Fe alternative-null sensitivity remains in alternative_null_results.csv;",
                   "this contract table does not replace that predeclared sensitivity analysis"),
    stringsAsFactors = FALSE)
  write_table(null_design, file.path(td, "final_null_design.csv"))

  cal <- lapply(seq_along(sp), function(i)
    calibrate_species(dat[[sp[i]]]$lines, sp[i], cfg,
                      args[["calibration-n"]], args[["null-n"]],
                      args$seed + 300000L * i, args$parallel))
  csum <- do.call(rbind, lapply(cal, `[[`, "summary"))
  cp <- do.call(rbind, lapply(cal, `[[`, "pvalues"))
  write_table(csum, file.path(td, "final_species_calibration.csv"))
  write_table(cp, file.path(td, "final_species_calibration_pvalues.csv"))

  inj <- do.call(rbind, lapply(seq_along(sp), function(i)
    inject_species(dat[[sp[i]]]$lines, sp[i], cfg, args[["injection-n"]],
                   args[["null-n"]], args$seed + 400000L * i, args$parallel)))
  write_table(inj, file.path(td, "final_species_injection_recovery.csv"))

  art <- lapply(sp, function(s) artifact_species(dat[[s]]$lines, dat[[s]]$flow, s, cfg))
  checks <- do.call(rbind, lapply(art, `[[`, "checks"))
  influence <- do.call(rbind, lapply(art, `[[`, "influence"))
  write_table(checks, file.path(td, "final_species_artifact_checks.csv"))
  write_table(influence, file.path(td, "final_species_influence.csv"))

  prov <- do.call(rbind, lapply(sp, function(s) {
    f <- dat[[s]]$flow
    data.frame(
      species = s, input_path = sprintf("data/%s_lines.csv", s),
      source_provider = "NIST spectral-line dataset",
      source_retrieval_version = NA_character_,
      source_retrieval_date = NA_character_,
      source_retrieval_metadata_complete = FALSE,
      ion = 2L, source_field = f$source_field, transform = "ell = ln(wavenumber / cm^-1)",
      duplicate_policy = "sort ascending wavenumber; retain first exact duplicate wavenumber",
      missing_policy = "exclude missing/non-numeric selected source before nonpositive filter",
      retained_lines = f$retained, cleaning_reconciles = f$reconciles,
      note = paste("Repository does not currently contain authoritative NIST retrieval/version",
                   "metadata; do not invent it. Populate from a primary retrieval record before",
                   "claiming provenance-complete freeze."),
      stringsAsFactors = FALSE)
  }))
  write_table(prov, file.path(td, "final_source_provenance.csv"))

  plot_contract(root, obs, bdraw, ndraw, csum, inj)

  cat(sprintf("[final_contract] Fe/Co observed=%d rows; bootstrap=%d/species; null=%d/bin/species; calibration=%d/species; injection=%d/cell\n",
              nrow(obs), args[["bootstrap-n"]], args[["null-n"]],
              args[["calibration-n"]], args[["injection-n"]]))
  invisible(list(observed = obs, bootstrap = bsum, null = nsum,
                 calibration = csum, injection = inj, artifacts = checks,
                 provenance = prov))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_final_contract_checks\\.R$", .invoked_file)) main()
