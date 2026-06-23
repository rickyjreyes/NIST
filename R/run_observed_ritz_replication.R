#!/usr/bin/env Rscript
# run_observed_ritz_replication.R
# ---------------------------------------------------------------------------
# Source-field comparison (NOT fully independent replication): rerun the Fe II
# scan on three measurement representations kept strictly separate:
#   observed   : wavenumbers from the observed wavelength column
#   ritz       : wavenumbers from the Ritz wavelength column
#   wavenumber : the direct NIST wn(cm-1) column
#
# Observed and Ritz values often describe the SAME physical transitions, so this
# is a "measurement-representation comparison", not independent replication.
#
# Writes:
#   tables_r/statistical_audit/observed_ritz_replication.csv
#   figures_r/statistical_audit/observed_ritz_comparison.png (+ fig08)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

source_overlap <- function(lines_a, lines_b, tol = 1e-6) {
  a <- lines_a$wavenumber_cm; b <- lines_b$wavenumber_cm
  # overlap = number of a-values within tol (relative) of any b-value
  # use rounding to a fine grid for speed
  ra <- round(a, 4); rb <- round(b, 4)
  sum(ra %in% rb)
}

run_source_replication <- function(cfg, null_n = cfg$null_n, parallel = FALSE,
                                   bootstrap_n = 0L) {
  root <- audit_repo_root()
  raw <- read_nist_csv(file.path(root, "data/Fe_lines.csv"))
  k_grid <- audit_k_grid(cfg)
  sources <- c("wavenumber", "observed", "ritz")
  cleaned <- lapply(sources, function(s) clean_lines_source(raw, "Fe", 2L, s)$lines)
  names(cleaned) <- sources

  # reference k from the canonical wavenumber source
  k_ref <- run_scan_analysis(cleaned[["wavenumber"]], cfg$bins_primary, k_grid,
                             cfg$degree, cfg$baseline_sigma)$best$k_best

  rows <- lapply(sources, function(s) {
    lines <- cleaned[[s]]
    setup_rng(cfg$seed)
    res <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)
    nd <- null_distribution(res$ell, res$y, res$baseline, res$mu0, k_grid, cfg$degree,
                            null_n, parallel = parallel)
    tail <- sum(nd$max_vals >= res$best$deltaD)
    kb <- res$best$k_best
    ci_lo <- NA_real_; ci_hi <- NA_real_
    if (bootstrap_n > 0L) {
      bs <- sample.int(.Machine$integer.max, bootstrap_n)
      kbs <- vapply(seq_len(bootstrap_n), function(b) {
        set.seed(bs[b]); idx <- sample.int(nrow(lines), nrow(lines), replace = TRUE)
        bl <- lines[idx, , drop = FALSE]; bl <- bl[order(bl$ell), , drop = FALSE]
        r <- tryCatch(run_scan_analysis(bl, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)$best$k_best,
                      error = function(e) NA_real_); r
      }, numeric(1))
      ci_lo <- stats::quantile(kbs, 0.025, na.rm = TRUE)
      ci_hi <- stats::quantile(kbs, 0.975, na.rm = TRUE)
    }
    data.frame(
      source = s, n_lines = nrow(lines),
      ell_min = min(res$ell), ell_max = max(res$ell),
      k_best = kb, delta_log_x = delta_log_x(kb), scale_ratio = scale_ratio(kb),
      amplitude = res$best$amplitude, deltaD = res$best$deltaD,
      global_p = emp_p(tail, null_n), null_n = null_n,
      ci_lo = ci_lo, ci_hi = ci_hi,
      same_reference_region = in_peak_region(kb, k_ref, "relative", cfg$tol_primary),
      abs_diff_from_wn = abs(kb - k_ref),
      rel_diff_from_wn = abs(kb - k_ref) / k_ref,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$overlap_with_wavenumber <- vapply(sources, function(s)
    source_overlap(cleaned[[s]], cleaned[["wavenumber"]]), integer(1))
  out
}

plot_observed_ritz <- function(rep_tab) {
  ggplot2::ggplot(rep_tab, ggplot2::aes(k_best, source, colour = source)) +
    { if (all(is.finite(rep_tab$ci_lo)))
        ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                               orientation = "y", width = 0.15) }+
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_colour_manual(values = c(wavenumber = species_colour("Fe"),
                                            observed = "#0072B2", ritz = "#009E73"),
                                 guide = "none") +
    ggplot2::labs(
      title = "Fe II peak location across measurement representations",
      subtitle = "Observed / Ritz / direct wavenumber. Not fully independent (often the same transitions).",
      x = "k_best (rad / ln cm^-1)", y = NULL) +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  null_n <- 1000L; bootstrap_n <- 0L
  i <- which(argv == "--null-n"); if (length(i) == 1L) null_n <- as.integer(argv[i + 1L])
  j <- which(argv == "--bootstrap-n"); if (length(j) == 1L) bootstrap_n <- as.integer(argv[j + 1L])
  parallel <- any(argv == "--parallel")
  cfg <- default_audit_config()
  if (parallel) configure_parallel(TRUE, "auto")
  root <- audit_repo_root()
  rep_tab <- run_source_replication(cfg, null_n, parallel, bootstrap_n)
  write_table(rep_tab, file.path(root, "tables_r/statistical_audit/observed_ritz_replication.csv"))
  save_fig(plot_observed_ritz(rep_tab),
           file.path(root, "figures_r/statistical_audit/observed_ritz_comparison.png"))
  save_fig(plot_observed_ritz(rep_tab),
           file.path(root, "figures_r/statistical_audit/fig08_observed_ritz_comparison.png"))
  cat("[observed_ritz] k_best by source:\n")
  print(rep_tab[, c("source", "n_lines", "k_best", "global_p", "same_reference_region")])
  invisible(rep_tab)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_observed_ritz_replication\\.R$", .invoked_file)) {
  main()
}
