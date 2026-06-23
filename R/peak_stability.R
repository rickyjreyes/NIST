#!/usr/bin/env Rscript
# peak_stability.R
# ---------------------------------------------------------------------------
# PREDEFINED peak-region stability from the data-bootstrap draws.
#
# Peak-region tolerances are defined BEFORE reading bootstrap results. Two
# definitions are supported:
#   1. absolute : |k_boot - k_ref| <= tol_absolute
#   2. relative : |k_boot - k_ref| / k_ref <= tol_relative   (PRIMARY)
# The primary confirmatory definition is relative tolerance = 2%. Sensitivity
# is also reported at 1%, 2%, 5%.
#
# Stability classes are DESCRIPTIVE, not universal statistical laws:
#   high >= 80% ; moderate 50-80% ; low < 50%.
#
# Writes:
#   tables_r/statistical_audit/peak_stability.csv
#   figures_r/statistical_audit/peak_selection_stability.png (+ fig04)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

# identify_competing_regions(): cluster bootstrap k that fall OUTSIDE the
# reference region into descriptive "competing regions" by rounding to a grid,
# and report the largest competitor's share.
competing_share <- function(k, k_ref, tol_rel) {
  outside <- k[is.finite(k) & abs(k - k_ref) / k_ref > tol_rel]
  if (length(outside) == 0L) return(list(share = 0, centre = NA_real_))
  # bin outside draws to ~0.5 rad clusters
  cl <- round(outside / 0.5) * 0.5
  tab <- sort(table(cl), decreasing = TRUE)
  list(share = as.numeric(tab[1]) / length(k), centre = as.numeric(names(tab)[1]))
}

peak_stability_table <- function(draws, k_ref, tol_grid = c(0.01, 0.02, 0.05),
                                 tol_absolute = NULL, primary = 0.02) {
  k <- draws$k_best[draws$converged & is.finite(draws$k_best)]
  n <- length(k)
  rows <- lapply(tol_grid, function(tol) {
    pct <- mean(abs(k - k_ref) / k_ref <= tol)
    comp <- competing_share(k, k_ref, tol)
    data.frame(
      tolerance_type = "relative", tolerance = tol,
      is_primary = isTRUE(abs(tol - primary) < 1e-12),
      k_reference = k_ref, n_valid = n,
      pct_in_reference = pct,
      pct_largest_competitor = comp$share,
      competitor_centre = comp$centre,
      median_k = stats::median(k), iqr_k = stats::IQR(k),
      ci_lo = unname(stats::quantile(k, 0.025)),
      ci_hi = unname(stats::quantile(k, 0.975)),
      stability_class = stability_class(pct),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (!is.null(tol_absolute)) {
    pct <- mean(abs(k - k_ref) <= tol_absolute)
    out <- rbind(out, data.frame(
      tolerance_type = "absolute", tolerance = tol_absolute, is_primary = FALSE,
      k_reference = k_ref, n_valid = n, pct_in_reference = pct,
      pct_largest_competitor = NA_real_, competitor_centre = NA_real_,
      median_k = stats::median(k), iqr_k = stats::IQR(k),
      ci_lo = unname(stats::quantile(k, 0.025)),
      ci_hi = unname(stats::quantile(k, 0.975)),
      stability_class = stability_class(pct), stringsAsFactors = FALSE))
  }
  out
}

plot_peak_selection <- function(draws, k_ref, stab, primary = 0.02) {
  k <- draws$k_best[draws$converged & is.finite(draws$k_best)]
  df <- data.frame(k = k)
  band_lo <- k_ref * (1 - primary); band_hi <- k_ref * (1 + primary)
  pct <- stab$pct_in_reference[stab$is_primary][1]
  ggplot2::ggplot(df, ggplot2::aes(k)) +
    ggplot2::annotate("rect", xmin = band_lo, xmax = band_hi, ymin = -Inf, ymax = Inf,
                      fill = species_colour("Fe"), alpha = 0.12) +
    ggplot2::geom_histogram(bins = 60, fill = species_colour("Fe"), colour = "white",
                            linewidth = 0.1) +
    ggplot2::geom_vline(xintercept = k_ref, linetype = 2, colour = "grey15") +
    ggplot2::labs(
      title = "Fe II bootstrap peak-selection stability",
      subtitle = sprintf("Shaded = +/-%.0f%% reference region; %.1f%% of resamples select it.",
                         100 * primary, 100 * pct),
      x = "bootstrap k_best (rad / ln cm^-1)", y = "count") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  cfg <- default_audit_config()
  root <- audit_repo_root()
  draws_path <- file.path(root, "tables_r/statistical_audit/bootstrap_peak_draws.csv")
  if (!file.exists(draws_path)) stop("Run bootstrap_peak_uncertainty.R first (no draws found).")
  draws <- utils::read.csv(draws_path)
  pe_path <- file.path(root, "tables_r/statistical_audit/peak_estimates.csv")
  k_ref <- if (file.exists(pe_path)) utils::read.csv(pe_path)$k_best[1] else stats::median(draws$k_best, na.rm = TRUE)

  # absolute tolerance set to one grid step for reference
  k_grid <- audit_k_grid(cfg)
  step <- diff(k_grid)[1]
  stab <- peak_stability_table(draws, k_ref, cfg$tol_grid, tol_absolute = 5 * step,
                               primary = cfg$tol_primary)
  write_table(stab, file.path(root, "tables_r/statistical_audit/peak_stability.csv"))

  p <- plot_peak_selection(draws, k_ref, stab, cfg$tol_primary)
  save_fig(p, file.path(root, "figures_r/statistical_audit/peak_selection_stability.png"))
  save_fig(p, file.path(root, "figures_r/statistical_audit/fig04_peak_stability.png"))

  prim <- stab[stab$is_primary, ][1, ]
  cat(sprintf("[peak_stability] primary(2%%): %.1f%% in-region, class=%s, median_k=%.4f\n",
              100 * prim$pct_in_reference, prim$stability_class, prim$median_k))
  invisible(stab)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("peak_stability\\.R$", .invoked_file)) {
  main()
}
