#!/usr/bin/env Rscript
# run_bin_grid.R
# ---------------------------------------------------------------------------
# Bin-stability grid for the canonical Fe II scan over the predefined grid
#   60, 80, 100, 120, 140, 160, 180, 200, 220, 240.
# For every bin count: retained lines, k_best, log-period, scale ratio,
# amplitude, phase, deltaD, scan-global p, null count, peak-region identity.
# A run is NOT called stable merely because some significant peak appears -
# we measure whether the SAME reference region is selected.
#
# Writes:
#   tables_r/statistical_audit/bin_grid_results.csv
#   tables_r/statistical_audit/bin_stability_summary.csv
#   figures_r/statistical_audit/frequency_by_bin_count.png (+ fig03)
#   figures_r/statistical_audit/bin_stability_heatmap.png
#   figures_r/statistical_audit/bin_scan_curves.png
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

run_bin_grid <- function(lines, cfg, k_ref, null_n = cfg$null_n, bins_grid = cfg$bins_grid) {
  k_grid <- audit_k_grid(cfg)
  rows <- list(); curves <- list()
  for (b in bins_grid) {
    setup_rng(cfg$seed)
    res <- run_scan_analysis(lines, b, k_grid, cfg$degree, cfg$baseline_sigma)
    nd <- null_distribution(res$ell, res$y, res$baseline, res$mu0, k_grid,
                            cfg$degree, null_n)
    tail <- sum(nd$max_vals >= res$best$deltaD)
    kb <- res$best$k_best
    rows[[as.character(b)]] <- data.frame(
      bins = b, n_lines = res$n_lines, k_best = kb,
      delta_log_x = delta_log_x(kb), scale_ratio = scale_ratio(kb),
      amplitude = res$best$amplitude, phase = res$best$phase,
      deltaD = res$best$deltaD, n_obs = res$n_obs,
      global_p = emp_p(tail, null_n), global_tail = tail, null_n = null_n,
      in_reference_region = in_peak_region(kb, k_ref, "relative", cfg$tol_primary),
      stringsAsFactors = FALSE)
    curves[[as.character(b)]] <- data.frame(bins = b, k = res$scan$k,
                                            deltaD = res$scan$deltaD)
    cat(sprintf("[bin] bins=%d k_best=%.4f deltaD=%.2f global_p=%.4g region=%s\n",
                b, kb, res$best$deltaD, emp_p(tail, null_n),
                rows[[as.character(b)]]$in_reference_region))
  }
  list(grid = do.call(rbind, rows), curves = do.call(rbind, curves))
}

# count mode switches: number of times the selected region toggles in/out of
# the reference region across the ordered bin grid.
count_mode_switches <- function(in_region) sum(diff(as.integer(in_region)) != 0)

bin_stability_summary <- function(grid, cfg) {
  k <- grid$k_best
  data.frame(
    mean_k = mean(k), median_k = stats::median(k), sd_k = stats::sd(k),
    iqr_k = stats::IQR(k), cv_k = cv(k),
    max_relative_drift = (max(k) - min(k)) / stats::median(k),
    n_mode_switches = count_mode_switches(grid$in_reference_region),
    pct_in_reference = mean(grid$in_reference_region),
    stability_class = stability_class(mean(grid$in_reference_region)),
    n_bins_tested = nrow(grid),
    tol_relative = cfg$tol_primary,
    stringsAsFactors = FALSE)
}

plot_frequency_by_bin <- function(grid) {
  ggplot2::ggplot(grid, ggplot2::aes(bins, k_best)) +
    ggplot2::geom_line(colour = species_colour("Fe"), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(shape = in_reference_region),
                        colour = species_colour("Fe"), size = 3) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                                name = "in reference region") +
    ggplot2::labs(title = "Fe II selected frequency across the bin grid",
                  subtitle = "Connected-dot plot; crosses mark bin settings outside the reference peak region.",
                  x = "bin count", y = "k_best (rad / ln cm^-1)") +
    theme_audit()
}

plot_bin_heatmap <- function(curves, grid) {
  ggplot2::ggplot(curves, ggplot2::aes(k, factor(bins), fill = deltaD)) +
    ggplot2::geom_raster() +
    ggplot2::geom_point(data = grid, ggplot2::aes(k_best, factor(bins)),
                        inherit.aes = FALSE, colour = "white", shape = 3, size = 1.5) +
    viridis::scale_fill_viridis(name = "deltaD", option = "magma") +
    ggplot2::labs(title = "Fe II deviance-improvement scan by bin count",
                  subtitle = "Crosses mark the selected k_best per bin count.",
                  x = "k (rad / ln cm^-1)", y = "bin count") +
    theme_audit()
}

plot_bin_scan_curves <- function(curves) {
  ggplot2::ggplot(curves, ggplot2::aes(k, deltaD, colour = factor(bins))) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.85) +
    viridis::scale_colour_viridis(discrete = TRUE, name = "bins") +
    ggplot2::labs(title = "Fe II scan curves overlaid across bin counts",
                  x = "k (rad / ln cm^-1)", y = "deltaD") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  null_n <- NULL
  i <- which(argv == "--null-n"); if (length(i) == 1L) null_n <- as.integer(argv[i + 1L])
  cfg <- default_audit_config(null_n = if (is.null(null_n)) 2000L else null_n)
  root <- audit_repo_root()
  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              "Fe", 2L, "wavenumber")$lines
  k_grid <- audit_k_grid(cfg)
  k_ref <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree,
                             cfg$baseline_sigma)$best$k_best

  bg <- run_bin_grid(lines, cfg, k_ref)
  write_table(bg$grid, file.path(root, "tables_r/statistical_audit/bin_grid_results.csv"))
  summ <- bin_stability_summary(bg$grid, cfg)
  write_table(summ, file.path(root, "tables_r/statistical_audit/bin_stability_summary.csv"))

  save_fig(plot_frequency_by_bin(bg$grid),
           file.path(root, "figures_r/statistical_audit/frequency_by_bin_count.png"))
  save_fig(plot_frequency_by_bin(bg$grid),
           file.path(root, "figures_r/statistical_audit/fig03_frequency_by_bin_count.png"))
  save_fig(plot_bin_heatmap(bg$curves, bg$grid),
           file.path(root, "figures_r/statistical_audit/bin_stability_heatmap.png"))
  save_fig(plot_bin_scan_curves(bg$curves),
           file.path(root, "figures_r/statistical_audit/bin_scan_curves.png"))

  cat(sprintf("[bin_grid] pct_in_reference=%.1f%% class=%s mode_switches=%d cv_k=%.4f\n",
              100 * summ$pct_in_reference, summ$stability_class,
              summ$n_mode_switches, summ$cv_k))
  invisible(bg)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_bin_grid\\.R$", .invoked_file)) {
  main()
}
