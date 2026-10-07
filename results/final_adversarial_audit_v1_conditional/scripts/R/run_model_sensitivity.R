#!/usr/bin/env Rscript
# run_model_sensitivity.R
# ---------------------------------------------------------------------------
# Specification (multiverse) analysis over a DECLARED grid of analysis choices:
#   baseline sigma in {4,5,6,7,8}
#   polynomial degree in {0,1,2}
#   bins in {120,160,200}
#   source field in {wavenumber, observed, ritz}
# Each specification records k_best, deltaD, amplitude, global p (modest null),
# peak-region assignment and convergence. The figure shows whether the claimed
# peak exists across a coherent region of specifications or only in isolation.
#
# Writes:
#   tables_r/statistical_audit/specification_results.csv
#   figures_r/statistical_audit/specification_curve.png   (+ fig10)
#   figures_r/statistical_audit/specification_heatmap.png
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

build_spec_grid <- function(cfg) {
  sigma_grid <- if (isTRUE(cfg$fast)) c(5, 6, 7) else cfg$sigma_grid
  degree_grid <- if (isTRUE(cfg$fast)) c(1L, 2L) else cfg$degree_grid
  base <- expand.grid(
    sigma = sigma_grid, degree = degree_grid,
    bins = c(120L, 160L, 200L), source = "wavenumber",
    stringsAsFactors = FALSE)
  src <- expand.grid(
    sigma = cfg$baseline_sigma, degree = cfg$degree,
    bins = cfg$bins_primary, source = c("observed", "ritz"),
    stringsAsFactors = FALSE)
  unique(rbind(base, src))
}

run_specifications <- function(cfg, k_ref, spec_null_n = 200L, parallel = FALSE) {
  root <- audit_repo_root()
  raw <- read_nist_csv(file.path(root, "data/Fe_lines.csv"))
  k_grid <- audit_k_grid(cfg)
  grid <- build_spec_grid(cfg)
  # pre-clean per source (cached)
  src_lines <- lapply(unique(grid$source), function(s)
    clean_lines_source(raw, "Fe", 2L, s)$lines)
  names(src_lines) <- unique(grid$source)

  one <- function(i) {
    g <- grid[i, ]
    lines <- src_lines[[g$source]]
    setup_rng(cfg$seed)
    res <- tryCatch(run_scan_analysis(lines, g$bins, k_grid, g$degree, g$sigma),
                    error = function(e) NULL)
    if (is.null(res)) {
      return(data.frame(g, k_best = NA_real_, deltaD = NA_real_, amplitude = NA_real_,
                        global_p = NA_real_, in_reference_region = NA, converged = FALSE,
                        n_lines = nrow(lines), stringsAsFactors = FALSE))
    }
    nd <- null_distribution(res$ell, res$y, res$baseline, res$mu0, k_grid,
                            g$degree, spec_null_n)
    tail <- sum(nd$max_vals >= res$best$deltaD)
    kb <- res$best$k_best
    data.frame(g, k_best = kb, deltaD = res$best$deltaD, amplitude = res$best$amplitude,
               global_p = emp_p(tail, spec_null_n),
               in_reference_region = in_peak_region(kb, k_ref, "relative", cfg$tol_primary),
               converged = TRUE, n_lines = nrow(lines), stringsAsFactors = FALSE)
  }
  out <- audit_lapply(seq_len(nrow(grid)), one, parallel = parallel)
  res <- do.call(rbind, out)
  res$spec_null_n <- spec_null_n
  res$spec_id <- sprintf("sig%g_deg%d_bin%d_%s", res$sigma, res$degree, res$bins, res$source)
  res[order(res$in_reference_region, res$deltaD), ]
}

plot_specification_curve <- function(spec) {
  s <- spec[is.finite(spec$k_best), ]
  s <- s[order(s$k_best), ]
  s$rank <- seq_len(nrow(s))
  ggplot2::ggplot(s, ggplot2::aes(rank, k_best, colour = in_reference_region)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = species_colour("Fe"), `FALSE` = "#999999"),
                                 name = "in reference region") +
    ggplot2::labs(
      title = "Specification curve: Fe II selected frequency across analysis choices",
      subtitle = "Each point is one (sigma, degree, bins, source) specification, ordered by k_best.",
      x = "specification (ordered by k_best)", y = "k_best (rad / ln cm^-1)") +
    theme_audit()
}

plot_specification_heatmap <- function(spec) {
  s <- spec[spec$source == "wavenumber" & is.finite(spec$deltaD), ]
  ggplot2::ggplot(s, ggplot2::aes(factor(sigma), factor(degree), fill = in_reference_region)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::facet_wrap(~bins, labeller = ggplot2::label_both) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = species_colour("Fe"), `FALSE` = "#DDDDDD"),
                               name = "reference region") +
    ggplot2::labs(title = "Specification region map (wavenumber source)",
                  subtitle = "Whether each (sigma, degree, bins) selects the reference peak region.",
                  x = "baseline sigma", y = "polynomial degree") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  spec_null_n <- 200L
  i <- which(argv == "--null-n"); if (length(i) == 1L) spec_null_n <- as.integer(argv[i + 1L])
  parallel <- any(argv == "--parallel")
  cfg <- default_audit_config(fast = is_fast(argv))
  if (parallel) configure_parallel(TRUE, "auto")
  root <- audit_repo_root()
  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              "Fe", 2L, "wavenumber")$lines
  k_grid <- audit_k_grid(cfg)
  k_ref <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)$best$k_best

  spec <- run_specifications(cfg, k_ref, spec_null_n, parallel)
  write_table(spec, file.path(root, "tables_r/statistical_audit/specification_results.csv"))
  save_fig(plot_specification_curve(spec),
           file.path(root, "figures_r/statistical_audit/specification_curve.png"))
  save_fig(plot_specification_curve(spec),
           file.path(root, "figures_r/statistical_audit/fig10_specification_curve.png"))
  save_fig(plot_specification_heatmap(spec),
           file.path(root, "figures_r/statistical_audit/specification_heatmap.png"))
  cat(sprintf("[model_sensitivity] %d specs; %.1f%% in reference region\n",
              nrow(spec), 100 * mean(spec$in_reference_region, na.rm = TRUE)))
  invisible(spec)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_model_sensitivity\\.R$", .invoked_file)) {
  main()
}
