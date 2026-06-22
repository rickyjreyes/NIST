#!/usr/bin/env Rscript
# build_effect_size_table.R
# ---------------------------------------------------------------------------
# Significance (pointwise vs scan-global) AND interpretable effect sizes for
# Fe II and the neighbouring ions. Both are produced here because they share
# one parametric-null run per species (cheaper and guaranteed consistent).
#
# Significance:
#   - pointwise p at the FIXED reference k (prespecified), null = deltaD at
#     that single k index over null draws;
#   - scan-global p using the maximum statistic across the WHOLE k grid in
#     every null realisation (look-elsewhere corrected).
#   - empirical p = (r+1)/(B+1); zero exceedances -> resolution floor 1/(B+1).
#
# Effect sizes (NOT a p-value converted to sigma):
#   - observed deltaD, excess over null-max mean, null-standardized z (with an
#     explicit non-Gaussian caveat), amplitude, deltaD per line, deltaD per
#     bin, and the percentile of the observed statistic in the null.
#
# Writes:
#   tables_r/statistical_audit/significance_results.csv
#   tables_r/statistical_audit/effect_sizes.csv
#   figures_r/statistical_audit/effect_size_forest.png (+ fig07)
#   figures_r/statistical_audit/null_ridgelines.png    (+ fig06)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

# analyze_with_null(): reference scan + parametric null (global max + pointwise
# at the reference k). Returns everything significance & effect size need.
analyze_with_null <- function(lines, bins, cfg, null_n = cfg$null_n) {
  k_grid <- audit_k_grid(cfg)
  res <- run_scan_analysis(lines, bins, k_grid, cfg$degree, cfg$baseline_sigma)
  k_idx <- nearest_k_index(k_grid, res$best$k_best)
  nd <- null_distribution(res$ell, res$y, res$baseline, res$mu0, k_grid,
                          cfg$degree, null_n, fixed_k_index = k_idx)
  list(res = res, null_max = nd$max_vals, null_point = nd$point_vals,
       k_idx = k_idx, n_lines = nrow(lines), bins = bins, null_n = null_n)
}

significance_row <- function(a, analysis_id) {
  obs <- a$res$best$deltaD
  obs_point <- a$res$scan$deltaD[a$k_idx]
  B <- a$null_n
  global_tail <- sum(a$null_max >= obs)
  point_tail <- sum(a$null_point >= obs_point)
  data.frame(
    analysis_id = analysis_id,
    pointwise_statistic = obs_point,
    pointwise_tail_count = point_tail,
    pointwise_B = B,
    pointwise_p = emp_p(point_tail, B),
    global_statistic = obs,
    global_tail_count = global_tail,
    global_B = B,
    global_p = emp_p(global_tail, B),
    empirical_resolution = resolution_floor(B),
    zero_exceedance = (global_tail == 0L),
    scan_k_min = a$res$scan$k[1],
    scan_k_max = a$res$scan$k[length(a$res$scan$k)],
    n_scanned_frequencies = length(a$res$scan$k),
    stringsAsFactors = FALSE
  )
}

effect_size_row <- function(a, species, analysis_id) {
  obs <- a$res$best$deltaD
  nm <- mean(a$null_max); ns <- stats::sd(a$null_max)
  pct <- mean(a$null_max < obs)
  data.frame(
    analysis_id = analysis_id,
    species = species,
    deltaD = obs,
    null_mean = nm,
    null_sd = ns,
    excess_deltaD = obs - nm,
    z_null = if (ns > 0) (obs - nm) / ns else NA_real_,
    amplitude = a$res$best$amplitude,
    deltaD_per_line = obs / a$n_lines,
    deltaD_per_bin = obs / a$bins,
    n_lines = a$n_lines,
    bins = a$bins,
    percentile_in_null = pct,
    null_n = a$null_n,
    stringsAsFactors = FALSE
  )
}

plot_effect_forest <- function(eff) {
  eff <- eff[order(eff$z_null), ]
  eff$species <- factor(eff$species, levels = eff$species)
  ggplot2::ggplot(eff, ggplot2::aes(x = z_null, y = species, colour = species)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = z_null, yend = species),
                          linewidth = 0.6) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = SPECIES_PALETTE, guide = "none") +
    ggplot2::labs(
      title = "Null-standardized effect size by ion (II)",
      subtitle = "z_null = (observed deltaD - null-max mean) / null-max sd. Descriptive distance; null is non-Gaussian.",
      x = "z_null (standardized distance from scan-global null)", y = NULL,
      caption = "Compare alongside line/bin counts; raw deltaD is not comparable across species.") +
    theme_audit()
}

plot_null_ridges <- function(null_list, obs_list) {
  long <- do.call(rbind, lapply(names(null_list), function(sp) {
    data.frame(species = sp, value = null_list[[sp]], stringsAsFactors = FALSE)
  }))
  obs_df <- data.frame(species = names(obs_list), value = unlist(obs_list))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = value, y = species, fill = species))
  if (requireNamespace("ggridges", quietly = TRUE)) {
    p <- p + ggridges::geom_density_ridges(scale = 1.5, alpha = 0.7, colour = "white")
  } else {
    p <- p + ggplot2::geom_violin(ggplot2::aes(group = species), alpha = 0.7)
  }
  p +
    ggplot2::geom_point(data = obs_df, ggplot2::aes(x = value, y = species),
                        colour = "black", shape = 18, size = 3,
                        inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = SPECIES_PALETTE, guide = "none") +
    ggplot2::labs(
      title = "Scan-global null distributions vs observed deltaD",
      subtitle = "Ridges = null max-deltaD per ion; diamonds = observed deltaD.",
      x = "max deltaD over scan (null)", y = NULL) +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  null_n <- NULL
  if (length(argv) > 0L) {
    i <- which(argv == "--null-n"); if (length(i) == 1L) null_n <- as.integer(argv[i + 1L])
  }
  cfg <- default_audit_config(null_n = if (is.null(null_n)) 5000L else null_n)
  setup_rng(cfg$seed)
  root <- audit_repo_root()

  species_list <- c("Fe", "Cr", "Mn", "Co", "Ni", "Ti")
  sig_rows <- list(); eff_rows <- list()
  null_list <- list(); obs_list <- list()

  # Fe II confirmatory significance at bins 120/160/200
  raw_fe <- read_nist_csv(file.path(root, "data/Fe_lines.csv"))
  fe_lines <- clean_lines_source(raw_fe, "Fe", 2L, "wavenumber")$lines
  for (b in c(120L, 160L, 200L)) {
    setup_rng(cfg$seed)
    a <- analyze_with_null(fe_lines, b, cfg)
    id <- sprintf("fe_ion2_wn_bin%d", b)
    sig_rows[[id]] <- significance_row(a, id)
    if (b == cfg$bins_primary) {
      eff_rows[["Fe"]] <- effect_size_row(a, "Fe", id)
      null_list[["Fe"]] <- a$null_max; obs_list[["Fe"]] <- a$res$best$deltaD
    }
    cat(sprintf("[sig] Fe bins=%d global_p=%.4g pointwise_p=%.4g tail=%d\n",
                b, sig_rows[[id]]$global_p, sig_rows[[id]]$pointwise_p,
                sig_rows[[id]]$global_tail_count))
  }

  # neighbouring ions at primary bins (significance + effect size)
  for (sp in setdiff(species_list, "Fe")) {
    path <- file.path(root, sprintf("data/%s_lines.csv", sp))
    if (!file.exists(path)) next
    lines <- clean_lines_source(read_nist_csv(path), sp, 2L, "wavenumber")$lines
    if (nrow(lines) < cfg$min_lines) {
      cat(sprintf("[skip] %s has only %d lines (< %d)\n", sp, nrow(lines), cfg$min_lines)); next
    }
    setup_rng(cfg$seed)
    a <- analyze_with_null(lines, cfg$bins_primary, cfg)
    id <- sprintf("%s_ion2_wn_bin%d", tolower(sp), cfg$bins_primary)
    sig_rows[[id]] <- significance_row(a, id)
    eff_rows[[sp]] <- effect_size_row(a, sp, id)
    null_list[[sp]] <- a$null_max; obs_list[[sp]] <- a$res$best$deltaD
    cat(sprintf("[sig] %s global_p=%.4g z_null=%.2f\n", sp,
                sig_rows[[id]]$global_p, eff_rows[[sp]]$z_null))
  }

  sig <- do.call(rbind, sig_rows); rownames(sig) <- NULL
  eff <- do.call(rbind, eff_rows); rownames(eff) <- NULL
  write_table(sig, file.path(root, "tables_r/statistical_audit/significance_results.csv"))
  write_table(eff, file.path(root, "tables_r/statistical_audit/effect_sizes.csv"))

  pf <- plot_effect_forest(eff)
  save_fig(pf, file.path(root, "figures_r/statistical_audit/effect_size_forest.png"))
  save_fig(pf, file.path(root, "figures_r/statistical_audit/fig07_effect_size_forest.png"))
  pr <- plot_null_ridges(null_list, obs_list)
  save_fig(pr, file.path(root, "figures_r/statistical_audit/null_ridgelines.png"))
  save_fig(pr, file.path(root, "figures_r/statistical_audit/fig06_null_ridgelines.png"))

  cat(sprintf("[effect/significance] %d significance rows, %d effect rows\n",
              nrow(sig), nrow(eff)))
  invisible(list(sig = sig, eff = eff))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("build_effect_size_table\\.R$", .invoked_file)) {
  main()
}
