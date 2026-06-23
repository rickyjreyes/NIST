#!/usr/bin/env Rscript
# global_multiple_testing.R
# ---------------------------------------------------------------------------
# Multiplicity correction across the DECLARED analysis family (read from the
# registry). Three p-value-based procedures plus a family-wise max-statistic
# procedure:
#   A. Benjamini-Hochberg FDR
#   B. Holm
#   C. Bonferroni
#   D. Family-wise maximum-statistic
#
# The family-max procedure simulates, for each family-level realisation, a null
# scan of every searched analysis and records the LARGEST statistic across the
# family (over species, ions, bins, sources, baseline settings, ...). Each
# observed statistic is then compared to that family-max null distribution.
#
# DEPENDENCE CAVEAT (written to output): the analyses share the same underlying
# NIST line list and are therefore NOT independent. Bonferroni/Holm treat them
# as independent and are conservative; the family-max procedure here combines
# per-analysis parametric nulls and is an approximation of the fully-joint null.
# Both facts are reported rather than hidden.
#
# Writes:
#   tables_r/statistical_audit/multiple_testing.csv
#   tables_r/statistical_audit/family_max_null.csv
#   figures_r/statistical_audit/multiplicity_comparison.png (+ fig12)
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

# assemble_family(): observed per-analysis global statistics for the searched
# family. Each row: species, ion, source, bins, sigma, degree, observed deltaD.
assemble_family <- function(cfg) {
  root <- audit_repo_root()
  raw_cache <- new.env()
  get_lines <- function(sp, src) {
    key <- paste(sp, src)
    if (is.null(raw_cache[[key]])) {
      path <- file.path(root, sprintf("data/%s_lines.csv", sp))
      raw_cache[[key]] <- clean_lines_source(read_nist_csv(path), sp, 2L, src)$lines
    }
    raw_cache[[key]]
  }
  fam <- list()
  add <- function(id, sp, src, bins, sigma, degree)
    fam[[length(fam) + 1L]] <<- data.frame(analysis_id = id, species = sp, source = src,
                                           bins = bins, sigma = sigma, degree = degree,
                                           stringsAsFactors = FALSE)
  # bin grid (Fe wn)
  for (b in cfg$bins_grid) add(sprintf("fe_ion2_wn_bingrid%d", b), "Fe", "wavenumber", b, cfg$baseline_sigma, cfg$degree)
  # sources (Fe)
  for (s in c("observed", "ritz")) add(sprintf("fe_ion2_%s_bin%d", s, cfg$bins_primary), "Fe", s, cfg$bins_primary, cfg$baseline_sigma, cfg$degree)
  # neighbours
  for (sp in c("Cr", "Mn", "Co", "Ni", "Ti")) {
    lines <- get_lines(sp, "wavenumber")
    if (nrow(lines) >= cfg$min_lines)
      add(sprintf("%s_ion2_wn_bin%d", tolower(sp), cfg$bins_primary), sp, "wavenumber", cfg$bins_primary, cfg$baseline_sigma, cfg$degree)
  }
  famdf <- do.call(rbind, fam)
  # compute observed deltaD + mu0 for each
  k_grid <- audit_k_grid(cfg)
  obs <- lapply(seq_len(nrow(famdf)), function(i) {
    g <- famdf[i, ]
    lines <- get_lines(g$species, g$source)
    res <- run_scan_analysis(lines, g$bins, k_grid, g$degree, g$sigma)
    list(deltaD = res$best$deltaD, ell = res$ell, y = res$y,
         baseline = res$baseline, mu0 = res$mu0)
  })
  famdf$observed_stat <- vapply(obs, function(o) o$deltaD, numeric(1))
  list(famdf = famdf, obs = obs, k_grid = k_grid)
}

# family_max_null(): per-analysis parametric null-max distributions, combined
# into the family-max distribution by independent-draw maximum across analyses.
family_max_null <- function(fam, cfg, family_n, parallel = FALSE) {
  k_grid <- fam$k_grid
  per_analysis <- lapply(seq_len(nrow(fam$famdf)), function(i) {
    o <- fam$obs[[i]]; g <- fam$famdf[i, ]
    nd <- null_distribution(o$ell, o$y, o$baseline, o$mu0, k_grid, g$degree,
                            family_n, parallel = parallel)
    nd$max_vals
  })
  M <- do.call(cbind, per_analysis)         # family_n x n_analyses
  fam_max <- apply(M, 1L, max)
  list(per_analysis = per_analysis, family_max = fam_max)
}

main <- function(argv = commandArgs(TRUE)) {
  family_n <- 200L
  i <- which(argv == "--family-n"); if (length(i) == 1L) family_n <- as.integer(argv[i + 1L])
  j <- which(argv == "--null-n"); if (length(j) == 1L) family_n <- as.integer(argv[j + 1L])
  parallel <- any(argv == "--parallel")
  cfg <- default_audit_config()
  if (parallel) configure_parallel(TRUE, "auto")
  setup_rng(cfg$seed)
  root <- audit_repo_root()

  fam <- assemble_family(cfg)
  fmn <- family_max_null(fam, cfg, family_n, parallel)
  fam_max <- fmn$family_max

  famdf <- fam$famdf
  # per-analysis scan-global p (own null) and family-adjusted p (family-max null)
  famdf$scan_global_p <- vapply(seq_len(nrow(famdf)), function(i) {
    emp_p(sum(fmn$per_analysis[[i]] >= famdf$observed_stat[i]), family_n)
  }, numeric(1))
  famdf$family_max_p <- vapply(famdf$observed_stat, function(s)
    emp_p(sum(fam_max >= s), family_n), numeric(1))

  p <- famdf$scan_global_p
  famdf$bonferroni_p <- pmin(1, p * length(p))
  famdf$holm_p <- stats::p.adjust(p, method = "holm")
  famdf$bh_fdr <- stats::p.adjust(p, method = "BH")
  famdf$family_size <- length(p)
  famdf$dependence_note <- "shared NIST line list; corrections approximate, see header"

  write_table(famdf, file.path(root, "tables_r/statistical_audit/multiple_testing.csv"))
  write_table(data.frame(family_max_deltaD = fam_max),
              file.path(root, "tables_r/statistical_audit/family_max_null.csv"))

  # comparison figure: raw vs scan-global vs family-adjusted for top analyses
  long <- rbind(
    data.frame(analysis_id = famdf$analysis_id, type = "scan-global", p = famdf$scan_global_p),
    data.frame(analysis_id = famdf$analysis_id, type = "BH-FDR", p = famdf$bh_fdr),
    data.frame(analysis_id = famdf$analysis_id, type = "Holm", p = famdf$holm_p),
    data.frame(analysis_id = famdf$analysis_id, type = "family-max", p = famdf$family_max_p))
  ord <- famdf$analysis_id[order(famdf$observed_stat)]
  long$analysis_id <- factor(long$analysis_id, levels = ord)
  pl <- ggplot2::ggplot(long, ggplot2::aes(p, analysis_id, colour = type, shape = type)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0.05, linetype = 3, colour = "grey40") +
    viridis::scale_colour_viridis(discrete = TRUE, end = 0.85, name = NULL) +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
    ggplot2::labs(
      title = "Multiplicity-corrected p-values across the declared family",
      subtitle = sprintf("family size = %d; family-max from %d realisations. Dotted line = 0.05.",
                         length(p), family_n),
      x = "p-value", y = NULL,
      caption = "Analyses share one line list; corrections are approximate (see table header).") +
    theme_audit()
  save_fig(pl, file.path(root, "figures_r/statistical_audit/multiplicity_comparison.png"))
  save_fig(pl, file.path(root, "figures_r/statistical_audit/fig12_multiplicity_comparison.png"))

  cat(sprintf("[multiple_testing] family=%d; min BH-FDR=%.4g; min family-max p=%.4g\n",
              length(p), min(famdf$bh_fdr), min(famdf$family_max_p)))
  invisible(famdf)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("global_multiple_testing\\.R$", .invoked_file)) {
  main()
}
