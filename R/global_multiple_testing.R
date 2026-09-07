#!/usr/bin/env Rscript
# global_multiple_testing.R
# ---------------------------------------------------------------------------
# Multiplicity correction across the FULL DECLARED searched family:
#   A. Bonferroni FWER (valid without independence)
#   B. Holm FWER (valid without independence)
#   C. Benjamini-Hochberg FDR (standard FDR reference; dependence-sensitive)
#   D. Benjamini-Yekutieli FDR (valid under arbitrary dependence)
#   E. Family-wise maximum-statistic calibration
#
# The family includes the complete predefined Fe II preprocessing/specification
# multiverse (bin counts, baseline sigma, polynomial degree, source field) plus
# neighbouring ion-II scans. Duplicate specifications are removed before the
# correction is computed.
#
# DEPENDENCE CAVEAT: many Fe analyses reuse the same transition list. Holm and
# Bonferroni do not require independence; BY is the dependence-robust FDR
# control. The family-max simulation below combines marginal per-analysis
# parametric null draws and is therefore an approximation to a fully joint
# generative null. That approximation is reported rather than hidden.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

assemble_family <- function(cfg) {
  root <- audit_repo_root()
  raw_cache <- new.env()
  get_lines <- function(sp, src) {
    key <- paste(sp, src, sep = "::")
    if (is.null(raw_cache[[key]])) {
      path <- file.path(root, sprintf("data/%s_lines.csv", sp))
      raw_cache[[key]] <- clean_lines_source(read_nist_csv(path), sp, 2L, src)$lines
    }
    raw_cache[[key]]
  }

  fam <- list()
  add <- function(sp, src, bins, sigma, degree, family_component) {
    fam[[length(fam) + 1L]] <<- data.frame(
      species = sp, source = src, bins = as.integer(bins), sigma = as.numeric(sigma),
      degree = as.integer(degree), family_component = family_component,
      stringsAsFactors = FALSE)
  }

  # Predeclared Fe bin grid.
  for (b in cfg$bins_grid)
    add("Fe", "wavenumber", b, cfg$baseline_sigma, cfg$degree, "bin_grid")

  # Full predeclared preprocessing/model multiverse used by run_model_sensitivity.R.
  for (b in c(120L, 160L, 200L))
    for (s in cfg$sigma_grid)
      for (d in cfg$degree_grid)
        add("Fe", "wavenumber", b, s, d, "specification_multiverse")

  # Alternate line-position representations at the primary canonical settings.
  for (src in c("observed", "ritz"))
    add("Fe", src, cfg$bins_primary, cfg$baseline_sigma, cfg$degree, "source_replication")

  # Neighbouring ion-II controls at canonical settings.
  for (sp in c("Cr", "Mn", "Co", "Ni", "Ti")) {
    lines <- get_lines(sp, "wavenumber")
    if (nrow(lines) >= cfg$min_lines)
      add(sp, "wavenumber", cfg$bins_primary, cfg$baseline_sigma, cfg$degree, "neighbour_control")
  }

  famdf <- unique(do.call(rbind, fam)[, c("species", "source", "bins", "sigma", "degree", "family_component")])
  # If the same numerical specification entered through two declared components,
  # collapse it to one hypothesis and retain all component labels.
  key <- with(famdf, paste(species, source, bins, sigma, degree, sep = "|"))
  groups <- split(seq_len(nrow(famdf)), key)
  famdf <- do.call(rbind, lapply(groups, function(ii) {
    x <- famdf[ii[1], , drop = FALSE]
    x$family_component <- paste(sort(unique(famdf$family_component[ii])), collapse = "+")
    x
  }))
  rownames(famdf) <- NULL

  src_tag <- function(x) ifelse(x == "wavenumber", "wn", x)
  famdf$analysis_id <- sprintf(
    "%s_ion2_%s_bin%d_sig%s_deg%d",
    tolower(famdf$species), src_tag(famdf$source), famdf$bins,
    format(famdf$sigma, trim = TRUE, scientific = FALSE), famdf$degree)

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

family_max_null <- function(fam, cfg, family_n, parallel = FALSE) {
  per_analysis <- lapply(seq_len(nrow(fam$famdf)), function(i) {
    o <- fam$obs[[i]]; g <- fam$famdf[i, ]
    nd <- null_distribution(o$ell, o$y, o$baseline, o$mu0, fam$k_grid, g$degree,
                            family_n, parallel = parallel)
    nd$max_vals
  })
  M <- do.call(cbind, per_analysis)
  list(per_analysis = per_analysis, family_max = apply(M, 1L, max))
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
  famdf <- fam$famdf
  fam_max <- fmn$family_max

  famdf$scan_global_p <- vapply(seq_len(nrow(famdf)), function(i)
    emp_p(sum(fmn$per_analysis[[i]] >= famdf$observed_stat[i]), family_n), numeric(1))
  famdf$family_max_p <- vapply(famdf$observed_stat, function(s)
    emp_p(sum(fam_max >= s), family_n), numeric(1))

  p <- famdf$scan_global_p
  famdf$bonferroni_p <- pmin(1, p * length(p))
  famdf$holm_p <- stats::p.adjust(p, method = "holm")
  famdf$bh_fdr <- stats::p.adjust(p, method = "BH")
  famdf$by_fdr <- stats::p.adjust(p, method = "BY")
  famdf$family_size <- length(p)
  famdf$dependence_note <- paste(
    "shared transition lists create dependence; Holm/Bonferroni valid without independence;",
    "BY controls FDR under arbitrary dependence; family-max uses approximate marginal-null coupling")

  # Put the canonical primary Fe II row first for stable downstream selection.
  is_primary <- with(famdf, species == "Fe" & source == "wavenumber" & bins == cfg$bins_primary &
                            sigma == cfg$baseline_sigma & degree == cfg$degree)
  famdf <- famdf[order(!is_primary, famdf$species, famdf$source, famdf$bins, famdf$sigma, famdf$degree), ]
  rownames(famdf) <- NULL

  write_table(famdf, file.path(root, "tables_r/statistical_audit/multiple_testing.csv"))
  write_table(data.frame(family_max_deltaD = fam_max),
              file.path(root, "tables_r/statistical_audit/family_max_null.csv"))

  long <- rbind(
    data.frame(analysis_id = famdf$analysis_id, type = "scan-global", p = famdf$scan_global_p),
    data.frame(analysis_id = famdf$analysis_id, type = "BH-FDR", p = famdf$bh_fdr),
    data.frame(analysis_id = famdf$analysis_id, type = "BY-FDR", p = famdf$by_fdr),
    data.frame(analysis_id = famdf$analysis_id, type = "Holm", p = famdf$holm_p),
    data.frame(analysis_id = famdf$analysis_id, type = "Bonferroni", p = famdf$bonferroni_p),
    data.frame(analysis_id = famdf$analysis_id, type = "family-max", p = famdf$family_max_p))
  ord <- famdf$analysis_id[order(famdf$observed_stat)]
  long$analysis_id <- factor(long$analysis_id, levels = ord)
  pl <- ggplot2::ggplot(long, ggplot2::aes(p, analysis_id, colour = type, shape = type)) +
    ggplot2::geom_point(size = 2.1, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0.05, linetype = 3, colour = "grey40") +
    viridis::scale_colour_viridis(discrete = TRUE, end = 0.85, name = NULL) +
    ggplot2::labs(
      title = "Multiplicity correction across the full declared analysis family",
      subtitle = sprintf("family size = %d; family-max from %d realisations; dotted line = 0.05",
                         length(p), family_n),
      x = "adjusted/global p-value", y = NULL,
      caption = "Holm/Bonferroni need no independence; BY is dependence-robust FDR; family-max coupling is approximate.") +
    theme_audit()
  save_fig(pl, file.path(root, "figures_r/statistical_audit/multiplicity_comparison.png"), height = 10)
  save_fig(pl, file.path(root, "figures_r/statistical_audit/fig12_multiplicity_comparison.png"), height = 10)

  cat(sprintf("[multiple_testing] family=%d; primary family-max p=%.4g; primary Holm=%.4g; primary BY=%.4g\n",
              length(p), famdf$family_max_p[1], famdf$holm_p[1], famdf$by_fdr[1]))
  invisible(famdf)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("global_multiple_testing\\.R$", .invoked_file)) {
  main()
}
