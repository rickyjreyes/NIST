#!/usr/bin/env Rscript
# run_alternative_nulls.R
# ---------------------------------------------------------------------------
# Adversarial scan-global null-model stress test for the primary Fe II scan.
#
# The canonical audit uses a fitted-smooth Poisson parametric bootstrap with the
# observed Gaussian baseline held fixed.  This module deliberately broadens the
# null family and, where indicated, re-estimates the smoothing baseline inside
# every synthetic replicate.  The goal is not to select the most favourable
# null; it is to report the WORST (largest) calibrated scan-global p-value across
# a declared set of reasonable alternatives.
#
# Declared nulls:
#   1. poisson_fixed_baseline      -- canonical fitted-Poisson bootstrap
#   2. poisson_refit_baseline      -- same generator, but refit smooth baseline
#   3. conditional_multinomial     -- fixes total line count, refits baseline
#   4. negative_binomial_refit     -- overdispersed count null, refits baseline
#   5. block_residual_refit        -- block-permuted Pearson residual structure,
#                                     converted to a rate perturbation, then
#                                     Poisson sampled and baseline-refit
#
# The block-residual construction is an adversarial diagnostic rather than an
# exact generative model.  It preserves short-range residual chunks while
# destroying their global ordering, so it tests sensitivity to local baseline
# misspecification/correlation without preserving a coherent long-range phase.
#
# Writes:
#   tables_r/statistical_audit/alternative_null_results.csv
#   tables_r/statistical_audit/alternative_null_maxima.csv
#   figures_r/statistical_audit/fig17_alternative_nulls.png
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) {
  source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this),
                   "audit_utils.R"))
}

estimate_nb_size <- function(y, mu0) {
  mu0 <- pmax(as.numeric(mu0), EPS)
  y <- as.numeric(y)
  excess <- sum((y - mu0)^2) - sum(mu0)
  if (!is.finite(excess) || excess <= EPS) return(Inf)
  size <- sum(mu0^2) / excess
  if (!is.finite(size) || size <= EPS) return(Inf)
  size
}

block_residual_mean <- function(y, mu0, block_len = 8L) {
  y <- as.numeric(y); mu0 <- pmax(as.numeric(mu0), EPS)
  block_len <- max(1L, as.integer(block_len))
  r <- (y - mu0) / sqrt(mu0)
  r <- r - mean(r, na.rm = TRUE)
  blocks <- split(seq_along(r), ceiling(seq_along(r) / block_len))
  ord <- sample.int(length(blocks), length(blocks), replace = FALSE)
  rp <- unlist(lapply(blocks[ord], function(ii) r[ii]), use.names = FALSE)
  rp <- rp[seq_along(r)]
  mu_star <- pmax(mu0 + rp * sqrt(mu0), EPS)
  # Preserve the fitted expected total count so only residual structure changes.
  if (sum(mu_star) > 0) mu_star <- mu_star * sum(mu0) / sum(mu_star)
  mu_star
}

simulate_counts <- function(model, y, mu0, nb_size, block_len) {
  if (model %in% c("poisson_fixed_baseline", "poisson_refit_baseline")) {
    return(stats::rpois(length(mu0), mu0))
  }
  if (model == "conditional_multinomial") {
    prob <- pmax(mu0, EPS); prob <- prob / sum(prob)
    return(as.numeric(stats::rmultinom(1L, size = as.integer(round(sum(y))), prob = prob)))
  }
  if (model == "negative_binomial_refit") {
    if (!is.finite(nb_size)) return(stats::rpois(length(mu0), mu0))
    return(stats::rnbinom(length(mu0), mu = mu0, size = nb_size))
  }
  if (model == "block_residual_refit") {
    mu_star <- block_residual_mean(y, mu0, block_len)
    return(stats::rpois(length(mu_star), mu_star))
  }
  stop("unknown null model: ", model)
}

summarize_null_model <- function(model, maxima, observed_stat, B, nb_size = NA_real_,
                                 block_len = NA_integer_) {
  tail <- sum(maxima >= observed_stat, na.rm = TRUE)
  ci <- binom_ci(tail, B)
  data.frame(
    null_model = model,
    B = as.integer(B),
    observed_deltaD = observed_stat,
    tail_count = as.integer(tail),
    scan_global_p = emp_p(tail, B),
    resolution_floor = resolution_floor(B),
    exceedance_rate_ci_lo = unname(ci["lower"]),
    exceedance_rate_ci_hi = unname(ci["upper"]),
    null_median = stats::median(maxima, na.rm = TRUE),
    null_q95 = unname(stats::quantile(maxima, 0.95, na.rm = TRUE)),
    null_q99 = unname(stats::quantile(maxima, 0.99, na.rm = TRUE)),
    null_max = max(maxima, na.rm = TRUE),
    nb_size = nb_size,
    block_len = block_len,
    stringsAsFactors = FALSE
  )
}

run_one_model <- function(model, res, k_grid, B, seed, block_len = 8L,
                          parallel = FALSE) {
  setup_rng(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  nb_size <- estimate_nb_size(res$y, res$mu0)
  refit <- model != "poisson_fixed_baseline"
  one <- function(i) {
    set.seed(seeds[i])
    y0 <- simulate_counts(model, res$y, res$mu0, nb_size, block_len)
    base0 <- if (refit) pmax(gaussian_filter_nearest(y0, res$baseline_sigma), EPS) else res$baseline
    sk <- scan_k(res$ell, y0, base0, k_grid, res$degree)
    sk$best$deltaD
  }
  vals <- unlist(audit_lapply(seq_len(B), one, parallel = parallel), use.names = FALSE)
  list(
    summary = summarize_null_model(model, vals, res$best$deltaD, B,
                                   nb_size = if (model == "negative_binomial_refit") nb_size else NA_real_,
                                   block_len = if (model == "block_residual_refit") block_len else NA_integer_),
    maxima = data.frame(null_model = model, replicate = seq_len(B),
                        max_deltaD = vals, stringsAsFactors = FALSE)
  )
}

model_metadata <- function() {
  data.frame(
    null_model = c("poisson_fixed_baseline", "poisson_refit_baseline",
                   "conditional_multinomial", "negative_binomial_refit",
                   "block_residual_refit"),
    description = c(
      "Fitted smooth Poisson null; observed smoothing baseline held fixed (canonical audit null).",
      "Fitted smooth Poisson null; Gaussian smoothing baseline re-estimated in every replicate.",
      "Conditional multinomial null with total line count fixed; smoothing baseline re-estimated.",
      "Negative-binomial overdispersed count null using method-of-moments size; smoothing baseline re-estimated.",
      "Block-permuted Pearson-residual rate perturbation followed by Poisson sampling; smoothing baseline re-estimated."
    ),
    interpretation = c(
      "reference",
      "tests baseline-estimation uncertainty",
      "conditions on observed total count",
      "tests extra-Poisson variance",
      "adversarial local-structure / baseline-misspecification diagnostic"
    ), stringsAsFactors = FALSE)
}

plot_alternative_nulls <- function(maxima, results) {
  lab <- setNames(sprintf("p=%.4g", results$scan_global_p), results$null_model)
  maxima$null_model <- factor(maxima$null_model, levels = rev(results$null_model))
  ggplot2::ggplot(maxima, ggplot2::aes(max_deltaD, null_model)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65, fill = "grey88") +
    ggplot2::geom_vline(xintercept = results$observed_deltaD[1], linetype = 2,
                        colour = species_colour("Fe")) +
    ggplot2::stat_summary(fun = stats::median, geom = "point", size = 2) +
    ggplot2::geom_text(data = results,
      ggplot2::aes(x = Inf, y = factor(null_model, levels = rev(results$null_model)),
                   label = lab[null_model]), inherit.aes = FALSE,
      hjust = 1.05, vjust = -0.5, size = 3) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Adversarial alternative-null scan-global calibration",
      subtitle = "Dashed line = observed Fe II maximum statistic; p-values use (r+1)/(B+1).",
      x = "null replicate maximum deltaD across full k scan", y = NULL,
      caption = "Largest p across declared nulls is the conservative summary; zero exceedances remain resolution-limited.") +
    theme_audit()
}

main <- function(argv = commandArgs(TRUE)) {
  B <- 5000L; block_len <- 8L; seed <- 20260517L
  i <- which(argv == "--null-n"); if (length(i) == 1L) B <- as.integer(argv[i + 1L])
  i <- which(argv == "--block-len"); if (length(i) == 1L) block_len <- as.integer(argv[i + 1L])
  i <- which(argv == "--seed"); if (length(i) == 1L) seed <- as.integer(argv[i + 1L])
  parallel <- any(argv == "--parallel")
  cfg <- default_audit_config(seed = seed, null_n = B, fast = is_fast(argv))
  if (isTRUE(cfg$fast)) B <- 20L
  if (parallel) configure_parallel(TRUE, "auto")
  setup_rng(cfg$seed)

  root <- audit_repo_root()
  lines <- clean_lines_source(read_nist_csv(file.path(root, "data/Fe_lines.csv")),
                              "Fe", 2L, "wavenumber")$lines
  k_grid <- audit_k_grid(cfg)
  res <- run_scan_analysis(lines, cfg$bins_primary, k_grid, cfg$degree, cfg$baseline_sigma)

  models <- model_metadata()$null_model
  model_seeds <- sample.int(.Machine$integer.max, length(models))
  runs <- lapply(seq_along(models), function(j)
    run_one_model(models[j], res, k_grid, B, model_seeds[j], block_len, parallel))

  tab <- do.call(rbind, lapply(runs, `[[`, "summary"))
  tab <- merge(tab, model_metadata(), by = "null_model", sort = FALSE)
  tab <- tab[match(models, tab$null_model), ]
  tab$worst_case <- tab$scan_global_p == max(tab$scan_global_p, na.rm = TRUE)
  maxima <- do.call(rbind, lapply(runs, `[[`, "maxima"))

  write_table(tab, file.path(root, "tables_r/statistical_audit/alternative_null_results.csv"))
  write_table(maxima, file.path(root, "tables_r/statistical_audit/alternative_null_maxima.csv"))
  save_fig(plot_alternative_nulls(maxima, tab),
           file.path(root, "figures_r/statistical_audit/fig17_alternative_nulls.png"),
           width = 10, height = 6.5)

  worst <- tab[which.max(tab$scan_global_p), ]
  cat(sprintf("[alternative_nulls] %d models x %d replicates; worst-case p=%.4g (%s)\n",
              nrow(tab), B, worst$scan_global_p, worst$null_model))
  invisible(tab)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_alternative_nulls\\.R$", .invoked_file)) {
  main()
}
