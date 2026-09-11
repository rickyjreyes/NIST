# test-final-adversarial-audit.R
# Fast unit tests for the final adversarial audit. These tests exercise pure
# helpers and small synthetic vectors only; they do not run Monte Carlo scans.

.src_final_module <- function(name) {
  p <- file.path(dirname(AUDIT_PATH), name)
  e <- new.env(parent = globalenv())
  sys.source(p, envir = e)
  e
}

test_that("final audit modules source and CLI runner parses without executing", {
  # Library-style modules are safe to source directly under testthat.
  for (name in c("global_multiple_testing.R", "run_alternative_nulls.R",
                 "run_holdout_replication.R", "build_final_adversarial_summary.R",
                 "run_resolution_mode_diagnostics.R", "model_comparison.R")) {
    e <- .src_final_module(name)
    expect_true(exists("main", envir = e, inherits = FALSE), info = name)
  }

  # The top-level CLI runner intentionally resolves repository paths at source
  # time. Under testthat the working directory differs from normal CLI use, so
  # syntax/parse safety is the appropriate no-side-effect CI check here.
  runner <- file.path(dirname(AUDIT_PATH), "run_final_adversarial_audit.R")
  expect_silent(parse(file = runner))
})

test_that("negative-binomial size reduces to Poisson when no overdispersion is estimated", {
  alt <- .src_final_module("run_alternative_nulls.R")
  mu <- rep(20, 50)
  y <- rep(20, 50)
  expect_true(is.infinite(alt$estimate_nb_size(y, mu)))
})

test_that("negative-binomial size is finite under clear extra-Poisson dispersion", {
  alt <- .src_final_module("run_alternative_nulls.R")
  mu <- rep(20, 100)
  y <- rep(c(2, 38), 50)
  size <- alt$estimate_nb_size(y, mu)
  expect_true(is.finite(size))
  expect_gt(size, 0)
})

test_that("block-residual adversarial mean remains positive and preserves expected total", {
  alt <- .src_final_module("run_alternative_nulls.R")
  setup_rng(123)
  mu <- seq(10, 30, length.out = 40)
  y <- round(mu + rep(c(-3, 2, 4, -2), 10))
  star <- alt$block_residual_mean(y, mu, block_len = 5L)
  expect_length(star, length(mu))
  expect_true(all(is.finite(star)))
  expect_true(all(star > 0))
  expect_equal(sum(star), sum(mu), tolerance = 1e-10)
})

test_that("alternative-null summary uses plus-one empirical p and never reports zero", {
  alt <- .src_final_module("run_alternative_nulls.R")
  maxima <- rep(10, 99)
  s <- alt$summarize_null_model("x", maxima, observed_stat = 100, B = 99)
  expect_equal(s$tail_count, 0L)
  expect_equal(s$scan_global_p, 1 / 100)
  expect_equal(s$resolution_floor, 1 / 100)
  expect_gt(s$scan_global_p, 0)
})

test_that("alternative-null model registry contains the canonical and adversarial generators", {
  alt <- .src_final_module("run_alternative_nulls.R")
  models <- alt$model_metadata()$null_model
  expect_setequal(models, c(
    "poisson_fixed_baseline",
    "poisson_refit_baseline",
    "conditional_multinomial",
    "negative_binomial_refit",
    "block_residual_refit"
  ))
})

test_that("strict holdout verdict cannot hide a partial failure", {
  s <- .src_final_module("build_final_adversarial_summary.R")
  expect_equal(s$holdout_verdict(c(0.01, 0.02), c(TRUE, TRUE)), "pass")
  expect_equal(s$holdout_verdict(c(0.01, 0.20), c(TRUE, TRUE)), "mixed")
  expect_equal(s$holdout_verdict(c(0.20, 0.30), c(TRUE, TRUE)), "fail")
  expect_equal(s$holdout_verdict(c(0.01, 0.02), c(TRUE, FALSE)), "mixed")
})

test_that("descriptive robustness thresholds are explicit pass/fail rather than softened", {
  s <- .src_final_module("build_final_adversarial_summary.R")
  expect_equal(s$threshold_verdict(0.80, 0.80, TRUE), "pass")
  expect_equal(s$threshold_verdict(0.70, 0.80, TRUE), "fail")
  expect_equal(s$threshold_verdict(0.049, 0.05, FALSE), "pass")
  expect_equal(s$threshold_verdict(0.051, 0.05, FALSE), "fail")
})

test_that("primary multiplicity row recognises canonical metadata", {
  s <- .src_final_module("build_final_adversarial_summary.R")
  x <- data.frame(
    analysis_id = c("fe_ion2_wn_bin160_sig5_deg1", "fe_ion2_wn_bin160_sig6_deg1"),
    species = c("Fe", "Fe"), source = c("wavenumber", "wavenumber"),
    bins = c(160, 160), sigma = c(5, 6), degree = c(1, 1),
    family_max_p = c(0.2, 0.01), stringsAsFactors = FALSE)
  row <- s$primary_row(x)
  expect_equal(row$analysis_id, "fe_ion2_wn_bin160_sig6_deg1")
  expect_equal(row$family_max_p, 0.01)
})

test_that("baseline-refit holdout test requires a smoothing sigma", {
  ho <- .src_final_module("run_holdout_replication.R")
  expect_error(
    ho$fixed_k_test_p(c(1, 2, 3), c(4, 5, 6), c(4, 5, 6),
                      degree = 1, k = 2, B = 2, refit_baseline = TRUE),
    "sigma is required"
  )
})

test_that("held-out polynomial prediction reuses the training transform", {
  mc <- .src_final_module("model_comparison.R")
  tr <- c(1, 2, 3, 4)
  te <- c(10, 11)
  tf <- mc$poly_transform(tr)

  # The explicit frozen transform must reproduce the canonical design on train.
  expect_equal(mc$design_poly_from_transform(tr, 2, tf), design_poly(tr, 2),
               tolerance = 1e-12)

  # Test points are mapped with TRAINING center/scale, not their own statistics.
  expected_z <- (te - mean(tr)) / sqrt(mean((tr - mean(tr))^2))
  x_te <- mc$design_poly_from_transform(te, 1, tf)
  expect_equal(x_te[, 2], expected_z, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(x_te[, 2], design_poly(te, 1)[, 2])))
})

test_that("matched smoothing holds approximate ell-space width fixed", {
  d <- .src_final_module("run_resolution_mode_diagnostics.R")
  expect_equal(d$matched_sigma(c(60, 80, 100, 160, 240)),
               c(2.25, 3.00, 3.75, 6.00, 9.00))
})

test_that("resolution diagnostic extracts distinct ranked local peaks", {
  d <- .src_final_module("run_resolution_mode_diagnostics.R")
  scan <- data.frame(k = 1:7, deltaD = c(0, 3, 1, 4, 1, 2, 0))
  peaks <- d$extract_top_peaks(scan, 3L)
  expect_equal(peaks$k, c(4, 2, 6))
  expect_equal(peaks$peak_rank, 1:3)
})

test_that("resolution branch labels distinguish frozen external and Fe primary modes", {
  d <- .src_final_module("run_resolution_mode_diagnostics.R")
  expect_equal(d$classify_peak_branch(9.63, d$GWTC_K_FROZEN, 31.3265, 0.02),
               "external_gwtc_9p602")
  expect_equal(d$classify_peak_branch(31.30, d$GWTC_K_FROZEN, 31.3265, 0.02),
               "primary_fe")
  expect_equal(d$classify_peak_branch(20, d$GWTC_K_FROZEN, 31.3265, 0.02),
               "other")
})
