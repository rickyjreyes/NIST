# test-final-adversarial-audit.R
# Fast unit tests for the final adversarial audit. These tests exercise pure
# helpers and small synthetic vectors only; they do not run Monte Carlo scans.

.src_final_module <- function(name) {
  p <- file.path(dirname(AUDIT_PATH), name)
  e <- new.env(parent = globalenv())
  sys.source(p, envir = e)
  e
}

test_that("final audit entry points source without executing their main routines", {
  for (name in c("global_multiple_testing.R", "run_alternative_nulls.R",
                 "run_holdout_replication.R", "build_final_adversarial_summary.R",
                 "run_final_adversarial_audit.R")) {
    e <- .src_final_module(name)
    expect_true(exists("main", envir = e, inherits = FALSE), info = name)
  }
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

test_that("primary multiplicity row recognises canonical metadata and legacy ids", {
  s <- .src_final_module("build_final_adversarial_summary.R")
  x <- data.frame(
    analysis_id = c("fe_ion2_wn_bin160_sig5_deg1", "fe_ion2_wn_bin160_sig6_deg1"),
    species = c("Fe", "Fe"), source = c("wavenumber", "wavenumber"),
    bins = c(160, 160), sigma = c(5, 6), degree = c(1, 1),
    family_max_p = c(0.2, 0.01), stringsAsFactors = FALSE)
  row <- s$primary_row(x)
  expect_equal(row$analysis_id, "fe_ion2_wn_bin160_sig6_deg1")
  expect_equal(row$family_max_p, 0.01)

  legacy <- data.frame(
    analysis_id = c("fe_ion2_wn_bingrid120", "fe_ion2_wn_bingrid160"),
    family_max_p = c(0.2, 0.01), stringsAsFactors = FALSE)
  old <- s$primary_row(legacy)
  expect_equal(old$analysis_id, "fe_ion2_wn_bingrid160")
})

test_that("baseline-refit holdout test requires a smoothing sigma", {
  ho <- .src_final_module("run_holdout_replication.R")
  expect_error(
    ho$fixed_k_test_p(c(1, 2, 3), c(4, 5, 6), c(4, 5, 6),
                      degree = 1, k = 2, B = 2, refit_baseline = TRUE),
    "sigma is required"
  )
})
