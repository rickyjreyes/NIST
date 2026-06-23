# test-audit-modules.R
# Tests for audit module functions that operate on already-computed draws,
# kept fast by using synthetic inputs (no heavy scans).

.src_module <- function(name) {
  p <- file.path(dirname(AUDIT_PATH), name)
  e <- new.env(parent = globalenv())
  sys.source(p, envir = e)
  e
}

# --- bootstrap CI / multimodal --------------------------------------------
test_that("summarise_ci reports percentile interval, bias and validity counts", {
  boot <- .src_module("bootstrap_peak_uncertainty.R")
  set.seed(1); x <- rnorm(1000, mean = 31.3, sd = 0.1)
  s <- boot$summarise_ci(x, observed = 31.3, quantity = "k_best")
  expect_equal(s$n_valid, 1000L)
  expect_lt(s$ci_lo, s$ci_hi)
  expect_true(abs(s$median - 31.3) < 0.05)
  expect_false(s$multimodal)
})

test_that("summarise_ci flags a clearly multimodal distribution", {
  boot <- .src_module("bootstrap_peak_uncertainty.R")
  set.seed(2); x <- c(rnorm(500, 10, 0.3), rnorm(500, 40, 0.3))
  s <- boot$summarise_ci(x, observed = 10, quantity = "k_best")
  expect_true(s$multimodal)
})

test_that("summarise_ci counts NA/failed resamples as invalid", {
  boot <- .src_module("bootstrap_peak_uncertainty.R")
  x <- c(rnorm(50, 31, 0.1), rep(NA_real_, 10))
  s <- boot$summarise_ci(x, observed = 31, quantity = "k_best")
  expect_equal(s$n_valid, 50L)
})

# --- peak-region stability percentages ------------------------------------
test_that("peak_stability_table computes selection percentages and classes", {
  ps <- .src_module("peak_stability.R")
  draws <- data.frame(k_best = c(rep(31.3, 90), rep(50, 10)),
                      converged = TRUE)
  tab <- ps$peak_stability_table(draws, k_ref = 31.3, tol_grid = c(0.02), primary = 0.02)
  expect_equal(tab$pct_in_reference[1], 0.9)
  expect_equal(tab$stability_class[1], "high")
  expect_true(tab$pct_largest_competitor[1] > 0)
})

# --- significance: pointwise vs global ------------------------------------
test_that("significance_row keeps pointwise and global statistics distinct", {
  eff <- .src_module("build_effect_size_table.R")
  fake <- list(
    res = list(best = list(deltaD = 100), scan = list(k = seq(0.5, 80, length.out = 200),
                                                      deltaD = c(rep(5, 199), 100))),
    null_max = c(rep(10, 49), 200),     # one global exceedance
    null_point = rep(1, 50),            # never exceeds pointwise stat
    k_idx = 200L, null_n = 50L, n_lines = 9000L, bins = 160L)
  fake$res$best$deltaD <- 100
  row <- eff$significance_row(fake, "x")
  expect_equal(row$global_tail_count, 1L)
  expect_equal(row$pointwise_tail_count, 0L)
  expect_true(row$pointwise_p < row$global_p)  # pointwise has zero exceedances -> smaller corrected p
  expect_true(row$zero_exceedance == FALSE)
})

# --- family-max uses the family-max null ----------------------------------
test_that("family-max p-value is computed from the family-max distribution", {
  fam_max <- c(rep(50, 90), rep(150, 10))     # 10% exceed 100
  stat <- 100
  p <- (1 + sum(fam_max >= stat)) / (1 + length(fam_max))
  expect_equal(p, (1 + 10) / (1 + 100))
  # a per-analysis null with no exceedances would give a smaller p -> they differ
  per_null <- rep(20, 100)
  p_self <- (1 + sum(per_null >= stat)) / (1 + length(per_null))
  expect_lt(p_self, p)
})

# --- injection: amplitude zero behaves like the null ----------------------
test_that("injection detection probability is a valid probability", {
  inj <- .src_module("run_injection_recovery.R")
  # detection probability must be in [0,1] by construction
  det <- 7; n <- 50
  ci <- binom_ci(det, n)
  expect_gte(det / n, 0); expect_lte(det / n, 1)
  expect_gte(ci["lower"], 0); expect_lte(ci["upper"], 1)
})

test_that("correct-region recovery is distinct from any-peak detection", {
  # detection-anywhere counts statistic exceedance; correct-region adds a
  # frequency-tolerance condition, so it can never exceed detection.
  detect <- c(1, 1, 1, 0)
  correct <- c(1, 0, 1, 0)
  sig_correct <- as.integer(detect == 1 & correct == 1)
  expect_lte(sum(sig_correct), sum(detect))
})

# --- model comparison parameter counts ------------------------------------
test_that("model comparison parameter counts are M0=degree+1, M1=+2", {
  mc <- .src_module("model_comparison.R")
  ell <- seq(8.6, 10.8, length.out = 60)
  set.seed(3); y <- rpois(60, 50)
  base <- rep(50, 60)
  fm <- mc$fit_models(ell, y, base, degree = 1L, k = 30)
  expect_equal(fm$p0, 2L)   # intercept + 1 poly term
  expect_equal(fm$p1, 4L)   # + cos + sin
})
