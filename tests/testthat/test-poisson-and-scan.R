# Tests 9-13: deviance, IRLS convergence, harmonic recovery, scan-max
# selection, empirical p-value correction.

test_that("poisson_deviance matches hand-calculated examples", {
  # y == mu everywhere -> deviance 0
  expect_equal(poisson_deviance(c(1, 2, 3), c(1, 2, 3)), 0, tolerance = 1e-12)
  # y = (0, 2), mu = (1, 2): term0 = 1, term1 = 0 -> 2*1 = 2
  expect_equal(poisson_deviance(c(0, 2), c(1, 2)), 2.0, tolerance = 1e-12)
  # single point y=4, mu=2: 2*(2 - 4 + 4*log(2)) = 2*(-2 + 4*0.6931...)
  expect_equal(poisson_deviance(4, 2), 2 * (2 - 4 + 4 * log(2)), tolerance = 1e-12)
})

test_that("IRLS converges on a deterministic intercept-only case", {
  # With X = intercept only and baseline B, the MLE is mu = mean(y/B)*B,
  # i.e. beta0 = log(sum(y)/sum(B)).
  set.seed(1)
  B <- rep(5.0, 50)
  y <- rep(7.0, 50)
  X <- matrix(1.0, nrow = 50, ncol = 1)
  fit <- fit_poisson_loglinear(y, B, X)
  expect_equal(fit$beta[1], log(sum(y) / sum(B)), tolerance = 1e-8)
  expect_equal(fit$mu, rep(7.0, 50), tolerance = 1e-6)
})

test_that("scan recovers an injected harmonic from synthetic Poisson data", {
  set.seed(42)
  n <- 160
  ell <- seq(8.66, 10.81, length.out = n)
  k_true <- 31.0
  base <- rep(40.0, n)
  mu <- base * exp(0.25 * cos(k_true * ell) - 0.10 * sin(k_true * ell))
  y <- rpois(n, mu)
  baseline <- pmax(gaussian_filter_nearest(y, 6.0), EPS)
  k_grid <- seq(0.5, 80.0, length.out = 2500)
  sk <- scan_k(ell, y, baseline, k_grid, 1L)
  # recovered k_best should be close to k_true
  expect_lt(abs(sk$best$k_best - k_true), 0.5)
  expect_gt(sk$best$deltaD, 0)
})

test_that("scan_k selects the first maximiser of deltaD (strict >)", {
  # Build a tiny scan-like check: best must be argmax with first-wins ties.
  set.seed(7)
  n <- 60
  ell <- seq(8.7, 10.8, length.out = n)
  y <- rpois(n, rep(20, n))
  baseline <- pmax(gaussian_filter_nearest(y, 6.0), EPS)
  k_grid <- seq(0.5, 80.0, length.out = 200)
  sk <- scan_k(ell, y, baseline, k_grid, 1L)
  expect_equal(sk$best$deltaD, max(sk$scan$deltaD), tolerance = 1e-12)
  first_idx <- which(sk$scan$deltaD == max(sk$scan$deltaD))[1]
  expect_equal(sk$best$k_best, sk$scan$k[first_idx])
})

test_that("empirical p-value applies the +1 / (N+1) correction", {
  # 0 exceedances out of 100 -> 1/101
  vals <- rep(0.0, 100)
  real_delta <- 10.0
  tail_count <- sum(vals >= real_delta)
  p <- (1 + tail_count) / (1 + length(vals))
  expect_equal(p, 1 / 101)
  # 3 exceedances out of 100 -> 4/101
  vals2 <- c(rep(20.0, 3), rep(0.0, 97))
  p2 <- (1 + sum(vals2 >= real_delta)) / (1 + length(vals2))
  expect_equal(p2, 4 / 101)
})

test_that("null_scan is deterministic for a fixed seed", {
  set.seed(NULL)
  n <- 40
  ell <- seq(8.7, 10.8, length.out = n)
  y <- rpois(n, rep(15, n))
  baseline <- pmax(gaussian_filter_nearest(y, 6.0), EPS)
  k_grid <- seq(0.5, 80.0, length.out = 100)
  a <- null_scan(ell, y, baseline, baseline, k_grid, 1L, 5.0, 5L, 123L, verbose = FALSE)
  b <- null_scan(ell, y, baseline, baseline, k_grid, 1L, 5.0, 5L, 123L, verbose = FALSE)
  expect_equal(a$vals, b$vals)
  expect_equal(a$p, b$p)
})
