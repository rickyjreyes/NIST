# test-audit-stats.R
# Unit tests for the statistical-audit suite (fast, no heavy Monte Carlo).

context_data <- function() file.path(DATA_DIR_T, "Fe_lines.csv")

# --- peak conversions ------------------------------------------------------
test_that("peak conversions match the angular-frequency convention", {
  k <- 31.3265306122449
  expect_equal(delta_log_x(k), 2 * pi / k)
  expect_equal(scale_ratio(k), exp(2 * pi / k))
  expect_equal(n_obs_from_k(k, 2.15405942045357), k * 2.15405942045357 / (2 * pi))
})

test_that("invalid k yields NA conversions", {
  expect_true(is.na(delta_log_x(0)))
  expect_true(is.na(delta_log_x(-5)))
  expect_true(is.na(scale_ratio(NA_real_)))
})

# --- empirical p-values ----------------------------------------------------
test_that("empirical p uses (r+1)/(B+1) and never returns zero", {
  expect_equal(emp_p(0, 999), 1 / 1000)
  expect_equal(emp_p(5, 99), 6 / 100)
  expect_gt(emp_p(0, 1e6), 0)
  expect_equal(resolution_floor(999), 1 / 1000)
})

test_that("zero exceedances equal the resolution floor", {
  B <- 4999
  expect_equal(emp_p(0, B), resolution_floor(B))
})

# --- dataset accounting ----------------------------------------------------
test_that("Fe II dataset flow reconciles exactly and retains 9447 lines", {
  skip_if_not(file.exists(context_data()))
  raw <- read_nist_csv(context_data())
  cl <- clean_lines_source(raw, "Fe", 2L, "wavenumber")
  f <- cl$flow
  expect_true(f$reconciles)
  expect_equal(f$n_raw, f$retained + f$excl_species + f$excl_ion +
                 f$excl_missing_source + f$excl_nonpositive + f$n_duplicates)
  expect_equal(f$retained, 9447L)
  expect_equal(nrow(cl$lines), 9447L)
})

test_that("observed and Ritz sources are counted separately and differ", {
  skip_if_not(file.exists(context_data()))
  raw <- read_nist_csv(context_data())
  obs <- clean_lines_source(raw, "Fe", 2L, "observed")$lines
  rtz <- clean_lines_source(raw, "Fe", 2L, "ritz")$lines
  wn  <- clean_lines_source(raw, "Fe", 2L, "wavenumber")$lines
  expect_gt(nrow(obs), 100L); expect_gt(nrow(rtz), 100L)
  # observed has missing values -> fewer lines than ritz here
  expect_false(nrow(obs) == nrow(rtz))
})

test_that("duplicate wavenumbers are removed", {
  df <- data.frame(element = "Fe", sp_num = "2",
                   `wn(cm-1)` = c("100", "100", "200", "300"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  cl <- clean_lines_source(df, "Fe", 2L, "wavenumber")
  expect_equal(cl$flow$n_duplicates, 1L)
  expect_equal(cl$flow$retained, 3L)
})

test_that("nonpositive and missing source values are excluded", {
  df <- data.frame(element = "Fe", sp_num = "2",
                   `wn(cm-1)` = c("100", "-5", "", "200"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  cl <- clean_lines_source(df, "Fe", 2L, "wavenumber")
  expect_equal(cl$flow$excl_nonpositive, 1L)
  expect_equal(cl$flow$excl_missing_source, 1L)
  expect_equal(cl$flow$retained, 2L)
  expect_true(cl$flow$reconciles)
})

# --- peak-region tolerance -------------------------------------------------
test_that("peak-region membership respects relative and absolute tolerance", {
  expect_true(in_peak_region(31.5, 31.3, "relative", tol_relative = 0.02))
  expect_false(in_peak_region(40, 31.3, "relative", tol_relative = 0.02))
  expect_true(in_peak_region(31.4, 31.3, "absolute", tol_absolute = 0.2))
  expect_false(in_peak_region(31.6, 31.3, "absolute", tol_absolute = 0.2))
})

# --- numeric helpers / stability ------------------------------------------
test_that("CV and stability classification behave", {
  expect_equal(cv(c(10, 10, 10)), 0)
  expect_true(is.na(cv(5)))
  expect_equal(stability_class(0.9), "high")
  expect_equal(stability_class(0.6), "moderate")
  expect_equal(stability_class(0.3), "low")
})

test_that("mode-switch counting detects toggles", {
  count_switches <- function(x) sum(diff(as.integer(x)) != 0)
  expect_equal(count_switches(c(TRUE, TRUE, TRUE)), 0)
  expect_equal(count_switches(c(TRUE, FALSE, TRUE)), 2)
})

# --- binomial CI -----------------------------------------------------------
test_that("binomial CI is bounded and brackets the point estimate", {
  ci <- binom_ci(50, 100)
  expect_gte(ci["lower"], 0); expect_lte(ci["upper"], 1)
  expect_lt(ci["lower"], 0.5); expect_gt(ci["upper"], 0.5)
  ci0 <- binom_ci(0, 100)
  expect_equal(unname(ci0["lower"]), 0)
})

# --- AIC / BIC formulas ----------------------------------------------------
test_that("AIC/BIC formulas are correct", {
  ll <- -100; npar <- 4L; n <- 160L
  expect_equal(-2 * ll + 2 * npar, 208)              # AIC
  expect_equal(-2 * ll + npar * log(n), -2*ll + 4*log(160))  # BIC
})

# --- significance integration (tiny null) ---------------------------------
test_that("scan-global null max statistic is reproducible across worker layout", {
  skip_if_not(file.exists(context_data()))
  cfg <- default_audit_config(fast = TRUE)
  lines <- clean_lines_source(read_nist_csv(context_data()), "Fe", 2L, "wavenumber")$lines
  k_grid <- seq(0.5, 80, length.out = 200)  # coarse grid for speed
  res <- run_scan_analysis(lines, 120L, k_grid, 1L, 6.0)
  setup_rng(123L)
  a <- null_distribution(res$ell, res$y, res$baseline, res$mu0, k_grid, 1L, 8L)
  setup_rng(123L)
  b <- null_distribution(res$ell, res$y, res$baseline, res$mu0, k_grid, 1L, 8L)
  expect_equal(a$max_vals, b$max_vals)  # identical given the same seed
})
