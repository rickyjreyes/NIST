# Test 14 (+ peak/ratio helpers): branch-report and rational approximation.

test_that("branch_report computes n_obs and target windings", {
  k_best <- 31.3265306122449
  delta_ell <- 2.1540594204535672
  br <- branch_report(k_best, delta_ell)
  expect_equal(br$n_obs, k_best * delta_ell / (2 * pi), tolerance = 1e-12)
  # target k for a given n is 2*pi*n/delta_ell
  row11 <- br$nearest[br$nearest$n == 11 & br$nearest$branch == "integer_1_to_40", ]
  expect_equal(row11$k_target[1], 2 * pi * 11 / delta_ell, tolerance = 1e-12)
  # nearest list is sorted ascending by abs_k_error and capped at 10
  expect_lte(nrow(br$nearest), 10)
  expect_true(!is.unsorted(br$nearest$abs_k_error))
})

test_that("limit_denominator reproduces simple rational approximations", {
  # 3.3421161825726142 -> 127/38 (from the canonical Fe II 160 run)
  fr <- limit_denominator(double_to_bigq(3.3421161825726142), 64L)
  expect_equal(as.character(gmp::numerator(fr)), "127")
  expect_equal(as.character(gmp::denominator(fr)), "38")
  # exact small fraction stays itself
  fr2 <- limit_denominator(double_to_bigq(0.5), 64L)
  expect_equal(as.character(gmp::numerator(fr2)), "1")
  expect_equal(as.character(gmp::denominator(fr2)), "2")
})

test_that("find_peaks_distance finds separated local maxima", {
  x <- c(0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0)
  peaks <- find_peaks_distance(x, distance = 1)
  # 0-based indices of the three maxima: 1, 6, 11
  expect_equal(sort(peaks), c(1, 6, 11))
})
