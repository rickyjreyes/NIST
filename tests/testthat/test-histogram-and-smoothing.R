# Tests 7-8: histogram construction and Gaussian smoothing edge behavior.

test_that("np_histogram_uniform matches NumPy uniform-bin semantics", {
  # 10 values across [0, 10] into 5 bins -> width 2 each.
  a <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10)
  h <- np_histogram_uniform(a, 5L, 0, 10)
  # bins: [0,2),[2,4),[4,6),[6,8),[8,10]; last bin closed includes 8 and 10
  expect_equal(h$counts, c(2, 2, 2, 2, 2))
  expect_equal(h$edges, seq(0, 10, length.out = 6))
})

test_that("histogram right edge is inclusive on the last bin", {
  a <- c(0.0, 1.0)   # 1.0 is the right edge
  h <- np_histogram_uniform(a, 1L, 0, 1)
  expect_equal(h$counts, 2L)
})

test_that("build_binned returns centers and edges", {
  lines <- data.frame(ell = c(0, 1, 2, 3, 4))
  b <- build_binned(lines, 4L, 0, 4)
  expect_equal(nrow(b), 4)
  expect_equal(b$edge_lo, c(0, 1, 2, 3))
  expect_equal(b$edge_hi, c(1, 2, 3, 4))
  expect_equal(b$ell, c(0.5, 1.5, 2.5, 3.5))
})

test_that("gaussian_filter_nearest uses nearest-edge extension (constant in -> constant out)", {
  # A constant signal must be returned unchanged under nearest-edge smoothing.
  y <- rep(7.0, 20)
  out <- gaussian_filter_nearest(y, sigma = 3.0)
  expect_equal(out, y, tolerance = 1e-12)
})

test_that("gaussian kernel has SciPy radius floor(truncate*sigma+0.5)", {
  w <- gaussian_kernel1d(6.0, truncate = 4.0)
  expect_equal(length(w), 2L * 24L + 1L)  # radius 24
  expect_equal(sum(w), 1.0, tolerance = 1e-15)
})

test_that("gaussian_filter_nearest preserves total mass less than reflect at a step edge", {
  # Compare to a hand-rolled nearest-mode reference for a small case.
  y <- c(0, 0, 0, 10, 0, 0, 0)
  sigma <- 1.0
  w <- gaussian_kernel1d(sigma)
  lw <- (length(w) - 1L) %/% 2L
  ref <- numeric(length(y))
  for (i in seq_along(y)) {
    acc <- 0
    for (j in seq.int(-lw, lw)) {
      idx <- min(max(i + j, 1L), length(y))
      acc <- acc + w[j + lw + 1L] * y[idx]
    }
    ref[i] <- acc
  }
  expect_equal(gaussian_filter_nearest(y, sigma), ref, tolerance = 1e-14)
})
