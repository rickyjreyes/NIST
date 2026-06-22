# Tests 1-2: Excel-style cell cleaning and numeric extraction from NIST cells.

test_that("clean_cell strips NIST/Excel quoting", {
  expect_equal(clean_cell('="123.45"'), "123.45")
  expect_equal(clean_cell('="49998.80"'), "49998.80")
  expect_equal(clean_cell("=123.45"), "123.45")
  expect_equal(clean_cell('  "999.0"  '), "999.0")
  expect_equal(clean_cell("plain"), "plain")
})

test_that("numeric_from_nist_cell extracts leading numeric token", {
  v <- numeric_from_nist_cell(c('="49998.80"', "2300d?", "250bl(Fe III)", "(0)", "", NA))
  expect_equal(v[1], 49998.80)
  expect_equal(v[2], 2300.0)
  expect_equal(v[3], 250.0)
  expect_equal(v[4], 0.0)
  expect_true(is.na(v[5]))
  expect_true(is.na(v[6]))
})

test_that("numeric_from_nist_cell handles scientific notation", {
  v <- numeric_from_nist_cell(c('="3.46e+08"', "2.0e+04"))
  expect_equal(v[1], 3.46e8)
  expect_equal(v[2], 2.0e4)
})
