# Tests 3-6: column detection, wavelength conversion, ion filtering, dedup.

test_that("find_column detects wavenumber column names", {
  df <- data.frame(check.names = FALSE,
    element = c("Fe", "Fe"), sp_num = c("2", "2"),
    `wn(cm-1)` = c('="40000.00"', '="20000.00"'), stringsAsFactors = FALSE)
  colnames(df)[3] <- "wn(cm-1)"
  ex <- extract_wavenumber_cm(df)
  expect_match(tolower(ex$source), "wn")
  expect_equal(ex$wn[1], 40000.0)
  expect_equal(ex$wn[2], 20000.0)
})

test_that("wavelength-to-wavenumber conversion uses 1e7/nm", {
  df <- data.frame(check.names = FALSE,
    `obs_wl_air(nm)` = c("500.0", "1000.0"), stringsAsFactors = FALSE)
  colnames(df)[1] <- "obs_wl_air(nm)"
  ex <- extract_wavenumber_cm(df)
  expect_match(ex$source, "converted_from")
  expect_equal(ex$wn[1], 1e7 / 500.0)   # 20000
  expect_equal(ex$wn[2], 1e7 / 1000.0)  # 10000
})

test_that("clean_lines filters by ion (sp_num)", {
  df <- data.frame(check.names = FALSE,
    element = rep("Fe", 4), sp_num = c("1", "2", "2", "3"),
    `wn(cm-1)` = c('="10000"', '="20000"', '="30000"', '="40000"'),
    stringsAsFactors = FALSE)
  colnames(df)[3] <- "wn(cm-1)"
  cl <- clean_lines(df, ion = 2)
  expect_equal(nrow(cl), 2)
  expect_setequal(cl$wavenumber_cm, c(20000.0, 30000.0))
})

test_that("clean_lines removes duplicate wavenumbers and sorts", {
  df <- data.frame(check.names = FALSE,
    element = rep("Fe", 4), sp_num = rep("2", 4),
    `wn(cm-1)` = c('="30000"', '="10000"', '="10000"', '="20000"'),
    stringsAsFactors = FALSE)
  colnames(df)[3] <- "wn(cm-1)"
  cl <- clean_lines(df, ion = 2)
  expect_equal(cl$wavenumber_cm, c(10000.0, 20000.0, 30000.0))  # sorted, deduped
  expect_equal(cl$ell, log(c(10000.0, 20000.0, 30000.0)))
})

test_that("clean_lines drops non-positive / non-finite wavenumbers", {
  df <- data.frame(check.names = FALSE,
    element = rep("Fe", 3), sp_num = rep("2", 3),
    `wn(cm-1)` = c('="0"', '="-5"', '="15000"'),
    stringsAsFactors = FALSE)
  colnames(df)[3] <- "wn(cm-1)"
  cl <- clean_lines(df, ion = 2)
  expect_equal(cl$wavenumber_cm, 15000.0)
})
