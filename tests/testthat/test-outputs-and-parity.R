# Tests 15-16: output-file creation via the CLI, and numerical parity against
# the committed Python reference outputs.

test_that("CLI produces the full output contract on a synthetic CSV", {
  scanner <- file.path(REPO_ROOT_T, "R", "nist_wct_log_spectral_scan.R")
  skip_if_not(file.exists(scanner))

  # Build a small synthetic NIST-style CSV with an injected harmonic.
  set.seed(99)
  n_lines <- 600
  ell <- runif(n_lines, log(6000), log(50000))
  wn <- exp(ell)
  tmp_csv <- tempfile(fileext = ".csv")
  df <- data.frame(check.names = FALSE,
    element = rep("Xx", n_lines),
    sp_num = rep("2", n_lines),
    `wn(cm-1)` = sprintf('="%.4f"', wn),
    stringsAsFactors = FALSE)
  colnames(df)[3] <- "wn(cm-1)"
  utils::write.csv(df, tmp_csv, row.names = FALSE, quote = TRUE)

  out_dir <- tempfile("out_")
  status <- system2("Rscript",
    c(scanner, "--csv", tmp_csv, "--ion", "2", "--bins", "40",
      "--n-k", "100", "--null-n", "0", "--min-lines", "100",
      "--out-dir", out_dir),
    stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0L)

  for (f in c("nist_lines_clean.csv", "nist_binned_spectrum.csv",
              "nist_scan_curve.csv", "nist_summary.json",
              "nist_spectrum_fit.png", "nist_scan_curve.png")) {
    expect_true(file.exists(file.path(out_dir, f)), info = f)
  }
  # null file should NOT exist when null-n == 0
  expect_false(file.exists(file.path(out_dir, "nist_null.csv")))

  s <- jsonlite::fromJSON(file.path(out_dir, "nist_summary.json"))
  expect_equal(s$n_unique_lines, n_lines)
  expect_true(is.null(s$scan_null_p))   # null serialized to JSON null
})

test_that("R pipeline matches committed Python reference outputs (Fe II, 160)", {
  fe_csv <- file.path(REPO_ROOT_T, "data", "Fe_lines.csv")
  py_dir <- file.path(REPO_ROOT_T, "outputs", "fe_ion2_160")
  skip_if_not(file.exists(fe_csv), "data/Fe_lines.csv not present")
  skip_if_not(file.exists(file.path(py_dir, "nist_summary.json")),
              "committed Python outputs not present")

  raw <- read_nist_csv(fe_csv)
  lines <- clean_lines(raw, ion = 2)
  binned <- build_binned(lines, 160L, NULL, NULL)
  y <- binned$count
  ell <- binned$ell
  baseline <- pmax(gaussian_filter_nearest(y, 6.0), EPS)
  k_grid <- seq(0.5, 80.0, length.out = 2500L)
  sk <- scan_k(ell, y, baseline, k_grid, 1L)
  delta_ell <- max(ell) - min(ell)
  br <- branch_report(sk$best$k_best, delta_ell)

  ps <- jsonlite::fromJSON(file.path(py_dir, "nist_summary.json"))

  # Deterministic parity tolerances from the task specification.
  expect_lte(abs(sk$best$k_best - ps$best$k_best), 1e-10)
  expect_lte(abs(sk$best$deltaD - ps$best$deltaD) / abs(ps$best$deltaD), 1e-5)
  expect_lte(abs(br$n_obs - ps$branch_report$n_obs), 1e-8)

  # Baseline column parity against committed binned spectrum.
  pb <- utils::read.csv(file.path(py_dir, "nist_binned_spectrum.csv"),
                        check.names = FALSE)$baseline
  expect_lte(max(abs(baseline - pb) / pmax(abs(pb), 1e-12)), 1e-5)

  # Scan-curve deltaD parity against committed scan curve.
  pd <- utils::read.csv(file.path(py_dir, "nist_scan_curve.csv"),
                        check.names = FALSE)$deltaD
  expect_lte(max(abs(sk$scan$deltaD - pd) / pmax(abs(pd), 1e-12)), 1e-5)
})
