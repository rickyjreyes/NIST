context("final adversarial acceptance contract")

ROOT <- normalizePath(file.path("..", ".."), mustWork = TRUE)
repo_path <- function(...) file.path(ROOT, ...)

test_that("frozen final protocol budgets and seed are exact", {
  cfg <- jsonlite::fromJSON(repo_path("config", "final_adversarial_audit.json"))
  expect_equal(cfg$seed, 20260517)
  expect_equal(cfg$budgets$bootstrap_n, 5000)
  expect_equal(cfg$budgets$null_n, 5000)
  expect_equal(cfg$budgets$calibration_n, 10000)
  expect_equal(cfg$budgets$injection_n, 2000)
  expect_equal(as.numeric(cfg$decision_rules$calibration_levels),
               c(0.10, 0.05, 0.01, 0.001))
})

test_that("strict public entry point routes to uncapped v2 runner", {
  txt <- paste(readLines(repo_path("R", "render_statistical_audit.R"), warn = FALSE), collapse = "\n")
  expect_match(txt, "run_final_adversarial_audit_v2\\.R")
  expect_match(txt, "--strict")
  expect_false(grepl("min\\(cfg\\$null_n, 500L\\)", txt))
  expect_false(grepl("min\\(cfg\\$bootstrap_n, 1000L\\)", txt))
})

test_that("contract runner keeps Fe and Co separate", {
  txt <- paste(readLines(repo_path("R", "run_final_contract_checks.R"), warn = FALSE), collapse = "\n")
  expect_match(txt, 'c\\("Fe", "Co"\\)')
  expect_match(txt, "fixed-frequency p is descriptive")
  expect_match(txt, "final_species_null_summary\\.csv")
  expect_match(txt, "final_species_calibration\\.csv")
  expect_match(txt, "final_species_injection_recovery\\.csv")
})

test_that("zero-exceedance Monte Carlo is never p zero", {
  source(repo_path("R", "audit_utils.R"), local = TRUE)
  expect_equal(emp_p(0, 5000), 1 / 5001)
  expect_gt(emp_p(0, 5000), 0)
})

test_that("freeze gate has exactly the declared classification labels", {
  txt <- paste(readLines(repo_path("R", "finalize_adversarial_audit.R"), warn = FALSE), collapse = "\n")
  expect_match(txt, "PASS — ADVERSARIAL AUDIT SURVIVED")
  expect_match(txt, "CONDITIONAL")
  expect_match(txt, "FAIL")
  expect_match(txt, "results/final_adversarial_audit_v1")
  expect_match(txt, "if \\(!identical\\(classification, PASS_LABEL\\)\\)")
})

test_that("corrected finalizer avoids classifier shadowing and platform path regex", {
  txt <- paste(readLines(repo_path("R", "finalize_adversarial_audit_v2.R"), warn = FALSE), collapse = "\n")
  expect_match(txt, "ctab <- cal")
  expect_match(txt, "base::c\\(0.10, 0.05, 0.01, 0.001\\)")
  expect_match(txt, "basename\\(json_path\\)")
})

test_that("source provenance cannot silently self-complete", {
  p <- read.csv(repo_path("config", "nist_source_provenance.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("Fe", "Co") %in% p$species))
  expect_true(all(is.na(p$retrieval_date) | !nzchar(trimws(p$retrieval_date))))
  expect_true(all(is.na(p$retrieval_version) | !nzchar(trimws(p$retrieval_version))))
})
