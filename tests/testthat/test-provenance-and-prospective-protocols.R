testthat::test_that("new release/prospective scripts parse", {
  files <- c(
    "verify_source_provenance.R",
    "freeze_release_snapshot.R",
    "run_resolution_invariant_protocol_v1.R",
    "run_independent_atomic_replication_v1.R",
    "finalize_adversarial_audit_v2.R"
  )
  for (f in files) {
    testthat::expect_error(parse(file = .find_R(f)), NA,
                           info = paste("R parse failed for", f))
  }
})

testthat::test_that("provenance config records unrecoverable fields explicitly", {
  root <- normalizePath(file.path(dirname(AUDIT_PATH), ".."))
  p <- file.path(root, "config", "nist_source_provenance.csv")
  x <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE,
                       fileEncoding = "UTF-8")
  if (length(names(x))) names(x)[1] <- sub("^\\ufeff", "", names(x)[1])
  testthat::expect_equal(sort(x$species), c("Co", "Fe"))
  testthat::expect_true(all(x$retrieval_date_status == "UNRECOVERABLE"))
  testthat::expect_true(all(x$query_export_settings_status == "UNRECOVERABLE"))
  testthat::expect_true(all(x$retrieval_version == "5.12"))
  testthat::expect_true(all(file.exists(file.path(root, x$provenance_evidence))))
})

testthat::test_that("resolution protocol is frozen before its outcome", {
  root <- normalizePath(file.path(dirname(AUDIT_PATH), ".."))
  p <- jsonlite::fromJSON(file.path(root, "config", "resolution_invariant_protocol_v1.json"))
  testthat::expect_identical(p$status, "FROZEN_BEFORE_FINAL_RUN")
  testthat::expect_equal(p$target_k, 31.3265306122449, tolerance = 1e-12)
  testthat::expect_equal(as.integer(p$bin_grid), seq(60L, 240L, by = 20L))
  testthat::expect_equal(p$smoothing$reference_bins, 160)
  testthat::expect_equal(p$smoothing$reference_sigma_bins, 6)
  testthat::expect_equal(p$null$replicates_per_bin, 5000)
  testthat::expect_equal(p$injection$replicates_per_bin_per_amplitude, 2000)
  testthat::expect_true(isTRUE(p$reporting$preserve_original_fixed_sigma_result))
})

testthat::test_that("independent catalog target is frozen and cannot be retuned", {
  root <- normalizePath(file.path(dirname(AUDIT_PATH), ".."))
  p <- jsonlite::fromJSON(file.path(root, "config", "independent_atomic_replication_v1.json"))
  testthat::expect_identical(p$status, "FROZEN_BEFORE_EXTERNAL_OUTCOME_ACCESS")
  testthat::expect_equal(p$selection$kurucz_element_code, 26.01)
  testthat::expect_equal(p$primary_test$target_k, 31.3265306122449, tolerance = 1e-12)
  testthat::expect_false(isTRUE(p$primary_test$retune_target_on_replication_data))
  testthat::expect_equal(p$primary_test$null_replicates, 5000)
  testthat::expect_true(isTRUE(p$exploratory_scan$cannot_change_primary_classification))
  testthat::expect_identical(p$source$raw_sha256, "")
})
