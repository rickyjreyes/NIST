#!/usr/bin/env Rscript
# run_final_adversarial_audit.R
# ---------------------------------------------------------------------------
# Strict one-command final adversarial statistical audit.
#
# Unlike render_statistical_audit.R's general-purpose development/full runner,
# this entry point has NO silent per-module Monte Carlo caps.  Every declared
# budget is explicit below and can be overridden on the command line.  Any
# module failure aborts the run by default, so a partially completed audit
# cannot be mistaken for a final result.
#
# Default final budgets:
#   bootstrap_n            5000
#   scan/null_n            5000
#   source_bootstrap_n     2000
#   sensitivity_null_n     1000
#   holdout_null_n         5000
#   family_n               5000
#   calibration_n         10000
#   injection_n            2000
#   injection_null_n       5000
#   alternative_null_n     5000 per declared null model
#
# Example:
#   Rscript R/run_final_adversarial_audit.R --parallel true --render-report true
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
source(file.path(R_DIR, "audit_utils.R"))

parse_args <- function(argv) {
  d <- list(
    `bootstrap-n` = 5000L,
    `null-n` = 5000L,
    `source-bootstrap-n` = 2000L,
    `sensitivity-null-n` = 1000L,
    `holdout-null-n` = 5000L,
    `family-n` = 5000L,
    `calibration-n` = 10000L,
    `injection-n` = 2000L,
    `injection-null-n` = 5000L,
    `alternative-null-n` = 5000L,
    `block-len` = 8L,
    parallel = FALSE,
    `render-report` = TRUE,
    `allow-partial` = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!key %in% names(d)) { i <- i + 1L; next }
    if (key %in% c("parallel", "render-report", "allow-partial")) {
      nxt <- if (i + 1L <= length(argv)) tolower(argv[i + 1L]) else NA_character_
      if (!is.na(nxt) && nxt %in% c("true", "false")) {
        d[[key]] <- nxt == "true"; i <- i + 2L
      } else {
        d[[key]] <- TRUE; i <- i + 1L
      }
    } else {
      if (i + 1L > length(argv)) stop("missing value for --", key)
      d[[key]] <- as.integer(argv[i + 1L]); i <- i + 2L
    }
  }
  d
}

run_checked <- function(name, argv, allow_partial = FALSE) {
  cat(sprintf("\n=== FINAL ADVERSARIAL [%s] ===\n", name))
  t0 <- Sys.time()
  err <- NULL
  ok <- tryCatch({
    e <- new.env(parent = globalenv())
    sys.source(file.path(R_DIR, name), envir = e)
    if (!exists("main", envir = e, inherits = FALSE)) stop("module has no main(): ", name)
    e$main(argv)
    TRUE
  }, error = function(e) { err <<- conditionMessage(e); FALSE })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!ok) {
    message(sprintf("[FAIL] %s: %s", name, err))
    if (!allow_partial) stop(sprintf("final adversarial audit aborted at %s: %s", name, err))
  }
  data.frame(module = name, ok = ok, seconds = secs,
             error = if (is.null(err)) "" else err, stringsAsFactors = FALSE)
}

require_outputs <- function(root, names) {
  td <- file.path(root, "tables_r/statistical_audit")
  miss <- names[!file.exists(file.path(td, names))]
  if (length(miss)) stop("final audit required outputs missing: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

main <- function(argv = commandArgs(TRUE)) {
  args <- parse_args(argv)
  root <- audit_repo_root()
  setup_rng(20260517L)
  pflag <- if (isTRUE(args$parallel)) "--parallel" else character(0)
  status <- list()
  run <- function(name, a = character(0)) {
    status[[length(status) + 1L]] <<- run_checked(name, a, args[["allow-partial"]])
  }

  # Core provenance/data accounting.
  run("build_analysis_registry.R")
  run("build_dataset_flow.R")

  # Bootstrap uncertainty and peak selection stability.
  run("bootstrap_peak_uncertainty.R",
      c("--bootstrap-n", as.character(args[["bootstrap-n"]]), pflag))
  run("peak_stability.R")

  # Main scan-global null evidence and declared preprocessing/specification grid.
  run("build_effect_size_table.R",
      c("--null-n", as.character(args[["null-n"]]), pflag))
  run("run_bin_grid.R",
      c("--null-n", as.character(args[["null-n"]]), pflag))
  run("run_model_sensitivity.R",
      c("--null-n", as.character(args[["sensitivity-null-n"]]), pflag))

  # In-sample vs held-out model behavior and representation replication.
  run("model_comparison.R")
  run("run_observed_ritz_replication.R",
      c("--null-n", as.character(args[["null-n"]]),
        "--bootstrap-n", as.character(args[["source-bootstrap-n"]]), pflag))
  run("run_holdout_replication.R",
      c("--null-n", as.character(args[["holdout-null-n"]]), pflag))

  # Full declared-family look-elsewhere/multiplicity correction: no hidden cap.
  run("global_multiple_testing.R",
      c("--family-n", as.character(args[["family-n"]]), pflag))

  # Type-I calibration and power.
  run("calibrate_false_positive_rate.R",
      c("--calibration-n", as.character(args[["calibration-n"]]),
        "--null-n", as.character(args[["null-n"]]), pflag))
  run("run_injection_recovery.R",
      c("--injection-n", as.character(args[["injection-n"]]),
        "--null-n", as.character(args[["injection-null-n"]]), pflag))

  # Adversarial alternative nulls, including baseline re-estimation.
  run("run_alternative_nulls.R",
      c("--null-n", as.character(args[["alternative-null-n"]]),
        "--block-len", as.character(args[["block-len"]]),
        "--seed", "20260517", pflag))

  # Reuse existing deterministic parity/dashboard artifacts without invoking
  # the general orchestrator (which has development-oriented caps).
  renderer <- new.env(parent = globalenv())
  sys.source(file.path(R_DIR, "render_statistical_audit.R"), envir = renderer)
  td <- file.path(root, "tables_r/statistical_audit")
  fd <- file.path(root, "figures_r/statistical_audit")
  tryCatch(renderer$build_parity(file.path(root, "outputs"), file.path(root, "outputs_r"), td, fd),
           error = function(e) {
             if (!args[["allow-partial"]]) stop("Python/R parity failed: ", conditionMessage(e))
             message("[parity] ", conditionMessage(e))
           })
  tryCatch(renderer$build_dashboard(fd, td), error = function(e) message("[dashboard] ", conditionMessage(e)))
  tryCatch(renderer$build_claim_matrix(td), error = function(e) message("[legacy claim matrix] ", conditionMessage(e)))

  # New cross-cutting calibrated evidence + explicit non-passing claim record.
  run("build_final_adversarial_summary.R")

  required <- c(
    "significance_results.csv", "multiple_testing.csv", "null_calibration.csv",
    "holdout_results.csv", "injection_recovery.csv", "specification_results.csv",
    "alternative_null_results.csv", "calibrated_global_evidence.csv",
    "final_adversarial_claim_matrix.csv", "failed_claims.csv"
  )
  if (!args[["allow-partial"]]) require_outputs(root, required)

  report_path <- NA_character_
  if (isTRUE(args[["render-report"]])) {
    report_path <- renderer$render_report(file.path(root, "reports/nist_final_adversarial_audit.qmd"),
                                           root, FALSE)
  }

  stat <- do.call(rbind, status)
  manifest <- data.frame(
    timestamp = audit_timestamp(),
    git_commit = git_commit(),
    canonical_seed = 20260517L,
    bootstrap_n = args[["bootstrap-n"]],
    null_n = args[["null-n"]],
    source_bootstrap_n = args[["source-bootstrap-n"]],
    sensitivity_null_n = args[["sensitivity-null-n"]],
    holdout_null_n = args[["holdout-null-n"]],
    family_n = args[["family-n"]],
    calibration_n = args[["calibration-n"]],
    injection_n = args[["injection-n"]],
    injection_null_n = args[["injection-null-n"]],
    alternative_null_n = args[["alternative-null-n"]],
    alternative_block_len = args[["block-len"]],
    parallel = args$parallel,
    allow_partial = args[["allow-partial"]],
    all_modules_ok = all(stat$ok),
    report = if (is.na(report_path)) "not rendered" else report_path,
    stringsAsFactors = FALSE
  )
  write_table(manifest, file.path(td, "final_adversarial_run_manifest.csv"))
  write_table(stat, file.path(td, "final_adversarial_module_status.csv"))

  ev <- utils::read.csv(file.path(td, "calibrated_global_evidence.csv"), stringsAsFactors = FALSE)
  fc <- utils::read.csv(file.path(td, "failed_claims.csv"), stringsAsFactors = FALSE)
  cat("\n=====================================================\n")
  cat("NIST FINAL ADVERSARIAL STATISTICAL AUDIT\n")
  cat("=====================================================\n")
  cat(sprintf("largest reported global p : %.6g\n", ev$largest_reported_global_p[1]))
  cat("  (diagnostic maximum across scan-global, family-max, and alternative-null p-values; not a combined p)\n")
  cat(sprintf("non-passing claims         : %d\n", nrow(fc)))
  cat(sprintf("modules                    : %d (%d ok, %d failed)\n",
              nrow(stat), sum(stat$ok), sum(!stat$ok)))
  cat(sprintf("report                     : %s\n", manifest$report[1]))
  cat("\nNo failed or mixed claim is suppressed from failed_claims.csv.\n")
  invisible(list(status = stat, manifest = manifest, evidence = ev, failed_claims = fc))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_final_adversarial_audit\\.R$", .invoked_file)) {
  main()
}
