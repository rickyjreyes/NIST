#!/usr/bin/env Rscript
# run_final_adversarial_audit.R
# ---------------------------------------------------------------------------
# Strict one-command final adversarial statistical audit.
#
# This entry point has NO silent Monte Carlo caps and NO partial-success mode.
# Any module failure aborts immediately so stale tables from a prior run cannot
# be summarized as current final evidence.
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

# Locate R/ robustly both when invoked directly and when sys.source()d from the
# test suite. commandArgs() describes the top-level process, not necessarily the
# currently sourced file, so it cannot be the sole locator.
locate_final_r_dir <- function() {
  all_args <- commandArgs(FALSE)
  file_args <- sub("^--file=", "", all_args[grep("^--file=", all_args)])
  invoked_dir <- if (length(file_args)) dirname(file_args[1]) else NA_character_
  helper_dir <- if (exists("AUDIT_PATH", inherits = TRUE)) dirname(get("AUDIT_PATH", inherits = TRUE)) else NA_character_
  candidates <- unique(c(invoked_dir, helper_dir, "R", file.path(getwd(), "R"), "."))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (d in candidates) {
    if (file.exists(file.path(d, "audit_utils.R")) &&
        file.exists(file.path(d, "render_statistical_audit.R"))) {
      return(normalizePath(d, mustWork = TRUE))
    }
  }
  stop("could not locate repository R/ directory for final adversarial audit")
}

R_DIR <- locate_final_r_dir()
if (!exists("audit_repo_root", mode = "function")) source(file.path(R_DIR, "audit_utils.R"))

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
    `render-report` = TRUE
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!key %in% names(d)) { i <- i + 1L; next }
    if (key %in% c("parallel", "render-report")) {
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
  nums <- setdiff(names(d), c("parallel", "render-report"))
  bad <- nums[!vapply(d[nums], function(x) is.finite(x) && x > 0, logical(1))]
  if (length(bad)) stop("all final audit budgets must be positive: ", paste(bad, collapse = ", "))
  d
}

run_checked <- function(name, argv = character(0)) {
  cat(sprintf("\n=== FINAL ADVERSARIAL [%s] ===\n", name))
  t0 <- Sys.time()
  e <- new.env(parent = globalenv())
  ok <- FALSE; err <- ""
  tryCatch({
    sys.source(file.path(R_DIR, name), envir = e)
    if (!exists("main", envir = e, inherits = FALSE)) stop("module has no main(): ", name)
    e$main(argv)
    ok <- TRUE
  }, error = function(x) { err <<- conditionMessage(x) })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  row <- data.frame(module = name, ok = ok, seconds = secs, error = err,
                    stringsAsFactors = FALSE)
  if (!ok) stop(sprintf("final adversarial audit aborted at %s: %s", name, err))
  row
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
  td <- file.path(root, "tables_r/statistical_audit")
  fd <- file.path(root, "figures_r/statistical_audit")
  setup_rng(20260517L)
  pflag <- if (isTRUE(args$parallel)) "--parallel" else character(0)
  status <- list()
  run <- function(name, a = character(0)) {
    status[[length(status) + 1L]] <<- run_checked(name, a)
  }

  run("build_analysis_registry.R")
  run("build_dataset_flow.R")

  run("bootstrap_peak_uncertainty.R",
      c("--bootstrap-n", as.character(args[["bootstrap-n"]]), pflag))
  run("peak_stability.R")

  run("build_effect_size_table.R",
      c("--null-n", as.character(args[["null-n"]]), pflag))
  run("run_bin_grid.R",
      c("--null-n", as.character(args[["null-n"]]), pflag))
  run("run_model_sensitivity.R",
      c("--null-n", as.character(args[["sensitivity-null-n"]]), pflag))

  run("model_comparison.R")
  run("run_observed_ritz_replication.R",
      c("--null-n", as.character(args[["null-n"]]),
        "--bootstrap-n", as.character(args[["source-bootstrap-n"]]), pflag))
  run("run_holdout_replication.R",
      c("--null-n", as.character(args[["holdout-null-n"]]), pflag))

  # Full family now includes the complete declared sigma/degree/bin multiverse.
  run("global_multiple_testing.R",
      c("--family-n", as.character(args[["family-n"]]), pflag))

  run("calibrate_false_positive_rate.R",
      c("--calibration-n", as.character(args[["calibration-n"]]),
        "--null-n", as.character(args[["null-n"]]), pflag))
  run("run_injection_recovery.R",
      c("--injection-n", as.character(args[["injection-n"]]),
        "--null-n", as.character(args[["injection-null-n"]]), pflag))

  run("run_alternative_nulls.R",
      c("--null-n", as.character(args[["alternative-null-n"]]),
        "--block-len", as.character(args[["block-len"]]),
        "--seed", "20260517", pflag))

  # Reuse deterministic parity/dashboard functions, but do not invoke the old
  # general orchestrator because it contains development-oriented caps.
  renderer <- new.env(parent = globalenv())
  sys.source(file.path(R_DIR, "render_statistical_audit.R"), envir = renderer)
  renderer$build_parity(file.path(root, "outputs"), file.path(root, "outputs_r"), td, fd)
  tryCatch(renderer$build_dashboard(fd, td), error = function(e) message("[dashboard] ", conditionMessage(e)))
  # Preserve the historical claim matrix too; the new one is stricter.
  tryCatch(renderer$build_claim_matrix(td), error = function(e) message("[legacy claim matrix] ", conditionMessage(e)))

  run("build_final_adversarial_summary.R")

  required <- c(
    "analysis_registry.csv", "dataset_flow_by_species.csv", "peak_confidence_intervals.csv",
    "significance_results.csv", "multiple_testing.csv", "null_calibration.csv",
    "holdout_results.csv", "model_comparison.csv", "injection_recovery.csv",
    "specification_results.csv", "alternative_null_results.csv",
    "calibrated_global_evidence.csv", "final_adversarial_claim_matrix.csv", "failed_claims.csv"
  )
  require_outputs(root, required)

  # Write CURRENT provenance before rendering so the report cannot embed an old
  # manifest from a previous run.
  stat <- do.call(rbind, status)
  expected_report <- file.path(root, "reports/rendered/nist_final_adversarial_audit.html")
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
    all_modules_ok = all(stat$ok),
    report_requested = args[["render-report"]],
    report_expected_path = expected_report,
    stringsAsFactors = FALSE
  )
  write_table(manifest, file.path(td, "final_adversarial_run_manifest.csv"))
  write_table(stat, file.path(td, "final_adversarial_module_status.csv"))

  report_path <- NA_character_
  if (isTRUE(args[["render-report"]])) {
    report_path <- renderer$render_report(file.path(root, "reports/nist_final_adversarial_audit.qmd"),
                                           root, FALSE)
    if (is.na(report_path)) message("[report] final tables are valid, but HTML was not rendered")
  }

  ev <- utils::read.csv(file.path(td, "calibrated_global_evidence.csv"), stringsAsFactors = FALSE)
  fc <- utils::read.csv(file.path(td, "failed_claims.csv"), stringsAsFactors = FALSE)
  cat("\n=====================================================\n")
  cat("NIST FINAL ADVERSARIAL STATISTICAL AUDIT\n")
  cat("=====================================================\n")
  cat(sprintf("largest reported global p : %.6g\n", ev$largest_reported_global_p[1]))
  cat("  (diagnostic maximum across reported global/multiplicity/alternative-null p-values; not a combined p)\n")
  cat(sprintf("non-passing claims         : %d\n", nrow(fc)))
  cat(sprintf("modules                    : %d (%d ok)\n", nrow(stat), sum(stat$ok)))
  cat(sprintf("report                     : %s\n",
              if (is.na(report_path)) "not rendered; use CSV outputs" else report_path))
  cat("\nNo failed or mixed claim is suppressed from failed_claims.csv.\n")
  invisible(list(status = stat, manifest = manifest, evidence = ev, failed_claims = fc))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_final_adversarial_audit\\.R$", .invoked_file)) main()
