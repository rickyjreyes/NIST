#!/usr/bin/env Rscript
# render_statistical_audit.R
# ---------------------------------------------------------------------------
# Public audit entry point.
#
# Development/non-strict invocations retain the historical orchestrator exactly
# in R/render_statistical_audit_legacy.R.
#
# --strict is different by design: it is routed to the uncapped, fail-fast
# final acceptance-contract runner. This prevents the former development caps
# (500/1000 simulations in selected modules) from silently weakening a command
# that requests 5000/10000/2000 final budgets.
# ---------------------------------------------------------------------------

.args_all <- commandArgs(FALSE)
.file_args <- sub("^--file=", "", .args_all[grep("^--file=", .args_all)])
.top_file <- if (length(.file_args)) basename(.file_args[1]) else ""
.direct <- identical(.top_file, "render_statistical_audit.R")
.argv <- commandArgs(TRUE)

if (.direct && any(.argv == "--strict")) {
  .self <- if (length(.file_args)) normalizePath(.file_args[1], mustWork = FALSE) else
    normalizePath(file.path("R", "render_statistical_audit.R"), mustWork = FALSE)
  .rdir <- dirname(.self)
  .runner <- file.path(.rdir, "run_final_adversarial_audit_v2.R")
  if (!file.exists(.runner)) stop("strict final runner missing: ", .runner)

  .rscript <- file.path(R.home("bin"), "Rscript")
  .status <- system2(.rscript, c(.runner, .argv))
  if (!identical(as.integer(.status), 0L))
    stop("strict final adversarial audit failed with exit status ", .status)
} else {
  # Sourcing the legacy file also exports build_parity(), build_dashboard(),
  # build_claim_matrix(), render_report(), and main() for the established final
  # runner. Its invocation guard only executes when this public renderer is the
  # top-level script, preserving historical non-strict CLI behaviour.
  .self <- if (length(.file_args)) normalizePath(.file_args[1], mustWork = FALSE) else
    normalizePath(file.path("R", "render_statistical_audit.R"), mustWork = FALSE)
  .legacy <- file.path(dirname(.self), "render_statistical_audit_legacy.R")
  if (!file.exists(.legacy)) stop("legacy audit orchestrator missing: ", .legacy)
  sys.source(.legacy, envir = environment())
}
