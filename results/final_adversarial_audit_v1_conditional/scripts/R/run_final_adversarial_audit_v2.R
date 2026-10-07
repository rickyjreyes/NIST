#!/usr/bin/env Rscript
# run_final_adversarial_audit_v2.R
# Canonical acceptance-contract runner for the final NIST audit.
#
# It deliberately reuses R/run_final_adversarial_audit.R for the established
# broad adversarial program, then adds the explicit Fe/Co contract checks and
# final evidence classifier/freeze gate. No detector retuning occurs here.

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)

parse_v2_args <- function(argv) {
  d <- list(`bootstrap-n` = 5000L, `null-n` = 5000L,
            `calibration-n` = 10000L, `injection-n` = 2000L,
            seed = 20260517L, parallel = FALSE, `render-report` = TRUE,
            force = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!key %in% names(d)) { i <- i + 1L; next }
    if (key %in% c("parallel", "render-report", "force")) {
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
  if (!identical(as.integer(d$seed), 20260517L))
    stop("strict final audit requires the frozen canonical seed 20260517")
  for (k in c("bootstrap-n", "null-n", "calibration-n", "injection-n")) {
    if (!is.finite(d[[k]]) || d[[k]] <= 0L) stop("invalid --", k)
  }
  d
}

run_main_from <- function(file, argv) {
  e <- new.env(parent = globalenv())
  sys.source(file.path(R_DIR, file), envir = e)
  if (!exists("main", envir = e, inherits = FALSE))
    stop("module has no main(): ", file)
  e$main(argv)
}

main <- function(argv = commandArgs(TRUE)) {
  a <- parse_v2_args(argv)
  pflag <- c("--parallel", if (isTRUE(a$parallel)) "true" else "false")
  rflag <- c("--render-report", if (isTRUE(a[["render-report"]])) "true" else "false")

  cat("\n=== PHASE 1: established full adversarial program ===\n")
  run_main_from("run_final_adversarial_audit.R", c(
    "--bootstrap-n", as.character(a[["bootstrap-n"]]),
    "--null-n", as.character(a[["null-n"]]),
    "--calibration-n", as.character(a[["calibration-n"]]),
    "--injection-n", as.character(a[["injection-n"]]),
    pflag, rflag
  ))

  cat("\n=== PHASE 2: Fe/Co acceptance-contract checks ===\n")
  run_main_from("run_final_contract_checks.R", c(
    "--bootstrap-n", as.character(a[["bootstrap-n"]]),
    "--null-n", as.character(a[["null-n"]]),
    "--calibration-n", as.character(a[["calibration-n"]]),
    "--injection-n", as.character(a[["injection-n"]]),
    "--seed", as.character(a$seed),
    pflag
  ))

  cat("\n=== PHASE 3: classify and freeze only if PASS ===\n")
  run_main_from("finalize_adversarial_audit_v2.R", character(0))
  invisible(TRUE)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_final_adversarial_audit_v2\\.R$", .invoked_file)) main()
