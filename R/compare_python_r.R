#!/usr/bin/env Rscript
# compare_python_r.R
# ---------------------------------------------------------------------------
# Compare a canonical Python output directory against an R output directory and
# report numerical parity on the deterministic quantities.
#
# Usage:
#   Rscript R/compare_python_r.R --py outputs/fe_ion2_160 --r outputs_r/fe_ion2_160
#
# Tolerances (deterministic quantities):
#   absolute k_best difference      <= 1e-10
#   relative deltaD difference      <= 1e-5
#   relative baseline differences   <= 1e-5
#   relative scan-deltaD / amplitude<= 1e-5
#   absolute n_obs difference       <= 1e-8
#
# The parametric Poisson bootstrap (nist_null.csv, scan_null_p) is NOT compared
# numerically: R's rpois and NumPy's Generator.poisson differ by construction.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages(library(jsonlite)))

parse_pair <- function(argv) {
  out <- list(py = NULL, r = NULL)
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (key %in% c("py", "r")) {
      out[[key]] <- argv[i + 1L]
      i <- i + 2L
    } else {
      stop(sprintf("Unknown option: %s", argv[i]))
    }
  }
  if (is.null(out$py) || is.null(out$r)) stop("Both --py and --r are required")
  out
}

rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-12)

read_col <- function(path, col) {
  df <- utils::read.csv(path, check.names = FALSE)
  as.numeric(df[[col]])
}

compare_dirs <- function(py_dir, r_dir) {
  checks <- list()
  add <- function(name, value, tol, kind) {
    pass <- is.finite(value) && value <= tol
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, kind = kind, value = value, tol = tol,
      pass = pass, stringsAsFactors = FALSE
    )
  }

  ps <- fromJSON(file.path(py_dir, "nist_summary.json"))
  rs <- fromJSON(file.path(r_dir, "nist_summary.json"))

  add("k_best (abs)", abs(ps$best$k_best - rs$best$k_best), 1e-10, "abs")
  add("deltaD (rel)", rel(rs$best$deltaD, ps$best$deltaD), 1e-5, "rel")
  add("D_base (rel)", rel(rs$best$D_base, ps$best$D_base), 1e-5, "rel")
  add("amplitude (rel)", rel(rs$best$amplitude, ps$best$amplitude), 1e-5, "rel")
  add("n_obs (abs)", abs(ps$branch_report$n_obs - rs$branch_report$n_obs), 1e-8, "abs")
  add("ell_min (abs)", abs(ps$ell_min - rs$ell_min), 1e-8, "abs")
  add("delta_ell (abs)", abs(ps$delta_ell - rs$delta_ell), 1e-8, "abs")

  pb <- read_col(file.path(py_dir, "nist_binned_spectrum.csv"), "baseline")
  rb <- read_col(file.path(r_dir, "nist_binned_spectrum.csv"), "baseline")
  add("baseline (max rel)", max(rel(rb, pb)), 1e-5, "rel")

  pc <- read_col(file.path(py_dir, "nist_binned_spectrum.csv"), "count")
  rc <- read_col(file.path(r_dir, "nist_binned_spectrum.csv"), "count")
  add("count (max abs)", max(abs(rc - pc)), 0, "abs")

  pd <- read_col(file.path(py_dir, "nist_scan_curve.csv"), "deltaD")
  rd <- read_col(file.path(r_dir, "nist_scan_curve.csv"), "deltaD")
  add("scan deltaD (max rel)", max(rel(rd, pd)), 1e-5, "rel")

  pa <- read_col(file.path(py_dir, "nist_scan_curve.csv"), "amplitude")
  ra <- read_col(file.path(r_dir, "nist_scan_curve.csv"), "amplitude")
  add("scan amplitude (max rel)", max(rel(ra, pa)), 1e-5, "rel")

  do.call(rbind, checks)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- parse_pair(argv)
  tab <- compare_dirs(args$py, args$r)
  cat(sprintf("Comparing PY=%s  R=%s\n\n", args$py, args$r))
  fmt <- "%-26s %-5s %14.3e %10.1e  %s\n"
  cat(sprintf("%-26s %-5s %14s %10s  %s\n", "check", "kind", "value", "tol", "result"))
  for (i in seq_len(nrow(tab))) {
    cat(sprintf(fmt, tab$check[i], tab$kind[i], tab$value[i], tab$tol[i],
                if (tab$pass[i]) "PASS" else "FAIL"))
  }
  cat("\n")
  if (all(tab$pass)) {
    cat("PARITY: PASS (all deterministic tolerances satisfied)\n")
    quit(status = 0)
  } else {
    cat("PARITY: FAIL\n")
    quit(status = 1)
  }
}

.invoked_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(.invoked_file) > 0L && grepl("compare_python_r\\.R$", .invoked_file)) {
  main()
}
