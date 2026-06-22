#!/usr/bin/env Rscript
# make_verdict.R
# ---------------------------------------------------------------------------
# R reproduction of scripts/make_verdict.py: build the PASS/FAIL verdict for
# the Fe II log-cosine bin-stability claim from an R master-results table.
#
# Reads:  tables_r/nist_master_results_r.csv
# Writes: tables_r/nist_verdict_r.csv
#
# The gating logic and thresholds are identical to make_verdict.py:
# PASS only if ALL hold:
#   - rows exist for Fe ion=2 bins in {120, 160, 200}
#   - n_unique_lines == 9447 for all three
#   - k_best identical across the three (tolerance 1e-6)
#   - scan_null_p <= 1 / (null_n + 1)  (zero exceedances)
#   - n_obs in [10.6, 10.9]
#   - deltaD > 0
# Otherwise FAIL. No PARTIAL verdict.
# ---------------------------------------------------------------------------

.invoked_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.invoked_file) > 0L) dirname(normalizePath(.invoked_file)) else "R"
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

MASTER_CSV <- file.path(REPO_ROOT, "tables_r", "nist_master_results_r.csv")
VERDICT_CSV <- file.path(REPO_ROOT, "tables_r", "nist_verdict_r.csv")

REQUIRED_BINS <- c(120L, 160L, 200L)
REQUIRED_LINES <- 9447L
K_TOL <- 1e-6
N_OBS_MIN <- 10.6
N_OBS_MAX <- 10.9

to_float <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) return(NA_real_)
  suppressWarnings(as.numeric(x))
}
to_int <- function(x) {
  f <- to_float(x)
  if (is.na(f)) return(NA_integer_)
  as.integer(round(f))
}

read_master <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.csv(path, colClasses = "character", check.names = FALSE,
                  stringsAsFactors = FALSE)
}

find_feii_rows <- function(rows) {
  found <- list()
  if (nrow(rows) == 0L) return(found)
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    sp <- trimws(as.character(r$species))
    ion <- to_int(r$ion)
    bins <- to_int(r$bins)
    if (!is.na(ion) && !is.na(bins) && tolower(sp) == "fe" &&
        ion == 2L && bins %in% REQUIRED_BINS) {
      key <- as.character(bins)
      if (is.null(found[[key]])) found[[key]] <- r
    }
  }
  found
}

# evaluate(): returns list(verdict, reasons, feii). Mirrors make_verdict.evaluate.
evaluate <- function(rows) {
  reasons <- character(0)
  feii <- find_feii_rows(rows)

  missing <- REQUIRED_BINS[!as.character(REQUIRED_BINS) %in% names(feii)]
  if (length(missing) > 0L) {
    reasons <- c(reasons, sprintf("missing Fe II rows for bins %s",
                                  paste0("[", paste(missing, collapse = ", "), "]")))
    return(list(verdict = "FAIL", reasons = reasons, feii = feii))
  }

  for (b in REQUIRED_BINS) {
    r <- feii[[as.character(b)]]
    status <- trimws(as.character(r$status))
    if (!status %in% c("ok", "reused")) {
      reasons <- c(reasons, sprintf("bins=%d has non-success status '%s'", b, status))
    }
  }

  lines_set <- unique(vapply(REQUIRED_BINS, function(b) to_int(feii[[as.character(b)]]$n_unique_lines), integer(1)))
  if (!identical(sort(lines_set), as.integer(REQUIRED_LINES))) {
    reasons <- c(reasons, sprintf("n_unique_lines mismatch: %s (expected %d)",
                                  paste(sort(lines_set), collapse = ","), REQUIRED_LINES))
  }

  k_vals <- vapply(REQUIRED_BINS, function(b) to_float(feii[[as.character(b)]]$k_best), double(1))
  if (any(is.na(k_vals))) {
    reasons <- c(reasons, sprintf("k_best missing or non-numeric: %s",
                                  paste(k_vals, collapse = ",")))
  } else if (max(k_vals) - min(k_vals) > K_TOL) {
    reasons <- c(reasons, sprintf("k_best not stable within tol %g: %s",
                                  K_TOL, paste(k_vals, collapse = ",")))
  }

  for (b in REQUIRED_BINS) {
    r <- feii[[as.character(b)]]
    null_n <- to_int(r$null_n)
    p <- to_float(r$scan_null_p)
    if (is.na(null_n) || null_n <= 0L) {
      reasons <- c(reasons, sprintf("bins=%d invalid null_n=%s", b, as.character(r$null_n)))
    } else {
      threshold <- 1.0 / (null_n + 1)
      if (is.na(p) || p > threshold + 1e-12) {
        reasons <- c(reasons, sprintf("bins=%d scan_null_p=%s > 1/(null_n+1)=%s (need zero exceedances)",
                                      b, as.character(p), as.character(threshold)))
      }
    }
    n_obs <- to_float(r$n_obs)
    if (is.na(n_obs) || !(n_obs >= N_OBS_MIN && n_obs <= N_OBS_MAX)) {
      reasons <- c(reasons, sprintf("bins=%d n_obs=%s outside [%g, %g]",
                                    b, as.character(n_obs), N_OBS_MIN, N_OBS_MAX))
    }
    dd <- to_float(r$deltaD)
    if (is.na(dd) || dd <= 0) {
      reasons <- c(reasons, sprintf("bins=%d deltaD=%s not > 0", b, as.character(dd)))
    }
  }

  verdict <- if (length(reasons) == 0L) "PASS" else "FAIL"
  list(verdict = verdict, reasons = reasons, feii = feii)
}

write_verdict_csv <- function(verdict, reasons, feii) {
  dir.create(dirname(VERDICT_CSV), showWarnings = FALSE, recursive = TRUE)
  cols <- c("claim", "verdict",
            "bins_120_k_best", "bins_160_k_best", "bins_200_k_best",
            "bins_120_n_obs", "bins_160_n_obs", "bins_200_n_obs",
            "bins_120_deltaD", "bins_160_deltaD", "bins_200_deltaD",
            "bins_120_p", "bins_160_p", "bins_200_p", "reasons")
  row <- setNames(as.list(rep("", length(cols))), cols)
  row$claim <- "Fe II log-cosine k_best bin-stability (bins 120/160/200)"
  row$verdict <- verdict
  getv <- function(b, field) {
    r <- feii[[as.character(b)]]
    if (is.null(r)) "" else as.character(r[[field]])
  }
  for (b in REQUIRED_BINS) {
    row[[sprintf("bins_%d_k_best", b)]] <- getv(b, "k_best")
    row[[sprintf("bins_%d_n_obs", b)]] <- getv(b, "n_obs")
    row[[sprintf("bins_%d_deltaD", b)]] <- getv(b, "deltaD")
    row[[sprintf("bins_%d_p", b)]] <- getv(b, "scan_null_p")
  }
  row$reasons <- paste(reasons, collapse = "; ")
  df <- as.data.frame(row[cols], stringsAsFactors = FALSE)
  utils::write.csv(df, VERDICT_CSV, row.names = FALSE)
  cat(sprintf("[save] %s\n", sub(paste0("^", REPO_ROOT, "/?"), "", VERDICT_CSV)))
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  rows <- read_master(MASTER_CSV)
  if (nrow(rows) == 0L) {
    cat(sprintf("[error] no rows in %s. Run R/run_feii_bin_stability.R first.\n", MASTER_CSV))
    quit(status = 2)
  }
  ev <- evaluate(rows)
  write_verdict_csv(ev$verdict, ev$reasons, ev$feii)
  cat(sprintf("[verdict] %s\n", ev$verdict))
  if (length(ev$reasons) > 0L) {
    for (r in ev$reasons) cat(sprintf("  - %s\n", r))
  }
  quit(status = 0)
}

if (length(.invoked_file) > 0L && grepl("make_verdict\\.R$", .invoked_file)) {
  main()
}
