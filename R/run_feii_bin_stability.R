#!/usr/bin/env Rscript
# run_feii_bin_stability.R
# ---------------------------------------------------------------------------
# Run the Fe II log-cosine bin-stability ladder (bins 120, 160, 200) with the
# R scanner and aggregate the results into an R-generated master table:
#
#     tables_r/nist_master_results_r.csv
#
# R outputs go to outputs_r/fe_ion2_{120,160,200}/ so that the canonical Python
# results under outputs/ are never touched. If a folder already contains
# nist_summary.json that run is reused unless --force is supplied.
#
# Canonical Python reference values (regression targets, NOT hard-coded here):
#   120 bins: k_best=31.3265306122449  n_obs~10.717134  deltaD~355.7443
#   160 bins: k_best=31.3265306122449  n_obs~10.739649  deltaD~315.7277
#   200 bins: k_best=31.3265306122449  n_obs~10.753158  deltaD~259.1788
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages(library(jsonlite)))

.invoked_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
SCRIPT_DIR <- if (length(.invoked_file) > 0L) dirname(normalizePath(.invoked_file)) else "R"
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

SCANNER <- file.path(SCRIPT_DIR, "nist_wct_log_spectral_scan.R")
FE_CSV <- file.path(REPO_ROOT, "data", "Fe_lines.csv")
OUT_ROOT <- file.path(REPO_ROOT, "outputs_r")
TABLES_DIR <- file.path(REPO_ROOT, "tables_r")
MASTER_CSV <- file.path(TABLES_DIR, "nist_master_results_r.csv")

BINS_LADDER <- c(120L, 160L, 200L)
MASTER_COLUMNS <- c(
  "species", "ion", "bins", "n_unique_lines", "k_best", "n_obs", "deltaD",
  "scan_null_p", "tail_count_ge", "null_n", "ell_min", "ell_max", "delta_ell",
  "out_dir", "status"
)

parse_args <- function(argv) {
  vals <- list(`null-n` = 500L, `min-lines` = 100L, `n-k` = 2500L,
               seed = 20260517L, force = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (key == "force") {
      vals$force <- TRUE
      i <- i + 1L
    } else if (key %in% names(vals)) {
      vals[[key]] <- if (key == "seed" || key %in% c("null-n", "min-lines", "n-k")) {
        as.integer(argv[i + 1L])
      } else {
        argv[i + 1L]
      }
      i <- i + 2L
    } else {
      stop(sprintf("Unknown option: %s", argv[i]))
    }
  }
  vals
}

load_summary <- function(out_dir) {
  p <- file.path(out_dir, "nist_summary.json")
  if (!file.exists(p)) return(NULL)
  jsonlite::fromJSON(p)
}

summary_to_row <- function(species, ion, bins, out_dir, status) {
  row <- setNames(as.list(rep("", length(MASTER_COLUMNS))), MASTER_COLUMNS)
  row$species <- species
  row$ion <- ion
  row$bins <- bins
  row$out_dir <- sub(paste0("^", REPO_ROOT, "/?"), "", out_dir)
  row$status <- status
  s <- load_summary(out_dir)
  if (!is.null(s)) {
    row$n_unique_lines <- s$n_unique_lines
    row$k_best <- s$best$k_best
    row$n_obs <- s$branch_report$n_obs
    row$deltaD <- s$best$deltaD
    row$scan_null_p <- if (is.null(s$scan_null_p)) "" else s$scan_null_p
    row$tail_count_ge <- if (is.null(s$tail_count_ge)) "" else s$tail_count_ge
    row$null_n <- s$null_n
    row$ell_min <- s$ell_min
    row$ell_max <- s$ell_max
    row$delta_ell <- s$delta_ell
  }
  row
}

run_scanner <- function(bins, null_n, out_dir, min_lines, seed, n_k) {
  cmd_args <- c(
    SCANNER,
    "--csv", FE_CSV,
    "--ion", "2",
    "--bins", as.character(bins),
    "--null-n", as.character(null_n),
    "--min-lines", as.character(min_lines),
    "--out-dir", out_dir,
    "--n-k", as.character(n_k),
    "--seed", as.character(seed)
  )
  cat("[run] Rscript", paste(cmd_args, collapse = " "), "\n")
  status <- system2("Rscript", cmd_args, stdout = "", stderr = "")
  status
}

write_master <- function(rows) {
  dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
  df <- do.call(rbind, lapply(rows, function(r) {
    as.data.frame(r[MASTER_COLUMNS], stringsAsFactors = FALSE)
  }))
  # write full numeric precision
  utils::write.csv(df, MASTER_CSV, row.names = FALSE)
  cat(sprintf("[save] %s\n", sub(paste0("^", REPO_ROOT, "/?"), "", MASTER_CSV)))
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- parse_args(argv)
  if (!file.exists(SCANNER)) stop(sprintf("scanner not found: %s", SCANNER))
  if (!file.exists(FE_CSV)) stop(sprintf("Fe CSV not found: %s", FE_CSV))

  dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
  rows <- list()
  for (bins in BINS_LADDER) {
    out_dir <- file.path(OUT_ROOT, sprintf("fe_ion2_%d", bins))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    existing <- load_summary(out_dir)
    if (!is.null(existing) && !args$force) {
      cat(sprintf("[reuse] %s (nist_summary.json present)\n",
                  sub(paste0("^", REPO_ROOT, "/?"), "", out_dir)))
      status <- "reused"
    } else {
      rc <- run_scanner(bins, args[["null-n"]], out_dir,
                        args[["min-lines"]], args$seed, args[["n-k"]])
      status <- if (rc == 0L) "ok" else sprintf("failed_rc%d", rc)
    }
    rows[[length(rows) + 1L]] <- summary_to_row("Fe", 2L, bins, out_dir, status)
  }
  write_master(rows)
  invisible(0L)
}

if (length(.invoked_file) > 0L && grepl("run_feii_bin_stability\\.R$", .invoked_file)) {
  main()
}
