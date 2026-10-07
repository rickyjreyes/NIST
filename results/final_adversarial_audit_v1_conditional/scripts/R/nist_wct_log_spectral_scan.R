#!/usr/bin/env Rscript
# nist_wct_log_spectral_scan.R
# ---------------------------------------------------------------------------
# Command-line driver for the NIST Atomic Spectra log-cosine spectral scan.
#
# Independent R computational reproduction of the canonical Python reference
# scripts/nist_wct_log_spectral_scan_FIXED.py. Same scientific definitions,
# scan grid, output meanings and PASS/FAIL inputs.
#
# Example:
#   Rscript R/nist_wct_log_spectral_scan.R \
#       --csv data/Fe_lines.csv --ion 2 --bins 160 \
#       --null-n 500 --min-lines 100 --out-dir outputs_r/fe_ion2_160
#
# No NIST endorsement, certification, or validation is claimed.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  # locate this script to source the shared library next to it
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
  if (length(file_arg) == 0L) {
    # interactive / sourced fallback
    script_dir <- "R"
  } else {
    script_dir <- dirname(normalizePath(file_arg))
  }
  source(file.path(script_dir, "nist_scan_lib.R"))
}))

OUT_DIR_DEFAULT <- "outputs_nist_wct"

# --- minimal "--flag value" argument parser --------------------------------
# Kept dependency-free; supports the documented flags only.
parse_args <- function(argv) {
  spec <- list(
    csv = list(type = "character", default = NULL),
    `out-dir` = list(type = "character", default = OUT_DIR_DEFAULT),
    ion = list(type = "integer", default = NA_integer_),
    bins = list(type = "integer", default = 160L),
    `ell-min` = list(type = "double", default = NA_real_),
    `ell-max` = list(type = "double", default = NA_real_),
    `k-min` = list(type = "double", default = 0.5),
    `k-max` = list(type = "double", default = 80.0),
    `n-k` = list(type = "integer", default = 2500L),
    degree = list(type = "integer", default = 1L),
    `baseline-sigma` = list(type = "double", default = 6.0),
    `null-n` = list(type = "integer", default = 500L),
    seed = list(type = "integer", default = 20260517L),
    `min-lines` = list(type = "integer", default = 100L)
  )
  vals <- lapply(spec, function(s) s$default)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("-h", "--help")) {
      cat("Usage: Rscript R/nist_wct_log_spectral_scan.R --csv FILE [options]\n")
      cat("Options: ", paste0("--", names(spec)), "\n")
      quit(status = 0)
    }
    if (!startsWith(a, "--")) {
      stop(sprintf("Unexpected argument: %s", a))
    }
    key <- sub("^--", "", a)
    if (!key %in% names(spec)) {
      stop(sprintf("Unknown option: --%s", key))
    }
    if (i + 1L > length(argv)) stop(sprintf("Missing value for --%s", key))
    raw <- argv[i + 1L]
    vals[[key]] <- switch(spec[[key]]$type,
      character = raw,
      integer = as.integer(raw),
      double = as.numeric(raw)
    )
    i <- i + 2L
  }
  if (is.null(vals$csv)) stop("--csv is required")
  vals
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- parse_args(argv)

  out_dir <- args[["out-dir"]]
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ion <- if (is.na(args$ion)) NULL else args$ion
  ell_min <- if (is.na(args[["ell-min"]])) NULL else args[["ell-min"]]
  ell_max <- if (is.na(args[["ell-max"]])) NULL else args[["ell-max"]]

  raw <- read_nist_csv(args$csv)
  cat(sprintf("[raw] rows=%d cols=%s\n", nrow(raw), paste(colnames(raw), collapse = ", ")))

  lines <- clean_lines(raw, ion = ion)
  cat(sprintf("[clean] usable unique lines=%d\n", nrow(lines)))
  if (nrow(lines) < args[["min-lines"]]) {
    stop(sprintf("Only %d usable lines. Need at least %d.", nrow(lines), args[["min-lines"]]))
  }

  utils::write.csv(lines, file.path(out_dir, "nist_lines_clean.csv"), row.names = FALSE)

  binned <- build_binned(lines, args$bins, ell_min, ell_max)
  y <- binned$count
  ell <- binned$ell
  baseline <- pmax(gaussian_filter_nearest(y, args[["baseline-sigma"]]), EPS)
  binned$baseline <- baseline
  utils::write.csv(binned, file.path(out_dir, "nist_binned_spectrum.csv"), row.names = FALSE)

  k_grid <- seq(args[["k-min"]], args[["k-max"]], length.out = args[["n-k"]])
  cat(sprintf("[scan] bins=%d k=%g..%g n_k=%d\n",
              args$bins, args[["k-min"]], args[["k-max"]], args[["n-k"]]))
  sk <- scan_k(ell, y, baseline, k_grid, args$degree)
  scan <- sk$scan
  best <- sk$best
  mu0 <- sk$mu0
  utils::write.csv(scan, file.path(out_dir, "nist_scan_curve.csv"), row.names = FALSE)

  null_vals <- NULL
  p <- NA_real_
  tail <- NA_integer_
  if (args[["null-n"]] > 0L) {
    ns <- null_scan(ell, y, baseline, mu0, k_grid, args$degree,
                    best$deltaD, args[["null-n"]], args$seed)
    null_vals <- ns$vals
    p <- ns$p
    tail <- ns$tail
    utils::write.csv(data.frame(null_max_deltaD = null_vals),
                     file.path(out_dir, "nist_null.csv"), row.names = FALSE)
  }

  delta_ell <- max(ell) - min(ell)
  br <- branch_report(best$k_best, delta_ell)
  ratios <- peak_ratios(scan)
  if (nrow(ratios) > 0L) {
    utils::write.csv(ratios, file.path(out_dir, "nist_peak_ratios.csv"), row.names = FALSE)
  }

  write_summary(out_dir, args, raw, lines, ell, delta_ell, best, p, tail, br, ratios)

  save_plots(out_dir, binned, baseline, mu0, scan, best, null_vals)

  cat("\n[done]\n")
  cat(jsonlite::toJSON(list(
    k_best = best$k_best, deltaD = best$deltaD,
    scan_null_p = p,
    tail_count_ge = tail,
    null_n = args[["null-n"]], n_obs = br$n_obs,
    n_unique_lines = nrow(lines)
  ), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"), "\n")
}

# write_summary(): emit nist_summary.json matching the Python key structure.
write_summary <- function(out_dir, args, raw, lines, ell, delta_ell, best,
                          p, tail, br, ratios) {
  nearest_list <- lapply(seq_len(nrow(br$nearest)), function(i) {
    r <- br$nearest[i, ]
    list(branch = r$branch, n = r$n, k_target = r$k_target,
         k_error = r$k_error, abs_k_error = r$abs_k_error,
         n_error = r$n_error, abs_n_error = r$abs_n_error)
  })
  top_ratios <- if (nrow(ratios) > 0L) {
    rr <- utils::head(ratios, 20L)
    lapply(seq_len(nrow(rr)), function(i) {
      r <- rr[i, ]
      list(k_low = r$k_low, k_high = r$k_high, ratio = r$ratio,
           rational = r$rational, rational_error = r$rational_error)
    })
  } else {
    list()
  }
  summary <- list(
    csv = args$csv,
    ion_filter = if (is.na(args$ion)) NA_integer_ else args$ion,
    n_raw_rows = nrow(raw),
    n_unique_lines = nrow(lines),
    wavenumber_min_cm = min(lines$wavenumber_cm),
    wavenumber_max_cm = max(lines$wavenumber_cm),
    ell_min = min(ell),
    ell_max = max(ell),
    delta_ell = delta_ell,
    best = best,
    scan_null_p = p,
    tail_count_ge = tail,
    null_n = args[["null-n"]],
    branch_report = list(n_obs = br$n_obs, nearest = nearest_list),
    top_peak_ratios = top_ratios
  )
  writeLines(
    jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"),
    file.path(out_dir, "nist_summary.json")
  )
}

# save_plots(): base-graphics equivalents of the Python matplotlib figures.
save_plots <- function(out_dir, binned, baseline, mu0, scan, best, null_vals) {
  grDevices::png(file.path(out_dir, "nist_spectrum_fit.png"),
                 width = 10, height = 5, units = "in", res = 180)
  plot(binned$ell, binned$count, type = "s", col = "grey40",
       xlab = "ell = ln(wavenumber cm^-1)", ylab = "counts/bin",
       main = "")
  graphics::lines(binned$ell, baseline, col = "blue")
  graphics::lines(binned$ell, mu0, col = "red")
  graphics::legend("topright", legend = c("line counts", "smooth baseline", "base fit"),
                   col = c("grey40", "blue", "red"), lty = 1, bty = "n")
  grDevices::dev.off()

  grDevices::png(file.path(out_dir, "nist_scan_curve.png"),
                 width = 10, height = 5, units = "in", res = 180)
  plot(scan$k, scan$deltaD, type = "l", xlab = "k", ylab = "DeltaD", main = "")
  graphics::abline(v = best$k_best, lty = 2)
  graphics::legend("topright", legend = sprintf("k_best=%.4g", best$k_best),
                   lty = 2, bty = "n")
  grDevices::dev.off()

  if (!is.null(null_vals)) {
    grDevices::png(file.path(out_dir, "nist_null.png"),
                   width = 10, height = 5, units = "in", res = 180)
    graphics::hist(null_vals, breaks = 50, main = "", xlab = "null max deltaD")
    graphics::abline(v = best$deltaD, lty = 2)
    graphics::legend("topright", legend = sprintf("obs=%.4g", best$deltaD),
                     lty = 2, bty = "n")
    grDevices::dev.off()
  }
}

# Auto-run only when executed directly via Rscript (the --file= argument names
# this script), not when sourced for unit tests.
.invoked_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(.invoked_file) > 0L && grepl("nist_wct_log_spectral_scan\\.R$", .invoked_file)) {
  main()
}
