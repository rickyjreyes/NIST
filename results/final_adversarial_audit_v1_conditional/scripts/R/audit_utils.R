# audit_utils.R
# ---------------------------------------------------------------------------
# Shared utilities for the NIST statistical-audit suite.
#
# This module sits *around* the canonical reference scanner. It sources the
# unchanged R/nist_scan_lib.R (the canonical computational reproduction of
# scripts/nist_wct_log_spectral_scan_FIXED.py) and adds reusable helpers for
# the statistical audit: deterministic RNG setup, empirical-p helpers, peak
# conversions, source-specific dataset construction, a colour-blind-safe
# species palette, a shared ggplot theme, figure/table writers, and small
# parallel/registry/git helpers.
#
# No numerical convention in nist_scan_lib.R is modified here. The audit
# computes additional statistics from the same deterministic pipeline.
#
# No NIST endorsement, certification, or validation is claimed. NIST is the
# data provider only.
# ---------------------------------------------------------------------------

# --- locate and source the canonical library -------------------------------

audit_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
  if (length(file_arg) == 0L) {
    # interactive / sourced fallback: try common locations
    cands <- c("R", file.path(getwd(), "R"), ".")
    for (d in cands) if (file.exists(file.path(d, "nist_scan_lib.R"))) return(normalizePath(d))
    return(normalizePath("R", mustWork = FALSE))
  }
  dirname(normalizePath(file_arg))
}

.AUDIT_DIR <- audit_script_dir()
if (!exists("scan_k", mode = "function")) {
  source(file.path(.AUDIT_DIR, "nist_scan_lib.R"))
}

# --- repository roots ------------------------------------------------------

audit_repo_root <- function() normalizePath(file.path(.AUDIT_DIR, ".."), mustWork = FALSE)

# --- deterministic RNG -----------------------------------------------------

# setup_rng(): use the parallel-safe L'Ecuyer-CMRG generator so that results
# do not depend on the number of workers. All audit entry points call this.
setup_rng <- function(seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  invisible(seed)
}

# --- empirical Monte Carlo p-values ----------------------------------------

# emp_p(): the standard +1 corrected empirical p-value, p = (r + 1)/(B + 1).
# Never returns 0.
emp_p <- function(tail_count, B) (tail_count + 1) / (B + 1)

# resolution_floor(): the empirical resolution floor 1/(B+1). When r = 0 the
# corrected estimate equals this value; it is a resolution-limited bound, NOT
# an exact p-value below it.
resolution_floor <- function(B) 1 / (B + 1)

# --- peak-derived conversions ----------------------------------------------
#
# Convention (matching nist_scan_lib.R): the harmonic model uses cos(k*ell) and
# sin(k*ell), so k is an ANGULAR frequency in radians per unit of ell =
# ln(wavenumber/cm^-1). The log-period in ell is therefore 2*pi/k and the
# multiplicative scale ratio between successive log-periodic features is
# exp(2*pi/k). We never silently switch to a cyclic frequency f = k/(2*pi).

delta_log_x <- function(k) {
  out <- 2 * pi / k
  out[!is.finite(k) | k <= 0] <- NA_real_
  out
}

scale_ratio <- function(k) {
  d <- delta_log_x(k)
  exp(d)
}

# n_obs convention from branch_report(): k * delta_ell / (2*pi).
n_obs_from_k <- function(k, delta_ell) k * delta_ell / (2 * pi)

# --- source-specific wavenumber extraction ---------------------------------
#
# The audit explicitly distinguishes three measurement representations:
#   "wavenumber" : the direct NIST wn(cm-1) column (canonical scanner default)
#   "observed"   : derived from the observed air/vac wavelength column
#   "ritz"       : derived from the Ritz air/vac wavelength column
# These are NOT silently mixed. extract_wavenumber_source() returns the raw
# per-row wavenumber vector (with NA where the chosen field is missing) plus a
# descriptive label, without dropping rows.

extract_wavenumber_source <- function(df, source = c("wavenumber", "observed", "ritz")) {
  source <- match.arg(source)
  if (source == "wavenumber") {
    wn_col <- find_column(df, c(
      "wn(cm-1)", "wn", "wavenumber", "wavenumber_cm", "wavenumber in cm-1"
    ))
    if (is.null(wn_col)) {
      return(list(wn = rep(NA_real_, nrow(df)), source = "wn_missing", field = NA_character_))
    }
    return(list(wn = numeric_from_nist_cell(df[[wn_col]]),
                source = wn_col, field = wn_col))
  }
  if (source == "observed") {
    col <- find_column(df, c("obs_wl_air(nm)", "obs_wl_vac(nm)", "obs_wl_air", "obs_wl"))
  } else {
    col <- find_column(df, c("ritz_wl_air(nm)", "ritz_wl_vac(nm)", "ritz_wl_air", "ritz_wl"))
  }
  if (is.null(col)) {
    return(list(wn = rep(NA_real_, nrow(df)),
                source = sprintf("%s_missing", source), field = NA_character_))
  }
  wl_nm <- numeric_from_nist_cell(df[[col]])
  wn <- 1.0e7 / wl_nm
  list(wn = wn, source = sprintf("converted_from_%s", col), field = col)
}

# clean_lines_source(): like clean_lines() but with an explicit source field
# and a fully reconciled accounting record. Returns
#   list(lines = data.frame, flow = named list of counts).
# The sequential flow reconciles exactly:
#   n_raw = retained + excl_species + excl_ion + excl_missing_source
#           + excl_nonpositive + n_duplicates
clean_lines_source <- function(df, species = NULL, ion = NULL,
                               source = c("wavenumber", "observed", "ritz")) {
  source <- match.arg(source)
  n_raw <- nrow(df)
  work <- df

  # 1. species filter (by element column when present and species supplied)
  excl_species <- 0L
  if (!is.null(species) && "element" %in% colnames(work)) {
    el <- trimws(as.character(work[["element"]]))
    keep <- !is.na(el) & tolower(el) == tolower(species)
    excl_species <- sum(!keep)
    work <- work[keep, , drop = FALSE]
  }

  # 2. ion filter (by sp_num)
  excl_ion <- 0L
  if (!is.null(ion) && "sp_num" %in% colnames(work)) {
    sp <- numeric_from_nist_cell(work[["sp_num"]])
    keep <- !is.na(sp) & sp == as.numeric(ion)
    excl_ion <- sum(!keep)
    work <- work[keep, , drop = FALSE]
  }
  n_valid_species_ion <- nrow(work)

  # independent diagnostic missing counts (on the species/ion subset)
  obs_wn  <- extract_wavenumber_source(work, "observed")$wn
  ritz_wn <- extract_wavenumber_source(work, "ritz")$wn
  wn_wn   <- extract_wavenumber_source(work, "wavenumber")$wn
  missing_observed   <- sum(!is.finite(obs_wn))
  missing_ritz       <- sum(!is.finite(ritz_wn))
  missing_wavenumber <- sum(!is.finite(wn_wn))

  # 3. source wavenumber + missing/nonnumeric exclusion
  ex <- extract_wavenumber_source(work, source)
  wn <- ex$wn
  na_source <- !is.finite(wn) & !(is.finite(wn) & wn <= 0)
  excl_missing_source <- sum(is.na(wn))
  nonnumeric_count <- excl_missing_source
  keep_present <- !is.na(wn)
  work <- work[keep_present, , drop = FALSE]
  wn <- wn[keep_present]

  # 4. nonpositive / nonfinite exclusion
  finite_pos <- is.finite(wn) & wn > 0
  excl_nonpositive <- sum(!finite_pos)
  work <- work[finite_pos, , drop = FALSE]
  wn <- wn[finite_pos]

  # attach derived columns
  work$wavenumber_cm <- wn
  work$ell <- log(wn)
  work$wavenumber_source <- ex$source

  # 5. sort and drop duplicate wavenumbers (keep first)
  ord <- order(work$wavenumber_cm)
  work <- work[ord, , drop = FALSE]
  dup <- duplicated(work$wavenumber_cm)
  n_duplicates <- sum(dup)
  work <- work[!dup, , drop = FALSE]
  rownames(work) <- NULL
  retained <- nrow(work)

  flow <- list(
    species = if (is.null(species)) NA_character_ else species,
    ion = if (is.null(ion)) NA_integer_ else as.integer(ion),
    source = source,
    source_field = ex$field,
    n_raw = n_raw,
    n_valid_species_ion = n_valid_species_ion,
    excl_species = as.integer(excl_species),
    excl_ion = as.integer(excl_ion),
    missing_observed = as.integer(missing_observed),
    missing_ritz = as.integer(missing_ritz),
    missing_wavenumber = as.integer(missing_wavenumber),
    nonnumeric_source = as.integer(nonnumeric_count),
    excl_missing_source = as.integer(excl_missing_source),
    excl_nonpositive = as.integer(excl_nonpositive),
    n_duplicates = as.integer(n_duplicates),
    retained = as.integer(retained),
    wavelength_min_nm = if (retained > 0) 1.0e7 / max(work$wavenumber_cm) else NA_real_,
    wavelength_max_nm = if (retained > 0) 1.0e7 / min(work$wavenumber_cm) else NA_real_,
    wavenumber_min_cm = if (retained > 0) min(work$wavenumber_cm) else NA_real_,
    wavenumber_max_cm = if (retained > 0) max(work$wavenumber_cm) else NA_real_,
    ell_min = if (retained > 0) min(work$ell) else NA_real_,
    ell_max = if (retained > 0) max(work$ell) else NA_real_
  )
  # exact reconciliation check
  flow$reconciles <- (flow$n_raw == flow$retained + flow$excl_species +
    flow$excl_ion + flow$excl_missing_source + flow$excl_nonpositive +
    flow$n_duplicates)
  list(lines = work, flow = flow)
}

# --- a single scan analysis ------------------------------------------------
#
# run_scan_analysis(): bin the lines, build the Gaussian baseline, run the
# canonical scan_k(), and return everything the audit needs. This is the one
# place that wires the canonical pipeline together so every audit module uses
# the exact same definitions.
run_scan_analysis <- function(lines, bins, k_grid, degree = 1L, baseline_sigma = 6.0,
                              ell_min = NULL, ell_max = NULL) {
  binned <- build_binned(lines, bins, ell_min, ell_max)
  y <- binned$count
  ell <- binned$ell
  baseline <- pmax(gaussian_filter_nearest(y, baseline_sigma), EPS)
  binned$baseline <- baseline
  sk <- scan_k(ell, y, baseline, k_grid, degree)
  # delta_ell uses the binned ell-range (bin centres), matching the canonical
  # scanner's branch_report()/n_obs convention exactly.
  delta_ell <- max(ell) - min(ell)
  list(
    binned = binned, ell = ell, y = y, baseline = baseline,
    scan = sk$scan, best = sk$best, mu0 = sk$mu0,
    delta_ell = delta_ell, n_obs = n_obs_from_k(sk$best$k_best, delta_ell),
    bins = bins, degree = degree, baseline_sigma = baseline_sigma,
    n_lines = nrow(lines)
  )
}

# nearest grid index to a target k
nearest_k_index <- function(k_grid, k_target) which.min(abs(k_grid - k_target))

# --- null distribution (global max + pointwise at fixed k) ------------------
#
# null_distribution(): parametric Poisson bootstrap under the fitted smooth
# null mu0. For each replicate it draws y0 ~ Poisson(mu0), reruns the full
# k-grid scan, and records BOTH the scan-global maximum deltaD (look-elsewhere
# corrected) and the deltaD at a prespecified fixed k index (pointwise). This
# keeps pointwise and global significance strictly separate.
#
# Determinism: a vector of per-replicate seeds is drawn ONCE from the current
# stream, then each replicate sets its own seed. Results are therefore
# identical whether run sequentially or in parallel, and independent of the
# number of workers.
null_distribution <- function(ell, y, baseline, mu0, k_grid, degree,
                              null_n, fixed_k_index = NULL, verbose = FALSE,
                              parallel = FALSE) {
  seeds <- sample.int(.Machine$integer.max, null_n)
  one <- function(i) {
    set.seed(seeds[i])
    y0 <- rpois(length(mu0), mu0)
    res <- scan_k(ell, y0, baseline, k_grid, degree)
    c(max = res$best$deltaD,
      point = if (is.null(fixed_k_index)) NA_real_ else res$scan$deltaD[fixed_k_index])
  }
  out <- audit_lapply(seq_len(null_n), one, parallel = parallel)
  m <- do.call(rbind, out)
  list(max_vals = m[, "max"],
       point_vals = if (is.null(fixed_k_index)) NULL else m[, "point"])
}

# --- parallel helper (graceful fallback) -----------------------------------
#
# audit_lapply(): map over a list, using future.apply::future_lapply when the
# package is available and parallel=TRUE, otherwise base lapply. Deterministic
# RNG is preserved via L'Ecuyer-CMRG streams (future.seed = seed).
audit_lapply <- function(X, FUN, parallel = FALSE, seed = NULL, ...) {
  use_future <- isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)
  if (use_future) {
    future.apply::future_lapply(X, FUN, future.seed = if (is.null(seed)) TRUE else seed, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

# configure_parallel(): set up a future plan when requested and available.
configure_parallel <- function(parallel = FALSE, workers = "auto") {
  if (!isTRUE(parallel)) return(invisible(FALSE))
  if (!requireNamespace("future", quietly = TRUE)) {
    message("[parallel] 'future' not installed; running sequentially.")
    return(invisible(FALSE))
  }
  nw <- if (identical(workers, "auto")) max(1L, future::availableCores() - 1L) else as.integer(workers)
  future::plan(future::multisession, workers = nw)
  message(sprintf("[parallel] future multisession with %d workers", nw))
  invisible(TRUE)
}

# --- colour-blind-safe species palette and theme ---------------------------
#
# One consistent species colour mapping throughout, using the Okabe-Ito
# colour-blind-safe palette. Fe is given the strongest (orange) so it stands
# out without hiding neighbours.
SPECIES_PALETTE <- c(
  Fe = "#D55E00",  # vermillion  (focus)
  Cr = "#0072B2",  # blue
  Mn = "#009E73",  # green
  Co = "#CC79A7",  # purple/pink
  Ni = "#56B4E9",  # sky blue
  Ti = "#E69F00",  # orange
  null = "#999999" # grey for null/reference
)

species_colour <- function(sp) {
  sp <- as.character(sp)
  out <- SPECIES_PALETTE[sp]
  out[is.na(out)] <- "#444444"
  unname(out)
}

# theme_audit(): a clean publication theme used across all figures.
theme_audit <- function(base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      plot.caption = ggplot2::element_text(colour = "grey45", size = base_size - 3),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# --- figure / table writers ------------------------------------------------

# save_fig(): write a ggplot (or base plotting closure) to a 300-dpi PNG and a
# vector companion (SVG via svglite, else PDF). Returns the PNG path.
save_fig <- function(plot, path, width = 9, height = 6, dpi = 300, vector = TRUE) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) {
    if (!is.null(png_dev)) {
      ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, device = ragg::agg_png)
    } else {
      ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
    }
    if (vector) {
      vec_path <- sub("\\.png$", ".svg", path)
      ok <- tryCatch({
        if (requireNamespace("svglite", quietly = TRUE)) {
          ggplot2::ggsave(vec_path, plot, width = width, height = height, device = svglite::svglite)
        } else {
          ggplot2::ggsave(sub("\\.png$", ".pdf", path), plot, width = width, height = height)
        }
        TRUE
      }, error = function(e) FALSE)
    }
  } else if (is.function(plot)) {
    grDevices::png(path, width = width, height = height, units = "in", res = dpi)
    plot(); grDevices::dev.off()
  }
  path
}

# write_table(): write a data.frame to CSV with full numeric precision.
write_table <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

# --- git / provenance ------------------------------------------------------

git_commit <- function() {
  out <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  if (length(out) == 0L || is.na(out[1])) return(NA_character_)
  trimws(out[1])
}

# audit_timestamp(): ISO-8601 UTC timestamp.
audit_timestamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

# --- region tolerance helpers ----------------------------------------------
#
# in_peak_region(): TRUE when k is within the chosen tolerance of the reference
# peak. Two definitions are supported and documented in R/peak_stability.R.
in_peak_region <- function(k, k_ref, type = c("relative", "absolute"),
                           tol_relative = 0.02, tol_absolute = NULL) {
  type <- match.arg(type)
  if (type == "relative") {
    abs(k - k_ref) / k_ref <= tol_relative
  } else {
    if (is.null(tol_absolute)) stop("tol_absolute required for absolute tolerance")
    abs(k - k_ref) <= tol_absolute
  }
}

# --- small numeric helpers -------------------------------------------------

cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || mean(x) == 0) return(NA_real_)
  stats::sd(x) / abs(mean(x))
}

iqr_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::IQR(x)
}

# --- shared audit configuration --------------------------------------------
#
# default_audit_config(): the single source of truth for the canonical Fe II
# analysis parameters and the predefined grids used throughout the audit. The
# scan grid and degree/baseline defaults match the canonical scanner so that
# the audit's reference peak coincides with the committed results.
default_audit_config <- function(seed = 20260517L,
                                 bootstrap_n = 2000L,
                                 null_n = 5000L,
                                 calibration_n = 2000L,
                                 injection_n = 500L,
                                 fast = FALSE) {
  if (isTRUE(fast)) {
    bootstrap_n <- 80L; null_n <- 25L; calibration_n <- 25L; injection_n <- 8L
  }
  list(
    seed = as.integer(seed),
    species = "Fe", ion = 2L,
    k_min = 0.5, k_max = 80.0, n_k = 2500L,
    degree = 1L, baseline_sigma = 6.0,
    min_lines = 100L,
    bins_primary = 160L,
    bins_grid = c(60L, 80L, 100L, 120L, 140L, 160L, 180L, 200L, 220L, 240L),
    sigma_grid = c(4, 5, 6, 7, 8),
    degree_grid = c(0L, 1L, 2L),
    bootstrap_n = as.integer(bootstrap_n),
    null_n = as.integer(null_n),
    calibration_n = as.integer(calibration_n),
    injection_n = as.integer(injection_n),
    # peak-region tolerances (relative is the primary confirmatory definition)
    tol_primary = 0.02,
    tol_grid = c(0.01, 0.02, 0.05),
    injection_amplitudes = if (isTRUE(fast)) c(0.00, 0.03, 0.05, 0.10)
                           else c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10),
    fast = isTRUE(fast)
  )
}

audit_k_grid <- function(cfg) seq(cfg$k_min, cfg$k_max, length.out = cfg$n_k)

# is_fast(): detect the --fast development flag in an argv vector.
is_fast <- function(argv) any(argv == "--fast")

# stability_class(): descriptive (not universal) stability labels.
stability_class <- function(pct) {
  ifelse(pct >= 0.80, "high",
         ifelse(pct >= 0.50, "moderate", "low"))
}

# binomial Wilson-ish CI via Clopper-Pearson (exact) for a rate.
binom_ci <- function(successes, n, conf = 0.95) {
  if (n == 0L) return(c(lower = NA_real_, upper = NA_real_))
  alpha <- 1 - conf
  lower <- if (successes == 0L) 0 else stats::qbeta(alpha / 2, successes, n - successes + 1)
  upper <- if (successes == n) 1 else stats::qbeta(1 - alpha / 2, successes + 1, n - successes)
  c(lower = lower, upper = upper)
}
