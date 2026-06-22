# nist_scan_lib.R
# ---------------------------------------------------------------------------
# Reusable functions for the NIST Atomic Spectra log-cosine spectral scan.
#
# This is an independent R *computational reproduction* of the canonical
# Python reference implementation:
#
#     scripts/nist_wct_log_spectral_scan_FIXED.py
#
# The numerical conventions below are chosen deliberately to reproduce the
# Python/NumPy/SciPy behaviour bit-for-bit on the deterministic pipeline
# (reading -> cleaning -> binning -> Gaussian baseline -> IRLS Poisson fit ->
# deviance scan). See R/README.md for the documented conventions and the
# known Python/R RNG difference in the parametric bootstrap.
#
# No NIST endorsement, certification, or validation is claimed. NIST is the
# data provider only.
# ---------------------------------------------------------------------------

# Epsilon floor, matching the Python reference EPS = 1e-12.
EPS <- 1e-12

# WCT reference winding numbers (kept identical to the Python reference).
WCT_REFERENCE_NS <- c(10.0, 15.0, 20.0)
WCT_FOLDED_NS <- c(20.0 / 3.0, 15.0, 40.0 / 3.0)


# --- cell cleaning ---------------------------------------------------------

# clean_cell(): strip NIST/Excel export decoration from a single cell.
# NIST often exports numeric cells as ="123.45". This mirrors the Python
# clean_cell() exactly: handle the ="..." form, a bare leading "=", and then
# strip surrounding double quotes and whitespace.
clean_cell <- function(x) {
  if (length(x) != 1L) {
    return(vapply(x, clean_cell, character(1), USE.NAMES = FALSE))
  }
  if (is.na(x)) {
    return(x)
  }
  s <- trimws(as.character(x))
  if (startsWith(s, '="') && endsWith(s, '"')) {
    s <- substr(s, 3L, nchar(s) - 1L)
  } else if (startsWith(s, "=")) {
    s <- substr(s, 2L, nchar(s))
  }
  # strip surrounding double quotes, then surrounding whitespace
  s <- gsub('^"+|"+$', "", s)
  s <- trimws(s)
  s
}

# numeric_from_nist_cell(): clean a vector of cells and extract the leading
# numeric token, mirroring Python to_float_series(). Flags such as "2300d?",
# "250bl(Fe III)" or "(0)" are reduced to their leading number; empty cells
# and pure-non-numeric cells become NA.
numeric_from_nist_cell <- function(x) {
  cleaned <- vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    clean_cell(v)
  }, character(1), USE.NAMES = FALSE)
  # Match the first signed integer/float with optional exponent, like the
  # Python regex ([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?).
  out <- rep(NA_real_, length(cleaned))
  matched <- regexpr("[-+]?[0-9]+(?:\\.[0-9]*)?(?:[eE][-+]?[0-9]+)?", cleaned)
  pos <- !is.na(matched) & matched > 0L
  tok <- regmatches(cleaned, matched)   # drops NA / non-matches
  out[pos] <- suppressWarnings(as.numeric(tok))
  out
}


# --- CSV reading -----------------------------------------------------------

# read_nist_csv(): read a NIST CSV robustly and apply clean_cell to every
# cell. Columns are returned as character, BOM stripped, names trimmed.
read_nist_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("CSV not found: %s", normalizePath(path, mustWork = FALSE)))
  }
  df <- utils::read.csv(
    path,
    colClasses = "character",
    check.names = FALSE,
    quote = "\"",
    na.strings = character(0),
    stringsAsFactors = FALSE
  )
  nm <- colnames(df)
  # strip BOM and trim, mirroring Python's str(c).strip().lstrip('﻿')
  nm <- sub("^﻿", "", nm)
  nm <- trimws(nm)
  colnames(df) <- nm
  # clean every cell
  for (j in seq_along(df)) {
    df[[j]] <- clean_cell(df[[j]])
  }
  df
}


# --- column detection ------------------------------------------------------

# norm_col(): normalise a column name to lowercase alnum tokens joined by "_".
norm_col <- function(c) {
  s <- tolower(as.character(c))
  s <- gsub("[^a-z0-9]+", "_", s)
  gsub("^_+|_+$", "", s)
}

# find_column(): locate a column by a list of candidate names, first by exact
# normalised match, then by substring containment (either direction).
find_column <- function(df, candidates) {
  cols <- colnames(df)
  norm <- setNames(cols, vapply(cols, norm_col, character(1)))
  for (cand in candidates) {
    nc <- norm_col(cand)
    if (nc %in% names(norm)) {
      return(unname(norm[[nc]]))
    }
  }
  for (cand in candidates) {
    nc <- norm_col(cand)
    for (i in seq_along(norm)) {
      k <- names(norm)[i]
      if (nzchar(nc) && nzchar(k) && (grepl(nc, k, fixed = TRUE) || grepl(k, nc, fixed = TRUE))) {
        return(unname(norm[[i]]))
      }
    }
  }
  NULL
}

# extract_wavenumber_cm(): return list(wn = numeric vector, source = label).
# Prefers an explicit wavenumber column; otherwise converts a wavelength-in-nm
# column via wavenumber_cm = 1e7 / wavelength_nm.
extract_wavenumber_cm <- function(df) {
  wn_col <- find_column(df, c(
    "wn(cm-1)", "wn", "wavenumber", "wavenumber_cm", "wavenumber in cm-1",
    "Ritz(cm-1)", "Observed(cm-1)"
  ))
  if (!is.null(wn_col)) {
    wn <- numeric_from_nist_cell(df[[wn_col]])
    if (sum(!is.na(wn)) > 0L) {
      return(list(wn = wn, source = wn_col))
    }
  }
  wl_nm_col <- find_column(df, c(
    "ritz_wl_air(nm)", "obs_wl_air(nm)", "ritz_wl_vac(nm)", "obs_wl_vac(nm)",
    "wavelength_nm", "wavelength nm"
  ))
  if (!is.null(wl_nm_col)) {
    wl_nm <- numeric_from_nist_cell(df[[wl_nm_col]])
    wn <- 1.0e7 / wl_nm
    return(list(wn = wn, source = sprintf("converted_from_%s", wl_nm_col)))
  }
  stop(sprintf(
    "Could not find wavenumber or wavelength column. Columns are: %s",
    paste(colnames(df), collapse = ", ")
  ))
}


# --- line cleaning ---------------------------------------------------------

# clean_lines(): optional ion filter via sp_num, extract wavenumber, keep
# finite positive wavenumbers, compute ell = log(wavenumber_cm), sort by
# wavenumber and drop duplicate wavenumbers. Returns a data.frame mirroring
# the Python columns we rely on downstream.
clean_lines <- function(df, ion = NULL) {
  out <- df
  if (!is.null(ion) && "sp_num" %in% colnames(out)) {
    sp <- numeric_from_nist_cell(out[["sp_num"]])
    keep <- !is.na(sp) & sp == as.numeric(ion)
    out <- out[keep, , drop = FALSE]
  }
  ex <- extract_wavenumber_cm(out)
  wn <- ex$wn
  finite_pos <- is.finite(wn) & wn > 0
  out <- out[finite_pos, , drop = FALSE]
  wn <- wn[finite_pos]
  out$wavenumber_cm <- wn
  out$ell <- log(wn)
  out$wavenumber_source <- ex$source
  # sort by wavenumber, then drop duplicate wavenumbers keeping first
  ord <- order(out$wavenumber_cm)
  out <- out[ord, , drop = FALSE]
  dup <- duplicated(out$wavenumber_cm)
  out <- out[!dup, , drop = FALSE]
  rownames(out) <- NULL
  out
}


# --- binning ---------------------------------------------------------------

# np_histogram_uniform(): replicate numpy.histogram(a, bins, range) for the
# uniform-bin fast path, including NumPy's floating-point edge corrections.
# Returns list(counts, edges).
np_histogram_uniform <- function(a, bins, first_edge, last_edge) {
  a <- as.numeric(a)
  bin_edges <- seq(first_edge, last_edge, length.out = bins + 1L)
  # keep only values within [first_edge, last_edge]
  keep <- a >= first_edge & a <= last_edge
  tmp <- a[keep]
  counts <- integer(bins)
  if (length(tmp) > 0L) {
    norm <- bins / (last_edge - first_edge)
    f_indices <- (tmp - first_edge) * norm
    indices <- as.integer(floor(f_indices))      # 0-based bin indices
    # values exactly at last_edge land in index == bins -> clamp to bins-1
    indices[indices == bins] <- bins - 1L
    # NumPy rounding corrections against the actual edges
    # decrement where value < left edge of its bin
    dec <- tmp < bin_edges[indices + 1L]          # bin_edges is 1-based
    indices[dec] <- indices[dec] - 1L
    # increment where value >= right edge of its bin (except last bin)
    inc <- (tmp >= bin_edges[indices + 2L]) & (indices != (bins - 1L))
    inc[is.na(inc)] <- FALSE
    indices[inc] <- indices[inc] + 1L
    tab <- tabulate(indices + 1L, nbins = bins)
    counts <- tab
  }
  list(counts = counts, edges = bin_edges)
}

# build_binned(): histogram ell into `bins` equal bins over [ell_min, ell_max]
# (defaulting to the data range), returning bin centers, counts, and edges.
build_binned <- function(lines, bins, ell_min = NULL, ell_max = NULL) {
  ell <- as.numeric(lines$ell)
  if (is.null(ell_min)) ell_min <- min(ell)
  if (is.null(ell_max)) ell_max <- max(ell)
  h <- np_histogram_uniform(ell, bins, ell_min, ell_max)
  edges <- h$edges
  centers <- 0.5 * (edges[-length(edges)] + edges[-1L])
  data.frame(
    ell = centers,
    count = as.numeric(h$counts),
    edge_lo = edges[-length(edges)],
    edge_hi = edges[-1L]
  )
}


# --- Gaussian smoothing ----------------------------------------------------

# gaussian_kernel1d(): build the order-0 Gaussian weights exactly as SciPy's
# _gaussian_kernel1d: phi_x = exp(-0.5 * x^2 / sigma^2) over x in [-r, r],
# normalised to sum 1. radius r = floor(truncate * sigma + 0.5).
gaussian_kernel1d <- function(sigma, truncate = 4.0) {
  lw <- as.integer(truncate * sigma + 0.5)
  x <- seq.int(-lw, lw)
  sigma2 <- sigma * sigma
  phi <- exp(-0.5 / sigma2 * (x * x))
  phi / sum(phi)
}

# gaussian_filter_nearest(): replicate scipy.ndimage.gaussian_filter1d with
# mode="nearest". The order-0 kernel is symmetric, so correlation equals
# convolution; edges are extended by replicating the nearest sample.
gaussian_filter_nearest <- function(y, sigma, truncate = 4.0) {
  y <- as.numeric(y)
  n <- length(y)
  w <- gaussian_kernel1d(sigma, truncate)
  lw <- (length(w) - 1L) %/% 2L
  out <- numeric(n)
  for (i in seq_len(n)) {
    acc <- 0.0
    for (j in seq.int(-lw, lw)) {
      idx <- i + j
      if (idx < 1L) idx <- 1L
      if (idx > n) idx <- n
      acc <- acc + w[j + lw + 1L] * y[idx]
    }
    out[i] <- acc
  }
  out
}


# --- Poisson model ---------------------------------------------------------

# poisson_deviance(): D(y || mu) = 2 sum(mu - y + y log(y/mu)), with zero-count
# terms (y == 0) contributing only mu, and mu floored at EPS.
poisson_deviance <- function(y, mu) {
  y <- as.numeric(y)
  mu <- pmax(as.numeric(mu), EPS)
  term <- mu - y
  nz <- y > 0
  term[nz] <- term[nz] + y[nz] * log(y[nz] / mu[nz])
  2.0 * sum(term)
}

# lstsq_minnorm(): minimum-norm least-squares solution via SVD, used as the
# singular-system fallback (mirrors numpy.linalg.lstsq with rcond=None).
lstsq_minnorm <- function(A, b) {
  sv <- svd(A)
  d <- sv$d
  tol <- max(dim(A)) * max(d) * .Machine$double.eps
  dinv <- ifelse(d > tol, 1.0 / d, 0.0)
  sv$v %*% (dinv * (t(sv$u) %*% b))
}

# design_poly(): centered/standardized polynomial design matrix of the given
# degree. z = (ell - mean(ell)); divide by population std (ddof = 0) when > 0;
# columns are z^0 .. z^degree.
design_poly <- function(ell, degree) {
  z <- ell - mean(ell)
  s <- sqrt(mean(z * z))   # population std, matching numpy.std default
  if (s > 0) z <- z / s
  do.call(cbind, lapply(0:degree, function(d) z^d))
}

# fit_poisson_loglinear(): IRLS for log(mu) = log(B) + X beta, replicating the
# Python reference step-for-step:
#   beta init: beta[1] = log(max(sum(y),EPS)/max(sum(B),EPS)); others 0
#   eta clipped to [-10, 10]; mu = max(B exp(eta), EPS)
#   working response z = eta + (y - mu)/mu; weights W = mu
#   A = X' (W X) + ridge I; b = X' (W z); solve, lstsq fallback
#   converge when max|beta2 - beta| < 1e-8; max 60 iterations
# Returns list(beta, deviance, mu).
fit_poisson_loglinear <- function(y, baseline, X, ridge = 1e-8, max_iter = 60L) {
  y <- as.numeric(y)
  B <- pmax(as.numeric(baseline), EPS)
  X <- as.matrix(X)
  p <- ncol(X)
  beta <- numeric(p)
  beta[1] <- log(max(sum(y), EPS) / max(sum(B), EPS))
  I <- diag(p)
  for (iter in seq_len(max_iter)) {
    eta <- pmin(pmax(as.numeric(X %*% beta), -10), 10)
    mu <- pmax(B * exp(eta), EPS)
    z <- eta + (y - mu) / mu
    W <- mu
    A <- t(X) %*% (W * X) + ridge * I
    bvec <- t(X) %*% (W * z)
    beta2 <- tryCatch(
      solve(A, bvec),
      error = function(e) lstsq_minnorm(A, bvec)
    )
    beta2 <- as.numeric(beta2)
    if (max(abs(beta2 - beta)) < 1e-8) {
      beta <- beta2
      break
    }
    beta <- beta2
  }
  eta <- pmin(pmax(as.numeric(X %*% beta), -10), 10)
  mu <- pmax(B * exp(eta), EPS)
  list(beta = beta, deviance = poisson_deviance(y, mu), mu = mu)
}


# --- scan ------------------------------------------------------------------

# scan_k(): for each k in k_grid fit the base + cos/sin harmonic model and
# record deltaD = D_base - D_harmonic. The best k is the first maximiser of
# deltaD (strict >), matching the Python reference. Returns list(scan, best,
# mu0).
scan_k <- function(ell, y, baseline, k_grid, degree) {
  X0 <- design_poly(ell, degree)
  base_fit <- fit_poisson_loglinear(y, baseline, X0)
  D0 <- base_fit$deviance
  mu0 <- base_fit$mu

  n <- length(k_grid)
  k_v <- numeric(n); dd_v <- numeric(n)
  amp_v <- numeric(n); ph_v <- numeric(n); Dh_v <- numeric(n)
  best <- NULL
  for (i in seq_len(n)) {
    k <- k_grid[i]
    X <- cbind(X0, cos(k * ell), sin(k * ell))
    fit <- fit_poisson_loglinear(y, baseline, X)
    D <- fit$deviance
    delta <- D0 - D
    b <- fit$beta
    a_coef <- b[length(b) - 1L]   # cos coefficient
    b_coef <- b[length(b)]        # sin coefficient
    amp <- sqrt(a_coef^2 + b_coef^2)
    phase <- atan2(-b_coef, a_coef)
    k_v[i] <- k; dd_v[i] <- delta; Dh_v[i] <- D
    amp_v[i] <- amp; ph_v[i] <- phase
    if (is.null(best) || delta > best$deltaD) {
      best <- list(k_best = k, deltaD = delta, D_base = D0,
                   D_harmonic = D, amplitude = amp, phase = phase)
    }
  }
  scan <- data.frame(
    k = k_v, deltaD = dd_v, D_base = D0,
    D_harmonic = Dh_v, amplitude = amp_v, phase = ph_v
  )
  list(scan = scan, best = best, mu0 = mu0)
}

# null_scan(): parametric Poisson bootstrap. For each replicate draw
# y_null ~ Poisson(mu0), run the full k-grid scan, and store the maximum
# deltaD. The empirical p-value applies the standard +1 correction.
#
# NOTE: R's rpois and NumPy's Generator.poisson use different algorithms, so
# individual null draws will NOT match Python for a given seed. Repeated R
# runs with the same seed DO match each other. See R/README.md.
null_scan <- function(ell, y, baseline, mu0, k_grid, degree, real_delta,
                      null_n, seed, verbose = TRUE) {
  set.seed(seed)
  vals <- numeric(null_n)
  for (i in seq_len(null_n)) {
    y0 <- rpois(length(mu0), mu0)
    res <- scan_k(ell, y0, baseline, k_grid, degree)
    vals[i] <- res$best$deltaD
    if (verbose && (i %% 100L == 0L || i == null_n)) {
      cat(sprintf("[null] %d/%d\n", i, null_n))
    }
  }
  tail_count <- sum(vals >= real_delta)
  p <- (1 + tail_count) / (1 + length(vals))
  list(vals = vals, p = p, tail = as.integer(tail_count))
}


# --- branch report ---------------------------------------------------------

# branch_report(): compute n_obs = k_best * delta_ell / (2 pi) and the
# nearest target windings across the three branches, sorted by |k error|,
# top 10. Returns list(n_obs, nearest = data.frame).
branch_report <- function(k_best, delta_ell) {
  n_obs <- k_best * delta_ell / (2.0 * pi)
  branches <- list(
    list(label = "koide_10_15_20", ns = WCT_REFERENCE_NS),
    list(label = "folded_4over9", ns = WCT_FOLDED_NS),
    list(label = "integer_1_to_40", ns = as.numeric(1:40))
  )
  rows <- list()
  for (br in branches) {
    for (n in br$ns) {
      k_t <- 2.0 * pi * n / delta_ell
      rows[[length(rows) + 1L]] <- data.frame(
        branch = br$label,
        n = n,
        k_target = k_t,
        k_error = k_best - k_t,
        abs_k_error = abs(k_best - k_t),
        n_error = n_obs - n,
        abs_n_error = abs(n_obs - n),
        stringsAsFactors = FALSE
      )
    }
  }
  allr <- do.call(rbind, rows)
  allr <- allr[order(allr$abs_k_error), , drop = FALSE]
  nearest <- utils::head(allr, 10L)
  rownames(nearest) <- NULL
  list(n_obs = n_obs, nearest = nearest)
}


# --- peak detection & ratios -----------------------------------------------

# local_maxima_1d(): replicate scipy.signal._local_maxima_1d, returning 0-based
# midpoint indices of local maxima (with plateau midpoints).
local_maxima_1d <- function(x) {
  n <- length(x)
  mids <- integer(0)
  i <- 2L            # 1-based; corresponds to python i=1
  i_max <- n - 1L    # python i_max = len-1 (0-based exclusive) -> stop before last
  # python loops while i < i_max with 0-based i. Translate to 1-based:
  # 0-based i in [1, i_max-1]; 1-based in [2, i_max].
  while (i <= i_max) {
    if (x[i - 1L] < x[i]) {
      i_ahead <- i + 1L
      while (i_ahead <= i_max && x[i_ahead] == x[i]) {
        i_ahead <- i_ahead + 1L
      }
      if (x[i_ahead] < x[i]) {
        left_edge <- i
        right_edge <- i_ahead - 1L
        midpoint <- (left_edge + right_edge) %/% 2L
        mids <- c(mids, midpoint - 1L)   # store 0-based
        i <- i_ahead
      }
    }
    i <- i + 1L
  }
  mids
}

# select_by_peak_distance(): replicate scipy.signal._select_by_peak_distance.
# peaks are 0-based positions, priority is the height used to break ties.
select_by_peak_distance <- function(peaks, priority, distance) {
  np <- length(peaks)
  keep <- rep(TRUE, np)
  distance_ <- ceiling(distance)
  order_idx <- order(priority)          # ascending priority -> 1-based
  for (i in seq.int(np, 1L)) {          # highest priority first
    j <- order_idx[i]
    if (!keep[j]) next
    k <- j - 1L
    while (k >= 1L && (peaks[j] - peaks[k]) < distance_) {
      keep[k] <- FALSE
      k <- k - 1L
    }
    k <- j + 1L
    while (k <= np && (peaks[k] - peaks[j]) < distance_) {
      keep[k] <- FALSE
      k <- k + 1L
    }
  }
  keep
}

# find_peaks_distance(): scipy.signal.find_peaks(x, distance=distance) for the
# distance-only case. Returns 0-based peak indices.
find_peaks_distance <- function(x, distance) {
  peaks <- local_maxima_1d(x)
  if (length(peaks) == 0L) return(integer(0))
  if (!is.null(distance)) {
    keep <- select_by_peak_distance(peaks, x[peaks + 1L], distance)
    peaks <- peaks[keep]
  }
  peaks
}

# limit_denominator(): exact reproduction of Python Fraction.limit_denominator
# using gmp big rationals. `frac` is a gmp bigq equal to the double's exact
# value; returns a bigq with denominator <= max_denominator.
limit_denominator <- function(frac, max_denominator = 64L) {
  num <- gmp::numerator(frac)
  den <- gmp::denominator(frac)
  if (den <= max_denominator) {
    return(frac)
  }
  p0 <- gmp::as.bigz(0); q0 <- gmp::as.bigz(1)
  p1 <- gmp::as.bigz(1); q1 <- gmp::as.bigz(0)
  n <- num; d <- den
  mdz <- gmp::as.bigz(max_denominator)
  repeat {
    a <- n %/% d
    q2 <- q0 + a * q1
    if (q2 > mdz) break
    tmp_p1 <- p0 + a * p1
    p0 <- p1; q0 <- q1
    p1 <- tmp_p1; q1 <- q2
    new_d <- n - a * d
    n <- d; d <- new_d
  }
  k <- (mdz - q0) %/% q1
  bound1 <- gmp::as.bigq(p0 + k * p1, q0 + k * q1)
  bound2 <- gmp::as.bigq(p1, q1)
  if (abs(bound2 - frac) <= abs(bound1 - frac)) {
    bound2
  } else {
    bound1
  }
}

# double_to_bigq(): exact gmp bigq representation of an R double, matching
# Python's Fraction(float) which is exact.
double_to_bigq <- function(x) {
  if (x == 0) return(gmp::as.bigq(0L, 1L))
  sign <- if (x < 0) -1L else 1L
  ax <- abs(x)
  # decompose ax = mant * 2^exp with mant in [0.5, 1); 53-bit mantissa
  e <- floor(log2(ax))
  # bring mantissa to an integer by scaling with 2^53
  # m = round(ax / 2^e * 2^52) gives 53-bit integer mantissa
  # Use a robust approach: find power of two P such that ax * 2^P is integer.
  P <- 0L
  m <- ax
  # increase P until m is integral (bounded by ~1100 for doubles)
  while (m != floor(m) && P < 1100L) {
    m <- m * 2
    P <- P + 1L
  }
  # m is now an exact integer value; represent via gmp from string to avoid
  # precision loss for large magnitudes
  mant_big <- gmp::as.bigz(formatC(m, format = "f", digits = 0))
  num <- sign * mant_big
  den <- gmp::as.bigz(2L)^P
  gmp::as.bigq(num, den)
}

# peak_ratios(): detect scan-curve peaks, take the top peaks by deltaD, and
# build the pairwise k-ratio table with rational approximations (denominator
# <= 64). Mirrors the Python reference output schema.
peak_ratios <- function(scan, top_n = 12L) {
  y <- as.numeric(scan$deltaD)
  distance <- max(1L, length(y) %/% 80L)
  peaks <- find_peaks_distance(y, distance)   # 0-based
  if (length(peaks) == 0L) {
    peaks <- which.max(y) - 1L
  }
  peak_rows <- peaks + 1L                      # 1-based into scan
  pk <- scan[peak_rows, , drop = FALSE]
  pk <- pk[order(-pk$deltaD), , drop = FALSE]
  pk <- utils::head(pk, top_n)
  kk <- as.numeric(pk$k)
  m <- length(kk)
  if (m < 2L) return(data.frame())
  rows <- list()
  for (i in seq_len(m - 1L)) {
    for (j in seq.int(i + 1L, m)) {
      pair <- sort(c(kk[i], kk[j]))
      lo <- pair[1]; hi <- pair[2]
      r <- hi / lo
      fr <- limit_denominator(double_to_bigq(r), 64L)
      num <- as.numeric(as.character(gmp::numerator(fr)))
      den <- as.numeric(as.character(gmp::denominator(fr)))
      rows[[length(rows) + 1L]] <- data.frame(
        k_low = lo, k_high = hi, ratio = r,
        rational = sprintf("%s/%s",
                           as.character(gmp::numerator(fr)),
                           as.character(gmp::denominator(fr))),
        rational_error = abs(r - num / den),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$rational_error), , drop = FALSE]
  rownames(out) <- NULL
  out
}
