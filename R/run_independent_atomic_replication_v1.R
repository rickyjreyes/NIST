#!/usr/bin/env Rscript
# run_independent_atomic_replication_v1.R
# ---------------------------------------------------------------------------
# Prospective independent-catalog replication runner.
#
# The primary target is frozen in config/independent_atomic_replication_v1.json
# before the external outcome is inspected. This script never retunes k for the
# primary classification. An exploratory full scan is generated only after the
# fixed-target primary result has been computed.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
source(file.path(R_DIR, "audit_utils.R"))

.prov <- new.env(parent = environment())
sys.source(file.path(R_DIR, "verify_source_provenance.R"), envir = .prov)

load_protocol <- function(root, simplify = TRUE) {
  p <- file.path(root, "config", "independent_atomic_replication_v1.json")
  if (!file.exists(p)) stop("missing frozen replication protocol: ", p)
  jsonlite::fromJSON(p, simplifyVector = simplify)
}

protocol_path <- function(root) file.path(root, "config", "independent_atomic_replication_v1.json")

pin_source_hash <- function(root) {
  p <- protocol_path(root)
  protocol <- load_protocol(root, simplify = FALSE)
  rel <- protocol$source$expected_local_path
  src <- file.path(root, rel)
  if (!file.exists(src)) stop("external source file not found: ", src)
  h <- .prov$sha256_file_exact(src)
  protocol$source$raw_sha256 <- h
  jsonlite::write_json(protocol, p, pretty = TRUE, auto_unbox = TRUE, na = "null")
  cat("Pinned external source SHA-256 without parsing scientific content:\n")
  cat(rel, "  ", h, "\n", sep = "")
  cat("Commit the protocol hash update before running the primary replication.\n")
  invisible(h)
}

verify_source_hash <- function(root, protocol) {
  rel <- as.character(protocol$source$expected_local_path)
  p <- file.path(root, rel)
  if (!file.exists(p)) stop("external source file not found: ", p)
  expected <- tolower(trimws(as.character(protocol$source$raw_sha256)))
  if (!.prov$is_hex_sha256(expected))
    stop("external source SHA-256 is not pinned; run --pin-source-hash, inspect, and commit before outcome access")
  actual <- .prov$sha256_file_exact(p)
  if (!identical(expected, actual)) stop("external source SHA-256 mismatch")
  list(path = p, relative_path = rel, sha256 = actual)
}

parse_num <- function(x) suppressWarnings(as.numeric(trimws(x)))

# Kurucz CD-ROM 23 / GFALL fixed-width fields used here:
#  1-11   wavelength
# 12-18   log(gf)
# 19-24   element/ion code (Z + charge/100; Fe II = 26.01)
# 25-36   first level energy (cm^-1)
# 37-41   first J
# 42      spacer
# 43-52   first level label
# 53-64   second level energy (cm^-1)
# 65-69   second J
# Only element code and the two level energies are used in the primary parser.
parse_gfall_feii <- function(path, protocol, chunk_n = 100000L) {
  code_target <- as.numeric(protocol$selection$kurucz_element_code)
  wmin <- as.numeric(protocol$selection[["wavenumber_min_cm-1"]])
  wmax <- as.numeric(protocol$selection[["wavenumber_max_cm-1"]])
  require_nonneg_lo <- isTRUE(protocol$selection$require_nonnegative_lower_energy)
  require_nonneg_hi <- isTRUE(protocol$selection$require_nonnegative_upper_energy)

  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  chunks <- list(); j <- 1L
  n_raw <- 0L; n_code <- 0L; n_energy_valid <- 0L; n_in_range <- 0L

  repeat {
    lines <- readLines(con, n = chunk_n, warn = FALSE)
    if (!length(lines)) break
    n_raw <- n_raw + length(lines)
    long_enough <- nchar(lines, type = "bytes") >= 64L
    if (!any(long_enough)) next
    z <- lines[long_enough]

    code <- parse_num(substr(z, 19L, 24L))
    hit <- is.finite(code) & abs(code - code_target) < 5e-4
    n_code <- n_code + sum(hit)
    if (!any(hit)) next
    z <- z[hit]

    e1 <- parse_num(substr(z, 25L, 36L))
    e2 <- parse_num(substr(z, 53L, 64L))
    valid <- is.finite(e1) & is.finite(e2)
    if (require_nonneg_lo) valid <- valid & e1 >= 0
    if (require_nonneg_hi) valid <- valid & e2 >= 0
    n_energy_valid <- n_energy_valid + sum(valid)
    if (!any(valid)) next
    e1 <- e1[valid]; e2 <- e2[valid]

    wn <- abs(e2 - e1)
    keep <- is.finite(wn) & wn > 0 & wn >= wmin & wn <= wmax
    n_in_range <- n_in_range + sum(keep)
    if (!any(keep)) next
    wn <- wn[keep]
    chunks[[j]] <- data.frame(wavenumber = wn, ell = log(wn), stringsAsFactors = FALSE)
    j <- j + 1L
  }

  selected <- if (length(chunks)) do.call(rbind, chunks) else
    data.frame(wavenumber = numeric(), ell = numeric())
  if (nrow(selected)) selected <- selected[order(selected$wavenumber), , drop = FALSE]
  dup <- duplicated(selected$wavenumber)
  n_duplicates <- sum(dup)
  selected <- selected[!dup, , drop = FALSE]
  rownames(selected) <- NULL

  summary <- data.frame(
    source_rows = n_raw,
    feii_code_rows = n_code,
    nonnegative_finite_energy_rows = n_energy_valid,
    in_frozen_wavenumber_range_rows = n_in_range,
    exact_duplicate_wavenumbers_removed = n_duplicates,
    retained_lines = nrow(selected),
    wavenumber_min = if (nrow(selected)) min(selected$wavenumber) else NA_real_,
    wavenumber_max = if (nrow(selected)) max(selected$wavenumber) else NA_real_,
    stringsAsFactors = FALSE)
  list(lines = selected, summary = summary)
}

fixed_target_analysis <- function(lines, protocol) {
  target <- as.numeric(protocol$primary_test$target_k)
  bins <- as.integer(protocol$primary_test$bins)
  sigma <- as.numeric(protocol$primary_test$baseline_sigma_bins)
  degree <- as.integer(protocol$primary_test$degree)
  # A one-element k grid ensures the primary calculation cannot look elsewhere.
  run_scan_analysis(lines, bins, c(target), degree, sigma)
}

fixed_target_null <- function(res, protocol, parallel = FALSE) {
  B <- as.integer(protocol$primary_test$null_replicates)
  seed <- as.integer(protocol$primary_test$seed)
  target <- as.numeric(protocol$primary_test$target_k)
  sigma <- as.numeric(protocol$primary_test$baseline_sigma_bins)
  degree <- as.integer(protocol$primary_test$degree)
  obs <- res$best$deltaD

  setup_rng(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(i) {
    set.seed(seeds[i])
    y0 <- stats::rpois(length(res$mu0), res$mu0)
    baseline0 <- pmax(gaussian_filter_nearest(y0, sigma), EPS)
    scan_k(res$ell, y0, baseline0, c(target), degree)$best$deltaD
  }
  sims <- unlist(audit_lapply(seq_len(B), one, parallel = parallel, seed = seed),
                 use.names = FALSE)
  tail <- sum(sims >= obs)
  sdev <- stats::sd(sims)
  data.frame(
    observed_deltaD = obs,
    null_n = B,
    exceedances = as.integer(tail),
    empirical_p = emp_p(tail, B),
    null_mean = mean(sims),
    null_sd = sdev,
    null_z = if (is.finite(sdev) && sdev > 0) (obs - mean(sims)) / sdev else NA_real_,
    null_q95 = unname(stats::quantile(sims, 0.95, type = 8)),
    stringsAsFactors = FALSE) -> summary
  list(summary = summary, draws = sims)
}

exploratory_scan <- function(lines, protocol) {
  ex <- protocol$exploratory_scan
  k_grid <- seq(as.numeric(ex$k_min), as.numeric(ex$k_max), length.out = as.integer(ex$n_k))
  run_scan_analysis(lines,
                    as.integer(protocol$primary_test$bins),
                    k_grid,
                    as.integer(protocol$primary_test$degree),
                    as.numeric(protocol$primary_test$baseline_sigma_bins))
}

classify_primary <- function(selection, primary, protocol) {
  minimum <- as.integer(protocol$selection$minimum_retained_lines_for_primary_test)
  if (selection$retained_lines[1] < minimum) return("INCONCLUSIVE")
  p <- primary$empirical_p[1]
  if (!is.finite(p)) return("INCONCLUSIVE")
  if (p <= as.numeric(protocol$primary_test$alpha)) "PASS" else "FAIL"
}

git_commit <- function(root) {
  z <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE,
                        stderr = FALSE), error = function(e) character())
  if (length(z)) trimws(z[1]) else NA_character_
}

main <- function(argv = commandArgs(TRUE)) {
  root <- audit_repo_root()
  if ("--pin-source-hash" %in% argv) return(pin_source_hash(root))

  protocol <- load_protocol(root, simplify = TRUE)
  if (!identical(protocol$status, "FROZEN_BEFORE_EXTERNAL_OUTCOME_ACCESS"))
    stop("replication protocol is not marked frozen")
  parallel <- any(argv == "--parallel")
  if (parallel) configure_parallel(TRUE, "auto")

  source_info <- verify_source_hash(root, protocol)
  parsed <- parse_gfall_feii(source_info$path, protocol)

  td <- file.path(root, "tables_r", "statistical_audit")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  write_table(parsed$summary, file.path(td, "independent_atomic_replication_v1_selection.csv"))

  minimum <- as.integer(protocol$selection$minimum_retained_lines_for_primary_test)
  if (parsed$summary$retained_lines[1] < minimum) {
    verdict <- "INCONCLUSIVE"
    machine <- list(
      schema_version = "independent_atomic_replication_v1_result",
      protocol = "config/independent_atomic_replication_v1.json",
      git_commit = git_commit(root),
      source_sha256 = source_info$sha256,
      selection = as.list(parsed$summary[1, , drop = FALSE]),
      classification = verdict,
      reason = paste0("retained lines below frozen minimum of ", minimum))
    jsonlite::write_json(machine, file.path(root, "independent_atomic_replication_v1_result.json"),
                         pretty = TRUE, auto_unbox = TRUE, digits = NA)
    cat("Independent catalog replication: INCONCLUSIVE (insufficient retained lines)\n")
    return(invisible(machine))
  }

  # PRIMARY OUTCOME: fixed target only. Classification is fixed before any full scan.
  res_fixed <- fixed_target_analysis(parsed$lines, protocol)
  null <- fixed_target_null(res_fixed, protocol, parallel)
  primary <- cbind(
    data.frame(
      target_k = as.numeric(protocol$primary_test$target_k),
      retained_lines = nrow(parsed$lines),
      bins = as.integer(protocol$primary_test$bins),
      baseline_sigma_bins = as.numeric(protocol$primary_test$baseline_sigma_bins),
      degree = as.integer(protocol$primary_test$degree),
      amplitude = res_fixed$best$amplitude,
      phase = res_fixed$best$phase,
      stringsAsFactors = FALSE),
    null$summary)
  verdict <- classify_primary(parsed$summary, primary, protocol)

  write_table(primary, file.path(td, "independent_atomic_replication_v1_primary.csv"))
  write_table(data.frame(null_fixed_deltaD = null$draws),
              file.path(td, "independent_atomic_replication_v1_null_draws.csv"))

  # Only now is a full exploratory scan allowed.
  ex <- exploratory_scan(parsed$lines, protocol)
  ex_table <- ex$scan
  write_table(ex_table, file.path(td, "independent_atomic_replication_v1_exploratory_scan.csv"))
  ex_summary <- data.frame(
    k_best = ex$best$k_best,
    deltaD_best = ex$best$deltaD,
    amplitude_best = ex$best$amplitude,
    phase_best = ex$best$phase,
    target_relative_offset = abs(ex$best$k_best - protocol$primary_test$target_k) /
      protocol$primary_test$target_k,
    primary_classification_unchanged = TRUE,
    stringsAsFactors = FALSE)
  write_table(ex_summary,
              file.path(td, "independent_atomic_replication_v1_exploratory_summary.csv"))

  machine <- list(
    schema_version = "independent_atomic_replication_v1_result",
    protocol = "config/independent_atomic_replication_v1.json",
    git_commit = git_commit(root),
    source = list(
      catalog = protocol$source$name,
      input = source_info$relative_path,
      sha256 = source_info$sha256),
    selection = as.list(parsed$summary[1, , drop = FALSE]),
    primary = as.list(primary[1, , drop = FALSE]),
    classification = verdict,
    exploratory = as.list(ex_summary[1, , drop = FALSE]),
    claim_boundary = paste(
      "catalog-level replication only; does not by itself constitute an independent laboratory",
      "experiment or establish WCT as the causal mechanism"))
  jsonlite::write_json(machine, file.path(root, "independent_atomic_replication_v1_result.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = NA)

  cat("\nINDEPENDENT ATOMIC CATALOG REPLICATION V1\n")
  cat("source SHA-256: ", source_info$sha256, "\n", sep = "")
  print(parsed$summary, row.names = FALSE)
  print(primary, row.names = FALSE)
  cat("classification: ", verdict, "\n", sep = "")
  cat("exploratory best k (cannot change classification): ", ex$best$k_best, "\n", sep = "")
  invisible(machine)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("run_independent_atomic_replication_v1\\.R$", .invoked_file)) main()
