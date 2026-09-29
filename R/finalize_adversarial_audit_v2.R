#!/usr/bin/env Rscript
# finalize_adversarial_audit_v2.R
# ---------------------------------------------------------------------------
# Corrected acceptance-layer wrapper around finalize_adversarial_audit.R.
#
# The underlying finalizer contains the complete report/freeze implementation.
# This wrapper fixes three acceptance-layer defects without changing any
# scientific statistic or decision threshold:
#   1. avoid shadowing base::c() with the calibration table;
#   2. treat NA/blank source metadata as incomplete provenance;
#   3. make the machine-result hash manifest path handling platform-neutral.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)

# Source the complete finalizer into this environment. Its direct-invocation
# guard does not fire when this wrapper is called from the v2 audit runner.
sys.source(file.path(R_DIR, "finalize_adversarial_audit.R"), envir = environment())
.legacy_finalize_main <- main

resolved_provenance <- function(root, contract_prov) {
  p <- file.path(root, "config/nist_source_provenance.csv")
  if (!file.exists(p)) return(contract_prov)
  cfg <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("species", "source_provider", "retrieval_date", "retrieval_version",
            "retrieval_url_or_query", "notes")
  if (!all(need %in% names(cfg))) return(contract_prov)

  out <- contract_prov
  for (i in seq_len(nrow(out))) {
    hit <- cfg[cfg$species == out$species[i], , drop = FALSE]
    if (nrow(hit) == 0L) next
    h <- hit[1, , drop = FALSE]
    out$source_provider[i] <- h$source_provider[1]
    out$source_retrieval_date[i] <- h$retrieval_date[1]
    out$source_retrieval_version[i] <- h$retrieval_version[1]
    vals <- as.character(c(h$source_provider[1], h$retrieval_date[1],
                           h$retrieval_version[1], h$retrieval_url_or_query[1]))
    out$source_retrieval_metadata_complete[i] <-
      all(!is.na(vals) & nzchar(trimws(vals)))
    out$note[i] <- h$notes[1]
  }
  out
}

species_gate <- function(species, obs, bs, ns, cal, inj, art, influence,
                         expected_k, expected_lines = NA_integer_) {
  o <- obs[obs$species == species, , drop = FALSE]
  b <- bs[bs$species == species, , drop = FALSE]
  n <- ns[ns$species == species, , drop = FALSE]
  ctab <- cal[cal$species == species, , drop = FALSE]
  j <- inj[inj$species == species, , drop = FALSE]
  a <- art[art$species == species, , drop = FALSE]
  po <- primary_row(obs, species)

  has_bins <- identical(sort(unique(as.integer(o$bins))), c(120L, 160L, 200L))
  line_ok <- if (is.finite(expected_lines)) all(o$n_lines == expected_lines) else TRUE
  reproduce_k <- is.finite(po$k_best[1]) && abs(po$k_best[1] - expected_k) <= 1e-9
  reproduce <- has_bins && line_ok && reproduce_k

  bootstrap_ok <- nrow(b) == 1L &&
    is.finite(b$candidate_retention_2pct[1]) &&
    b$candidate_retention_2pct[1] >= 0.80 &&
    is.finite(b$failure_rate[1]) && b$failure_rate[1] <= 0.01

  null_ok <- nrow(n) == 3L && all(is.finite(n$global_p)) && all(n$global_p <= 0.05)

  required_alpha <- base::c(0.10, 0.05, 0.01, 0.001)
  alpha_ok <- all(vapply(required_alpha, function(x)
    any(abs(ctab$nominal_alpha - x) < 1e-12 & truthy(ctab$compatible)), logical(1)))

  j0 <- j[abs(j$amplitude_ratio) < 1e-12, , drop = FALSE]
  j1 <- j[abs(j$amplitude_ratio - 1) < 1e-12, , drop = FALSE]
  injection_ok <- nrow(j0) == 1L && nrow(j1) == 1L &&
    is.finite(j0$detection_prob[1]) && j0$detection_prob[1] <= 0.05 &&
    is.finite(j1$sig_correct_prob[1]) && j1$sig_correct_prob[1] >= 0.80

  bin_ok <- nrow(o) == 3L && all(truthy(o$compatible_2pct))
  implementation_ok <- nrow(a) == 1L && truthy(a$cleaning_reconciles[1]) &&
    truthy(a$sorting_invariant[1])
  influence_ok <- nrow(a) == 1L && is.finite(a$influence_reference_fraction[1]) &&
    a$influence_reference_fraction[1] >= 0.80

  list(
    reproduce = reproduce,
    bootstrap = bootstrap_ok,
    canonical_null = null_ok,
    calibration = alpha_ok,
    injection = injection_ok,
    cross_binning = bin_ok,
    implementation = implementation_ok,
    influence = influence_ok,
    primary = po,
    fixed = primary_row(ns, species),
    details = list(has_bins = has_bins, line_ok = line_ok, reproduce_k = reproduce_k)
  )
}

hash_manifest <- function(root, json_path) {
  td <- file.path(root, "tables_r/statistical_audit")
  key <- c(
    "data/Fe_lines.csv", "data/Co_lines.csv",
    "config/final_adversarial_audit.json", "config/nist_source_provenance.csv",
    "R/nist_scan_lib.R", "R/audit_utils.R", "R/render_statistical_audit.R",
    "R/run_final_adversarial_audit.R", "R/run_final_adversarial_audit_v2.R",
    "R/run_final_contract_checks.R", "R/build_final_adversarial_summary.R",
    "R/finalize_adversarial_audit.R", "R/finalize_adversarial_audit_v2.R",
    "tables_r/statistical_audit/final_species_observed_results.csv",
    "tables_r/statistical_audit/final_species_bootstrap_summary.csv",
    "tables_r/statistical_audit/final_species_null_summary.csv",
    "tables_r/statistical_audit/final_species_calibration.csv",
    "tables_r/statistical_audit/final_species_injection_recovery.csv",
    "tables_r/statistical_audit/multiple_testing.csv",
    "tables_r/statistical_audit/alternative_null_results.csv",
    "tables_r/statistical_audit/final_adversarial_run_manifest.csv"
  )
  if (!is.null(json_path) && file.exists(json_path)) key <- c(key, basename(json_path))
  key <- unique(key[file.exists(file.path(root, key))])
  out <- data.frame(
    path = key,
    bytes = vapply(file.path(root, key), function(p) file.info(p)$size, numeric(1)),
    sha256 = vapply(file.path(root, key), sha256_file, character(1)),
    stringsAsFactors = FALSE)
  write_table(out, file.path(td, "provenance_hash_manifest.csv"))
  out
}

main <- function(argv = commandArgs(TRUE)) .legacy_finalize_main(argv)

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("finalize_adversarial_audit_v2\\.R$", .invoked_file)) main()
