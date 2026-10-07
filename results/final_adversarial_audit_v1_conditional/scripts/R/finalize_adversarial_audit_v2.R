#!/usr/bin/env Rscript
# finalize_adversarial_audit_v2.R
# ---------------------------------------------------------------------------
# Corrected acceptance-layer wrapper around finalize_adversarial_audit.R.
#
# This wrapper does not change any detector, statistic, threshold, null, or
# Monte-Carlo result. It strengthens the release/provenance layer by:
#   1. avoiding the historical calibration-table/base::c() name collision;
#   2. distinguishing a complete provenance reconstruction record from exact
#      historical retrieval/query recovery;
#   3. requiring real SHA-256 (never an MD5 value in a sha256 column);
#   4. making the hash-manifest path handling platform-neutral; and
#   5. correcting final report wording when historical fields are explicitly
#      documented as UNRECOVERABLE rather than merely left blank.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)

# Source the complete finalizer into this environment. Its direct-invocation
# guard does not fire when this wrapper is called from the v2 audit runner.
sys.source(file.path(R_DIR, "finalize_adversarial_audit.R"), envir = environment())
.legacy_finalize_main <- main

# Reuse the exact SHA-256 implementation from the standalone provenance checker
# without importing its main() into this environment.
.prov_env <- new.env(parent = environment())
sys.source(file.path(R_DIR, "verify_source_provenance.R"), envir = .prov_env)
sha256_file <- .prov_env$sha256_file_exact

.nonempty <- function(x) !is.na(x) && nzchar(trimws(as.character(x)))
.status <- function(x) if (.nonempty(x)) toupper(trimws(as.character(x))) else "MISSING"

.read_provenance_config <- function(root) {
  p <- file.path(root, "config/nist_source_provenance.csv")
  if (!file.exists(p)) return(NULL)
  cfg <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE,
                         fileEncoding = "UTF-8")
  if (length(names(cfg))) names(cfg)[1] <- sub("^\\ufeff", "", names(cfg)[1])
  cfg
}

resolved_provenance <- function(root, contract_prov) {
  cfg <- .read_provenance_config(root)
  if (is.null(cfg)) return(contract_prov)

  base_need <- c("species", "source_provider", "retrieval_date",
                 "retrieval_version", "retrieval_url_or_query", "notes")
  if (!all(base_need %in% names(cfg))) return(contract_prov)

  ext_need <- c("retrieval_date_status", "query_export_settings_status",
                "provenance_evidence", "input_file", "input_sha256")
  extended <- all(ext_need %in% names(cfg))

  out <- contract_prov
  extra_cols <- list(
    source_retrieval_url_or_query = rep(NA_character_, nrow(out)),
    source_retrieval_date_status = rep("MISSING", nrow(out)),
    source_query_export_settings_status = rep("MISSING", nrow(out)),
    source_provenance_record_complete = rep(FALSE, nrow(out)),
    source_exact_historical_provenance_recovered = rep(FALSE, nrow(out)),
    source_input_sha256 = rep(NA_character_, nrow(out)),
    source_input_sha256_verified = rep(FALSE, nrow(out)),
    source_provenance_evidence = rep(NA_character_, nrow(out)))
  for (nm in names(extra_cols)) if (!nm %in% names(out)) out[[nm]] <- extra_cols[[nm]]

  for (i in seq_len(nrow(out))) {
    hit <- cfg[cfg$species == out$species[i], , drop = FALSE]
    if (nrow(hit) == 0L) next
    h <- hit[1, , drop = FALSE]

    out$source_provider[i] <- h$source_provider[1]
    out$source_retrieval_date[i] <- h$retrieval_date[1]
    out$source_retrieval_version[i] <- h$retrieval_version[1]
    out$source_retrieval_url_or_query[i] <- h$retrieval_url_or_query[1]
    out$note[i] <- h$notes[1]

    if (!extended) {
      vals <- as.character(c(h$source_provider[1], h$retrieval_date[1],
                             h$retrieval_version[1], h$retrieval_url_or_query[1]))
      exact <- all(!is.na(vals) & nzchar(trimws(vals)))
      out$source_retrieval_metadata_complete[i] <- exact
      out$source_exact_historical_provenance_recovered[i] <- exact
      next
    }

    dstatus <- .status(h$retrieval_date_status[1])
    qstatus <- .status(h$query_export_settings_status[1])
    evidence <- as.character(h$provenance_evidence[1])
    input_rel <- as.character(h$input_file[1])
    expected_hash <- tolower(trimws(as.character(h$input_sha256[1])))

    allowed <- c("RECOVERED", "UNRECOVERABLE")
    status_ok <- dstatus %in% allowed && qstatus %in% allowed
    evidence_ok <- .nonempty(evidence) && file.exists(file.path(root, evidence))
    input_ok <- .nonempty(input_rel) && file.exists(file.path(root, input_rel))
    hash_pinned <- .prov_env$is_hex_sha256(expected_hash)
    actual_hash <- if (input_ok) sha256_file(file.path(root, input_rel)) else NA_character_
    hash_ok <- hash_pinned && !is.na(actual_hash) && identical(expected_hash, actual_hash)

    exact <- status_ok && dstatus == "RECOVERED" && qstatus == "RECOVERED" &&
      .nonempty(h$retrieval_date[1]) && .nonempty(h$retrieval_url_or_query[1]) && hash_ok
    record_complete <- status_ok && .nonempty(h$source_provider[1]) &&
      .nonempty(h$retrieval_version[1]) && evidence_ok && input_ok && hash_ok

    out$source_retrieval_date_status[i] <- dstatus
    out$source_query_export_settings_status[i] <- qstatus
    out$source_provenance_record_complete[i] <- record_complete
    out$source_exact_historical_provenance_recovered[i] <- exact
    out$source_input_sha256[i] <- if (hash_pinned) expected_hash else NA_character_
    out$source_input_sha256_verified[i] <- hash_ok
    out$source_provenance_evidence[i] <- evidence

    # Keep the legacy classifier conservative. A documented UNRECOVERABLE field
    # makes the reconstruction record complete but does not transform missing
    # historical facts into recovered facts or a PASS.
    out$source_retrieval_metadata_complete[i] <- exact
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
    "NIST_SOURCE_PROVENANCE_RECOVERY.md",
    "FINAL_ADVERSARIAL_AUDIT.md", "REPRODUCE_FINAL_ADVERSARIAL_AUDIT.txt",
    "R/nist_scan_lib.R", "R/audit_utils.R", "R/render_statistical_audit.R",
    "R/run_final_adversarial_audit.R", "R/run_final_adversarial_audit_v2.R",
    "R/run_final_contract_checks.R", "R/build_final_adversarial_summary.R",
    "R/finalize_adversarial_audit.R", "R/finalize_adversarial_audit_v2.R",
    "R/verify_source_provenance.R",
    "tables_r/statistical_audit/final_species_observed_results.csv",
    "tables_r/statistical_audit/final_species_bootstrap_summary.csv",
    "tables_r/statistical_audit/final_species_null_summary.csv",
    "tables_r/statistical_audit/final_species_calibration.csv",
    "tables_r/statistical_audit/final_species_injection_recovery.csv",
    "tables_r/statistical_audit/multiple_testing.csv",
    "tables_r/statistical_audit/alternative_null_results.csv",
    "tables_r/statistical_audit/final_adversarial_run_manifest.csv",
    "tables_r/statistical_audit/final_source_provenance_resolved.csv"
  )
  if (!is.null(json_path) && file.exists(json_path)) key <- c(key, basename(json_path))
  key <- unique(key[file.exists(file.path(root, key))])
  out <- data.frame(
    path = key,
    bytes = vapply(file.path(root, key), function(p) file.info(p)$size, numeric(1)),
    sha256 = vapply(file.path(root, key), sha256_file, character(1)),
    stringsAsFactors = FALSE)
  if (any(!vapply(out$sha256, .prov_env$is_hex_sha256, logical(1))))
    stop("provenance hash manifest contains a non-SHA-256 value")
  write_table(out, file.path(td, "provenance_hash_manifest.csv"))
  out
}

.postprocess_provenance_outputs <- function(root) {
  prov_path <- file.path(root, "tables_r/statistical_audit/final_source_provenance_resolved.csv")
  json_path <- file.path(root, "final_adversarial_audit.json")
  report_path <- file.path(root, "FINAL_ADVERSARIAL_AUDIT.md")
  if (!file.exists(prov_path) || !file.exists(json_path)) return(invisible(FALSE))

  prov <- utils::read.csv(prov_path, stringsAsFactors = FALSE, check.names = FALSE)
  species <- c("Fe", "Co")
  record_complete <- setNames(vapply(species, function(s) {
    z <- prov[prov$species == s, , drop = FALSE]
    nrow(z) == 1L && "source_provenance_record_complete" %in% names(z) &&
      truthy(z$source_provenance_record_complete[1])
  }, logical(1)), species)
  exact_recovered <- setNames(vapply(species, function(s) {
    z <- prov[prov$species == s, , drop = FALSE]
    nrow(z) == 1L && "source_exact_historical_provenance_recovered" %in% names(z) &&
      truthy(z$source_exact_historical_provenance_recovered[1])
  }, logical(1)), species)

  machine <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  machine$provenance_record_complete <- as.list(record_complete)
  machine$exact_historical_provenance_recovered <- as.list(exact_recovered)

  old <- unlist(machine$unresolved, use.names = FALSE)
  old <- old[!grepl("^Exact NIST (Fe|Co) II retrieval date/version/query provenance", old)]
  replacement <- character()
  for (s in species) {
    label <- paste(s, "II")
    if (!record_complete[[s]]) {
      replacement <- c(replacement, paste0(
        "NIST ", label, " provenance reconstruction is not yet machine-complete; ",
        "the frozen-input SHA-256 must be pinned and verified against the committed recovery record."))
    } else if (!exact_recovered[[s]]) {
      replacement <- c(replacement, paste0(
        "NIST ", label, " source identity and provenance reconstruction are documented, but the exact ",
        "historical retrieval/export date and exact historical query/export parameters are explicitly ",
        "UNRECOVERABLE after the recorded search."))
    }
  }
  machine$unresolved <- c(old, replacement)
  jsonlite::write_json(machine, json_path, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", digits = NA)

  if (file.exists(report_path)) {
    lines <- readLines(report_path, warn = FALSE, encoding = "UTF-8")
    lines <- gsub("| Exact source retrieval/version provenance |",
                  "| Exact historical retrieval/query provenance recovered |",
                  lines, fixed = TRUE)
    lines <- lines[!grepl("^- Exact NIST (Fe|Co) II retrieval date/version/query provenance", lines)]
    marker <- which(lines == "## Remaining uncertainties")
    if (length(marker) == 1L && length(replacement)) {
      insert_at <- marker + 1L
      while (insert_at <= length(lines) && lines[insert_at] == "") insert_at <- insert_at + 1L
      before <- if (insert_at > 1L) lines[seq_len(insert_at - 1L)] else character()
      after <- if (insert_at <= length(lines)) lines[insert_at:length(lines)] else character()
      lines <- c(before, paste0("- ", replacement), after)
    }
    writeLines(lines, report_path, useBytes = TRUE)
  }
  invisible(TRUE)
}

main <- function(argv = commandArgs(TRUE)) {
  res <- .legacy_finalize_main(argv)
  root <- audit_repo_root()
  .postprocess_provenance_outputs(root)
  # Re-hash after post-processing so the machine result/report hashes correspond
  # to the final bytes rather than the pre-correction intermediate files.
  hash_manifest(root, file.path(root, "final_adversarial_audit.json"))
  invisible(res)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("finalize_adversarial_audit_v2\\.R$", .invoked_file)) main()
