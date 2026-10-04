#!/usr/bin/env Rscript
# verify_source_provenance.R
# ---------------------------------------------------------------------------
# Verify the frozen Fe/Co input identity and provenance-disposition record.
#
# Historical fields may be explicitly UNRECOVERABLE only when accompanied by
# the committed recovery record. This script never infers a retrieval date or
# query from filesystem/publication timestamps.
#
# Usage:
#   Rscript R/verify_source_provenance.R --check
#   Rscript R/verify_source_provenance.R --write-config
#
# --write-config computes SHA-256 for the frozen inputs and pins those hashes in
# config/nist_source_provenance.csv. Commit that one-time update before the
# clean canonical release run. --check is read-only and fails on any mismatch.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
if (!exists("audit_repo_root", mode = "function")) source(file.path(R_DIR, "audit_utils.R"))

is_hex_sha256 <- function(x) {
  length(x) == 1L && !is.na(x) && grepl("^[0-9A-Fa-f]{64}$", trimws(x))
}

sha256_file_exact <- function(path) {
  if (!file.exists(path) || dir.exists(path)) stop("cannot hash missing/non-file path: ", path)

  normalize <- function(x) {
    x <- trimws(x)
    hit <- x[grepl("^[0-9A-Fa-f]{64}$", x)]
    if (length(hit)) tolower(hit[1]) else NA_character_
  }

  bin <- Sys.which("sha256sum")
  if (nzchar(bin)) {
    z <- tryCatch(system2(bin, path, stdout = TRUE, stderr = TRUE), error = function(e) character())
    if (length(z)) {
      token <- strsplit(trimws(z[1]), "[[:space:]]+")[[1]][1]
      if (is_hex_sha256(token)) return(tolower(token))
    }
  }

  bin <- Sys.which("shasum")
  if (nzchar(bin)) {
    z <- tryCatch(system2(bin, c("-a", "256", path), stdout = TRUE, stderr = TRUE),
                  error = function(e) character())
    if (length(z)) {
      token <- strsplit(trimws(z[1]), "[[:space:]]+")[[1]][1]
      if (is_hex_sha256(token)) return(tolower(token))
    }
  }

  # certutil is present on normal supported Windows installations and avoids
  # silently degrading the hash algorithm when Unix sha256sum is unavailable.
  bin <- Sys.which("certutil")
  if (nzchar(bin)) {
    z <- tryCatch(system2(bin, c("-hashfile", path, "SHA256"), stdout = TRUE, stderr = TRUE),
                  error = function(e) character())
    h <- normalize(gsub("[[:space:]]", "", z))
    if (is_hex_sha256(h)) return(h)
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    h <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
    if (is_hex_sha256(h)) return(tolower(h))
  }

  if (requireNamespace("openssl", quietly = TRUE)) {
    raw <- readBin(path, what = "raw", n = file.info(path)$size)
    h <- as.character(openssl::sha256(raw))
    if (is_hex_sha256(h)) return(tolower(h))
  }

  py <- Sys.which("python")
  if (!nzchar(py)) py <- Sys.which("python3")
  if (nzchar(py)) {
    code <- paste0(
      "import hashlib,sys; h=hashlib.sha256(); ",
      "f=open(sys.argv[1],'rb'); ",
      "[h.update(b) for b in iter(lambda:f.read(1048576),b'')]; ",
      "f.close(); print(h.hexdigest())")
    z <- tryCatch(system2(py, c("-c", shQuote(code), shQuote(path)), stdout = TRUE, stderr = TRUE),
                  error = function(e) character())
    h <- normalize(z)
    if (is_hex_sha256(h)) return(h)
  }

  stop("no working SHA-256 implementation found; refusing MD5 fallback")
}

parse_mode <- function(argv) {
  if ("--write-config" %in% argv) return("write")
  "check"
}

read_provenance_config <- function(path) {
  cfg <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                         fileEncoding = "UTF-8")
  # Be tolerant of a legacy UTF-8 BOM without requiring one.
  if (length(names(cfg))) names(cfg)[1] <- sub("^\\ufeff", "", names(cfg)[1])
  cfg
}

main <- function(argv = commandArgs(TRUE)) {
  root <- audit_repo_root()
  cfg_path <- file.path(root, "config", "nist_source_provenance.csv")
  if (!file.exists(cfg_path)) stop("missing provenance config: ", cfg_path)

  cfg <- read_provenance_config(cfg_path)
  required <- c(
    "species", "source_provider", "retrieval_date", "retrieval_date_status",
    "retrieval_version", "retrieval_url_or_query", "query_export_settings_status",
    "source_access_date", "local_file_first_seen", "provenance_evidence",
    "input_file", "input_sha256", "notes")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop("provenance config missing columns: ", paste(missing, collapse = ", "))

  if (!identical(sort(unique(cfg$species)), c("Co", "Fe")))
    stop("provenance config must contain exactly one Fe and one Co row")
  if (nrow(cfg) != 2L || anyDuplicated(cfg$species))
    stop("provenance config must contain exactly one row per Fe/Co species")

  allowed_status <- c("RECOVERED", "UNRECOVERABLE")
  if (any(!toupper(trimws(cfg$retrieval_date_status)) %in% allowed_status))
    stop("retrieval_date_status must be RECOVERED or UNRECOVERABLE")
  if (any(!toupper(trimws(cfg$query_export_settings_status)) %in% allowed_status))
    stop("query_export_settings_status must be RECOVERED or UNRECOVERABLE")

  evidence_ok <- vapply(cfg$provenance_evidence, function(rel) {
    !is.na(rel) && nzchar(trimws(rel)) && file.exists(file.path(root, rel))
  }, logical(1))
  if (any(!evidence_ok)) stop("provenance evidence record missing for one or more species")

  mode <- parse_mode(argv)
  actual <- character(nrow(cfg))
  for (i in seq_len(nrow(cfg))) {
    rel <- cfg$input_file[i]
    if (is.na(rel) || !nzchar(trimws(rel))) stop("blank input_file for ", cfg$species[i])
    actual[i] <- sha256_file_exact(file.path(root, rel))
  }

  if (mode == "write") {
    cfg$input_sha256 <- actual
    utils::write.csv(cfg, cfg_path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    cat("Pinned SHA-256 hashes in ", cfg_path, "\n", sep = "")
  } else {
    expected <- tolower(trimws(as.character(cfg$input_sha256)))
    if (any(!vapply(expected, is_hex_sha256, logical(1))))
      stop("input_sha256 is not pinned for every species; run once with --write-config, inspect, and commit")
    if (any(expected != actual)) {
      bad <- cfg$species[expected != actual]
      stop("frozen source SHA-256 mismatch for: ", paste(bad, collapse = ", "))
    }
  }

  date_status <- toupper(trimws(cfg$retrieval_date_status))
  query_status <- toupper(trimws(cfg$query_export_settings_status))
  exact_recovered <- date_status == "RECOVERED" & query_status == "RECOVERED"
  record_complete <- evidence_ok &
    !is.na(cfg$source_provider) & nzchar(trimws(cfg$source_provider)) &
    !is.na(cfg$retrieval_version) & nzchar(trimws(cfg$retrieval_version)) &
    vapply(actual, is_hex_sha256, logical(1))

  out <- data.frame(
    species = cfg$species,
    input_file = cfg$input_file,
    sha256 = actual,
    provenance_record_complete = record_complete,
    exact_historical_provenance_recovered = exact_recovered,
    retrieval_date_status = date_status,
    query_export_settings_status = query_status,
    stringsAsFactors = FALSE)

  print(out, row.names = FALSE)
  if (any(!record_complete)) stop("provenance reconstruction record is incomplete")
  invisible(out)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("verify_source_provenance\\.R$", .invoked_file)) main()
