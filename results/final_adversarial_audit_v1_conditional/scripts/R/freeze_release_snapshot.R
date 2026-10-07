#!/usr/bin/env Rscript
# freeze_release_snapshot.R
# ---------------------------------------------------------------------------
# Create an immutable-style release snapshot for the exact result produced by a
# clean pinned canonical run. Unlike the legacy PASS-only freezer, this script
# may preserve a CONDITIONAL or FAIL verdict in a verdict-labelled directory.
# It never promotes or changes the classification.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
source(file.path(R_DIR, "audit_utils.R"))
.prov <- new.env(parent = environment())
sys.source(file.path(R_DIR, "verify_source_provenance.R"), envir = .prov)

arg_value <- function(argv, flag, default = NA_character_) {
  i <- which(argv == flag)
  if (!length(i) || i[1] == length(argv)) return(default)
  argv[i[1] + 1L]
}

git_head <- function(root) {
  z <- system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE)
  trimws(z[1])
}

copy_file <- function(root, dest, rel, subdir = NULL) {
  src <- file.path(root, rel)
  if (!file.exists(src)) return(FALSE)
  dd <- if (is.null(subdir)) dest else file.path(dest, subdir)
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(dd, basename(src))
  if (file.exists(target)) return(TRUE)
  ok <- file.copy(src, target, overwrite = FALSE, copy.date = TRUE)
  if (!ok) stop("failed to copy release artifact: ", rel)
  TRUE
}

copy_dir_files <- function(src_dir, dest_dir, pattern = NULL) {
  if (!dir.exists(src_dir)) return(invisible(0L))
  files <- list.files(src_dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
  files <- files[!dir.exists(files)]
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(files)) return(invisible(0L))
  copied <- 0L
  for (src in files) {
    target <- file.path(dest_dir, basename(src))
    if (file.exists(target)) next
    ok <- file.copy(src, target, overwrite = FALSE, copy.date = TRUE)
    if (!ok) stop("failed copying release file: ", src)
    copied <- copied + 1L
  }
  invisible(copied)
}

safe_label <- function(x) {
  y <- tolower(trimws(x))
  if (grepl("^pass", y)) return("pass")
  if (grepl("^conditional", y)) return("conditional")
  if (grepl("^fail", y)) return("fail")
  stop("unknown classification: ", x)
}

main <- function(argv = commandArgs(TRUE)) {
  root <- audit_repo_root()
  input_commit <- arg_value(argv, "--input-commit", Sys.getenv("NIST_AUDIT_INPUT_COMMIT", ""))
  clean_start <- toupper(Sys.getenv("NIST_AUDIT_CLEAN_START", "")) == "TRUE"
  if (!clean_start) stop("snapshot refused: clean-start attestation is missing")
  if (!nzchar(input_commit) || !grepl("^[0-9A-Fa-f]{40}$", input_commit))
    stop("snapshot refused: --input-commit must be a full 40-hex commit SHA")
  head <- git_head(root)
  if (!identical(tolower(head), tolower(input_commit)))
    stop("snapshot refused: HEAD changed since the canonical run began")

  json_path <- file.path(root, "final_adversarial_audit.json")
  if (!file.exists(json_path)) stop("missing final_adversarial_audit.json")
  machine <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  if (!identical(tolower(as.character(machine$git_commit)), tolower(input_commit)))
    stop("final machine result is not tied to the clean input commit")

  # The source provenance hashes must already be pinned and verified before a
  # release snapshot is allowed, even when the historical date/query are
  # explicitly UNRECOVERABLE and the scientific verdict remains CONDITIONAL.
  .prov$main("--check")

  fe_class <- as.character(machine$classifications$Fe)
  label <- safe_label(fe_class)
  dest_name <- if (label == "pass") "final_adversarial_audit_v1" else
    paste0("final_adversarial_audit_v1_", label)
  dest <- file.path(root, "results", dest_name)
  if (dir.exists(dest)) stop("immutable release snapshot already exists: ", dest)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  # Core inputs/configuration/evidence.
  copy_file(root, dest, "data/Fe_lines.csv", "inputs")
  copy_file(root, dest, "data/Co_lines.csv", "inputs")
  copy_file(root, dest, "config/final_adversarial_audit.json", "config")
  copy_file(root, dest, "config/nist_source_provenance.csv", "config")
  copy_file(root, dest, "NIST_SOURCE_PROVENANCE_RECOVERY.md", "provenance")
  copy_file(root, dest, "FINAL_ADVERSARIAL_AUDIT.md")
  copy_file(root, dest, "final_adversarial_audit.json")
  copy_file(root, dest, "REPRODUCE_FINAL_ADVERSARIAL_AUDIT.txt")
  copy_file(root, dest, "reports/rendered/nist_final_adversarial_audit.html", "report")

  # Exact implementation/configuration present at the pinned commit. Existing
  # core config copies are skipped rather than treated as copy failures.
  copy_dir_files(file.path(root, "R"), file.path(dest, "scripts", "R"), "\\.R$")
  copy_dir_files(file.path(root, "config"), file.path(dest, "config"), "\\.(json|csv)$")
  copy_dir_files(file.path(root, "tables_r", "statistical_audit"),
                 file.path(dest, "tables"), "\\.csv$")
  copy_dir_files(file.path(root, "figures_r", "statistical_audit"),
                 file.path(dest, "figures"), "\\.(png|pdf|svg)$")

  context <- list(
    schema_version = "nist_clean_release_snapshot_v1",
    analysis_input_commit = input_commit,
    current_head = head,
    clean_working_tree_at_canonical_run_start = TRUE,
    classification = machine$classifications,
    canonical_seed = machine$canonical_seed,
    budgets = machine$budgets,
    provenance_record_complete = machine$provenance_record_complete,
    exact_historical_provenance_recovered = machine$exact_historical_provenance_recovered,
    snapshot_policy = paste(
      "immutable-style; directory name preserves PASS/CONDITIONAL/FAIL;",
      "a conditional snapshot must never be cited as a PASS release"))
  jsonlite::write_json(context, file.path(dest, "RUN_CONTEXT.json"), pretty = TRUE,
                       auto_unbox = TRUE, na = "null", digits = NA)

  # State metadata is written before HASHES.csv so it is itself integrity-pinned.
  writeLines(c(
    paste("analysis_input_commit:", input_commit),
    paste("classification_Fe:", machine$classifications$Fe),
    paste("classification_Co:", machine$classifications$Co),
    "clean_working_tree_at_canonical_run_start: true",
    paste("snapshot_directory:", sub("\\\\", "/", dest)),
    "policy: immutable-style; do not overwrite; CONDITIONAL is not PASS"
  ), file.path(dest, "GIT_STATE.txt"))

  files <- list.files(dest, recursive = TRUE, full.names = TRUE)
  files <- files[!dir.exists(files)]
  # HASHES.csv is itself excluded from its own manifest to avoid recursion.
  files <- files[basename(files) != "HASHES.csv"]
  rel <- substring(files, nchar(dest) + 2L)
  hashes <- data.frame(
    path = rel,
    bytes = vapply(files, function(p) file.info(p)$size, numeric(1)),
    sha256 = vapply(files, .prov$sha256_file_exact, character(1)),
    stringsAsFactors = FALSE)
  if (any(!vapply(hashes$sha256, .prov$is_hex_sha256, logical(1))))
    stop("release snapshot contains a non-SHA-256 digest")
  utils::write.csv(hashes, file.path(dest, "HASHES.csv"), row.names = FALSE)

  cat("Created release snapshot: ", dest, "\n", sep = "")
  cat("Preserved Fe classification: ", fe_class, "\n", sep = "")
  invisible(dest)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("freeze_release_snapshot\\.R$", .invoked_file)) main()
