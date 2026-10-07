#!/usr/bin/env Rscript
# finalize_adversarial_audit.R
# ---------------------------------------------------------------------------
# Evidence classification + immutable-style release finalizer.
#
# This script never changes the detector, target, null, or simulation outputs.
# It only reads completed final-audit artifacts, applies the declared decision
# rules, writes the required human/machine summaries, and freezes an immutable
# release directory ONLY when Fe II receives:
#   PASS — ADVERSARIAL AUDIT SURVIVED
#
# A CONDITIONAL or FAIL result is preserved as-is and is never "rescued" by
# redefining a target or rerunning a more favourable analysis.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
if (!exists("audit_repo_root", mode = "function")) source(file.path(R_DIR, "audit_utils.R"))

PASS_LABEL <- "PASS — ADVERSARIAL AUDIT SURVIVED"
CONDITIONAL_LABEL <- "CONDITIONAL"
FAIL_LABEL <- "FAIL"

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "PASS")
}

read_csv_required <- function(root, name) {
  p <- file.path(root, "tables_r/statistical_audit", name)
  if (!file.exists(p)) stop("required final-audit artifact missing: ", p)
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
}

read_csv_optional <- function(root, name) {
  p <- file.path(root, "tables_r/statistical_audit", name)
  if (!file.exists(p)) return(NULL)
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
}

git_full_commit <- function(root) {
  x <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"),
                        stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (length(x)) trimws(x[1]) else NA_character_
}

git_status <- function(root) {
  x <- tryCatch(system2("git", c("-C", root, "status", "--porcelain"),
                        stdout = TRUE, stderr = FALSE), error = function(e) character())
  paste(x, collapse = "\n")
}

sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  bin <- Sys.which("sha256sum")
  if (nzchar(bin)) {
    x <- tryCatch(system2(bin, path, stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
    if (length(x)) return(strsplit(trimws(x[1]), "[[:space:]]+")[[1]][1])
  }
  bin <- Sys.which("shasum")
  if (nzchar(bin)) {
    x <- tryCatch(system2(bin, c("-a", "256", path), stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
    if (length(x)) return(strsplit(trimws(x[1]), "[[:space:]]+")[[1]][1])
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    raw <- readBin(path, what = "raw", n = file.info(path)$size)
    return(as.character(openssl::sha256(raw)))
  }
  paste0("md5:", unname(tools::md5sum(path)))
}

load_protocol_config <- function(root) {
  p <- file.path(root, "config/final_adversarial_audit.json")
  if (!file.exists(p)) stop("missing protocol configuration: ", p)
  jsonlite::fromJSON(p, simplifyVector = TRUE)
}

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
    h <- hit[1, ]
    out$source_provider[i] <- h$source_provider
    out$source_retrieval_date[i] <- h$retrieval_date
    out$source_retrieval_version[i] <- h$retrieval_version
    out$source_retrieval_metadata_complete[i] <-
      all(nzchar(trimws(c(h$source_provider, h$retrieval_date,
                         h$retrieval_version, h$retrieval_url_or_query))))
    out$note[i] <- h$notes
  }
  out
}

primary_row <- function(df, species) {
  z <- df[df$species == species & df$bins == 160, , drop = FALSE]
  if (nrow(z) == 0L) stop("missing primary row for ", species)
  z[1, , drop = FALSE]
}

species_gate <- function(species, obs, bs, ns, cal, inj, art, influence,
                         expected_k, expected_lines = NA_integer_) {
  o <- obs[obs$species == species, , drop = FALSE]
  b <- bs[bs$species == species, , drop = FALSE]
  n <- ns[ns$species == species, , drop = FALSE]
  c <- cal[cal$species == species, , drop = FALSE]
  j <- inj[inj$species == species, , drop = FALSE]
  a <- art[art$species == species, , drop = FALSE]
  inf <- influence[influence$species == species, , drop = FALSE]
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

  required_alpha <- c(0.10, 0.05, 0.01, 0.001)
  alpha_ok <- all(vapply(required_alpha, function(x)
    any(abs(c$nominal_alpha - x) < 1e-12 & truthy(c$compatible)), logical(1)))

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

classify_species <- function(g, extra_hard = TRUE, provenance_complete = TRUE,
                             global_accounting_complete = TRUE) {
  hard <- c(g$reproduce, g$bootstrap, g$canonical_null, g$calibration,
            g$injection, g$cross_binning, g$implementation, extra_hard)
  if (any(!hard)) return(FAIL_LABEL)
  material <- c(g$influence, provenance_complete, global_accounting_complete)
  if (any(!material)) return(CONDITIONAL_LABEL)
  PASS_LABEL
}

format_p <- function(p, r = NA_integer_, B = NA_integer_) {
  if (!is.finite(p)) return("not established")
  if (is.finite(r) && is.finite(B) && r == 0L)
    return(sprintf("%.6g (= 1/(%d+1), zero exceedances; Monte-Carlo floor, not p=0)", p, B))
  sprintf("%.6g", p)
}

build_multiple_testing <- function(mult, species) {
  z <- mult[mult$species == species & mult$source == "wavenumber" &
              mult$bins == 160 & abs(mult$sigma - 6) < 1e-12 &
              mult$degree == 1, , drop = FALSE]
  if (nrow(z) == 0L) return(list(ok = FALSE, row = NULL, corrected = NA_real_))
  z <- z[1, , drop = FALSE]
  vals <- c(z$family_max_p[1], z$holm_p[1], z$bonferroni_p[1])
  ok <- all(is.finite(vals)) && all(vals <= 0.05)
  list(ok = ok, row = z, corrected = max(vals, na.rm = TRUE))
}

copy_named_outputs <- function(root, obs, bs, ns, cal, inj, mult) {
  td <- file.path(root, "tables_r/statistical_audit")
  write_table(obs[obs$species == "Fe", , drop = FALSE],
              file.path(td, "final_fe_ii_results.csv"))
  write_table(obs[obs$species == "Co", , drop = FALSE],
              file.path(td, "final_co_ii_results.csv"))
  write_table(bs, file.path(td, "bootstrap_summary.csv"))
  write_table(ns, file.path(td, "null_summary.csv"))
  write_table(cal, file.path(td, "calibration_summary.csv"))
  write_table(inj, file.path(td, "injection_recovery_summary.csv"))
  write_table(mult, file.path(td, "multiple_testing_accounting.csv"))
  write_table(obs[, c("species", "bins", "k_best", "delta_log_x", "scale_ratio",
                      "n_obs", "amplitude", "phase", "deltaD", "primary_k",
                      "relative_k_error", "compatible_2pct")],
              file.path(td, "cross_binning_robustness.csv"))
}

hash_manifest <- function(root, json_path) {
  td <- file.path(root, "tables_r/statistical_audit")
  key <- c(
    "data/Fe_lines.csv", "data/Co_lines.csv",
    "config/final_adversarial_audit.json", "config/nist_source_provenance.csv",
    "R/nist_scan_lib.R", "R/audit_utils.R", "R/render_statistical_audit.R",
    "R/run_final_adversarial_audit.R", "R/run_final_contract_checks.R",
    "R/build_final_adversarial_summary.R", "R/finalize_adversarial_audit.R",
    "tables_r/statistical_audit/final_species_observed_results.csv",
    "tables_r/statistical_audit/final_species_bootstrap_summary.csv",
    "tables_r/statistical_audit/final_species_null_summary.csv",
    "tables_r/statistical_audit/final_species_calibration.csv",
    "tables_r/statistical_audit/final_species_injection_recovery.csv",
    "tables_r/statistical_audit/multiple_testing.csv",
    "tables_r/statistical_audit/alternative_null_results.csv",
    "tables_r/statistical_audit/final_adversarial_run_manifest.csv"
  )
  if (!is.null(json_path)) key <- c(key, sub(paste0("^", root, "/?"), "", json_path))
  key <- unique(key[file.exists(file.path(root, key))])
  out <- data.frame(
    path = key,
    bytes = vapply(file.path(root, key), function(p) file.info(p)$size, numeric(1)),
    sha256 = vapply(file.path(root, key), sha256_file, character(1)),
    stringsAsFactors = FALSE)
  write_table(out, file.path(td, "provenance_hash_manifest.csv"))
  out
}

write_reproduce <- function(root) {
  txt <- paste(
    "# Canonical final adversarial audit",
    "Rscript R/render_statistical_audit.R \\",
    "  --bootstrap-n 5000 \\",
    "  --null-n 5000 \\",
    "  --calibration-n 10000 \\",
    "  --injection-n 2000 \\",
    "  --seed 20260517 \\",
    "  --parallel true \\",
    "  --force \\",
    "  --strict \\",
    "  --render-report true",
    "",
    "# Sequential fallback is statistically identical when future/future.apply are unavailable:",
    "# use --parallel false; the canonical seed and per-replicate RNG streams are unchanged.",
    sep = "\n")
  writeLines(txt, file.path(root, "REPRODUCE_FINAL_ADVERSARIAL_AUDIT.txt"))
  txt
}

freeze_release <- function(root, classification, hashes) {
  if (!identical(classification, PASS_LABEL)) return(NA_character_)
  dest <- file.path(root, "results/final_adversarial_audit_v1")
  if (dir.exists(dest)) {
    message("[freeze] immutable release already exists; preserving it unchanged: ", dest)
    return(dest)
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  copy_one <- function(rel, subdir = NULL) {
    src <- file.path(root, rel)
    if (!file.exists(src)) return(invisible(FALSE))
    dd <- if (is.null(subdir)) dest else file.path(dest, subdir)
    dir.create(dd, recursive = TRUE, showWarnings = FALSE)
    file.copy(src, file.path(dd, basename(src)), overwrite = FALSE,
              recursive = dir.exists(src), copy.date = TRUE)
  }

  copy_one("data/Fe_lines.csv", "inputs")
  copy_one("data/Co_lines.csv", "inputs")
  copy_one("config/final_adversarial_audit.json", "config")
  copy_one("config/nist_source_provenance.csv", "config")
  copy_one("FINAL_ADVERSARIAL_AUDIT.md")
  copy_one("final_adversarial_audit.json")
  copy_one("REPRODUCE_FINAL_ADVERSARIAL_AUDIT.txt")
  copy_one("reports/rendered/nist_final_adversarial_audit.html", "report")

  rfiles <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  dir.create(file.path(dest, "scripts/R"), recursive = TRUE, showWarnings = FALSE)
  file.copy(rfiles, file.path(dest, "scripts/R", basename(rfiles)),
            overwrite = FALSE, copy.date = TRUE)

  tables <- list.files(file.path(root, "tables_r/statistical_audit"),
                       pattern = "\\.csv$", full.names = TRUE)
  dir.create(file.path(dest, "tables"), recursive = TRUE, showWarnings = FALSE)
  file.copy(tables, file.path(dest, "tables", basename(tables)),
            overwrite = FALSE, copy.date = TRUE)

  figs <- list.files(file.path(root, "figures_r/statistical_audit"),
                     pattern = "\\.(png|svg|pdf)$", full.names = TRUE)
  dir.create(file.path(dest, "figures"), recursive = TRUE, showWarnings = FALSE)
  file.copy(figs, file.path(dest, "figures", basename(figs)),
            overwrite = FALSE, copy.date = TRUE)

  files <- list.files(dest, recursive = TRUE, full.names = TRUE)
  files <- files[!dir.exists(files)]
  rel <- substring(files, nchar(dest) + 2L)
  fh <- data.frame(path = rel,
                   bytes = vapply(files, function(p) file.info(p)$size, numeric(1)),
                   sha256 = vapply(files, sha256_file, character(1)),
                   stringsAsFactors = FALSE)
  utils::write.csv(fh, file.path(dest, "HASHES.csv"), row.names = FALSE)
  writeLines(c(
    paste("git_commit:", git_full_commit(root)),
    paste("working_tree_clean_at_finalize:", !nzchar(git_status(root))),
    "freeze_policy: immutable; do not overwrite this directory in later exploratory work"
  ), file.path(dest, "GIT_STATE.txt"))
  dest
}

write_report <- function(root, verdict, fe, co, prov_complete, mult_fe, mult_co,
                         alt_fe_ok, alt_fe_worst, ready, unresolved) {
  ffix <- fe$fixed
  cfix <- co$fixed
  q <- c(
    "# FINAL ADVERSARIAL AUDIT",
    "",
    "## Verdict at a glance",
    "",
    sprintf("1. **Does the Fe II result reproduce from the frozen raw inputs?** %s.",
            if (fe$reproduce) "Yes" else "No"),
    sprintf("2. **Does it survive bootstrap resampling?** %s.",
            if (fe$bootstrap) "Yes under the predeclared 80% peak-region stability criterion" else "No"),
    sprintf("3. **Do the declared null models produce comparable features?** %s.",
            if (fe$canonical_null && alt_fe_ok) "No at the declared 0.05 decision level" else "At least one required null check remains comparable or unresolved"),
    sprintf("4. **Is the reported statistical significance empirically calibrated?** %s.",
            if (fe$calibration) "Yes at 10%, 5%, 1%, and 0.1% within Monte-Carlo confidence intervals" else "No or unresolved"),
    sprintf("5. **Can injected signals of the observed magnitude be reliably recovered?** %s.",
            if (fe$injection) "Yes under the predeclared 80% significant-correct-location recovery criterion" else "No"),
    sprintf("6. **Is the result robust across the predeclared binning choices?** %s.",
            if (fe$cross_binning) "Yes across 120/160/200 bins" else "No"),
    sprintf("7. **Does Co II provide independent compatible evidence?** %s.",
            if (co$reproduce && co$bootstrap && co$canonical_null)
              "Co II survives its own species-level checks, but it is not pooled with Fe II and is not described as an independent same-frequency experimental replication"
            else "No strong independent species-level support is established"),
    sprintf("8. **What is the fixed-frequency significance?** Fe II primary descriptive fixed-k p = %s. The frequency was historically scan-selected, so this is **not** a global discovery p-value.",
            format_p(ffix$fixed_frequency_p[1], ffix$fixed_exceedances[1], ffix$B[1])),
    sprintf("9. **What is the appropriately corrected/global significance, if established?** Fe II scan-global p = %s; conservative reported FWER-adjusted value = %s.",
            format_p(ffix$global_p[1], ffix$global_exceedances[1], ffix$B[1]),
            if (is.finite(mult_fe$corrected)) format_p(mult_fe$corrected) else "not established"),
    sprintf("10. **What is the final Fe II classification?** **%s**.", verdict$Fe),
    sprintf("11. **What is the final Co II classification?** **%s**.", verdict$Co),
    sprintf("12. **Is the result ready to freeze?** %s.",
            if (ready) "Yes; the immutable-style release directory was created or already exists unchanged" else "No; unresolved items are listed below"),
    "",
    "## Scientific scope",
    "",
    "This audit can establish **a statistically robust empirical spectral regularity in the frozen NIST catalog**. It does **not**, by itself, validate Wave Confinement Theory as a physical theory, establish a causal mechanism, or constitute NIST endorsement. Those are separate claims requiring separate physical evidence.",
    "",
    "Fe II is primary. Co II is analyzed separately and is never combined post hoc with Fe II to increase significance. A fixed-frequency p-value at a frequency selected by a scan of the same data is labelled descriptive; the scan-global and declared-family corrections carry the discovery interpretation.",
    "",
    "## Core decision matrix",
    "",
    "| Check | Fe II | Co II |",
    "|---|---:|---:|",
    sprintf("| Frozen observed statistic reproduces | %s | %s |", fe$reproduce, co$reproduce),
    sprintf("| Bootstrap peak stability | %s | %s |", fe$bootstrap, co$bootstrap),
    sprintf("| Canonical scan-global null | %s | %s |", fe$canonical_null, co$canonical_null),
    sprintf("| Empirical p-value calibration | %s | %s |", fe$calibration, co$calibration),
    sprintf("| 1.0x observed injection recovery | %s | %s |", fe$injection, co$injection),
    sprintf("| 120/160/200 bin robustness | %s | %s |", fe$cross_binning, co$cross_binning),
    sprintf("| Sorting/cleaning implementation checks | %s | %s |", fe$implementation, co$implementation),
    sprintf("| Leave-decile-out influence stability | %s | %s |", fe$influence, co$influence),
    sprintf("| Declared-family/global accounting | %s | %s |", mult_fe$ok, mult_co$ok),
    sprintf("| Exact source retrieval/version provenance | %s | %s |", prov_complete["Fe"], prov_complete["Co"]),
    "",
    "## Fe II adversarial null sensitivity",
    "",
    sprintf("Worst declared Fe II alternative-null scan-global p: **%s** under **%s**. The original canonical null result is preserved regardless of this sensitivity result.",
            if (is.finite(alt_fe_worst$p)) format_p(alt_fe_worst$p) else "not available",
            alt_fe_worst$model),
    "",
    "## Remaining uncertainties",
    ""
  )
  if (length(unresolved) == 0L) q <- c(q, "- None within the declared statistical audit contract.")
  else q <- c(q, paste0("- ", unresolved))
  q <- c(q,
         "",
         "Finite Monte-Carlo results use the +1 correction, `(r + 1)/(B + 1)`. Zero exceedances are never reported as `p = 0`; they are reported at the simulation-resolution floor.",
         "",
         "See `final_adversarial_audit.json`, `tables_r/statistical_audit/`, `REPRODUCE_FINAL_ADVERSARIAL_AUDIT.txt`, and (only after a PASS freeze) `results/final_adversarial_audit_v1/` for machine-readable evidence and reproduction artifacts.")
  writeLines(q, file.path(root, "FINAL_ADVERSARIAL_AUDIT.md"))
  q
}

main <- function(argv = commandArgs(TRUE)) {
  root <- audit_repo_root()
  cfg <- load_protocol_config(root)
  manifest <- read_csv_required(root, "final_adversarial_run_manifest.csv")
  obs <- read_csv_required(root, "final_species_observed_results.csv")
  bs <- read_csv_required(root, "final_species_bootstrap_summary.csv")
  ns <- read_csv_required(root, "final_species_null_summary.csv")
  cal <- read_csv_required(root, "final_species_calibration.csv")
  inj <- read_csv_required(root, "final_species_injection_recovery.csv")
  art <- read_csv_required(root, "final_species_artifact_checks.csv")
  influence <- read_csv_required(root, "final_species_influence.csv")
  contract_prov <- read_csv_required(root, "final_source_provenance.csv")
  mult <- read_csv_required(root, "multiple_testing.csv")
  alt <- read_csv_required(root, "alternative_null_results.csv")

  budget_ok <- manifest$bootstrap_n[1] >= cfg$budgets$bootstrap_n &&
    manifest$null_n[1] >= cfg$budgets$null_n &&
    manifest$calibration_n[1] >= cfg$budgets$calibration_n &&
    manifest$injection_n[1] >= cfg$budgets$injection_n &&
    manifest$canonical_seed[1] == cfg$seed &&
    truthy(manifest$all_modules_ok[1])
  if (!budget_ok) stop("finalization refused: run manifest does not satisfy frozen final budgets/seed")

  prov <- resolved_provenance(root, contract_prov)
  write_table(prov, file.path(root, "tables_r/statistical_audit/final_source_provenance_resolved.csv"))
  pc <- setNames(vapply(c("Fe", "Co"), function(s) {
    z <- prov[prov$species == s, , drop = FALSE]
    nrow(z) == 1L && truthy(z$source_retrieval_metadata_complete[1])
  }, logical(1)), c("Fe", "Co"))

  fe <- species_gate("Fe", obs, bs, ns, cal, inj, art, influence,
                     expected_k = cfg$expected$Fe$k_primary,
                     expected_lines = cfg$expected$Fe$retained_lines)
  co <- species_gate("Co", obs, bs, ns, cal, inj, art, influence,
                     expected_k = cfg$expected$Co$k_primary,
                     expected_lines = cfg$expected$Co$retained_lines)

  mfe <- build_multiple_testing(mult, "Fe")
  mco <- build_multiple_testing(mult, "Co")

  alt_bad <- alt[is.finite(alt$scan_global_p) & alt$scan_global_p > 0.05, , drop = FALSE]
  alt_fe_ok <- nrow(alt) > 0L && nrow(alt_bad) == 0L
  wa <- if (nrow(alt)) alt[which.max(alt$scan_global_p), , drop = FALSE] else NULL
  alt_worst <- list(p = if (is.null(wa)) NA_real_ else wa$scan_global_p[1],
                    model = if (is.null(wa)) "not available" else wa$null_model[1])

  fe_class <- classify_species(fe, extra_hard = alt_fe_ok,
                               provenance_complete = pc["Fe"],
                               global_accounting_complete = mfe$ok)
  co_class <- classify_species(co, extra_hard = TRUE,
                               provenance_complete = pc["Co"],
                               global_accounting_complete = mco$ok)
  verdict <- list(Fe = fe_class, Co = co_class)

  unresolved <- character()
  if (!pc["Fe"]) unresolved <- c(unresolved,
    "Exact NIST Fe II retrieval date/version/query provenance is not yet populated in `config/nist_source_provenance.csv`; repository hashes alone must not be presented as a source-version record.")
  if (!pc["Co"]) unresolved <- c(unresolved,
    "Exact NIST Co II retrieval date/version/query provenance is not yet populated in `config/nist_source_provenance.csv`.")
  if (!fe$influence) unresolved <- c(unresolved,
    "Fe II leave-decile-out influence stability is below the audit-declared 80% robustness convention.")
  if (!co$influence) unresolved <- c(unresolved,
    "Co II leave-decile-out influence stability is below the audit-declared 80% robustness convention.")
  if (!mfe$ok) unresolved <- c(unresolved,
    "Fe II declared-family/global multiple-testing accounting does not meet the 0.05 criterion.")
  if (!mco$ok) unresolved <- c(unresolved,
    "Co II declared-family/global multiple-testing accounting does not meet the 0.05 criterion.")
  if (!alt_fe_ok) unresolved <- c(unresolved,
    "At least one declared Fe II alternative null produces a scan-global p-value above 0.05.")
  if (fe_class == FAIL_LABEL) unresolved <- c(unresolved,
    "At least one hard Fe II reproduction/bootstrap/null/calibration/injection/binning/implementation gate failed; the target must not be redefined to rescue it.")
  if (co_class == FAIL_LABEL) unresolved <- c(unresolved,
    "At least one hard Co II species-level gate failed; Co II must not be described as an automatic replication of Fe II.")

  copy_named_outputs(root, obs, bs, ns, cal, inj, mult)
  reproduction_command <- write_reproduce(root)

  machine <- list(
    schema_version = "final_adversarial_audit_v1",
    generated_utc = audit_timestamp(),
    git_commit = git_full_commit(root),
    working_tree_clean_at_finalize = !nzchar(git_status(root)),
    canonical_seed = cfg$seed,
    budgets = cfg$budgets,
    classifications = verdict,
    fe = list(
      reproduces = fe$reproduce, bootstrap_survives = fe$bootstrap,
      canonical_null_survives = fe$canonical_null,
      alternative_nulls_survive = alt_fe_ok,
      calibrated = fe$calibration, injection_recovery = fe$injection,
      cross_binning = fe$cross_binning, implementation_checks = fe$implementation,
      influence_stable = fe$influence,
      fixed_frequency = list(
        p = fe$fixed$fixed_frequency_p[1],
        exceedances = fe$fixed$fixed_exceedances[1],
        B = fe$fixed$B[1],
        interpretation = "descriptive only because k was scan-selected from the same data"),
      scan_global = list(
        p = fe$fixed$global_p[1],
        exceedances = fe$fixed$global_exceedances[1],
        B = fe$fixed$B[1]),
      corrected_global = list(
        conservative_fwer_p = mfe$corrected,
        family_size = if (is.null(mfe$row)) NA_integer_ else mfe$row$family_size[1])),
    co = list(
      reproduces = co$reproduce, bootstrap_survives = co$bootstrap,
      canonical_null_survives = co$canonical_null,
      calibrated = co$calibration, injection_recovery = co$injection,
      cross_binning = co$cross_binning, implementation_checks = co$implementation,
      influence_stable = co$influence,
      fixed_frequency = list(
        p = co$fixed$fixed_frequency_p[1],
        exceedances = co$fixed$fixed_exceedances[1],
        B = co$fixed$B[1],
        interpretation = "descriptive only because Co k was scan-selected from Co data"),
      scan_global = list(
        p = co$fixed$global_p[1],
        exceedances = co$fixed$global_exceedances[1],
        B = co$fixed$B[1]),
      corrected_global = list(
        conservative_fwer_p = mco$corrected,
        family_size = if (is.null(mco$row)) NA_integer_ else mco$row$family_size[1]),
      evidence_relationship = paste(
        "independent species-level statistical audit; not pooled with Fe II;",
        "not claimed as an independent same-frequency experiment")),
    provenance_complete = as.list(pc),
    source_scope = paste(
      "Statistically robust empirical spectral regularity in the frozen NIST catalog,",
      "if the declared gates pass. This artifact alone does not validate Wave Confinement",
      "Theory as a physical theory or establish a causal mechanism."),
    unresolved = unresolved,
    reproduction_command = reproduction_command
  )

  json_path <- file.path(root, "final_adversarial_audit.json")
  jsonlite::write_json(machine, json_path, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", digits = NA)

  ready <- identical(fe_class, PASS_LABEL)
  write_report(root, verdict, fe, co, pc, mfe, mco, alt_fe_ok, alt_worst,
               ready, unresolved)

  hashes <- hash_manifest(root, json_path)
  freeze <- freeze_release(root, fe_class, hashes)

  matrix <- data.frame(
    species = c("Fe", "Co"),
    reproduction = c(fe$reproduce, co$reproduce),
    bootstrap = c(fe$bootstrap, co$bootstrap),
    canonical_null = c(fe$canonical_null, co$canonical_null),
    calibration = c(fe$calibration, co$calibration),
    injection = c(fe$injection, co$injection),
    cross_binning = c(fe$cross_binning, co$cross_binning),
    implementation = c(fe$implementation, co$implementation),
    influence = c(fe$influence, co$influence),
    multiple_testing = c(mfe$ok, mco$ok),
    provenance_complete = unname(pc[c("Fe", "Co")]),
    classification = c(fe_class, co_class),
    stringsAsFactors = FALSE)
  write_table(matrix, file.path(root, "tables_r/statistical_audit/final_pass_conditional_fail_matrix.csv"))

  cat("\n=====================================================\n")
  cat("FINAL ADVERSARIAL AUDIT CLASSIFICATION\n")
  cat("=====================================================\n")
  cat("Fe II: ", fe_class, "\n", sep = "")
  cat("Co II: ", co_class, "\n", sep = "")
  cat("freeze: ", if (is.na(freeze)) "not created" else freeze, "\n", sep = "")
  if (length(unresolved)) {
    cat("unresolved:\n")
    for (u in unresolved) cat(" - ", u, "\n", sep = "")
  }
  invisible(list(verdict = verdict, machine = machine, matrix = matrix, freeze = freeze))
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("finalize_adversarial_audit\\.R$", .invoked_file)) main()
