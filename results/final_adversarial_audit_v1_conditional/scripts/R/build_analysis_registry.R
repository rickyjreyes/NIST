#!/usr/bin/env Rscript
# build_analysis_registry.R
# ---------------------------------------------------------------------------
# Build the DECLARED analysis registry: every analysis the audit runs, with a
# unique id and the full specification needed for reproducibility and for the
# family-wise multiple-testing correction. The multiplicity analysis operates
# on THIS registry, not on an undocumented subset of convenient results.
#
# Writes:
#   tables_r/statistical_audit/analysis_registry.csv
#   outputs_r/statistical_audit/analysis_registry.json
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("emp_p", mode = "function")) source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this), "audit_utils.R"))

# registry_row(): one fully specified analysis record.
registry_row <- function(id, species, ion, source, bins, cfg,
                         purpose, status,
                         baseline_sigma = cfg$baseline_sigma,
                         degree = cfg$degree,
                         k_min = cfg$k_min, k_max = cfg$k_max,
                         ell_min = NA_real_, ell_max = NA_real_,
                         searched = TRUE,
                         out_root = "outputs_r/statistical_audit") {
  data.frame(
    analysis_id = id,
    species = species,
    ion = ion,
    source_field = source,
    bins = bins,
    spectral_range = if (is.na(ell_min)) "full" else sprintf("[%.4f,%.4f]", ell_min, ell_max),
    transform = "log_wavenumber",
    baseline_sigma = baseline_sigma,
    poly_degree = degree,
    k_min = k_min,
    k_max = k_max,
    n_k = cfg$n_k,
    null_n = cfg$null_n,
    bootstrap_n = cfg$bootstrap_n,
    seed = cfg$seed,
    filter_settings = sprintf("element==%s & sp_num==%d", species, ion),
    duplicate_rule = "drop duplicate wavenumber, keep first after sort",
    min_lines_rule = sprintf(">= %d unique lines", cfg$min_lines),
    run_purpose = purpose,
    confirmatory = status,
    output_path = file.path(out_root, id),
    software = "R/nist_scan_lib.R (canonical reproduction)",
    timestamp = audit_timestamp(),
    git_commit = git_commit(),
    searched_family = searched,
    stringsAsFactors = FALSE
  )
}

# build_registry(): enumerate the declared analysis family from the config.
build_registry <- function(cfg = default_audit_config(), neighbors = c("Cr","Mn","Co","Ni","Ti")) {
  rows <- list()
  add <- function(r) rows[[length(rows) + 1L]] <<- r

  # --- confirmatory: canonical Fe II bin-stability ladder (120/160/200) -----
  for (b in c(120L, 160L, 200L)) {
    add(registry_row(sprintf("fe_ion2_wn_bin%d_confirm", b),
                     "Fe", 2L, "wavenumber", b, cfg,
                     purpose = "canonical Fe II bin-stability claim",
                     status = "confirmatory", searched = FALSE))
  }

  # --- exploratory: full predefined bin grid (Fe II, wavenumber) ------------
  for (b in cfg$bins_grid) {
    add(registry_row(sprintf("fe_ion2_wn_bingrid%d", b),
                     "Fe", 2L, "wavenumber", b, cfg,
                     purpose = "bin-stability grid", status = "exploratory"))
  }

  # --- exploratory: model-sensitivity grid (sigma x degree at primary bins) -
  for (s in cfg$sigma_grid) for (d in cfg$degree_grid) {
    add(registry_row(sprintf("fe_ion2_wn_sig%d_deg%d_bin%d", s, d, cfg$bins_primary),
                     "Fe", 2L, "wavenumber", cfg$bins_primary, cfg,
                     purpose = "model sensitivity (baseline sigma x degree)",
                     status = "exploratory",
                     baseline_sigma = s, degree = d))
  }

  # --- exploratory: source-field replication (observed / ritz / wavenumber) -
  for (src in c("observed", "ritz", "wavenumber")) {
    add(registry_row(sprintf("fe_ion2_%s_bin%d_src", src, cfg$bins_primary),
                     "Fe", 2L, src, cfg$bins_primary, cfg,
                     purpose = "observed/Ritz source-field replication",
                     status = "exploratory"))
  }

  # --- exploratory: neighbouring ions (ion 2, primary bins) -----------------
  for (sp in neighbors) {
    add(registry_row(sprintf("%s_ion2_wn_bin%d", tolower(sp), cfg$bins_primary),
                     sp, 2L, "wavenumber", cfg$bins_primary, cfg,
                     purpose = "neighbouring-ion comparison", status = "exploratory"))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

main <- function(argv = commandArgs(TRUE)) {
  cfg <- default_audit_config()
  reg <- build_registry(cfg)
  root <- audit_repo_root()
  csv_path <- file.path(root, "tables_r/statistical_audit/analysis_registry.csv")
  json_path <- file.path(root, "outputs_r/statistical_audit/analysis_registry.json")
  write_table(reg, csv_path)
  dir.create(dirname(json_path), showWarnings = FALSE, recursive = TRUE)
  writeLines(jsonlite::toJSON(reg, dataframe = "rows", pretty = TRUE, auto_unbox = TRUE,
                              na = "null"), json_path)
  cat(sprintf("[registry] %d analyses declared\n", nrow(reg)))
  cat(sprintf("[registry] confirmatory=%d exploratory=%d\n",
              sum(reg$confirmatory == "confirmatory"),
              sum(reg$confirmatory == "exploratory")))
  cat(sprintf("[save] %s\n", csv_path))
  cat(sprintf("[save] %s\n", json_path))
  invisible(reg)
}

.invoked_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(.invoked_file) > 0L && grepl("build_analysis_registry\\.R$", .invoked_file)) {
  main()
}
