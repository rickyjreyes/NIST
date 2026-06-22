#!/usr/bin/env Rscript
# render_statistical_audit.R
# ---------------------------------------------------------------------------
# Command-line orchestrator for the complete NIST statistical audit. Runs every
# audit module in dependency order, builds Python/R parity and the final claim
# matrix, optionally renders the Quarto report, and prints a completion summary.
#
# Each module is sourced into its OWN environment so the per-module main()
# functions do not clobber one another. Deterministic parallel RNG
# (L'Ecuyer-CMRG) is used so results do not depend on the number of workers.
#
# Example (development):
#   Rscript R/render_statistical_audit.R --fast --force --render-report true
#
# Example (final, long-running):
#   Rscript R/render_statistical_audit.R --bootstrap-n 5000 --null-n 5000 \
#       --calibration-n 10000 --injection-n 2000 --seed 20260517 \
#       --parallel true --force --strict --render-report true
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
AUDIT_R_DIR <- if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this)
source(file.path(AUDIT_R_DIR, "audit_utils.R"))

# --- argument parsing ------------------------------------------------------
parse_render_args <- function(argv) {
  d <- list(
    `data-root` = "data", `python-root` = "outputs", `r-root` = "outputs_r",
    `out-root` = "outputs_r/statistical_audit",
    `table-dir` = "tables_r/statistical_audit",
    `figure-dir` = "figures_r/statistical_audit",
    report = "reports/nist_statistical_audit.qmd",
    `bootstrap-n` = 2000L, `null-n` = 5000L, `calibration-n` = 2000L,
    `injection-n` = 500L, seed = 20260517L,
    parallel = FALSE, workers = "auto", `render-report` = TRUE,
    force = FALSE, strict = FALSE, fast = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!key %in% names(d)) { i <- i + 1L; next }
    # boolean flags may appear with or without an explicit value
    if (key %in% c("force", "strict", "fast")) {
      nxt <- if (i + 1L <= length(argv)) argv[i + 1L] else NA
      if (!is.na(nxt) && nxt %in% c("true", "false")) { d[[key]] <- (nxt == "true"); i <- i + 2L }
      else { d[[key]] <- TRUE; i <- i + 1L }
    } else {
      d[[key]] <- argv[i + 1L]; i <- i + 2L
    }
  }
  num <- c("bootstrap-n", "null-n", "calibration-n", "injection-n", "seed")
  for (k in num) d[[k]] <- as.integer(d[[k]])
  for (k in c("parallel", "render-report")) if (is.character(d[[k]])) d[[k]] <- (d[[k]] == "true")
  d
}

# run_module(): source a module into a fresh env and call its main(argv),
# capturing success/failure without aborting the whole pipeline.
run_module <- function(name, argv, results) {
  file <- file.path(AUDIT_R_DIR, name)
  cat(sprintf("\n=== [%s] ===\n", name))
  t0 <- Sys.time()
  ok <- tryCatch({
    e <- new.env(parent = globalenv())
    sys.source(file, envir = e)
    e$main(argv)
    TRUE
  }, error = function(err) { message(sprintf("[FAIL] %s: %s", name, conditionMessage(err))); FALSE })
  results[[name]] <- list(ok = ok, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  results
}

# --- python/R parity -------------------------------------------------------
build_parity <- function(py_root, r_root, table_dir, fig_dir) {
  dirs <- c("fe_ion2_120", "fe_ion2_160", "fe_ion2_200")
  rel <- function(a, b) abs(a - b) / max(abs(b), 1e-12)
  rows <- list()
  for (d in dirs) {
    pj <- file.path(py_root, d, "nist_summary.json")
    rj <- file.path(r_root, d, "nist_summary.json")
    if (!file.exists(pj) || !file.exists(rj)) next
    ps <- jsonlite::fromJSON(pj); rs <- jsonlite::fromJSON(rj)
    rows[[d]] <- data.frame(
      dir = d,
      k_best_abs_diff = abs(ps$best$k_best - rs$best$k_best),
      deltaD_rel_diff = rel(rs$best$deltaD, ps$best$deltaD),
      amplitude_rel_diff = rel(rs$best$amplitude, ps$best$amplitude),
      n_obs_abs_diff = abs(ps$branch_report$n_obs - rs$branch_report$n_obs),
      parity_pass = (abs(ps$best$k_best - rs$best$k_best) <= 1e-10 &&
                     rel(rs$best$deltaD, ps$best$deltaD) <= 1e-5),
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0L) return(NULL)
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  write_table(tab, file.path(table_dir, "python_r_parity.csv"))
  p <- ggplot2::ggplot(tab, ggplot2::aes(dir, k_best_abs_diff)) +
    ggplot2::geom_col(fill = species_colour("Fe"), width = 0.6) +
    ggplot2::geom_hline(yintercept = 1e-10, linetype = 2, colour = "grey40") +
    ggplot2::labs(title = "Python/R deterministic parity (|k_best| difference)",
                  subtitle = "Dashed line = 1e-10 tolerance. Bootstrap RNG differs by construction and is not compared.",
                  x = NULL, y = "|k_best_python - k_best_R|") +
    theme_audit()
  save_fig(p, file.path(fig_dir, "fig15_python_r_parity.png"))
  tab
}

# --- evidence dashboard ----------------------------------------------------
build_dashboard <- function(fig_dir, table_dir) {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  rd <- function(f) { p <- file.path(table_dir, f); if (file.exists(p)) utils::read.csv(p) else NULL }
  bg <- rd("bin_grid_results.csv"); cal <- rd("null_calibration.csv")
  inj <- rd("injection_recovery.csv"); eff <- rd("effect_sizes.csv")
  panels <- list()
  if (!is.null(bg)) panels[[length(panels)+1L]] <- ggplot2::ggplot(bg, ggplot2::aes(bins, k_best)) +
    ggplot2::geom_line(colour = species_colour("Fe")) + ggplot2::geom_point(colour = species_colour("Fe")) +
    ggplot2::labs(title = "Bin stability", x = "bins", y = "k_best") + theme_audit(9)
  if (!is.null(eff)) panels[[length(panels)+1L]] <- ggplot2::ggplot(eff, ggplot2::aes(stats::reorder(species, z_null), z_null, fill = species)) +
    ggplot2::geom_col() + ggplot2::scale_fill_manual(values = SPECIES_PALETTE, guide = "none") +
    ggplot2::coord_flip() + ggplot2::labs(title = "Effect size (z_null)", x = NULL, y = "z_null") + theme_audit(9)
  if (!is.null(cal)) panels[[length(panels)+1L]] <- ggplot2::ggplot(cal, ggplot2::aes(nominal_alpha, observed_fpr)) +
    ggplot2::geom_abline(linetype = 2, colour = "grey50") + ggplot2::geom_point(colour = species_colour("Fe"), size = 2) +
    ggplot2::labs(title = "Null calibration", x = "alpha", y = "FPR") + theme_audit(9)
  if (!is.null(inj)) { fe <- inj[inj$freq_name == "fe_reference", ]
    panels[[length(panels)+1L]] <- ggplot2::ggplot(fe, ggplot2::aes(amplitude, sig_correct_prob)) +
      ggplot2::geom_line(colour = species_colour("Fe")) + ggplot2::geom_point(colour = species_colour("Fe")) +
      ggplot2::labs(title = "Injection power", x = "amplitude", y = "power") + theme_audit(9) }
  if (length(panels) == 0L) return(NULL)
  comb <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(title = "Fe II statistical-audit evidence dashboard",
                               caption = "Computational/statistical evidence only; not independent physical validation.")
  save_fig(comb, file.path(fig_dir, "fig16_evidence_dashboard.png"), width = 11, height = 8)
  comb
}

# --- final claim matrix ----------------------------------------------------
build_claim_matrix <- function(table_dir) {
  rd <- function(f) { p <- file.path(table_dir, f); if (file.exists(p)) utils::read.csv(p) else NULL }
  sig <- rd("significance_results.csv"); bss <- rd("bin_stability_summary.csv")
  ci <- rd("peak_confidence_intervals.csv"); stab <- rd("peak_stability.csv")
  par <- rd("python_r_parity.csv"); orr <- rd("observed_ritz_replication.csv")
  cal <- rd("null_calibration.csv"); ho <- rd("holdout_results.csv")
  mc <- rd("model_comparison.csv")
  rows <- list()
  add <- function(claim, analysis, result, threshold, verdict, limitations, ref)
    rows[[length(rows)+1L]] <<- data.frame(claim = claim, analysis = analysis,
      result = result, threshold = threshold, verdict = verdict,
      limitations = limitations, output_reference = ref, stringsAsFactors = FALSE)

  if (!is.null(par)) add("Python/R computational agreement", "python_r_parity",
    sprintf("max |k_best| diff = %.2e", max(par$k_best_abs_diff)), "<= 1e-10",
    if (all(par$parity_pass)) "pass" else "fail",
    "deterministic pipeline only; RNG bootstrap excluded", "python_r_parity.csv")
  if (!is.null(sig)) { g <- sig[sig$analysis_id == "fe_ion2_wn_bin160", ][1, ]
    add("Fe II scan-global significance", "significance_results",
      sprintf("global p = %.4g (tail=%d/%d)", g$global_p, g$global_tail_count, g$global_B),
      "report (resolution-limited)",
      if (g$zero_exceedance) "pass (zero exceedances)" else "inconclusive",
      "p bounded by 1/(B+1); not an exact p below the floor", "significance_results.csv") }
  if (!is.null(bss)) add("Fe II bin stability", "bin_stability_summary",
    sprintf("%.0f%% of bins in reference region; CV(k)=%.4f", 100*bss$pct_in_reference, bss$cv_k),
    ">= 80% (descriptive)", if (bss$pct_in_reference >= 0.8) "pass" else "inconclusive",
    "descriptive threshold, not a universal law", "bin_stability_summary.csv")
  if (!is.null(stab)) { pr <- stab[stab$is_primary == "TRUE" | stab$is_primary == TRUE, ][1, ]
    add("Fe II bootstrap peak stability", "peak_stability",
      sprintf("%.0f%% of resamples in 2%% region", 100*pr$pct_in_reference),
      ">= 80% (descriptive)", if (pr$pct_in_reference >= 0.8) "pass" else "inconclusive",
      "data bootstrap; not independent data", "peak_stability.csv") }
  if (!is.null(orr)) add("Observed/Ritz consistency", "observed_ritz_replication",
    sprintf("all sources same region: %s", all(orr$same_reference_region)),
    "same predefined region", if (all(orr$same_reference_region)) "pass" else "inconclusive",
    "not fully independent (same transitions)", "observed_ritz_replication.csv")
  if (!is.null(ho)) add("Held-out predictive behaviour", "holdout_results",
    sprintf("%d/%d designs direction-consistent; median fixed-k p=%.3f",
            sum(ho$direction_consistent == TRUE | ho$direction_consistent == "TRUE"),
            nrow(ho), stats::median(ho$fixed_k_test_p)),
    "consistent on blocked holdout", "see result",
    "k locked from training; blocks not independent experiments", "holdout_results.csv")
  if (!is.null(cal)) { a05 <- cal[cal$nominal_alpha == 0.05, ][1, ]
    add("False-positive calibration", "null_calibration",
      sprintf("FPR@0.05 = %.3f [%.3f, %.3f]", a05$observed_fpr, a05$ci_lo, a05$ci_hi),
      "alpha within CI", if (isTRUE(a05$compatible) || a05$compatible == "TRUE") "pass" else "inconclusive",
      "wide CI at small calibration_n", "null_calibration.csv") }
  if (!is.null(mc)) { ho_gain <- mc$heldout_loglik_test[2] - mc$heldout_loglik_test[1]
    add("Model comparison (in-sample vs held-out)", "model_comparison",
    sprintf("in-sample AIC favours M1 by %.0f; held-out loglik gain=%.1f (%s)",
            max(mc$deltaAIC), ho_gain,
            if (ho_gain > 0) "transfers" else "does NOT transfer to held-out block"),
    "report both", "see result",
    "k from scan; AIC/BIC do not fully correct look-elsewhere; M0 not a physical model",
    "model_comparison.csv") }

  add("Independent experimental confirmation", "(none)", "not attempted", "n/a",
      "not established", "code-language change is not independent replication", "limitations")
  add("WCT physical mechanism / universal law", "(none)", "not attempted", "n/a",
      "not established", "out of scope; statistical pattern only", "limitations")

  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  write_table(tab, file.path(table_dir, "final_claim_matrix.csv"))
  tab
}

main <- function(argv = commandArgs(TRUE)) {
  args <- parse_render_args(argv)
  cfg <- default_audit_config(seed = args$seed, bootstrap_n = args[["bootstrap-n"]],
                              null_n = args[["null-n"]], calibration_n = args[["calibration-n"]],
                              injection_n = args[["injection-n"]], fast = args$fast)
  root <- audit_repo_root()
  table_dir <- file.path(root, args[["table-dir"]])
  fig_dir <- file.path(root, args[["figure-dir"]])
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  if (isTRUE(args$parallel)) configure_parallel(TRUE, args$workers)

  setup_rng(cfg$seed)
  # shared flags appended to every module's argv
  pflag <- c(if (isTRUE(args$parallel)) "--parallel" else character(0),
             if (isTRUE(args$fast)) "--fast" else character(0))
  bn <- as.character(cfg$bootstrap_n); nn <- as.character(cfg$null_n)
  cn <- as.character(cfg$calibration_n); injn <- as.character(cfg$injection_n)

  results <- list()
  results <- run_module("build_analysis_registry.R", character(0), results)
  results <- run_module("build_dataset_flow.R", character(0), results)
  results <- run_module("bootstrap_peak_uncertainty.R", c("--bootstrap-n", bn, pflag), results)
  results <- run_module("peak_stability.R", character(0), results)
  results <- run_module("build_effect_size_table.R", c("--null-n", nn, pflag), results)
  results <- run_module("run_bin_grid.R", c("--null-n", nn, pflag), results)
  sens_n <- if (isTRUE(args$fast)) 12L else min(cfg$null_n, 500L)
  results <- run_module("run_model_sensitivity.R", c("--null-n", as.character(sens_n), pflag), results)
  results <- run_module("model_comparison.R", character(0), results)
  results <- run_module("run_observed_ritz_replication.R", c("--null-n", nn, "--bootstrap-n", as.character(min(cfg$bootstrap_n, 1000L)), pflag), results)
  results <- run_module("run_holdout_replication.R", c("--null-n", as.character(min(cfg$null_n, 1000L)), pflag), results)
  fam_n <- if (isTRUE(args$fast)) 15L else min(cfg$null_n, 500L)
  results <- run_module("global_multiple_testing.R", c("--family-n", as.character(fam_n), pflag), results)
  results <- run_module("calibrate_false_positive_rate.R", c("--calibration-n", cn, "--null-n", nn, pflag), results)
  results <- run_module("run_injection_recovery.R", c("--injection-n", injn, "--null-n", as.character(min(cfg$null_n, 1000L)), pflag), results)

  # derived cross-cutting artifacts
  parity <- tryCatch(build_parity(file.path(root, args[["python-root"]]),
                                  file.path(root, args[["r-root"]]), table_dir, fig_dir),
                     error = function(e) { message("[parity] ", conditionMessage(e)); NULL })
  tryCatch(build_dashboard(fig_dir, table_dir), error = function(e) message("[dashboard] ", conditionMessage(e)))
  claim <- tryCatch(build_claim_matrix(table_dir), error = function(e) { message("[claims] ", conditionMessage(e)); NULL })

  # strict-mode gating
  if (isTRUE(args$strict)) {
    req_in <- vapply(c("fe_ion2_120", "fe_ion2_160", "fe_ion2_200"), function(d)
      file.exists(file.path(root, args[["r-root"]], d, "nist_summary.json")), logical(1))
    if (!all(req_in)) stop("[strict] required Fe II 120/160/200 R inputs missing")
    req_out <- c("significance_results.csv", "bin_stability_summary.csv",
                 "peak_confidence_intervals.csv", "final_claim_matrix.csv")
    miss <- req_out[!file.exists(file.path(table_dir, req_out))]
    if (length(miss) > 0L) stop(sprintf("[strict] required outputs missing: %s", paste(miss, collapse = ", ")))
  }

  # optional Quarto render
  report_path <- NA_character_
  if (isTRUE(args[["render-report"]])) {
    report_path <- render_report(file.path(root, args$report), root, args$fast)
  }

  print_summary(results, table_dir, fig_dir, report_path, args$fast)
  invisible(results)
}

render_report <- function(qmd, root, fast) {
  rendered_dir <- file.path(root, "reports/rendered")
  dir.create(rendered_dir, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(qmd)) { message("[report] qmd not found: ", qmd); return(NA_character_) }
  if (nzchar(Sys.which("quarto")) && requireNamespace("quarto", quietly = TRUE)) {
    out <- tryCatch({
      quarto::quarto_render(qmd, quiet = TRUE)
      html <- sub("\\.qmd$", ".html", qmd)
      dest <- file.path(rendered_dir, basename(html))
      if (file.exists(html) && normalizePath(html) != normalizePath(dest)) file.rename(html, dest)
      dest
    }, error = function(e) { message("[report] quarto render failed: ", conditionMessage(e)); NA_character_ })
    return(out)
  }
  message("[report] Quarto CLI not available; skipping HTML render (qmd is still committed).")
  NA_character_
}

print_summary <- function(results, table_dir, fig_dir, report_path, fast) {
  tabs <- list.files(table_dir, pattern = "\\.csv$")
  figs <- list.files(fig_dir, pattern = "\\.png$")
  ok <- vapply(results, function(r) r$ok, logical(1))
  cat("\n=====================================================\n")
  cat("NIST statistical audit - completion summary\n")
  cat("=====================================================\n")
  cat(sprintf("mode               : %s\n", if (fast) "FAST (development-only)" else "full-resolution"))
  cat(sprintf("modules run        : %d (%d ok, %d failed)\n", length(ok), sum(ok), sum(!ok)))
  if (any(!ok)) cat(sprintf("failed modules     : %s\n", paste(names(ok)[!ok], collapse = ", ")))
  cat(sprintf("tables generated   : %d\n", length(tabs)))
  cat(sprintf("figures generated  : %d\n", length(figs)))
  cat(sprintf("report             : %s\n", if (is.na(report_path)) "not rendered (see message above)" else report_path))
  cat(sprintf("git commit         : %s\n", git_commit()))
  for (n in names(results)) cat(sprintf("  %-38s %-4s %6.1fs\n", n,
                                          if (results[[n]]$ok) "ok" else "FAIL", results[[n]]$secs))
  if (fast) cat("\nNOTE: FAST mode outputs are development-only. Re-run at full resolution for final results.\n")
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("render_statistical_audit\\.R$", .invoked_file)) {
  main()
}
