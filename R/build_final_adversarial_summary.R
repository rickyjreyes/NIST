#!/usr/bin/env Rscript
# build_final_adversarial_summary.R
# ---------------------------------------------------------------------------
# Cross-cutting final evidence summary for the adversarial NIST audit.
#
# This script does NOT recompute scans.  It consumes the machine-readable
# outputs from the audit modules and produces:
#   - calibrated_global_evidence.csv : one-row conservative evidence dashboard
#   - final_adversarial_claim_matrix.csv : explicit pass/fail/mixed/not-established
#   - failed_claims.csv : every non-passing claim retained verbatim
#
# A failed or mixed claim is never converted to a neutral "see result" label.
# The largest reported global p-value is exposed as a conservative diagnostic
# summary; it is NOT represented as a new mathematically-combined p-value.
# ---------------------------------------------------------------------------

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
if (!exists("write_table", mode = "function")) {
  source(file.path(if (length(.this) == 0L || is.na(.this)) "R" else dirname(.this),
                   "audit_utils.R"))
}

is_true <- function(x) isTRUE(x) || identical(toupper(as.character(x)), "TRUE")

read_audit_table <- function(root, name, required = TRUE) {
  p <- file.path(root, "tables_r/statistical_audit", name)
  if (!file.exists(p)) {
    if (required) stop("required audit table missing: ", p)
    return(NULL)
  }
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
}

primary_row <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  if ("analysis_id" %in% names(df)) {
    hit <- which(df$analysis_id == "fe_ion2_wn_bin160")
    if (length(hit) == 0L) hit <- grep("^fe_ion2_wn.*160$", df$analysis_id)
    if (length(hit) > 0L) return(df[hit[1], , drop = FALSE])
  }
  if (all(c("species", "source", "bins") %in% names(df))) {
    hit <- which(df$species == "Fe" & df$source == "wavenumber" & df$bins == 160)
    if (length(hit) > 0L) return(df[hit[1], , drop = FALSE])
  }
  df[1, , drop = FALSE]
}

claim_row <- function(claim, analysis, result, threshold, verdict, limitations, ref) {
  data.frame(claim = claim, analysis = analysis, result = result,
             threshold = threshold, verdict = verdict,
             limitations = limitations, output_reference = ref,
             stringsAsFactors = FALSE)
}

holdout_verdict <- function(p, direction, alpha = 0.05) {
  ok <- is.finite(p) & p <= alpha & direction
  if (length(ok) == 0L || sum(ok) == 0L) return("fail")
  if (all(ok)) return("pass")
  "mixed"
}

threshold_verdict <- function(value, threshold, higher_is_better = TRUE) {
  if (!is.finite(value)) return("inconclusive")
  good <- if (higher_is_better) value >= threshold else value <= threshold
  if (good) "pass" else "fail"
}

build_summary <- function(root) {
  sig <- read_audit_table(root, "significance_results.csv")
  mult <- read_audit_table(root, "multiple_testing.csv")
  alt <- read_audit_table(root, "alternative_null_results.csv")
  cal <- read_audit_table(root, "null_calibration.csv")
  ho <- read_audit_table(root, "holdout_results.csv")
  mc <- read_audit_table(root, "model_comparison.csv")
  bins <- read_audit_table(root, "bin_stability_summary.csv")
  peak <- read_audit_table(root, "peak_stability.csv")
  spec <- read_audit_table(root, "specification_results.csv")
  inj <- read_audit_table(root, "injection_recovery.csv")

  s <- primary_row(sig)
  m <- primary_row(mult)
  a <- alt[which.max(alt$scan_global_p), , drop = FALSE]
  c05 <- cal[which.min(abs(cal$nominal_alpha - 0.05)), , drop = FALSE]
  mc0 <- mc[mc$model == "M0_smooth_null", , drop = FALSE]
  mc1 <- mc[mc$model == "M1_smooth_plus_logperiodic", , drop = FALSE]
  if (nrow(mc0) == 0L || nrow(mc1) == 0L) stop("model comparison rows M0/M1 missing")
  heldout_gain <- mc1$heldout_loglik_test[1] - mc0$heldout_loglik_test[1]
  primary_peak <- peak[peak$is_primary %in% c(TRUE, "TRUE"), , drop = FALSE]
  if (nrow(primary_peak) == 0L) primary_peak <- peak[1, , drop = FALSE]
  spec_frac <- mean(spec$in_reference_region %in% c(TRUE, "TRUE"), na.rm = TRUE)
  inj0 <- inj[inj$freq_name == "fe_reference" & abs(inj$amplitude) < 1e-12, , drop = FALSE]
  if (nrow(inj0) == 0L) inj0 <- inj[which.min(abs(inj$amplitude)), , drop = FALSE]

  p_candidates <- c(s$global_p[1], m$family_max_p[1], a$scan_global_p[1])
  largest_p <- max(p_candidates, na.rm = TRUE)

  evidence <- data.frame(
    analysis_id = "fe_ion2_wn_bin160",
    observed_deltaD = s$global_statistic[1],
    scan_global_p = s$global_p[1],
    scan_global_tail_count = s$global_tail_count[1],
    scan_global_B = s$global_B[1],
    family_max_p = m$family_max_p[1],
    bh_fdr = m$bh_fdr[1],
    holm_p = m$holm_p[1],
    bonferroni_p = m$bonferroni_p[1],
    family_size = m$family_size[1],
    worst_alternative_null = a$null_model[1],
    worst_alternative_null_p = a$scan_global_p[1],
    worst_alternative_null_B = a$B[1],
    largest_reported_global_p = largest_p,
    largest_reported_global_p_note = "max(scan-global, family-max, worst alternative-null); conservative diagnostic, not a newly combined formal p-value",
    calibration_fpr_alpha_0_05 = c05$observed_fpr[1],
    calibration_ci_lo = c05$ci_lo[1],
    calibration_ci_hi = c05$ci_hi[1],
    calibration_compatible = is_true(c05$compatible[1]),
    holdout_designs = nrow(ho),
    holdout_sig_direction_consistent = sum(is.finite(ho$fixed_k_test_p) & ho$fixed_k_test_p <= 0.05 & ho$direction_consistent %in% c(TRUE, "TRUE")),
    holdout_median_fixed_k_p = stats::median(ho$fixed_k_test_p, na.rm = TRUE),
    heldout_model_loglik_gain_M1_minus_M0 = heldout_gain,
    bin_reference_fraction = bins$pct_in_reference[1],
    bootstrap_peak_reference_fraction = primary_peak$pct_in_reference[1],
    specification_reference_fraction = spec_frac,
    injection_null_detection_prob = inj0$detection_prob[1],
    injection_null_n = inj0$n_sims[1],
    stringsAsFactors = FALSE
  )

  claims <- list()
  add <- function(...) claims[[length(claims) + 1L]] <<- claim_row(...)

  add("Fe II full-scan look-elsewhere evidence", "significance_results",
      sprintf("scan-global p=%.4g (%d/%d exceedances)", s$global_p[1], s$global_tail_count[1], s$global_B[1]),
      "p <= 0.05", threshold_verdict(s$global_p[1], 0.05, FALSE),
      "empirical p is bounded by Monte Carlo resolution; zero exceedances are not an exact smaller p",
      "significance_results.csv")

  add("Fe II declared-family multiplicity correction", "global_multiple_testing",
      sprintf("family-max p=%.4g; BH=%.4g; Holm=%.4g; Bonferroni=%.4g; family=%d",
              m$family_max_p[1], m$bh_fdr[1], m$holm_p[1], m$bonferroni_p[1], m$family_size[1]),
      "family-max p <= 0.05", threshold_verdict(m$family_max_p[1], 0.05, FALSE),
      "analyses share underlying NIST line lists; family-max construction is still an approximation to a fully joint generative null",
      "multiple_testing.csv")

  bad_alt <- alt$null_model[alt$scan_global_p > 0.05]
  add("Fe II robustness to declared alternative null models", "alternative_null_results",
      sprintf("worst-case p=%.4g under %s%s", a$scan_global_p[1], a$null_model[1],
              if (length(bad_alt)) paste0("; p>0.05 under: ", paste(bad_alt, collapse = ", ")) else ""),
      "all declared null-model scan-global p <= 0.05",
      if (length(bad_alt) == 0L) "pass" else "fail",
      "alternative nulls are stress tests, not a claim that one generator is the unique physical null",
      "alternative_null_results.csv")

  add("Synthetic-null false-positive calibration", "null_calibration",
      sprintf("FPR@0.05=%.4f [%.4f, %.4f]", c05$observed_fpr[1], c05$ci_lo[1], c05$ci_hi[1]),
      "nominal 0.05 lies within calibration CI",
      if (is_true(c05$compatible[1])) "pass" else "fail",
      "calibration tests the declared synthetic-null mechanism, not every possible misspecification",
      "null_calibration.csv")

  hverd <- holdout_verdict(ho$fixed_k_test_p, ho$direction_consistent %in% c(TRUE, "TRUE"))
  add("Frozen blocked holdout replication", "run_holdout_replication",
      sprintf("%d/%d designs have p<=0.05 and positive locked-k direction; median p=%.4g",
              sum(is.finite(ho$fixed_k_test_p) & ho$fixed_k_test_p <= 0.05 & ho$direction_consistent %in% c(TRUE, "TRUE")),
              nrow(ho), stats::median(ho$fixed_k_test_p, na.rm = TRUE)),
      "strict: every declared blocked design p<=0.05 with positive direction",
      hverd,
      "k is frozen from each training block; blocks still come from one line list and are not independent experiments",
      "holdout_results.csv")

  add("Held-out predictive transfer of smooth-plus-periodic model", "model_comparison",
      sprintf("held-out log-likelihood gain M1-M0 = %.4f", heldout_gain),
      "gain > 0", if (heldout_gain > 0) "pass" else "fail",
      "this explicitly preserves the previously observed failure when the in-sample model gain does not transfer",
      "model_comparison.csv")

  add("Declared Fe II bin stability", "bin_stability_summary",
      sprintf("%.1f%% of predefined bins select the reference region", 100 * bins$pct_in_reference[1]),
      ">= 80% (predeclared descriptive threshold)",
      threshold_verdict(bins$pct_in_reference[1], 0.80, TRUE),
      "descriptive robustness threshold, not a universal physical law",
      "bin_stability_summary.csv")

  add("Bootstrap peak-region stability", "peak_stability",
      sprintf("%.1f%% of primary bootstrap resamples select the 2%% reference region", 100 * primary_peak$pct_in_reference[1]),
      ">= 80% (predeclared descriptive threshold)",
      threshold_verdict(primary_peak$pct_in_reference[1], 0.80, TRUE),
      "resampling the same dataset does not constitute independent replication",
      "peak_stability.csv")

  add("Specification multiverse stability", "run_model_sensitivity",
      sprintf("%.1f%% of declared specifications select the reference region", 100 * spec_frac),
      ">= 80% (adversarial descriptive threshold)",
      threshold_verdict(spec_frac, 0.80, TRUE),
      "the threshold is a robustness convention; full table remains primary evidence",
      "specification_results.csv")

  add("Injection null type-I behaviour", "run_injection_recovery",
      sprintf("A=0 detection probability=%.4f (%d simulations)", inj0$detection_prob[1], inj0$n_sims[1]),
      "<= 0.05", threshold_verdict(inj0$detection_prob[1], 0.05, FALSE),
      "finite simulation uncertainty applies; calibration table is the dedicated type-I assessment",
      "injection_recovery.csv")

  add("Independent experimental confirmation", "none",
      "not attempted", "independent dataset/experiment required", "not established",
      "programming-language parity, alternate preprocessing, and same-line-list holdouts are not independent experiments",
      "limitations")

  add("WCT physical mechanism or universal atomic law", "none",
      "not tested by this statistical audit", "physical mechanism evidence required", "not established",
      "statistical structure alone cannot identify a causal WCT mechanism",
      "limitations")

  claim_tab <- do.call(rbind, claims)
  nonpassing <- claim_tab[!(claim_tab$verdict %in% c("pass")), , drop = FALSE]
  list(evidence = evidence, claims = claim_tab, nonpassing = nonpassing)
}

main <- function(argv = commandArgs(TRUE)) {
  root <- audit_repo_root()
  out <- build_summary(root)
  td <- file.path(root, "tables_r/statistical_audit")
  write_table(out$evidence, file.path(td, "calibrated_global_evidence.csv"))
  write_table(out$claims, file.path(td, "final_adversarial_claim_matrix.csv"))
  write_table(out$nonpassing, file.path(td, "failed_claims.csv"))
  cat(sprintf("[final_summary] largest reported global p=%.4g; %d/%d claims non-passing\n",
              out$evidence$largest_reported_global_p[1], nrow(out$nonpassing), nrow(out$claims)))
  invisible(out)
}

.invoked_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
if (length(.invoked_file) > 0L && grepl("build_final_adversarial_summary\\.R$", .invoked_file)) {
  main()
}
