# Final adversarial statistical audit

The repository contains a broad statistical-audit suite in `R/`.  Development runs remain useful for checking code paths quickly, but their small Monte Carlo budgets are not final evidence.

`R/run_final_adversarial_audit.R` is the strict final entry point.  It was added specifically to prevent development-resolution outputs or hidden per-module caps from being mistaken for the final audit.

## Run

```powershell
Rscript R/run_final_adversarial_audit.R --parallel true --render-report true
```

Default final budgets are explicit:

| Stage | Replicates |
|---|---:|
| dataset bootstrap | 5,000 |
| primary scan-global null | 5,000 |
| observed/Ritz bootstrap | 2,000 |
| specification-grid null per specification | 1,000 |
| frozen blocked-holdout null per design | 5,000 |
| declared-family maximum-statistic null | 5,000 |
| outer null-calibration datasets | 10,000 |
| injection/recovery per amplitude/frequency | 2,000 |
| injection decision null | 5,000 |
| alternative null per model | 5,000 |

The runner aborts on a module failure by default.  The full budgets can be overridden explicitly with command-line options, but they are never silently reduced by the final runner.

## Adversarial requirements

The final runner executes and reports all of the following:

- **Look-elsewhere correction:** the maximum statistic over the complete declared `k` scan.
- **Multiplicity:** family maximum-statistic, BH-FDR, Holm, and Bonferroni over the declared analysis family.
- **Bootstrap uncertainty:** peak frequency, amplitude, phase, period/scale-ratio uncertainty and peak-region selection stability.
- **Null calibration:** synthetic-null false-positive calibration with Monte Carlo intervals.
- **Frozen holdouts:** `k` is selected only on the training block and locked before confirmatory test-block evaluation. Exploratory test rescans are kept separate.
- **Injection/recovery:** detection and correct-frequency recovery are reported separately.
- **Preprocessing/specification sensitivity:** predefined bin counts, Gaussian baseline widths, polynomial degrees, and wavenumber/observed/Ritz representations.
- **Alternative nulls:** fixed-baseline Poisson, baseline-refit Poisson, conditional multinomial, overdispersed negative binomial, and block-residual baseline-misspecification stress tests.
- **Failure preservation:** every verdict other than `pass` is copied into `failed_claims.csv`. A failed holdout or predictive comparison cannot be rescued by a strong in-sample p-value.

## Primary outputs

All tables live under `tables_r/statistical_audit/`.

- `calibrated_global_evidence.csv` — compact primary evidence dashboard.
- `alternative_null_results.csv` — scan-global p-value under each declared null.
- `multiple_testing.csv` — declared-family corrections.
- `null_calibration.csv` — observed false-positive rate versus nominal alpha.
- `holdout_results.csv` — locked-frequency blocked validation.
- `injection_recovery.csv` — power and localisation recovery.
- `specification_results.csv` — preprocessing/model multiverse.
- `final_adversarial_claim_matrix.csv` — explicit claim-by-claim verdicts.
- `failed_claims.csv` — all failed, mixed, inconclusive, and not-established claims.
- `final_adversarial_run_manifest.csv` — exact budgets, commit, timestamp, and run status.
- `final_adversarial_module_status.csv` — per-module success/failure and runtime.

The rendered report is `reports/rendered/nist_final_adversarial_audit.html` when Quarto is available.

## Interpreting global evidence

`calibrated_global_evidence.csv` exposes three distinct global quantities for the primary Fe II analysis:

1. the full-scan maximum-statistic p-value;
2. the declared-family maximum-statistic p-value;
3. the largest p-value found among the declared alternative null models.

It also reports `largest_reported_global_p`, defined as the maximum of those three.  This is deliberately conservative, but it is only a **diagnostic summary**.  It is not a mathematically new combined p-value and must not be described as one.

Every empirical Monte Carlo p-value uses

```text
p = (r + 1) / (B + 1)
```

so zero exceedances are reported at the resolution floor `1/(B+1)`, never as `p=0` or as an unsupported extrapolation beyond the simulation budget.

## Claim boundaries

This audit can evaluate the robustness of statistical structure in the declared NIST analysis.  It does **not** establish independent experimental confirmation, a unique physical null, a WCT mechanism, a universal atomic law, or NIST endorsement.
