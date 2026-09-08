# Final adversarial statistical audit

The repository already contains a broad statistical-audit suite in `R/`. Development runs remain useful for checking code paths quickly, but their small Monte Carlo budgets are not final evidence.

`R/run_final_adversarial_audit.R` is the strict final entry point. It prevents development-resolution outputs, hidden per-module caps, or partial failed runs from being mistaken for the final audit.

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

The runner aborts on any module failure. Budgets can be overridden explicitly, but the final runner never silently reduces them.

## Adversarial requirements

The final runner executes and reports all of the following:

- **Look-elsewhere correction:** maximum statistic over the complete declared `k` scan.
- **Full-family multiplicity:** the family now includes the Fe bin grid, every declared sigma/degree/bin preprocessing specification, observed/Ritz representations, and neighbouring ion-II controls. It reports family maximum-statistic, Holm, Bonferroni, BH-FDR, and dependence-robust BY-FDR corrections.
- **Bootstrap uncertainty:** peak frequency, amplitude, phase, period/scale-ratio uncertainty, and peak-region selection stability.
- **Null calibration:** synthetic-null false-positive calibration with Monte Carlo intervals.
- **Frozen holdouts:** `k` is selected only on the training block and locked before confirmatory test-block evaluation. Exploratory test rescans remain separate.
- **Injection/recovery:** any-peak detection and correct-frequency recovery are reported separately.
- **Preprocessing/specification sensitivity:** bin counts, Gaussian baseline widths, polynomial degrees, and wavenumber/observed/Ritz representations.
- **Alternative nulls:** fixed-baseline Poisson, baseline-refit Poisson, conditional multinomial, overdispersed negative binomial, and block-residual baseline-misspecification stress tests.
- **Failure preservation:** every verdict other than `pass` is copied to `failed_claims.csv`. A failed holdout or predictive comparison cannot be rescued by a strong in-sample result.

The alternative-null family is a **post-signal adversarial stress-test family**, not historical preregistration. The 80% stability threshold used in the final claim matrix is likewise an audit-declared descriptive robustness criterion, not a universal physical threshold.

## Dependence handling

Many Fe analyses reuse the same transition list, so they are dependent. The audit therefore reports several complementary corrections:

- **Holm and Bonferroni:** FWER control that does not require independence.
- **Benjamini–Yekutieli:** FDR control valid under arbitrary dependence.
- **Benjamini–Hochberg:** standard FDR reference, retained for comparison.
- **Family maximum statistic:** useful but explicitly approximate here because it couples marginal per-analysis parametric-null draws rather than simulating one fully joint synthetic line-list experiment.

This distinction is kept in the output rather than collapsing every correction into one number.

## Primary outputs

All tables live under `tables_r/statistical_audit/`.

- `calibrated_global_evidence.csv` — compact primary evidence dashboard.
- `alternative_null_results.csv` — scan-global p-value under each declared null.
- `multiple_testing.csv` — full declared-family corrections.
- `null_calibration.csv` — observed false-positive rate versus nominal alpha.
- `holdout_results.csv` — locked-frequency blocked validation.
- `injection_recovery.csv` — power and localisation recovery.
- `specification_results.csv` — preprocessing/model multiverse.
- `final_adversarial_claim_matrix.csv` — explicit claim-by-claim verdicts.
- `failed_claims.csv` — all failed, mixed, inconclusive, and not-established claims.
- `final_adversarial_run_manifest.csv` — exact budgets, commit, timestamp, and run settings.
- `final_adversarial_module_status.csv` — per-module success and runtime.

The rendered report is `reports/rendered/nist_final_adversarial_audit.html` when Quarto is available.

## Interpreting global evidence

`calibrated_global_evidence.csv` exposes distinct quantities for the primary Fe II analysis: the full-scan maximum-statistic p-value; family-max, Holm, Bonferroni, BH, and BY multiplicity adjustments; and the worst p-value among the declared alternative null models.

It also reports `largest_reported_global_p`, defined as the maximum across those reported global/multiplicity/alternative-null values. This is deliberately conservative but is only a **diagnostic summary**. It is not a mathematically new combined p-value and must not be described as one.

Every empirical Monte Carlo p-value uses

```text
p = (r + 1) / (B + 1)
```

so zero exceedances are reported at the resolution floor `1/(B+1)`, never as `p=0` or as an unsupported extrapolation beyond the simulation budget.

## Claim boundaries

This audit can evaluate the robustness of statistical structure in the declared NIST analysis. It does **not** establish independent experimental confirmation, a unique physical null, a WCT mechanism, a universal atomic law, or NIST endorsement.

## Resolution-mode follow-up is separate

A resolution-dependent winner discovered by the frozen bin-grid audit is investigated by the standalone `R/run_resolution_mode_diagnostics.R` module documented in `RESOLUTION_MODE_DIAGNOSTICS.md`.

That follow-up deliberately does **not** run inside `run_final_adversarial_audit.R` and cannot overwrite its claim matrix. It compares the original fixed-sigma-bin pipeline with an approximately fixed `ell`-space smoothing width, tracks both the Fe primary mode and the previously frozen GWTC `k = 9.602325620315224`, stores the top five local peaks across a dense resolution grid, and performs conservative fixed-frequency null tests across Fe/neighboring-ion sources.

Run it only after preserving the final audit outputs:

```powershell
Rscript R/run_resolution_mode_diagnostics.R --parallel true --fixed-null-n 5000
```

The NIST/GWTC frequency comparison is explicitly labeled post-hoc exploratory follow-up even though the GWTC target itself was frozen previously.
