# Resolution-mode follow-up diagnostics

This follow-up exists because the frozen adversarial bin-grid test can reveal a resolution-dependent winner without explaining why it changes. It is intentionally separate from the final audit and cannot overwrite or rescue a failed robustness claim.

The original `R/run_bin_grid.R` result remains authoritative for the declared bin-stability verdict. If fewer than 80% of the predefined bin counts select the Fe II reference region, the final claim matrix reports that robustness claim as `fail` even if a later diagnostic explains the transition.

## Why this follow-up is needed

The canonical pipeline bins `ell = log(wavenumber)` and then applies a Gaussian baseline with `baseline_sigma = 6` measured in **bins**. Changing the number of histogram bins therefore changes both resolution and the approximate smoothing width in physical `ell` space.

The follow-up separates those effects by comparing:

- `fixed_sigma_bins`: the original convention, sigma = 6 at every bin count;
- `matched_sigma_ell`: sigma scaled as

```text
sigma_bins(B) = 6 * B / 160
```

so the approximate smoothing width in `ell` is held fixed relative to the 160-bin primary analysis.

For reference, this gives sigma values 2.25, 3.00, 3.75, 4.50, 5.25, 6.00, 6.75, 7.50, 8.25, and 9.00 for 60 through 240 bins in steps of 20.

## Two tracked frequencies

The diagnostic records the full winning mode and also the strength of two explicit targets at every resolution:

1. the canonical Fe II primary frequency selected by the 160-bin analysis;
2. `k = 9.602325620315224`, frozen previously in the independent GWTC repository.

The GWTC value is taken from `rickyjreyes/GWTC/tables/gwtc_frozen_mode.json` and the corresponding frozen holdout result. The GWTC frequency predates this NIST resolution follow-up. However, noticing that the coarse-bin NIST branch lies near that value is **post-hoc**, so the NIST/GWTC comparison is reported only as exploratory cross-domain follow-up evidence.

The dense map stores the top five distinct local maxima rather than only the global winner. This distinguishes a true disappearance of one mode from a winner swap between two simultaneously strong modes.

## Fixed external-frequency test

The script also tests the exact frozen GWTC frequency without re-optimizing `k` in NIST. At 160 bins it evaluates:

- Fe II wavenumber source;
- Fe II observed-wavelength-derived source;
- Fe II Ritz-wavelength-derived source;
- Cr II, Mn II, Co II, Ni II, and Ti II wavenumber sources.

For each dataset it reports the fixed-frequency deviance improvement and two empirical Poisson-null calibrations:

- observed smoothing baseline held fixed;
- smoothing baseline re-estimated in every null replicate.

The reported conservative fixed-frequency p-value is the larger of those two p-values. Every empirical p-value uses the plus-one correction `(r + 1)/(B + 1)`.

Because this test was motivated after seeing the NIST resolution behavior, it is **not** promoted to the original confirmatory family and is not used to repair the final audit verdict.

## Run

Run the final audit first and preserve its outputs. Then pull the latest branch and execute:

```powershell
Rscript R/run_resolution_mode_diagnostics.R --parallel true --fixed-null-n 5000
```

The dense observed-data resolution map defaults to 40 through 300 bins in steps of 5. The bounds can be changed explicitly:

```powershell
Rscript R/run_resolution_mode_diagnostics.R --parallel true --min-bin 40 --max-bin 300 --step-bin 5 --fixed-null-n 5000
```

A development-only smoke run is available with:

```powershell
Rscript R/run_resolution_mode_diagnostics.R --fast
```

`--fast` is not final evidence.

## Outputs

All outputs are written alongside the statistical-audit products:

- `tables_r/statistical_audit/resolution_mode_dense.csv`
  - winner at each bin count under both smoothing conventions;
  - approximate physical smoothing width;
  - deviance improvement at the frozen GWTC target and Fe primary target;
  - each target's strength relative to the winning peak.
- `tables_r/statistical_audit/resolution_mode_top_peaks.csv`
  - top five distinct local maxima at every resolution;
  - branch labels for the frozen GWTC target region, Fe primary region, and other modes.
- `tables_r/statistical_audit/external_gwtc_k_fixed_tests.csv`
  - exact fixed-frequency cross-ion/source tests;
  - fixed-baseline, refit-baseline, and conservative empirical p-values.
- `figures_r/statistical_audit/resolution_mode_transition.png`
  - selected frequency versus bin count under both smoothing conventions.
- `figures_r/statistical_audit/resolution_mode_strength.png`
  - strength of the two tracked modes relative to the global winner.

## Interpretation rules

The follow-up is designed to discriminate among several possibilities without protecting a preferred result:

- If the `k ~ 31.3` branch returns at coarse bin counts after matching physical smoothing width, the original coarse-bin failure is likely entangled with the smoothing convention.
- If `k ~ 9.6` remains dominant under matched smoothing, it is a genuine resolution-robust competing mode within this pipeline.
- If both tracked modes remain strong while the winner switches, the data support a multimodal spectrum more strongly than a single universally dominant frequency.
- If the frozen GWTC target fails under baseline refitting or across neighboring ions/sources, the cross-domain coincidence weakens.
- If the frozen GWTC target survives conservative fixed-frequency tests, that is noteworthy follow-up evidence but still does not establish WCT, a universal law, or independent experimental confirmation.
