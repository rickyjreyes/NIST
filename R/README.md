# R implementation — internals

This directory is an **independent R computational reproduction** of the
canonical Python reference scanner
`scripts/nist_wct_log_spectral_scan_FIXED.py`. It is a *reproduction of the
computation*, not an independent dataset and not an independent scientific
confirmation. No NIST endorsement, certification, or validation is claimed —
NIST is the data provider only. Python remains the currently established
canonical implementation; the R port tracks it and is validated by the parity
checks described below.

## Files

| File | Purpose |
|---|---|
| `nist_scan_lib.R` | Reusable functions (the numerical core). |
| `nist_wct_log_spectral_scan.R` | Command-line scanner (mirrors the Python CLI). |
| `run_feii_bin_stability.R` | Fe II 120/160/200 ladder → `tables_r/nist_master_results_r.csv`. |
| `make_verdict.R` | PASS/FAIL gate → `tables_r/nist_verdict_r.csv`. |
| `compare_python_r.R` | Numerical parity report between a Python and an R output dir. |

## Functions (`nist_scan_lib.R`)

- `clean_cell(x)` — strip NIST/Excel decoration (`="123.45"`, leading `=`,
  surrounding quotes/whitespace).
- `numeric_from_nist_cell(x)` — clean then extract the leading numeric token
  (`([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?`), mirroring the Python
  `to_float_series`.
- `read_nist_csv(path)` — robust CSV read (character columns, BOM-stripped
  names) with `clean_cell` applied to every cell.
- `find_column(df, candidates)` — normalized exact match, then substring
  containment, like the Python `find_column`/`norm_col`.
- `extract_wavenumber_cm(df)` — prefer a wavenumber column; otherwise convert a
  wavelength-in-nm column via `wavenumber_cm = 1e7 / wavelength_nm`.
- `clean_lines(df, ion)` — optional `sp_num` ion filter, keep finite positive
  wavenumbers, `ell = log(wavenumber_cm)`, sort, drop duplicate wavenumbers.
- `np_histogram_uniform(a, bins, lo, hi)` — replicate NumPy's uniform-bin
  histogram fast path, including the floating-point edge corrections and the
  closed right edge on the last bin.
- `build_binned(lines, bins, ell_min, ell_max)` — bin counts on `ell`, return
  centers/counts/edges.
- `gaussian_kernel1d(sigma, truncate)` / `gaussian_filter_nearest(y, sigma)` —
  reproduce `scipy.ndimage.gaussian_filter1d(..., mode="nearest")`.
- `poisson_deviance(y, mu)` — `2 sum(mu - y + y log(y/mu))`, zero-count safe,
  `mu` floored at `EPS = 1e-12`.
- `design_poly(ell, degree)` — centered, **population-std** (ddof = 0)
  standardized polynomial design matrix.
- `fit_poisson_loglinear(y, baseline, X)` — IRLS for
  `log(mu) = log(B) + X beta`.
- `scan_k(ell, y, baseline, k_grid, degree)` — base fit and per-`k`
  cos/sin harmonic fits; returns the scan table, the best row, and `mu0`.
- `null_scan(...)` — parametric Poisson bootstrap with the full look-elsewhere
  scan per replicate.
- `branch_report(k_best, delta_ell)` — `n_obs` and nearest target windings.
- `peak_ratios(scan)` — scan-curve peaks + pairwise rational approximations
  (denominator ≤ 64) via exact `gmp` arithmetic.

## Numerical conventions (chosen for bit-level parity with Python)

1. **Histogram** — NumPy uniform-bin fast path including the `tmp_a < edge`
   decrement / `tmp_a >= next_edge` increment corrections and the inclusive
   right edge on the final bin. `range` defaults to the data `[min, max]`.
2. **Gaussian baseline** — SciPy `gaussian_filter1d`, `mode="nearest"`,
   `truncate = 4.0`, radius `lw = floor(truncate*sigma + 0.5)` (= 24 at
   `sigma = 6`). The order-0 kernel is symmetric (correlation = convolution);
   edges are extended by replicating the nearest sample. Baseline floored at
   `EPS`.
3. **Design matrix** — `z = (ell - mean(ell))`, divided by the **population**
   standard deviation `sqrt(mean(z^2))` (NumPy `np.std` default `ddof = 0`),
   *not* R's sample `sd()`.
4. **IRLS** — coefficient init `beta[1] = log(sum(y)/sum(B))`, `eta` clipped to
   `[-10, 10]`, ridge `1e-8`, max 60 iterations, convergence
   `max|Δbeta| < 1e-8`, SVD minimum-norm least-squares fallback for singular
   systems, `EPS = 1e-12` floors.
5. **Best `k`** — first maximiser of `deltaD` (strict `>`), matching the Python
   tie-break.
6. **`delta_ell`** — range of **bin centers** (`max(ell) - min(ell)`), as in the
   Python `main`, used for `n_obs = k_best * delta_ell / (2π)`.
7. **Peak ratios** — `Fraction(float).limit_denominator(64)` reproduced exactly
   with `gmp` big rationals (`double_to_bigq` gives the double's exact value).

## Output schema

Each output directory contains the same files and column/key names as the
Python reference:

- `nist_lines_clean.csv` — cleaned line list plus `wavenumber_cm`, `ell`,
  `wavenumber_source`.
- `nist_binned_spectrum.csv` — `ell, count, edge_lo, edge_hi, baseline`.
- `nist_scan_curve.csv` — `k, deltaD, D_base, D_harmonic, amplitude, phase`.
- `nist_null.csv` — `null_max_deltaD` (only when `--null-n > 0`).
- `nist_peak_ratios.csv` — `k_low, k_high, ratio, rational, rational_error`.
- `nist_summary.json` — same key structure (`csv, ion_filter, n_raw_rows,
  n_unique_lines, wavenumber_min_cm, wavenumber_max_cm, ell_min, ell_max,
  delta_ell, best{...}, scan_null_p, tail_count_ge, null_n, branch_report{...},
  top_peak_ratios[...]`). The Python `engine` block (CuPy-only) is intentionally
  absent.
- `nist_spectrum_fit.png`, `nist_scan_curve.png`, `nist_null.png`.

## Known Python/R RNG difference (bootstrap only)

The deterministic pipeline (read → clean → bin → baseline → IRLS → deviance
scan) is reproduced to machine precision (relative differences ≤ ~1e-11; see
`R/compare_python_r.R`). The **parametric Poisson bootstrap is the one place R
and Python cannot match draw-for-draw**: R's `rpois` and NumPy's
`Generator.poisson(seed=...)` use different algorithms and bit streams. Hence:

- Deterministic, non-bootstrap quantities match Python closely.
- Repeated R runs with the same `--seed` are byte-identical to each other.
- Bootstrap **conclusions** are statistically consistent (for Fe II the null
  produces zero exceedances of the observed `deltaD`, so the empirical p-value
  is `1/(null_n+1)` in both implementations), but the individual
  `null_max_deltaD` values differ between R and Python.

Because of this, `compare_python_r.R` deliberately does **not** compare
`nist_null.csv` or `scan_null_p` numerically.

---

## Complete statistical audit

The scripts above are the **canonical reference scanner** and are unchanged. The
files below add a full, reproducible **statistical audit** *around* that scanner
(GitHub issue #3). The audit never modifies the numerical meaning of
`nist_scan_lib.R`, `nist_wct_log_spectral_scan.R`, or the Python reference; it
re-uses the canonical pipeline (`run_scan_analysis()` wraps the unchanged
`build_binned` → Gaussian baseline → IRLS Poisson fit → deviance `scan_k`).

### Installation

```powershell
Rscript .\R\check_audit_dependencies.R
```

This lists required vs optional packages, prints a single `install.packages(...)`
command for anything missing, and returns non-zero if a **required** package is
absent. Nothing is installed silently. Optional packages (`ggridges`, `ggrepel`,
`future`/`future.apply`, `gt`, `quarto`) degrade gracefully when missing.
Rendering the Quarto report additionally needs the **Quarto CLI** on `PATH`.

### Development run (small resample counts, fast)

```powershell
Rscript .\R\render_statistical_audit.R `
  --fast `
  --force `
  --render-report true
```

`--fast` uses development-only resample counts and labels outputs accordingly.
Use it to verify wiring, tables and figures before committing to a long run.

### Final run (full resolution — long-running)

```powershell
Rscript .\R\render_statistical_audit.R `
  --bootstrap-n 5000 `
  --null-n 5000 `
  --calibration-n 10000 `
  --injection-n 2000 `
  --seed 20260517 `
  --parallel true `
  --force `
  --strict `
  --render-report true
```

**Runtime implications (no promises):** the pure-R scan over the 2500-point
`k`-grid costs ~1 s per scan. Nulls, bootstraps, calibration and
injection–recovery each launch thousands of scans, so the full-resolution run is
expensive (hours, dominated by `--calibration-n` and `--injection-n`). Use
`--parallel true` with deterministic L'Ecuyer-CMRG streams (results are
independent of worker count). Run the final command yourself on adequate
hardware; do not assume it has already been executed.

### Output structure

```
tables_r/statistical_audit/    machine-readable CSV results
outputs_r/statistical_audit/   analysis_registry.json + per-analysis outputs
figures_r/statistical_audit/   300-dpi PNG + SVG/PDF figures (fig01..fig16)
reports/nist_statistical_audit.qmd   Quarto source
reports/rendered/              rendered HTML (when Quarto CLI is present)
```

Module → output map:

| Module | Key outputs |
|---|---|
| `check_audit_dependencies.R` | dependency report (stdout) |
| `build_analysis_registry.R` | `analysis_registry.csv/json` |
| `build_dataset_flow.R` | `dataset_flow*.csv`, `fig01` |
| `bootstrap_peak_uncertainty.R` | `peak_estimates.csv`, `bootstrap_peak_draws.csv`, `peak_confidence_intervals.csv`, `fig05` |
| `peak_stability.R` | `peak_stability.csv`, `fig04` |
| `build_effect_size_table.R` | `significance_results.csv`, `effect_sizes.csv`, `fig06`, `fig07` |
| `run_bin_grid.R` | `bin_grid_results.csv`, `bin_stability_summary.csv`, `fig03` + heatmaps |
| `run_model_sensitivity.R` | `specification_results.csv`, `fig10` + heatmap |
| `model_comparison.R` | `model_comparison.csv`, `fig11` |
| `run_observed_ritz_replication.R` | `observed_ritz_replication.csv`, `fig08` |
| `run_holdout_replication.R` | `holdout_results.csv`, `fig09` |
| `global_multiple_testing.R` | `multiple_testing.csv`, `family_max_null.csv`, `fig12` |
| `calibrate_false_positive_rate.R` | `null_calibration.csv`, `fig13` |
| `run_injection_recovery.R` | `injection_recovery.csv`, `fig14` + bias |
| `render_statistical_audit.R` | `python_r_parity.csv`, `final_claim_matrix.csv`, `fig15`, `fig16`, report |

### Statistical definitions and conventions

- **Frequency convention.** The harmonic model uses `cos(k·ell)`, `sin(k·ell)`
  with `ell = ln(wavenumber/cm⁻¹)`, so `k` is an **angular** frequency. The
  log-period is `Δlog x = 2π/k` and the multiplicative scale ratio is
  `exp(2π/k)`. `n_obs = k·Δell/(2π)`. Cyclic and angular frequency are never
  silently interchanged.
- **Empirical p-value.** `p = (r + 1)/(B + 1)`; never 0. With zero exceedances
  the corrected estimate equals the **resolution floor** `1/(B + 1)`; this is a
  resolution-limited bound, **not** an exact `p < 1/(B+1)`.
- **Pointwise vs scan-global p.** Pointwise p is at the fixed, prespecified `k`;
  scan-global p uses the **maximum** statistic across the whole `k`-grid in every
  null realisation (look-elsewhere correction). Both are reported.
- **Multiplicity family.** Read from the declared registry (bin grid, baseline
  σ × degree, source fields, neighbouring ions). Corrections: Benjamini–Hochberg
  FDR, Holm, Bonferroni, and a family-wise max-statistic. Because the analyses
  share one NIST line list they are **not** independent; this is stated in the
  output and the corrections are therefore approximate.
- **Peak-region tolerance.** Defined **before** reading bootstrap results.
  Primary confirmatory definition: relative tolerance **2 %**; sensitivity also
  at 1 %, 2 %, 5 %, plus an absolute-tolerance variant. Stability classes are
  descriptive: high ≥ 80 %, moderate 50–80 %, low < 50 %.
- **Observed/Ritz distinction.** Observed-, Ritz- and direct-wavenumber datasets
  are kept strictly separate; this is a **measurement-representation comparison**,
  not fully independent replication (observed and Ritz often describe the same
  transitions).
- **Holdout protocol.** Blocked splits in log-wavenumber space (lower/upper,
  alternating blocks, repeated blocked K-fold). `k` is estimated on training and
  **locked** before the test block is evaluated; the confirmatory test never
  rescans `k`. Exploratory rescans are reported separately and clearly labelled.
- **Calibration protocol.** Synthetic datasets are generated under the fitted
  smooth null; the full scan-global selection rule is applied; the observed
  false-positive rate is compared to nominal α ∈ {0.10, 0.05, 0.01} with exact
  binomial CIs. Exact calibration is not claimed when the Monte Carlo CI is wide.
- **Injection protocol.** A known log-periodic mode is injected on the fitted
  baseline over an amplitude grid {0, …, 0.10} at several frequencies (Fe
  reference, low/high controls, off-grid). Detection-anywhere, correct-region and
  globally-significant correct-region are distinguished.

### Interpretation limitations

This audit reports **computational and statistical** evidence only. It does
**not** establish independent experimental replication, an independent dataset, a
WCT physical mechanism, a universal atomic law, causal interpretation, or NIST
endorsement, and it cannot resolve significance below `1/(B+1)`. A change of
programming language is **not** independent replication. See
`tables_r/statistical_audit/final_claim_matrix.csv` and the report's executive
summary for the explicit supported-vs-not-established split.
