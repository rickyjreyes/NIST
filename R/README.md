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
