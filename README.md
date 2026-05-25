# NIST Atomic Spectra Log-Cosine Spectral Scan

This repository contains a reproducible log-cosine spectral-density scan over
**public NIST Atomic Spectra Database** line lists for transition metals
(Fe, Ni, Co, Cr, Mn, Ti). The primary result is a stable Fe II log-cosine
spectral-density mode in logarithmic wavenumber coordinates.

> This project uses public CSV exports from the NIST Atomic Spectra
> Database. It is **not** "NIST certified", endorsed, or validated by NIST.
> NIST is the data provider only.

## What the scanner does

For each line list (e.g. `data/Fe_lines.csv`) the scanner reads
wavenumber `wn(cm-1)` (handling NIST/Excel-style cells such as `="49998.80"`),
optionally filters by ionization stage `sp_num`, converts to the logarithmic
coordinate

```
ell = ln(wavenumber_cm^-1)
```

bins line counts on `ell`, and fits a one-mode Poisson log-linear model:

```
log mu(ell; k, a, b, c) = log mu0(ell) + c + a cos(k ell) + b sin(k ell)
```

against a Gaussian-smoothed baseline `mu0`. The scan statistic is the
deviance improvement

```
DeltaChi^2_star = max_k ( D[y || mu0] - D[y || mu(k)] )
```

over a `k` grid. The active-domain winding is

```
n = k * Delta_ell / (2 pi)
```

where `Delta_ell = ell_max - ell_min`. The null is a parametric Poisson
bootstrap from `mu0`, scanned with the same `k` grid; the empirical p-value is

```
p_hat = (1 + tail_count_ge) / (1 + null_n)
```

## Primary result (Fe II, ion = 2, n_unique_lines = 9447)

| Bins | k_best              | n_obs      | DeltaD       | null (tail/N) | p_hat                |
|---:|---|---:|---:|---|---|
| 120  | 31.3265306122449    | 10.717134  | 355.7443     | 0 / 5000      | 1/5001 ≈ 0.00019996  |
| 160  | 31.3265306122449    | 10.739649  | 315.7277     | 0 / 5000      | 1/5001 ≈ 0.00019996  |
| 200  | 31.3265306122449    | 10.753158  | 259.1788     | 0 / 5000      | 1/5001 ≈ 0.00019996  |

The headline observation is **bin-stability of `k_best = 31.3265306122449`
across 120 / 160 / 200 histogram bins**, with `n_obs` stable near `10.7` and
zero null exceedances at 5000 nulls per bin setting. The PASS/FAIL gate
encoded in `scripts/make_verdict.py` enforces this.

## Repository layout

```
README.md
RESULTS.md
requirements.txt
data/
  Fe_lines.csv  Ni_lines.csv  Co_lines.csv
  Cr_lines.csv  Mn_lines.csv  Ti_lines.csv
scripts/
  nist_wct_log_spectral_scan_FIXED.py    # canonical CPU scanner
  nist_wct_log_spectral_scan.py          # legacy one-shot Fe splitter
  nist_wct_log_spectral_scan_CUPY.py     # optional GPU mirror (CuPy)
  nist_batch_run_neighbors.py            # neighbor-ion batch runner
  run_feii_bin_stability.py              # Fe II 120/160/200 ladder + master CSV
  make_verdict.py                        # PASS/FAIL from master CSV
outputs/
  fe_ion2_120/  fe_ion2_160/  fe_ion2_200/   # canonical Fe II runs
  batch_neighbors/                            # canonical neighbor batch (written on demand)
tables/
  nist_master_results.csv
  nist_verdict.csv
tests/
  test_nist_scanner.py
```

### CPU canonical script vs optional CuPy script

- `scripts/nist_wct_log_spectral_scan_FIXED.py` is the **canonical** scanner.
  It is pure NumPy/SciPy and produces the values in the table above. The
  PASS/FAIL verdict is built off this script's output.
- `scripts/nist_wct_log_spectral_scan_CUPY.py` is an **optional** GPU mirror
  using CuPy. It is faster for large `null_n` but is not required for
  reproduction. If CuPy is not installed, ignore this script.

## Commands

Quick single-bin run (Fe II, bins=160, 500-null preview):

```bash
python scripts/nist_wct_log_spectral_scan_FIXED.py \
    --csv data/Fe_lines.csv --ion 2 --bins 160 \
    --null-n 500 --min-lines 100 \
    --out-dir outputs/fe_ion2_160
```

Bin-stability ladder (quick, 500 nulls per bin):

```bash
python scripts/run_feii_bin_stability.py --null-n 500
```

Full run (5000 nulls per bin) plus verdict:

```bash
python scripts/run_feii_bin_stability.py --null-n 5000
python scripts/make_verdict.py
```

Tests:

```bash
pytest -q
```

Neighbor batch (Ni / Co / Cr / Mn / Ti, ion=2) — exploratory:

```bash
python scripts/nist_batch_run_neighbors.py --preview-null 500 --no-promote
```

## Canonical vs legacy outputs

- The canonical Fe II outputs live in `outputs/fe_ion2_{120,160,200}/`.
  Only these (or whatever rows are listed in `tables/nist_master_results.csv`)
  feed the PASS/FAIL verdict.
- Older / exploratory runs at the repo root (for example
  `outputs_fe_ion2_120_cupy_full/`, `outputs_ni_ion2_160_preview/`,
  `outputs_nist_wct/`, etc.) are **preserved as provenance**. They are
  legacy, non-canonical, and are intentionally **not used** by the verdict
  unless they are explicitly referenced in `tables/nist_master_results.csv`.
  Do not delete them.

## Limitations

- This is a **line-density** scan, not a flux-spectrum scan.
- The data are public NIST Atomic Spectra Database CSV exports; line
  completeness varies by species and is not uniform across `ell`.
- The null is a parametric Poisson bootstrap from a Gaussian-smoothed
  baseline. Baseline-sigma sensitivity is documented and should be audited
  for any extension beyond Fe II.
- This repository makes no claim of NIST endorsement, certification, or
  validation.
