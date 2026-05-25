# NIST Log-Cosine Atomic Spectral Scan

This repository reproduces the NIST atomic-spectroscopy component of the paper:

**Reproducible Log-Cosine Spectral Structure in Logarithmic Coordinates Across Four Independent Detector Datasets**  
Richard J. Reyes, 2026

The primary result is a stable log-cosine spectral-density mode in **NIST Fe II atomic spectral lines**.

The Fe II scan uses public NIST Atomic Spectra Database line data, converts wavenumber \(\sigma\) in cm\(^{-1}\) to the logarithmic coordinate

\[
\ell = \ln(\sigma),
\]

and fits a one-mode log-cosine residual model against a smooth Poisson baseline:

\[
\log \mu(\ell;k,a,b,c)
=
\log \mu_0(\ell)
+
c
+
a\cos(k\ell)
+
b\sin(k\ell).
\]

The scan statistic is the scan-max deviance improvement:

\[
\Delta \chi^2_\star = \max_k \left(D[y||\mu_0] - D[y||\mu(k)]\right).
\]

The active-domain winding coordinate is

\[
n = \frac{k \Delta \ell}{2\pi}.
\]

## Primary Fe II Result

The Fe II result is bin-stable across 120, 160, and 200 histogram bins:

| Species | Ion | Bins | Lines | \(k_\star\) | \(n_\star\) | \(\Delta\chi^2_\star\) | Null Result |
|---|---:|---:|---:|---:|---:|---:|---|
| Fe | II | 120 | 9,447 | 31.3265 | 10.717 | 355.7 | 0 / 5000 exceedances |
| Fe | II | 160 | 9,447 | 31.3265 | 10.740 | 315.7 | 0 / 5000 exceedances |
| Fe | II | 200 | 9,447 | 31.3265 | 10.753 | 259.2 | 0 / 5000 exceedances |

The key observation is that \(k_\star = 31.3265\) is unchanged across the binning ladder, while the active-domain coordinate remains stable near \(n \approx 10.7\).

This is the main NIST detection.

## Repository Structure

```text
nist-logcos/
  README.md
  data/
    Fe_lines.csv
    Ni_lines.csv
    Co_lines.csv
    Cr_lines.csv
    Mn_lines.csv
    Ti_lines.csv
    data_source_notes.md
  scripts/
    nist_wct_log_spectral_scan_FIXED.py
  outputs/
    fe_ion2_120/
    fe_ion2_160/
    fe_ion2_200/
    ni_ion2_160/
    co_ion2_160/
    cr_ion2_160/
    mn_ion2_160/
    ti_ion2_160/
  synthetic/
    make_synthetic_feii_like_lines.py
    run_synthetic_challenge.ps1
  tables/
    nist_master_results.csv