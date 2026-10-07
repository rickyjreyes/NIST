# Independent atomic-catalog replication protocol v1

## Status

**FROZEN BEFORE EXTERNAL OUTCOME ACCESS.**

The primary purpose is to test whether the NIST Fe II target transfers to a separate atomic-line catalog **without selecting a new target frequency from the replication data**.

## External catalog

The frozen source family is the Kurucz CD-ROM 23 / GFALL atomic-line catalog maintained through the Harvard-Smithsonian/CfA Kurucz resources. This is a separate catalog from the frozen NIST ASD CSV inputs.

Before any scientific parsing or test is run:

1. obtain the raw `gfall.dat` source file;
2. place it at `external/kurucz/gfall.dat`;
3. compute its SHA-256 without inspecting the Fe II outcome;
4. pin that SHA-256 in `config/independent_atomic_replication_v1.json`;
5. commit the hash before running the primary test.

A source whose bytes are not pinned first is not eligible for the primary replication classification.

## Frozen Fe II selection

Kurucz species code `26.01` identifies Fe II.

The primary extraction keeps only transitions satisfying all of the following:

- element code exactly `26.01` within parser precision;
- finite lower and upper level energies;
- lower energy >= 0;
- upper energy >= 0;
- transition wavenumber defined as `abs(E_upper - E_lower)` in cm^-1;
- wavenumber in the frozen NIST support interval `5721.4978 ... 49990.73 cm^-1`;
- exact duplicate wavenumbers removed after ascending sort, retaining the first occurrence.

Using the level-energy difference avoids introducing an air/vacuum wavelength-conversion convention into the primary coordinate.

The primary test is automatically `INCONCLUSIVE` if fewer than 500 lines remain after the frozen selection.

## Frozen primary target

The primary target is

```text
k = 31.3265306122449
```

This value comes from the already frozen NIST Fe II analysis. The external data may not move, widen, or redefine it.

The replication coordinate and model are:

- `ell = ln(wavenumber / cm^-1)`;
- 160 bins;
- Gaussian baseline sigma = 6 bins;
- degree-1 smooth model;
- primary statistic = Poisson deviance improvement at the **fixed** target k.

## External-catalog null

The source-specific null is generated from the smoothed external-catalog baseline.

For each of 5,000 null replicates:

1. draw Poisson counts from the external smooth baseline;
2. re-estimate the same sigma=6 smooth baseline on that replicate;
3. evaluate the deviance improvement at the fixed target only;
4. do not scan or reselect k for the primary p-value.

The empirical primary p-value is `(r+1)/(5000+1)`. Zero exceedances are never reported as `p=0`.

## Effect size

Report, in addition to the empirical p-value:

- observed fixed-target deviance improvement;
- target fitted amplitude;
- target fitted phase;
- null mean and standard deviation;
- null-standardized z-score `(observed - mean(null)) / sd(null)`;
- null 95th percentile.

## Classification

- `PASS`: source/hash/provenance and parser gates pass, at least 500 selected Fe II lines remain, and fixed-target empirical `p <= 0.05`.
- `FAIL`: source/hash/provenance and parser gates pass, at least 500 selected lines remain, and fixed-target empirical `p > 0.05`.
- `INCONCLUSIVE`: the source identity/hash/parser/selection is invalid or fewer than 500 lines remain.

A `FAIL` must be preserved. Do not switch catalogs, change the range, retune k, change the binning/smoothing, or widen a tolerance after seeing the result.

## Exploratory scan

Only after the fixed-target statistic/classification has been computed, an exploratory full scan from `k=0.5 ... 80` on 2,500 points may be generated. It must be labelled exploratory and cannot alter the primary classification.

## Claim boundary

A successful external catalog replication would materially strengthen evidence that the NIST feature is not peculiar to the frozen NIST CSVs. It would still be catalog-level replication, not by itself an independent laboratory experiment or validation of Wave Confinement Theory as the causal mechanism.

The machine-readable source of truth is `config/independent_atomic_replication_v1.json`.
