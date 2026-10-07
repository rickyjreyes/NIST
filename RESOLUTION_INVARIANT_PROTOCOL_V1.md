# Resolution-invariant Fe II robustness protocol v1

## Status

**FROZEN BEFORE THE FINAL PROSPECTIVE RUN.**

This protocol is a prospective follow-up to the already observed fixed-`sigma=6` resolution sensitivity. It is not allowed to overwrite, rescue, or reinterpret the canonical 120/160/200-bin audit result or the historical fixed-sigma 60–240-bin result.

The historical fixed-bin-sigma dense-grid result remains reported as approximately 70% reference-region retention. The new test asks a different, predeclared question: whether the Fe II branch remains stable when the smoothing scale is held approximately fixed in the analyzed `ell = ln(wavenumber)` coordinate rather than fixed in histogram-bin units.

## Frozen target and data transformation

- Species: Fe II
- Source: NIST frozen `data/Fe_lines.csv`
- Coordinate: `ell = ln(wavenumber / cm^-1)`
- Primary target: `k = 31.3265306122449`
- Target tolerance: ±2% relative
- Scan range: `k = 0.5 ... 80`
- Grid: 2,500 k values
- Smooth model degree: 1

The target is inherited from the previously frozen final-audit configuration; it is not re-estimated from the prospective dense-grid result.

## Frozen bin grid

The declared resolutions are:

`60, 80, 100, 120, 140, 160, 180, 200, 220, 240`

Every declared resolution must be reported, including failures and branch swaps.

## Frozen smoothing convention

The primary prospective convention is

```text
sigma_bins(B) = 6 * B / 160
```

This is frozen because a Gaussian width expressed in bin units changes its approximate width in `ell` when the number of bins changes. Scaling sigma linearly with bin count keeps the smoothing width approximately invariant in `ell` relative to the canonical 160-bin, sigma=6 analysis.

This convention is justified by coordinate-scale invariance, not by whichever convention gives the preferred peak after the run.

For comparison only, the historical `sigma_bins = 6` result must remain in the report. It is not replaced.

## Peak tracking

At each resolution the runner must record:

1. the global winning peak;
2. the top five distinct local maxima;
3. the rank of the frozen Fe target branch;
4. `deltaD` at the frozen target;
5. `deltaD(target) / deltaD(best)`;
6. whether the winner is within the ±2% target region.

This prevents a winner swap between persistent branches from being mislabeled as disappearance of the target branch.

## Null calibration

For each declared bin count:

- generate 5,000 Poisson null replicates;
- use the same frozen coordinate, k grid, degree, and resolution-invariant smoothing rule;
- re-apply the smoothing/null construction inside each replicate;
- use the scan maximum as the global statistic;
- report `(r+1)/(B+1)` rather than `p=0`;
- quantify the frequency with which null data produce a target-region locked peak of comparable strength.

The smoothing rule fails the artifact-control requirement if the zero-signal target-lock rate exceeds 5%.

## Injection/recovery

At each declared bin count, run 2,000 replicates at amplitude ratios:

- 0.0 × observed target amplitude;
- 0.5 × observed target amplitude;
- 1.0 × observed target amplitude.

The injected frequency is always the frozen Fe target. Do not reselect it from injected outcomes.

Report detection probability, correct-location probability, significant-correct-location probability, recovered k bias, and amplitude bias for every declared bin count/amplitude cell.

## Decision rule

The prospective resolution-invariant test is classified as `PASS` only if all of the following hold:

1. at least 80% of declared bin counts select the frozen target region as the global winner;
2. the frozen target appears among the top five local maxima in at least 80% of declared bin counts;
3. the zero-signal false target-lock rate is ≤5%;
4. the 1.0× observed-amplitude significant-correct-location recovery rate is ≥80%;
5. no implementation or null-calibration failure invalidates the comparison.

Otherwise classify the prospective test as `FAIL` or `INCONCLUSIVE` according to the declared failure mode. Do not alter the bin grid, tolerance, target, smoothing rule, or thresholds after the run.

## Required outputs

- full per-bin observed results;
- top-five local-peak table;
- per-bin null summary;
- per-bin injection/recovery summary;
- explicit comparison with the historical fixed-sigma 70% result;
- final `PASS / FAIL / INCONCLUSIVE` classification;
- exact command, seed, commit, source hashes, and protocol hash.

The machine-readable source of truth is `config/resolution_invariant_protocol_v1.json`.
