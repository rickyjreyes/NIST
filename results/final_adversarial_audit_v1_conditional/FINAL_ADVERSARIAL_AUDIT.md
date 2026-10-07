# FINAL ADVERSARIAL AUDIT

## Verdict at a glance

1. **Does the Fe II result reproduce from the frozen raw inputs?** Yes.
2. **Does it survive bootstrap resampling?** Yes under the predeclared 80% peak-region stability criterion.
3. **Do the declared null models produce comparable features?** No at the declared 0.05 decision level.
4. **Is the reported statistical significance empirically calibrated?** Yes at 10%, 5%, 1%, and 0.1% within Monte-Carlo confidence intervals.
5. **Can injected signals of the observed magnitude be reliably recovered?** Yes under the predeclared 80% significant-correct-location recovery criterion.
6. **Is the result robust across the predeclared binning choices?** Yes across 120/160/200 bins.
7. **Does Co II provide independent compatible evidence?** Co II survives its own species-level checks, but it is not pooled with Fe II and is not described as an independent same-frequency experimental replication.
8. **What is the fixed-frequency significance?** Fe II primary descriptive fixed-k p = 0.00019996 (= 1/(5000+1), zero exceedances; Monte-Carlo floor, not p=0). The frequency was historically scan-selected, so this is **not** a global discovery p-value.
9. **What is the appropriately corrected/global significance, if established?** Fe II scan-global p = 0.00019996 (= 1/(5000+1), zero exceedances; Monte-Carlo floor, not p=0); conservative reported FWER-adjusted value = 0.0117976.
10. **What is the final Fe II classification?** **CONDITIONAL**.
11. **What is the final Co II classification?** **CONDITIONAL**.
12. **Is the result ready to freeze?** No; unresolved items are listed below.

## Scientific scope

This audit can establish **a statistically robust empirical spectral regularity in the frozen NIST catalog**. It does **not**, by itself, validate Wave Confinement Theory as a physical theory, establish a causal mechanism, or constitute NIST endorsement. Those are separate claims requiring separate physical evidence.

Fe II is primary. Co II is analyzed separately and is never combined post hoc with Fe II to increase significance. A fixed-frequency p-value at a frequency selected by a scan of the same data is labelled descriptive; the scan-global and declared-family corrections carry the discovery interpretation.

## Core decision matrix

| Check | Fe II | Co II |
|---|---:|---:|
| Frozen observed statistic reproduces | TRUE | TRUE |
| Bootstrap peak stability | TRUE | TRUE |
| Canonical scan-global null | TRUE | TRUE |
| Empirical p-value calibration | TRUE | TRUE |
| 1.0x observed injection recovery | TRUE | TRUE |
| 120/160/200 bin robustness | TRUE | TRUE |
| Sorting/cleaning implementation checks | TRUE | TRUE |
| Leave-decile-out influence stability | TRUE | TRUE |
| Declared-family/global accounting | TRUE | TRUE |
| Exact historical retrieval/query provenance recovered | FALSE | FALSE |

## Fe II adversarial null sensitivity

Worst declared Fe II alternative-null scan-global p: **0.00059988** under **block_residual_refit**. The original canonical null result is preserved regardless of this sensitivity result.

## Remaining uncertainties


- NIST Fe II source identity and provenance reconstruction are documented, but the exact historical retrieval/export date and exact historical query/export parameters are explicitly UNRECOVERABLE after the recorded search.
- NIST Co II source identity and provenance reconstruction are documented, but the exact historical retrieval/export date and exact historical query/export parameters are explicitly UNRECOVERABLE after the recorded search.
Finite Monte-Carlo results use the +1 correction, `(r + 1)/(B + 1)`. Zero exceedances are never reported as `p = 0`; they are reported at the simulation-resolution floor.

See `final_adversarial_audit.json`, `tables_r/statistical_audit/`, `REPRODUCE_FINAL_ADVERSARIAL_AUDIT.txt`, and (only after a PASS freeze) `results/final_adversarial_audit_v1/` for machine-readable evidence and reproduction artifacts.
