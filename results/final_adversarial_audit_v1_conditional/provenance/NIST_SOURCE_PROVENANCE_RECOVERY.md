# NIST source provenance recovery record

## Purpose

This record documents the historical provenance reconstruction for the frozen Fe II and Co II CSV inputs used by the NIST adversarial audit. It is deliberately conservative: an exact retrieval date or exact ASD export query is not inferred from later publication metadata, repository timestamps, or local filesystem timestamps.

## Frozen inputs

- `data/Fe_lines.csv`
- `data/Co_lines.csv`

The final release tooling must record and verify SHA-256 hashes for these exact bytes. MD5 is not an acceptable substitute for the provenance gate.

## Recovered source identity

For both species, the recoverable source identity is:

- Provider: NIST Atomic Spectra Database (ASD), Standard Reference Database 78
- ASD version documented by the May 28, 2026 publication: 5.12
- Public interface documented by that publication: NIST ASD Lines Form
- Publication-recorded source access date: May 27, 2026

The May 27 date is **not** promoted to the historical CSV retrieval date because the frozen local CSV files already existed by May 25, 2026.

## Historical search performed

The following sources were checked for the original export/retrieval record:

1. NTFS alternate data streams (`Zone.Identifier`) for both frozen CSV files. No download/referrer metadata was present.
2. Local file metadata. Both frozen CSV files had local creation/last-write timestamps of May 25, 2026 00:27:16 on the audited workstation.
3. Tracked repository contents and repository history for `physics.nist.gov`, `nist.gov`, ASD URLs, and related query fragments. No authoritative original export URL/query was recovered.
4. PowerShell command history for NIST, `physics.nist`, `Fe_lines`, and `Co_lines`. Later analysis commands were present; the original export/download command was not.
5. Chrome and Edge URL history across the available profiles, including the May 20–27, 2026 interval. No original NIST ASD export visit was recovered.
6. Chrome and Edge download databases (`downloads` and `downloads_url_chains`) across the available profiles. No Fe/Co ASD export record was recovered. The only later NIST download found was an unrelated periodic-table image.
7. The May 28, 2026 publication record. It establishes NIST ASD SRD 78, ASD v5.12, the Lines Form, and a May 27 source-access date, but not the earlier exact export time or exact query that produced the frozen CSV bytes.

## Disposition of unrecoverable fields

The exact historical retrieval/export date and exact historical ASD query/export parameters are classified as `UNRECOVERABLE` unless an authoritative contemporaneous record is later found.

`UNRECOVERABLE` means the historical field was searched for and not recovered. It does **not** mean a value may be inferred from the nearest known timestamp. In particular:

- May 25, 2026 is evidence that the frozen files existed locally by that time; it is not proof of the retrieval time.
- May 27, 2026 is a publication-recorded ASD access date; it is not proof that the frozen files were exported on that date.

A future newly generated ASD export may be used only as a documented comparison/reference export. It must not silently replace the historical frozen inputs.

## Provenance gate semantics

The release tooling distinguishes two questions:

1. `provenance_record_complete`: Are the source identity, recovered facts, unrecoverable-field dispositions, evidence record, and frozen-input SHA-256 checks documented and machine-checkable?
2. `exact_historical_provenance_recovered`: Were the exact original retrieval date and exact original query/export settings actually recovered?

A complete reconstruction record can therefore coexist with `exact_historical_provenance_recovered = FALSE`. The latter remains a scientific/release limitation and must not be converted into a recovered fact merely to obtain a PASS label.
