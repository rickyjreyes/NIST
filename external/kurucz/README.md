# Kurucz GFALL external replication input

This directory is a staging location for the raw external catalog used by `independent_atomic_replication_v1`.

The raw catalog is intentionally **not committed to this repository**. Its identity is pinned by SHA-256 in `config/independent_atomic_replication_v1.json` before any scientific parsing/outcome access.

## Frozen source family

- Kurucz CD-ROM 23 / GFALL atomic line catalog
- CfA catalog interface / provenance page: `https://lweb.cfa.harvard.edu/amp/ampdata/kurucz23/sekur.html`
- Upstream GFALL family: `kurucz.harvard.edu/linelists/gfall`

The CfA page documents that its application uses `gfall.dat`; it reports a source file dated 2012-09-19 and a download/update on 2017-02-28. The replication protocol uses the raw file actually obtained for the replication and pins **those exact bytes** by SHA-256; do not assume a current upstream file is byte-identical to the CfA snapshot.

## Staging

Place the chosen raw source at:

```text
external/kurucz/gfall.dat
```

Before parsing it, run:

```powershell
& $Rscript R/run_independent_atomic_replication_v1.R --pin-source-hash
```

That command computes SHA-256 only and writes it to the frozen protocol. It does not parse Fe II or expose the replication outcome.

Then inspect and commit the protocol change:

```powershell
git diff -- config/independent_atomic_replication_v1.json
git add config/independent_atomic_replication_v1.json
git commit -m "Pin independent Kurucz GFALL source hash"
git push origin agent/final-adversarial-stat-audit
```

Only **after that commit** run the primary replication:

```powershell
& $Rscript R/run_independent_atomic_replication_v1.R --parallel
```

The primary target is fixed before the external outcome and may not be retuned after this step.
