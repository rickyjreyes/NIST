#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
NIST neighboring-ion batch runner for log-cosine scans
-----------------------------------------------------

This script runs your existing scanner on each neighboring transition-metal
NIST CSV, summarizes results, and automatically promotes hits to full 5000-null
and bin-stability tests.

It DOES NOT change the statistic. It calls:

    nist_wct_log_spectral_scan_FIXED.py

so the outputs remain comparable to your Fe II anchor result.

Expected folder:
    nist_wct_log_spectral_scan_FIXED.py
    Ni_lines.csv
    Co_lines.csv
    Cr_lines.csv
    Mn_lines.csv
    Ti_lines.csv

Fast preview:
    python nist_batch_run_neighbors.py --preview-null 500 --full-null 5000

Only preview, no full promotion:
    python nist_batch_run_neighbors.py --preview-null 500 --no-promote

Outputs:
    outputs_nist_batch/
      batch_results.csv
      batch_results.json
      <element>_ion2_160_preview/
      <element>_ion2_160_full/
      <element>_ion2_120_full_bin120/
      <element>_ion2_200_full_bin200/
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

import numpy as np

DEFAULT_FILES = [
    "Ni_lines.csv",
    "Co_lines.csv",
    "Cr_lines.csv",
    "Mn_lines.csv",
    "Ti_lines.csv",
]

SUMMARY_COLUMNS = [
    "label", "csv", "ion", "bins", "null_n", "mode", "out_dir", "status",
    "n_unique_lines", "k_best", "deltaD", "scan_null_p", "tail_count_ge",
    "n_obs", "ell_min", "ell_max", "delta_ell",
    "wavenumber_min_cm", "wavenumber_max_cm", "error",
]


def label_from_csv(csv_path: str | Path) -> str:
    stem = Path(csv_path).stem
    return stem.replace("_lines", "").replace("-lines", "").lower()


def load_summary(out_dir: Path) -> Optional[Dict[str, Any]]:
    p = out_dir / "nist_summary.json"
    if not p.exists():
        return None
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def summarize(label: str, csv_path: str, ion: int, bins: int, null_n: int, mode: str, out_dir: Path, status: str, error: str = "") -> Dict[str, Any]:
    s = load_summary(out_dir)
    row = {
        "label": label,
        "csv": csv_path,
        "ion": ion,
        "bins": bins,
        "null_n": null_n,
        "mode": mode,
        "out_dir": str(out_dir),
        "status": status,
        "n_unique_lines": "",
        "k_best": "",
        "deltaD": "",
        "scan_null_p": "",
        "tail_count_ge": "",
        "n_obs": "",
        "ell_min": "",
        "ell_max": "",
        "delta_ell": "",
        "wavenumber_min_cm": "",
        "wavenumber_max_cm": "",
        "error": error,
    }
    if s:
        row.update({
            "n_unique_lines": s.get("n_unique_lines", ""),
            "k_best": s.get("best", {}).get("k_best", ""),
            "deltaD": s.get("best", {}).get("deltaD", ""),
            "scan_null_p": s.get("scan_null_p", ""),
            "tail_count_ge": s.get("tail_count_ge", ""),
            "n_obs": s.get("branch_report", {}).get("n_obs", ""),
            "ell_min": s.get("ell_min", ""),
            "ell_max": s.get("ell_max", ""),
            "delta_ell": s.get("delta_ell", ""),
            "wavenumber_min_cm": s.get("wavenumber_min_cm", ""),
            "wavenumber_max_cm": s.get("wavenumber_max_cm", ""),
        })
    return row


def should_skip(out_dir: Path, force: bool) -> bool:
    return (not force) and (out_dir / "nist_summary.json").exists()


def run_one(
    scanner: Path,
    csv_path: str,
    ion: int,
    bins: int,
    null_n: int,
    out_dir: Path,
    min_lines: int,
    n_k: int,
    k_min: float,
    k_max: float,
    baseline_sigma: float,
    degree: int,
    seed: int,
    force: bool,
    mode: str,
) -> Dict[str, Any]:
    label = label_from_csv(csv_path)

    if should_skip(out_dir, force):
        print(f"[skip] existing {out_dir}")
        return summarize(label, csv_path, ion, bins, null_n, mode, out_dir, "skipped_existing")

    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(scanner),
        "--csv", str(csv_path),
        "--ion", str(ion),
        "--bins", str(bins),
        "--null-n", str(null_n),
        "--min-lines", str(min_lines),
        "--out-dir", str(out_dir),
        "--n-k", str(n_k),
        "--k-min", str(k_min),
        "--k-max", str(k_max),
        "--baseline-sigma", str(baseline_sigma),
        "--degree", str(degree),
        "--seed", str(seed),
    ]

    print("\n" + "=" * 100)
    print("[run]", " ".join(cmd))
    print("=" * 100)

    try:
        completed = subprocess.run(cmd, check=False)
        if completed.returncode != 0:
            err = f"scanner exit code {completed.returncode}"
            print(f"[error] {err}")
            return summarize(label, csv_path, ion, bins, null_n, mode, out_dir, "failed", err)
        return summarize(label, csv_path, ion, bins, null_n, mode, out_dir, "ok")
    except Exception as exc:
        err = repr(exc)
        print(f"[error] {err}")
        return summarize(label, csv_path, ion, bins, null_n, mode, out_dir, "failed", err)


def p_value(row: Dict[str, Any]) -> float:
    try:
        return float(row.get("scan_null_p", "nan"))
    except Exception:
        return float("nan")


def write_results(rows: List[Dict[str, Any]], out_root: Path) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    csv_path = out_root / "batch_results.csv"
    json_path = out_root / "batch_results.json"

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in SUMMARY_COLUMNS})

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)

    print("\n[save]", csv_path)
    print("[save]", json_path)


def print_table(rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    cols = ["label", "mode", "bins", "n_unique_lines", "k_best", "deltaD", "scan_null_p", "tail_count_ge", "n_obs", "status"]
    print("\n" + "#" * 100)
    print("[batch summary]")
    print("#" * 100)
    widths = {c: max(len(c), max(len(str(r.get(c, ""))) for r in rows)) for c in cols}
    header = "  ".join(c.ljust(widths[c]) for c in cols)
    print(header)
    print("-" * len(header))
    for r in rows:
        print("  ".join(str(r.get(c, "")).ljust(widths[c]) for c in cols))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scanner", default="nist_wct_log_spectral_scan_FIXED.py")
    ap.add_argument("--csv", nargs="*", default=DEFAULT_FILES)
    ap.add_argument("--ion", type=int, default=2)
    ap.add_argument("--bins", type=int, default=160)
    ap.add_argument("--preview-null", type=int, default=500)
    ap.add_argument("--full-null", type=int, default=5000)
    ap.add_argument("--hit-p", type=float, default=0.01)
    ap.add_argument("--out-root", default="outputs_nist_batch")
    ap.add_argument("--min-lines", type=int, default=100)
    ap.add_argument("--n-k", type=int, default=2500)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=80.0)
    ap.add_argument("--baseline-sigma", type=float, default=6.0)
    ap.add_argument("--degree", type=int, default=1)
    ap.add_argument("--seed", type=int, default=20260517)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-promote", action="store_true", help="Do not run 5000-null/bin stability for preview hits.")
    args = ap.parse_args()

    scanner = Path(args.scanner)
    if not scanner.exists():
        raise FileNotFoundError(f"Scanner not found: {scanner.resolve()}")

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    all_rows: List[Dict[str, Any]] = []

    print("[batch] files:", args.csv)
    print("[batch] ion:", args.ion)
    print("[batch] preview null:", args.preview_null)
    print("[batch] full null:", args.full_null)
    print("[batch] promotion p threshold:", args.hit_p)

    preview_rows = []
    for csv_path in args.csv:
        if not Path(csv_path).exists():
            row = {
                "label": label_from_csv(csv_path),
                "csv": csv_path,
                "ion": args.ion,
                "bins": args.bins,
                "null_n": args.preview_null,
                "mode": "preview",
                "out_dir": "",
                "status": "missing_csv",
                "error": f"CSV not found: {Path(csv_path).resolve()}",
            }
            all_rows.append(row)
            preview_rows.append(row)
            continue

        label = label_from_csv(csv_path)
        out_dir = out_root / f"{label}_ion{args.ion}_{args.bins}_preview"
        row = run_one(
            scanner=scanner,
            csv_path=csv_path,
            ion=args.ion,
            bins=args.bins,
            null_n=args.preview_null,
            out_dir=out_dir,
            min_lines=args.min_lines,
            n_k=args.n_k,
            k_min=args.k_min,
            k_max=args.k_max,
            baseline_sigma=args.baseline_sigma,
            degree=args.degree,
            seed=args.seed,
            force=args.force,
            mode="preview",
        )
        all_rows.append(row)
        preview_rows.append(row)
        write_results(all_rows, out_root)

    if not args.no_promote:
        for row in preview_rows:
            p = p_value(row)
            if not np.isfinite(p) or p >= args.hit_p or row.get("status") not in ("ok", "skipped_existing"):
                continue

            csv_path = row["csv"]
            label = row["label"]
            print(f"\n[promote] {label}: preview p={p} < {args.hit_p}")

            for bins, mode in [(args.bins, "full"), (120, "full_bin120"), (200, "full_bin200")]:
                out_dir = out_root / f"{label}_ion{args.ion}_{bins}_{mode}"
                full_row = run_one(
                    scanner=scanner,
                    csv_path=csv_path,
                    ion=args.ion,
                    bins=bins,
                    null_n=args.full_null,
                    out_dir=out_dir,
                    min_lines=args.min_lines,
                    n_k=args.n_k,
                    k_min=args.k_min,
                    k_max=args.k_max,
                    baseline_sigma=args.baseline_sigma,
                    degree=args.degree,
                    seed=args.seed,
                    force=args.force,
                    mode=mode,
                )
                all_rows.append(full_row)
                write_results(all_rows, out_root)

    write_results(all_rows, out_root)
    print_table(all_rows)

    hits = [r for r in all_rows if r.get("status") in ("ok", "skipped_existing") and np.isfinite(p_value(r)) and p_value(r) < args.hit_p]
    hits = sorted(hits, key=lambda r: (p_value(r), -float(r.get("deltaD") or 0)))
    if hits:
        print("\n[hit list]")
        for r in hits:
            print(
                f"{r['label']} {r['mode']} bins={r['bins']} "
                f"k={r['k_best']} n={r['n_obs']} p={r['scan_null_p']} "
                f"tail={r['tail_count_ge']}/{r['null_n']} lines={r['n_unique_lines']}"
            )
    else:
        print("\n[hit list] no p<threshold hits")


if __name__ == "__main__":
    main()
