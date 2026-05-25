#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run the canonical Fe II log-cosine bin-stability ladder.

Bins: 120, 160, 200
Default null-n: 500 (quick reproducibility). Use --null-n 5000 for full.

Outputs go to:
    outputs/fe_ion2_120/
    outputs/fe_ion2_160/
    outputs/fe_ion2_200/

Master summary is written to:
    tables/nist_master_results.csv

If a canonical output folder already contains nist_summary.json, that run is
reused unless --force is supplied. No old run is deleted.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
SCANNER = REPO_ROOT / "scripts" / "nist_wct_log_spectral_scan_FIXED.py"
FE_CSV = REPO_ROOT / "data" / "Fe_lines.csv"
OUT_ROOT = REPO_ROOT / "outputs"
TABLES_DIR = REPO_ROOT / "tables"
MASTER_CSV = TABLES_DIR / "nist_master_results.csv"

BINS_LADDER = [120, 160, 200]

MASTER_COLUMNS = [
    "species",
    "ion",
    "bins",
    "n_unique_lines",
    "k_best",
    "n_obs",
    "deltaD",
    "scan_null_p",
    "tail_count_ge",
    "null_n",
    "ell_min",
    "ell_max",
    "delta_ell",
    "out_dir",
    "status",
]


def load_summary(out_dir: Path) -> Optional[Dict[str, Any]]:
    p = out_dir / "nist_summary.json"
    if not p.exists():
        return None
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def summary_to_row(species: str, ion: int, bins: int, out_dir: Path, status: str) -> Dict[str, Any]:
    row = {c: "" for c in MASTER_COLUMNS}
    row.update({
        "species": species,
        "ion": ion,
        "bins": bins,
        "out_dir": str(out_dir.relative_to(REPO_ROOT)) if out_dir.is_absolute() else str(out_dir),
        "status": status,
    })
    s = load_summary(out_dir)
    if s:
        best = s.get("best", {}) or {}
        br = s.get("branch_report", {}) or {}
        row["n_unique_lines"] = s.get("n_unique_lines", "")
        row["k_best"] = best.get("k_best", "")
        row["n_obs"] = br.get("n_obs", "")
        row["deltaD"] = best.get("deltaD", "")
        row["scan_null_p"] = s.get("scan_null_p", "")
        row["tail_count_ge"] = s.get("tail_count_ge", "")
        row["null_n"] = s.get("null_n", "")
        row["ell_min"] = s.get("ell_min", "")
        row["ell_max"] = s.get("ell_max", "")
        row["delta_ell"] = s.get("delta_ell", "")
    return row


def run_scanner(bins: int, null_n: int, out_dir: Path, min_lines: int, seed: int, n_k: int) -> int:
    cmd = [
        sys.executable, str(SCANNER),
        "--csv", str(FE_CSV),
        "--ion", "2",
        "--bins", str(bins),
        "--null-n", str(null_n),
        "--min-lines", str(min_lines),
        "--out-dir", str(out_dir),
        "--n-k", str(n_k),
        "--seed", str(seed),
    ]
    print("[run]", " ".join(cmd))
    completed = subprocess.run(cmd, check=False)
    return completed.returncode


def write_master(rows: List[Dict[str, Any]]) -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    with open(MASTER_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=MASTER_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in MASTER_COLUMNS})
    print(f"[save] {MASTER_CSV.relative_to(REPO_ROOT)}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-n", type=int, default=500)
    ap.add_argument("--min-lines", type=int, default=100)
    ap.add_argument("--n-k", type=int, default=2500)
    ap.add_argument("--seed", type=int, default=20260517)
    ap.add_argument("--force", action="store_true",
                    help="Re-run scanner even if canonical output already exists. "
                         "Existing files are overwritten by the scanner inside that "
                         "canonical folder; non-canonical legacy folders are never touched.")
    args = ap.parse_args()

    if not SCANNER.exists():
        print(f"[error] scanner not found: {SCANNER}", file=sys.stderr)
        return 2
    if not FE_CSV.exists():
        print(f"[error] Fe CSV not found: {FE_CSV}", file=sys.stderr)
        return 2

    rows: List[Dict[str, Any]] = []
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    for bins in BINS_LADDER:
        out_dir = OUT_ROOT / f"fe_ion2_{bins}"
        out_dir.mkdir(parents=True, exist_ok=True)
        existing = load_summary(out_dir)
        if existing and not args.force:
            print(f"[reuse] {out_dir.relative_to(REPO_ROOT)} (nist_summary.json present)")
            status = "reused"
        else:
            rc = run_scanner(bins, args.null_n, out_dir, args.min_lines, args.seed, args.n_k)
            status = "ok" if rc == 0 else f"failed_rc{rc}"
        rows.append(summary_to_row("Fe", 2, bins, out_dir, status))

    write_master(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
