#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Build the PASS/FAIL verdict for the Fe II log-cosine bin-stability claim.

Reads:
    tables/nist_master_results.csv

Writes:
    tables/nist_verdict.csv

Optionally updates RESULTS.md between the markers:
    <!-- BEGIN_VERDICT -->
    <!-- END_VERDICT -->

Verdict logic (PASS only if ALL hold):
    - rows exist for Fe ion=2 bins in {120, 160, 200}
    - n_unique_lines == 9447 for all three
    - k_best identical across the three (tolerance 1e-6)
    - scan_null_p <= 1 / (null_n + 1)  (i.e., zero exceedances)
    - n_obs in [10.6, 10.9]
    - deltaD > 0
Otherwise FAIL. No PARTIAL verdict.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
MASTER_CSV = REPO_ROOT / "tables" / "nist_master_results.csv"
VERDICT_CSV = REPO_ROOT / "tables" / "nist_verdict.csv"
RESULTS_MD = REPO_ROOT / "RESULTS.md"

REQUIRED_BINS = [120, 160, 200]
REQUIRED_LINES = 9447
K_TOL = 1e-6
N_OBS_MIN = 10.6
N_OBS_MAX = 10.9

BEGIN_MARK = "<!-- BEGIN_VERDICT -->"
END_MARK = "<!-- END_VERDICT -->"


def _to_float(x: Any) -> Optional[float]:
    try:
        if x is None or x == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def _to_int(x: Any) -> Optional[int]:
    f = _to_float(x)
    if f is None:
        return None
    return int(round(f))


def read_master(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def find_feii_rows(rows: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
    found: Dict[int, Dict[str, Any]] = {}
    for r in rows:
        sp = (r.get("species") or "").strip()
        ion = _to_int(r.get("ion"))
        bins = _to_int(r.get("bins"))
        if sp.lower() == "fe" and ion == 2 and bins in REQUIRED_BINS:
            if bins in found:
                continue
            found[bins] = r
    return found


def evaluate(rows: List[Dict[str, Any]]) -> Tuple[str, List[str], Dict[int, Dict[str, Any]]]:
    reasons: List[str] = []
    feii = find_feii_rows(rows)

    missing = [b for b in REQUIRED_BINS if b not in feii]
    if missing:
        reasons.append(f"missing Fe II rows for bins {missing}")
        return "FAIL", reasons, feii

    for b, r in feii.items():
        status = (r.get("status") or "").strip()
        if status not in ("ok", "reused"):
            reasons.append(f"bins={b} has non-success status '{status}'")

    lines_set = {_to_int(feii[b].get("n_unique_lines")) for b in REQUIRED_BINS}
    if lines_set != {REQUIRED_LINES}:
        reasons.append(f"n_unique_lines mismatch: {sorted(x for x in lines_set if x is not None)} (expected {REQUIRED_LINES})")

    k_vals = [_to_float(feii[b].get("k_best")) for b in REQUIRED_BINS]
    if any(k is None for k in k_vals):
        reasons.append(f"k_best missing or non-numeric: {k_vals}")
    else:
        if max(k_vals) - min(k_vals) > K_TOL:
            reasons.append(f"k_best not stable within tol {K_TOL}: {k_vals}")

    for b in REQUIRED_BINS:
        r = feii[b]
        null_n = _to_int(r.get("null_n"))
        p = _to_float(r.get("scan_null_p"))
        if null_n is None or null_n <= 0:
            reasons.append(f"bins={b} invalid null_n={r.get('null_n')}")
            continue
        threshold = 1.0 / (null_n + 1)
        if p is None or p > threshold + 1e-12:
            reasons.append(f"bins={b} scan_null_p={p} > 1/(null_n+1)={threshold} (need zero exceedances)")

        n_obs = _to_float(r.get("n_obs"))
        if n_obs is None or not (N_OBS_MIN <= n_obs <= N_OBS_MAX):
            reasons.append(f"bins={b} n_obs={n_obs} outside [{N_OBS_MIN}, {N_OBS_MAX}]")

        dd = _to_float(r.get("deltaD"))
        if dd is None or dd <= 0:
            reasons.append(f"bins={b} deltaD={dd} not > 0")

    verdict = "PASS" if not reasons else "FAIL"
    return verdict, reasons, feii


def write_verdict_csv(verdict: str, reasons: List[str], feii: Dict[int, Dict[str, Any]]) -> None:
    VERDICT_CSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ["claim", "verdict", "bins_120_k_best", "bins_160_k_best", "bins_200_k_best",
            "bins_120_n_obs", "bins_160_n_obs", "bins_200_n_obs",
            "bins_120_deltaD", "bins_160_deltaD", "bins_200_deltaD",
            "bins_120_p", "bins_160_p", "bins_200_p",
            "reasons"]
    row = {c: "" for c in cols}
    row["claim"] = "Fe II log-cosine k_best bin-stability (bins 120/160/200)"
    row["verdict"] = verdict
    for b in REQUIRED_BINS:
        r = feii.get(b) or {}
        row[f"bins_{b}_k_best"] = r.get("k_best", "")
        row[f"bins_{b}_n_obs"] = r.get("n_obs", "")
        row[f"bins_{b}_deltaD"] = r.get("deltaD", "")
        row[f"bins_{b}_p"] = r.get("scan_null_p", "")
    row["reasons"] = "; ".join(reasons)
    with open(VERDICT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerow(row)
    print(f"[save] {VERDICT_CSV.relative_to(REPO_ROOT)}")


def render_block(verdict: str, reasons: List[str], feii: Dict[int, Dict[str, Any]]) -> str:
    lines = [BEGIN_MARK, "", f"## Fe II Bin-Stability Verdict: **{verdict}**", ""]
    lines.append("| Bins | n_unique_lines | k_best | n_obs | deltaD | scan_null_p | null_n | status |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---|")
    for b in REQUIRED_BINS:
        r = feii.get(b) or {}
        lines.append(
            f"| {b} | {r.get('n_unique_lines','')} | {r.get('k_best','')} | "
            f"{r.get('n_obs','')} | {r.get('deltaD','')} | {r.get('scan_null_p','')} | "
            f"{r.get('null_n','')} | {r.get('status','')} |"
        )
    lines.append("")
    if reasons:
        lines.append("Failure reasons:")
        for r in reasons:
            lines.append(f"- {r}")
        lines.append("")
    else:
        lines.append("All gating criteria satisfied.")
        lines.append("")
    lines.append(END_MARK)
    return "\n".join(lines)


def update_results_md(block: str) -> None:
    if not RESULTS_MD.exists():
        RESULTS_MD.write_text(block + "\n", encoding="utf-8")
        return
    text = RESULTS_MD.read_text(encoding="utf-8")
    if BEGIN_MARK in text and END_MARK in text:
        start = text.index(BEGIN_MARK)
        end = text.index(END_MARK) + len(END_MARK)
        new_text = text[:start] + block + text[end:]
    else:
        sep = "" if text.endswith("\n") else "\n"
        new_text = text + sep + "\n" + block + "\n"
    RESULTS_MD.write_text(new_text, encoding="utf-8")
    print(f"[update] {RESULTS_MD.relative_to(REPO_ROOT)} (verdict block)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-results-md", action="store_true",
                    help="Do not update RESULTS.md, just write tables/nist_verdict.csv.")
    args = ap.parse_args()

    rows = read_master(MASTER_CSV)
    if not rows:
        print(f"[error] no rows in {MASTER_CSV}. Run scripts/run_feii_bin_stability.py first.",
              file=sys.stderr)
        return 2

    verdict, reasons, feii = evaluate(rows)
    write_verdict_csv(verdict, reasons, feii)
    block = render_block(verdict, reasons, feii)
    if not args.no_results_md:
        update_results_md(block)
    print(f"[verdict] {verdict}")
    if reasons:
        for r in reasons:
            print(f"  - {r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
