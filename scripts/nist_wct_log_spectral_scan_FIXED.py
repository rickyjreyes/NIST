#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
NIST Atomic Spectra WCT Log-Spectral Scanner
--------------------------------------------

Run:
    python nist_wct_log_spectral_scan.py --csv Fe_lines.csv --null-n 500 --min-lines 100

Full:
    python nist_wct_log_spectral_scan.py --csv Fe_lines.csv --null-n 5000 --min-lines 100

This script reads NIST CSV output and scans line density in:
    ell = ln(wavenumber_cm^-1)

It handles NIST/Excel-style cells like:
    ="49998.80"
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
except Exception as exc:
    raise RuntimeError("Missing scipy. Install with: pip install scipy") from exc

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None


EPS = 1e-12
OUT_DIR_DEFAULT = "outputs_nist_wct"

WCT_REFERENCE_NS = [10.0, 15.0, 20.0]
WCT_FOLDED_NS = [20.0 / 3.0, 15.0, 40.0 / 3.0]


def clean_cell(x: Any) -> Any:
    if pd.isna(x):
        return x
    s = str(x).strip()
    # NIST often exports cells as ="123.45"
    if s.startswith('="') and s.endswith('"'):
        s = s[2:-1]
    elif s.startswith("="):
        s = s[1:]
    s = s.strip('"').strip()
    return s


def to_float_series(s: pd.Series) -> pd.Series:
    cleaned = s.map(clean_cell).astype(str)
    # Strip non-numeric flags like "2300d?", "250bl(Fe III)", "(0)".
    cleaned = cleaned.str.replace(r"^\s*$", "", regex=True)
    cleaned = cleaned.str.extract(r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", expand=False)
    return pd.to_numeric(cleaned, errors="coerce")


def read_nist_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path.resolve()}")

    # NIST CSV rows may contain commas inside term/config fields.
    # pandas with quotechar usually works. Use python engine and skip bad only if needed.
    try:
        df = pd.read_csv(path, dtype=str, engine="python", quotechar='"', on_bad_lines="warn")
    except Exception:
        df = pd.read_csv(path, dtype=str, engine="python", quotechar='"', escapechar="\\", on_bad_lines="skip")

    df.columns = [str(c).strip().lstrip("\ufeff") for c in df.columns]
    # pandas 2.1+ replaces DataFrame.applymap with DataFrame.map.
    df_map = getattr(df, "map", None)
    if callable(df_map):
        try:
            return df_map(clean_cell)
        except TypeError:
            pass
    return df.applymap(clean_cell)


def norm_col(c: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(c).lower()).strip("_")


def find_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    norm = {norm_col(c): c for c in df.columns}
    for cand in candidates:
        nc = norm_col(cand)
        if nc in norm:
            return norm[nc]
    for cand in candidates:
        nc = norm_col(cand)
        for k, orig in norm.items():
            if nc in k or k in nc:
                return orig
    return None


def extract_wavenumber_cm(df: pd.DataFrame) -> Tuple[pd.Series, str]:
    wn_col = find_column(df, [
        "wn(cm-1)", "wn", "wavenumber", "wavenumber_cm", "wavenumber in cm-1",
        "Ritz(cm-1)", "Observed(cm-1)"
    ])
    if wn_col is not None:
        wn = to_float_series(df[wn_col])
        if wn.notna().sum() > 0:
            return wn, wn_col

    wl_nm_col = find_column(df, [
        "ritz_wl_air(nm)", "obs_wl_air(nm)", "ritz_wl_vac(nm)", "obs_wl_vac(nm)",
        "wavelength_nm", "wavelength nm"
    ])
    if wl_nm_col is not None:
        wl_nm = to_float_series(df[wl_nm_col])
        wn = 1.0e7 / wl_nm
        return wn, f"converted_from_{wl_nm_col}"

    raise KeyError(f"Could not find wavenumber or wavelength column. Columns are: {list(df.columns)}")


def clean_lines(df: pd.DataFrame, ion: Optional[int]) -> pd.DataFrame:
    out = df.copy()

    if ion is not None and "sp_num" in out.columns:
        sp = to_float_series(out["sp_num"])
        out = out[sp == float(ion)].copy()

    wn, source_col = extract_wavenumber_cm(out)
    out["wavenumber_cm"] = pd.to_numeric(wn, errors="coerce")
    out = out[np.isfinite(out["wavenumber_cm"]) & (out["wavenumber_cm"] > 0)].copy()
    out["ell"] = np.log(out["wavenumber_cm"].to_numpy(float))
    out["wavenumber_source"] = source_col
    out = out.sort_values("wavenumber_cm").drop_duplicates(subset=["wavenumber_cm"]).reset_index(drop=True)
    return out


def build_binned(lines: pd.DataFrame, bins: int, ell_min: Optional[float], ell_max: Optional[float]) -> pd.DataFrame:
    ell = lines["ell"].to_numpy(float)
    if ell_min is None:
        ell_min = float(np.min(ell))
    if ell_max is None:
        ell_max = float(np.max(ell))
    counts, edges = np.histogram(ell, bins=bins, range=(ell_min, ell_max))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return pd.DataFrame({"ell": centers, "count": counts.astype(float), "edge_lo": edges[:-1], "edge_hi": edges[1:]})


def poisson_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    y = np.asarray(y, float)
    mu = np.maximum(np.asarray(mu, float), EPS)
    term = mu - y
    nz = y > 0
    term[nz] += y[nz] * np.log(y[nz] / mu[nz])
    return float(2.0 * np.sum(term))


def fit_poisson_loglinear(y: np.ndarray, baseline: np.ndarray, X: np.ndarray, ridge: float = 1e-8, max_iter: int = 60):
    y = np.asarray(y, float)
    B = np.maximum(np.asarray(baseline, float), EPS)
    X = np.asarray(X, float)

    beta = np.zeros(X.shape[1], float)
    beta[0] = np.log(max(np.sum(y), EPS) / max(np.sum(B), EPS))

    for _ in range(max_iter):
        eta = np.clip(X @ beta, -10, 10)
        mu = np.maximum(B * np.exp(eta), EPS)
        z = eta + (y - mu) / mu
        W = mu
        A = X.T @ (W[:, None] * X) + ridge * np.eye(X.shape[1])
        b = X.T @ (W * z)
        try:
            beta2 = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            beta2 = np.linalg.lstsq(A, b, rcond=None)[0]
        if np.max(np.abs(beta2 - beta)) < 1e-8:
            beta = beta2
            break
        beta = beta2

    eta = np.clip(X @ beta, -10, 10)
    mu = np.maximum(B * np.exp(eta), EPS)
    return beta, poisson_deviance(y, mu), mu


def design_poly(ell: np.ndarray, degree: int) -> np.ndarray:
    z = ell - np.mean(ell)
    s = np.std(z)
    if s > 0:
        z = z / s
    return np.column_stack([z ** d for d in range(degree + 1)])


def scan_k(ell: np.ndarray, y: np.ndarray, baseline: np.ndarray, k_grid: np.ndarray, degree: int):
    X0 = design_poly(ell, degree)
    beta0, D0, mu0 = fit_poisson_loglinear(y, baseline, X0)

    rows = []
    best = None
    for k in k_grid:
        X = np.column_stack([X0, np.cos(k * ell), np.sin(k * ell)])
        beta, D, mu = fit_poisson_loglinear(y, baseline, X)
        delta = D0 - D
        amp = float(math.sqrt(beta[-2] ** 2 + beta[-1] ** 2))
        phase = float(math.atan2(-beta[-1], beta[-2]))
        row = {
            "k": float(k), "deltaD": float(delta), "D_base": float(D0),
            "D_harmonic": float(D), "amplitude": amp, "phase": phase
        }
        rows.append(row)
        if best is None or delta > best["deltaD"]:
            best = {
                "k_best": float(k), "deltaD": float(delta), "D_base": float(D0),
                "D_harmonic": float(D), "amplitude": amp, "phase": phase
            }
    return pd.DataFrame(rows), best, mu0


def null_scan(ell: np.ndarray, y: np.ndarray, baseline: np.ndarray, mu0: np.ndarray, k_grid: np.ndarray, degree: int, real_delta: float, null_n: int, seed: int):
    rng = np.random.default_rng(seed)
    vals = np.empty(null_n, float)
    for i in range(null_n):
        y0 = rng.poisson(mu0)
        _, b, _ = scan_k(ell, y0, baseline, k_grid, degree)
        vals[i] = b["deltaD"]
        if (i + 1) % 100 == 0 or i + 1 == null_n:
            print(f"[null] {i+1}/{null_n}")
    p = float((1 + np.sum(vals >= real_delta)) / (1 + len(vals)))
    return vals, p, int(np.sum(vals >= real_delta))


def branch_report(k_best: float, delta_ell: float) -> Dict[str, Any]:
    n_obs = k_best * delta_ell / (2.0 * math.pi)
    targets = []
    for label, ns in [
        ("koide_10_15_20", WCT_REFERENCE_NS),
        ("folded_4over9", WCT_FOLDED_NS),
        ("integer_1_to_40", list(map(float, range(1, 41)))),
    ]:
        for n in ns:
            k_t = 2.0 * math.pi * n / delta_ell
            targets.append({
                "branch": label,
                "n": float(n),
                "k_target": float(k_t),
                "k_error": float(k_best - k_t),
                "abs_k_error": float(abs(k_best - k_t)),
                "n_error": float(n_obs - n),
                "abs_n_error": float(abs(n_obs - n)),
            })
    return {"n_obs": float(n_obs), "nearest": sorted(targets, key=lambda r: r["abs_k_error"])[:10]}


def peak_ratios(scan: pd.DataFrame, top_n: int = 12) -> pd.DataFrame:
    y = scan["deltaD"].to_numpy(float)
    peaks, _ = find_peaks(y, distance=max(1, len(y) // 80))
    if len(peaks) == 0:
        peaks = np.array([int(np.argmax(y))])
    pk = scan.iloc[peaks].sort_values("deltaD", ascending=False).head(top_n).reset_index(drop=True)
    rows = []
    for i in range(len(pk)):
        for j in range(i + 1, len(pk)):
            lo, hi = sorted([float(pk.loc[i, "k"]), float(pk.loc[j, "k"])])
            r = hi / lo
            f = Fraction(r).limit_denominator(64)
            rows.append({
                "k_low": lo, "k_high": hi, "ratio": r,
                "rational": f"{f.numerator}/{f.denominator}",
                "rational_error": abs(r - f.numerator / f.denominator),
            })
    return pd.DataFrame(rows).sort_values("rational_error").reset_index(drop=True) if rows else pd.DataFrame()


def save_plots(out_dir: Path, binned: pd.DataFrame, baseline: np.ndarray, mu0: np.ndarray, scan: pd.DataFrame, best: Dict[str, Any], null_vals: Optional[np.ndarray]):
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.step(binned["ell"], binned["count"], where="mid", label="line counts")
    ax.plot(binned["ell"], baseline, label="smooth baseline")
    ax.plot(binned["ell"], mu0, label="base fit")
    ax.set_xlabel("ell = ln(wavenumber cm^-1)")
    ax.set_ylabel("counts/bin")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "nist_spectrum_fit.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(scan["k"], scan["deltaD"])
    ax.axvline(best["k_best"], linestyle="--", label=f"k_best={best['k_best']:.4g}")
    ax.set_xlabel("k")
    ax.set_ylabel("DeltaD")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "nist_scan_curve.png", dpi=180)
    plt.close(fig)

    if null_vals is not None:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(null_vals, bins=50)
        ax.axvline(best["deltaD"], linestyle="--", label=f"obs={best['deltaD']:.4g}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "nist_null.png", dpi=180)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out-dir", default=OUT_DIR_DEFAULT)
    ap.add_argument("--ion", type=int, default=None, help="Optional sp_num filter: 1=Fe I, 2=Fe II, 3=Fe III")
    ap.add_argument("--bins", type=int, default=160)
    ap.add_argument("--ell-min", type=float, default=None)
    ap.add_argument("--ell-max", type=float, default=None)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=80.0)
    ap.add_argument("--n-k", type=int, default=2500)
    ap.add_argument("--degree", type=int, default=1)
    ap.add_argument("--baseline-sigma", type=float, default=6.0)
    ap.add_argument("--null-n", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260517)
    ap.add_argument("--min-lines", type=int, default=100)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    raw = read_nist_csv(args.csv)
    print(f"[raw] rows={len(raw)} cols={list(raw.columns)}")

    lines = clean_lines(raw, ion=args.ion)
    print(f"[clean] usable unique lines={len(lines)}")
    if len(lines) < args.min_lines:
        raise RuntimeError(f"Only {len(lines)} usable lines. Need at least {args.min_lines}.")

    lines.to_csv(out_dir / "nist_lines_clean.csv", index=False)

    binned = build_binned(lines, args.bins, args.ell_min, args.ell_max)
    y = binned["count"].to_numpy(float)
    ell = binned["ell"].to_numpy(float)
    baseline = np.maximum(gaussian_filter1d(y, sigma=args.baseline_sigma, mode="nearest"), EPS)
    binned["baseline"] = baseline
    binned.to_csv(out_dir / "nist_binned_spectrum.csv", index=False)

    k_grid = np.linspace(args.k_min, args.k_max, args.n_k)
    print(f"[scan] bins={args.bins} k={args.k_min}..{args.k_max} n_k={args.n_k}")
    scan, best, mu0 = scan_k(ell, y, baseline, k_grid, args.degree)
    scan.to_csv(out_dir / "nist_scan_curve.csv", index=False)

    null_vals = None
    p = None
    tail = None
    if args.null_n > 0:
        null_vals, p, tail = null_scan(ell, y, baseline, mu0, k_grid, args.degree, best["deltaD"], args.null_n, args.seed)
        pd.DataFrame({"null_max_deltaD": null_vals}).to_csv(out_dir / "nist_null.csv", index=False)

    delta_ell = float(np.max(ell) - np.min(ell))
    br = branch_report(best["k_best"], delta_ell)
    ratios = peak_ratios(scan)
    if not ratios.empty:
        ratios.to_csv(out_dir / "nist_peak_ratios.csv", index=False)

    summary = {
        "csv": args.csv,
        "ion_filter": args.ion,
        "n_raw_rows": int(len(raw)),
        "n_unique_lines": int(len(lines)),
        "wavenumber_min_cm": float(lines["wavenumber_cm"].min()),
        "wavenumber_max_cm": float(lines["wavenumber_cm"].max()),
        "ell_min": float(np.min(ell)),
        "ell_max": float(np.max(ell)),
        "delta_ell": delta_ell,
        "best": best,
        "scan_null_p": p,
        "tail_count_ge": tail,
        "null_n": int(args.null_n),
        "branch_report": br,
        "top_peak_ratios": ratios.head(20).to_dict(orient="records") if not ratios.empty else [],
    }
    with open(out_dir / "nist_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    save_plots(out_dir, binned, baseline, mu0, scan, best, null_vals)

    print("\n[done]")
    print(json.dumps({
        "k_best": best["k_best"],
        "deltaD": best["deltaD"],
        "scan_null_p": p,
        "tail_count_ge": tail,
        "null_n": args.null_n,
        "n_obs": br["n_obs"],
        "n_unique_lines": len(lines),
    }, indent=2))
    print("\n[nearest branch]")
    print(pd.DataFrame(br["nearest"]).to_string(index=False))


if __name__ == "__main__":
    main()
