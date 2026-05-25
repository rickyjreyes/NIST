#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
NIST Atomic Spectra WCT Log-Spectral Scanner - CuPy accelerated
---------------------------------------------------------------

This is the CuPy version of the NIST scanner.

It keeps the same data model as the CPU script:
    ell = ln(wavenumber_cm^-1)
    counts are binned in ell
    baseline is Gaussian-smoothed counts
    model is Poisson log-linear:
        lambda_i = baseline_i * exp(poly_drift + a cos(k ell_i) + b sin(k ell_i))

CuPy accelerates:
    - full k scan
    - Poisson null scan-max ensemble

Run:
    python nist_wct_log_spectral_scan_CUPY.py --csv Fe_lines.csv --ion 2 --bins 160 --null-n 5000 --out-dir outputs_feii_cupy

Neighbor examples:
    python nist_wct_log_spectral_scan_CUPY.py --csv Ni_lines.csv --ion 2 --bins 160 --null-n 500 --out-dir outputs_ni_ion2_160_cupy

Notes:
    - If CuPy import fails, this script stops. It does not silently fall back.
    - Use --null-batch and --k-batch if GPU memory is tight.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import cupy as cp
except Exception as exc:
    raise RuntimeError(
        "CuPy is required for this script. Install the CUDA-matched build, e.g.:\n"
        "  pip install cupy-cuda12x\n"
        "or use your existing cupy-env."
    ) from exc

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
OUT_DIR_DEFAULT = "outputs_nist_wct_cupy"

WCT_REFERENCE_NS = [10.0, 15.0, 20.0]
WCT_FOLDED_NS = [20.0 / 3.0, 15.0, 40.0 / 3.0]


# =============================================================================
# CSV loading / cleaning
# =============================================================================

def clean_cell(x: Any) -> Any:
    if pd.isna(x):
        return x
    s = str(x).strip()
    if s.startswith('="') and s.endswith('"'):
        s = s[2:-1]
    elif s.startswith("="):
        s = s[1:]
    return s.strip('"').strip()


def to_float_series(s: pd.Series) -> pd.Series:
    cleaned = s.map(clean_cell).astype(str)
    cleaned = cleaned.str.replace(r"^\s*$", "", regex=True)
    cleaned = cleaned.str.extract(r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", expand=False)
    return pd.to_numeric(cleaned, errors="coerce")


def read_nist_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path.resolve()}")

    try:
        df = pd.read_csv(path, dtype=str, engine="python", quotechar='"', on_bad_lines="warn")
    except Exception:
        df = pd.read_csv(path, dtype=str, engine="python", quotechar='"', escapechar="\\", on_bad_lines="skip")

    df.columns = [str(c).strip().lstrip("\ufeff") for c in df.columns]

    # Avoid deprecated DataFrame.applymap on newer pandas.
    try:
        return df.map(clean_cell)
    except AttributeError:
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
        return 1.0e7 / wl_nm, f"converted_from_{wl_nm_col}"

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


# =============================================================================
# GPU math
# =============================================================================

def poisson_deviance_np(y: np.ndarray, mu: np.ndarray) -> float:
    y = np.asarray(y, float)
    mu = np.maximum(np.asarray(mu, float), EPS)
    term = mu - y
    nz = y > 0
    term[nz] += y[nz] * np.log(y[nz] / mu[nz])
    return float(2.0 * np.sum(term))


def design_poly_np(ell: np.ndarray, degree: int) -> np.ndarray:
    z = ell - np.mean(ell)
    s = np.std(z)
    if s > 0:
        z = z / s
    return np.column_stack([z ** d for d in range(degree + 1)]).astype(np.float64)


def build_X0_gpu(ell: np.ndarray, degree: int):
    return cp.asarray(design_poly_np(ell, degree), dtype=cp.float64)


def build_Xk_gpu(ell: np.ndarray, k_grid: np.ndarray, degree: int):
    X0 = design_poly_np(ell, degree)
    ell_g = cp.asarray(ell, dtype=cp.float64)
    k_g = cp.asarray(k_grid, dtype=cp.float64)
    phase = k_g[:, None] * ell_g[None, :]
    c = cp.cos(phase)
    s = cp.sin(phase)
    X0g = cp.asarray(X0, dtype=cp.float64)
    X0_rep = cp.broadcast_to(X0g[None, :, :], (len(k_grid), len(ell), X0.shape[1]))
    return cp.concatenate([X0_rep, c[:, :, None], s[:, :, None]], axis=2)


def batched_poisson_irls_single_y(y_np: np.ndarray, baseline_np: np.ndarray, X_g, max_iter: int, ridge: float):
    """
    Fit many designs for one y.

    X_g:
        shape (K,N,P) or (N,P)

    Returns:
        beta, deviance, mu
    """
    y_g = cp.asarray(y_np, dtype=cp.float64)
    B_g = cp.maximum(cp.asarray(baseline_np, dtype=cp.float64), EPS)

    if X_g.ndim == 2:
        X = X_g[None, :, :]
        squeeze = True
    else:
        X = X_g
        squeeze = False

    K, N, P = X.shape
    beta = cp.zeros((K, P), dtype=cp.float64)

    c0 = cp.log(cp.maximum(cp.sum(y_g), EPS) / cp.maximum(cp.sum(B_g), EPS))
    beta[:, 0] = c0

    I = cp.eye(P, dtype=cp.float64)[None, :, :]

    for _ in range(max_iter):
        eta = cp.einsum("knp,kp->kn", X, beta)
        eta = cp.clip(eta, -10.0, 10.0)
        mu = cp.maximum(B_g[None, :] * cp.exp(eta), EPS)
        z = eta + (y_g[None, :] - mu) / mu
        WX = mu[:, :, None] * X
        A = cp.einsum("knp,knq->kpq", X, WX) + ridge * I
        b = cp.einsum("knp,kn->kp", X, mu * z)
        beta2 = cp.linalg.solve(A, b)
        if float(cp.max(cp.abs(beta2 - beta)).get()) < 1e-8:
            beta = beta2
            break
        beta = beta2

    eta = cp.einsum("knp,kp->kn", X, beta)
    eta = cp.clip(eta, -10.0, 10.0)
    mu = cp.maximum(B_g[None, :] * cp.exp(eta), EPS)

    term = mu - y_g[None, :]
    nz = y_g > 0
    if bool(cp.any(nz).get()):
        term[:, nz] += y_g[None, nz] * cp.log(y_g[None, nz] / mu[:, nz])
    D = 2.0 * cp.sum(term, axis=1)

    if squeeze:
        return beta[0], D[0], mu[0]
    return beta, D, mu


def scan_k_gpu(ell: np.ndarray, y: np.ndarray, baseline: np.ndarray, k_grid: np.ndarray, degree: int, irls_iter: int, ridge: float, k_batch: int):
    X0_g = build_X0_gpu(ell, degree)
    beta0, D0_g, mu0_g = batched_poisson_irls_single_y(y, baseline, X0_g, irls_iter, ridge)
    D0 = float(D0_g.get())
    mu0 = cp.asnumpy(mu0_g)

    rows = []
    best = None

    for start in range(0, len(k_grid), k_batch):
        stop = min(start + k_batch, len(k_grid))
        kg = k_grid[start:stop]
        Xk = build_Xk_gpu(ell, kg, degree)
        beta, D_g, mu = batched_poisson_irls_single_y(y, baseline, Xk, irls_iter, ridge)
        D = cp.asnumpy(D_g)
        beta_np = cp.asnumpy(beta)

        delta = D0 - D
        amp = np.sqrt(beta_np[:, -2] ** 2 + beta_np[:, -1] ** 2)
        phase = np.arctan2(-beta_np[:, -1], beta_np[:, -2])

        for i, k in enumerate(kg):
            row = {
                "k": float(k),
                "deltaD": float(delta[i]),
                "D_base": D0,
                "D_harmonic": float(D[i]),
                "amplitude": float(amp[i]),
                "phase": float(phase[i]),
            }
            rows.append(row)
            if best is None or row["deltaD"] > best["deltaD"]:
                best = {
                    "k_best": row["k"],
                    "deltaD": row["deltaD"],
                    "D_base": row["D_base"],
                    "D_harmonic": row["D_harmonic"],
                    "amplitude": row["amplitude"],
                    "phase": row["phase"],
                }

        # free chunk
        del Xk, beta, D_g, mu
        cp.get_default_memory_pool().free_all_blocks()

    return pd.DataFrame(rows), best, mu0


def fit_base_many_y_gpu(Y_g, baseline_g, X0_g, max_iter: int, ridge: float):
    """
    Fit base model to many y vectors.

    Y_g: (B,N)
    X0_g: (N,P)
    Returns:
        D0: (B,)
        mu: (B,N)
    """
    Y = cp.asarray(Y_g, dtype=cp.float64)
    B, N = Y.shape
    X = cp.asarray(X0_g, dtype=cp.float64)
    P = X.shape[1]
    base = cp.maximum(cp.asarray(baseline_g, dtype=cp.float64), EPS)

    beta = cp.zeros((B, P), dtype=cp.float64)
    beta[:, 0] = cp.log(cp.maximum(cp.sum(Y, axis=1), EPS) / cp.maximum(cp.sum(base), EPS))
    I = cp.eye(P, dtype=cp.float64)[None, :, :]

    for _ in range(max_iter):
        eta = cp.einsum("np,bp->bn", X, beta)
        eta = cp.clip(eta, -10.0, 10.0)
        mu = cp.maximum(base[None, :] * cp.exp(eta), EPS)
        z = eta + (Y - mu) / mu
        WX = mu[:, :, None] * X[None, :, :]
        A = cp.einsum("np,bnq->bpq", X, WX) + ridge * I
        b = cp.einsum("np,bn->bp", X, mu * z)
        beta2 = cp.linalg.solve(A, b)
        if float(cp.max(cp.abs(beta2 - beta)).get()) < 1e-8:
            beta = beta2
            break
        beta = beta2

    eta = cp.einsum("np,bp->bn", X, beta)
    eta = cp.clip(eta, -10.0, 10.0)
    mu = cp.maximum(base[None, :] * cp.exp(eta), EPS)
    D = poisson_deviance_many_gpu(Y, mu)
    return D, mu


def poisson_deviance_many_gpu(Y, MU):
    term = MU - Y
    nz = Y > 0
    term = cp.where(nz, term + Y * cp.log(cp.maximum(Y, EPS) / MU), term)
    return 2.0 * cp.sum(term, axis=1)


def harmonic_deviance_many_y_for_kchunk(Y_g, baseline_g, X_g, max_iter: int, ridge: float):
    """
    Fit harmonic model for many y vectors and a chunk of k designs.

    Y_g: (B,N)
    X_g: (K,N,P)

    Returns D: (B,K)
    """
    Y = cp.asarray(Y_g, dtype=cp.float64)
    base = cp.maximum(cp.asarray(baseline_g, dtype=cp.float64), EPS)
    X = cp.asarray(X_g, dtype=cp.float64)

    Bn, N = Y.shape
    K, _, P = X.shape

    beta = cp.zeros((Bn, K, P), dtype=cp.float64)
    c0 = cp.log(cp.maximum(cp.sum(Y, axis=1), EPS) / cp.maximum(cp.sum(base), EPS))
    beta[:, :, 0] = c0[:, None]

    I = cp.eye(P, dtype=cp.float64)[None, None, :, :]

    for _ in range(max_iter):
        eta = cp.einsum("knp,bkp->bkn", X, beta)
        eta = cp.clip(eta, -10.0, 10.0)
        mu = cp.maximum(base[None, None, :] * cp.exp(eta), EPS)
        z = eta + (Y[:, None, :] - mu) / mu

        WX = mu[:, :, :, None] * X[None, :, :, :]
        A = cp.einsum("knp,bknq->bkpq", X, WX) + ridge * I
        b = cp.einsum("knp,bkn->bkp", X, mu * z)
        beta2 = cp.linalg.solve(A, b)

        if float(cp.max(cp.abs(beta2 - beta)).get()) < 1e-8:
            beta = beta2
            break
        beta = beta2

    eta = cp.einsum("knp,bkp->bkn", X, beta)
    eta = cp.clip(eta, -10.0, 10.0)
    mu = cp.maximum(base[None, None, :] * cp.exp(eta), EPS)

    Ybk = Y[:, None, :]
    term = mu - Ybk
    nz = Ybk > 0
    term = cp.where(nz, term + Ybk * cp.log(cp.maximum(Ybk, EPS) / mu), term)
    D = 2.0 * cp.sum(term, axis=2)
    return D


def null_scan_gpu(ell: np.ndarray, baseline: np.ndarray, mu0: np.ndarray, k_grid: np.ndarray, real_delta: float, null_n: int, seed: int, degree: int, irls_iter: int, ridge: float, null_batch: int, k_batch: int):
    rng = cp.random.default_rng(seed)

    baseline_g = cp.asarray(baseline, dtype=cp.float64)
    mu0_g = cp.asarray(mu0, dtype=cp.float64)
    X0_g = build_X0_gpu(ell, degree)

    vals = np.empty(null_n, dtype=np.float64)
    done = 0

    while done < null_n:
        bsz = min(null_batch, null_n - done)

        Y = rng.poisson(mu0_g[None, :], size=(bsz, len(mu0))).astype(cp.float64)
        D0, mu_base = fit_base_many_y_gpu(Y, baseline_g, X0_g, irls_iter, ridge)

        best_delta = cp.full((bsz,), -cp.inf, dtype=cp.float64)

        for start in range(0, len(k_grid), k_batch):
            stop = min(start + k_batch, len(k_grid))
            kg = k_grid[start:stop]
            Xk = build_Xk_gpu(ell, kg, degree)
            D1 = harmonic_deviance_many_y_for_kchunk(Y, baseline_g, Xk, irls_iter, ridge)
            delta = D0[:, None] - D1
            best_delta = cp.maximum(best_delta, cp.max(delta, axis=1))
            del Xk, D1, delta
            cp.get_default_memory_pool().free_all_blocks()

        vals[done:done + bsz] = cp.asnumpy(best_delta)
        done += bsz
        print(f"[null] {done}/{null_n}")

        del Y, D0, mu_base, best_delta
        cp.get_default_memory_pool().free_all_blocks()

    p = float((1 + np.sum(vals >= real_delta)) / (1 + len(vals)))
    return vals, p, int(np.sum(vals >= real_delta))


# =============================================================================
# Reports / plots
# =============================================================================

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


# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out-dir", default=OUT_DIR_DEFAULT)
    ap.add_argument("--ion", type=int, default=None, help="Optional sp_num filter: 1=neutral, 2=singly ionized, 3=doubly ionized")
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

    ap.add_argument("--irls-iter", type=int, default=35)
    ap.add_argument("--ridge", type=float, default=1e-8)
    ap.add_argument("--k-batch", type=int, default=256)
    ap.add_argument("--null-batch", type=int, default=64)

    args = ap.parse_args()

    print("[gpu] CuPy version:", cp.__version__)
    print("[gpu] CUDA device count:", cp.cuda.runtime.getDeviceCount())

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
    scan, best, mu0 = scan_k_gpu(
        ell=ell,
        y=y,
        baseline=baseline,
        k_grid=k_grid,
        degree=args.degree,
        irls_iter=args.irls_iter,
        ridge=args.ridge,
        k_batch=args.k_batch,
    )
    scan.to_csv(out_dir / "nist_scan_curve.csv", index=False)

    null_vals = None
    p = None
    tail = None
    if args.null_n > 0:
        null_vals, p, tail = null_scan_gpu(
            ell=ell,
            baseline=baseline,
            mu0=mu0,
            k_grid=k_grid,
            real_delta=best["deltaD"],
            null_n=args.null_n,
            seed=args.seed,
            degree=args.degree,
            irls_iter=args.irls_iter,
            ridge=args.ridge,
            null_batch=args.null_batch,
            k_batch=args.k_batch,
        )
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
        "engine": {
            "cupy": True,
            "cupy_version": cp.__version__,
            "k_batch": args.k_batch,
            "null_batch": args.null_batch,
            "irls_iter": args.irls_iter,
        },
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
