#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Lightweight tests for the canonical NIST log-cosine pipeline."""

from __future__ import annotations

import csv
import importlib.util
import io
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
DATA_DIR = REPO_ROOT / "data"
TABLES_DIR = REPO_ROOT / "tables"
SCANNER = SCRIPTS_DIR / "nist_wct_log_spectral_scan_FIXED.py"


def _load_scanner_module():
    spec = importlib.util.spec_from_file_location("scanner_fixed", SCANNER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_scanner_help_runs():
    result = subprocess.run(
        [sys.executable, str(SCANNER), "--help"],
        capture_output=True, text=True, cwd=str(REPO_ROOT),
    )
    assert result.returncode == 0, result.stderr
    assert "--csv" in result.stdout
    assert "--ion" in result.stdout
    assert "--bins" in result.stdout


def test_clean_cell_strips_nist_excel_quotes():
    scanner = _load_scanner_module()
    assert scanner.clean_cell('="123.45"') == "123.45"
    assert scanner.clean_cell('="49998.80"') == "49998.80"
    assert scanner.clean_cell("=123.45") == "123.45"
    assert scanner.clean_cell('  "999.0"  ') == "999.0"
    assert scanner.clean_cell("plain") == "plain"


def test_to_float_series_handles_flag_suffixes():
    scanner = _load_scanner_module()
    s = pd.Series(['="49998.80"', "2300d?", "250bl(Fe III)", "(0)", "", None])
    out = scanner.to_float_series(s)
    assert out.iloc[0] == pytest.approx(49998.80)
    assert out.iloc[1] == pytest.approx(2300.0)
    assert out.iloc[2] == pytest.approx(250.0)
    assert out.iloc[3] == pytest.approx(0.0)


def test_wavenumber_column_detection():
    scanner = _load_scanner_module()
    df = pd.DataFrame({
        "element": ["Fe", "Fe"],
        "sp_num": ["2", "2"],
        "wn(cm-1)": ['="40000.00"', '="20000.00"'],
    })
    wn, src = scanner.extract_wavenumber_cm(df)
    assert "wn" in src.lower()
    assert float(wn.iloc[0]) == pytest.approx(40000.0)
    assert float(wn.iloc[1]) == pytest.approx(20000.0)


def test_ion_filter_keeps_only_requested_sp_num():
    scanner = _load_scanner_module()
    df = pd.DataFrame({
        "element": ["Fe"] * 4,
        "sp_num": ["1", "2", "2", "3"],
        "wn(cm-1)": ['="10000"', '="20000"', '="30000"', '="40000"'],
    })
    cleaned = scanner.clean_lines(df, ion=2)
    assert len(cleaned) == 2
    assert set(cleaned["wavenumber_cm"].astype(float)) == {20000.0, 30000.0}


def test_read_nist_csv_handles_real_fe_csv():
    if not (DATA_DIR / "Fe_lines.csv").exists():
        pytest.skip("data/Fe_lines.csv not present in this checkout")
    scanner = _load_scanner_module()
    df = scanner.read_nist_csv(DATA_DIR / "Fe_lines.csv")
    assert "sp_num" in df.columns
    assert "wn(cm-1)" in df.columns or any("wn" in c.lower() for c in df.columns)


def test_master_results_has_required_columns_if_present():
    p = TABLES_DIR / "nist_master_results.csv"
    if not p.exists():
        pytest.skip("master results not generated yet")
    required = {
        "species", "ion", "bins", "n_unique_lines", "k_best", "n_obs",
        "deltaD", "scan_null_p", "tail_count_ge", "null_n",
        "ell_min", "ell_max", "delta_ell", "out_dir", "status",
    }
    with open(p, "r", encoding="utf-8") as f:
        header = next(csv.reader(f))
    assert required.issubset(set(header)), f"missing cols: {required - set(header)}"


def test_verdict_requires_all_three_bins():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        import make_verdict
    finally:
        sys.path.pop(0)

    # Only bins 120 and 160 present; bins 200 missing -> must FAIL.
    rows = [
        {"species": "Fe", "ion": "2", "bins": "120", "n_unique_lines": "9447",
         "k_best": "31.3265306122449", "n_obs": "10.7171", "deltaD": "355.74",
         "scan_null_p": "0.0002", "null_n": "5000", "status": "ok"},
        {"species": "Fe", "ion": "2", "bins": "160", "n_unique_lines": "9447",
         "k_best": "31.3265306122449", "n_obs": "10.7396", "deltaD": "315.73",
         "scan_null_p": "0.0002", "null_n": "5000", "status": "ok"},
    ]
    verdict, reasons, _ = make_verdict.evaluate(rows)
    assert verdict == "FAIL"
    assert any("missing" in r.lower() for r in reasons)


def test_verdict_passes_on_canonical_three_bin_table():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        import make_verdict
    finally:
        sys.path.pop(0)
    rows = []
    for b, dd, n_obs in [(120, 355.74, 10.7171), (160, 315.73, 10.7396), (200, 259.18, 10.7531)]:
        rows.append({
            "species": "Fe", "ion": "2", "bins": str(b), "n_unique_lines": "9447",
            "k_best": "31.3265306122449", "n_obs": str(n_obs), "deltaD": str(dd),
            "scan_null_p": "0.0001999600079984003", "null_n": "5000", "status": "ok",
        })
    verdict, reasons, _ = make_verdict.evaluate(rows)
    assert verdict == "PASS", reasons


def test_verdict_fails_when_k_best_not_stable():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        import make_verdict
    finally:
        sys.path.pop(0)
    rows = []
    k_vals = ["31.3265306122449", "31.3265306122449", "31.4"]
    for b, k, dd, n_obs in zip([120, 160, 200], k_vals, [355.74, 315.73, 259.18], [10.717, 10.740, 10.753]):
        rows.append({
            "species": "Fe", "ion": "2", "bins": str(b), "n_unique_lines": "9447",
            "k_best": k, "n_obs": str(n_obs), "deltaD": str(dd),
            "scan_null_p": "0.0002", "null_n": "5000", "status": "ok",
        })
    verdict, reasons, _ = make_verdict.evaluate(rows)
    assert verdict == "FAIL"
    assert any("k_best" in r for r in reasons)


def test_no_destructive_calls_in_scripts():
    """Preservation guard: no script may delete output directories."""
    forbidden = [
        re.compile(r"\bshutil\.rmtree\b"),
        re.compile(r"\bos\.removedirs\b"),
        re.compile(r"\bos\.rmdir\b"),
        re.compile(r"\bsubprocess\.[a-z_]+\([^)]*['\"]rm['\"][^)]*-rf"),
        re.compile(r"['\"]rm\s+-rf['\"]"),
        re.compile(r"['\"]rmdir\s+/s"),
    ]
    offenders = []
    for py in SCRIPTS_DIR.glob("*.py"):
        text = py.read_text(encoding="utf-8")
        for pat in forbidden:
            if pat.search(text):
                offenders.append(f"{py.name}: pattern {pat.pattern}")
    assert not offenders, "destructive call(s) found in scripts/: " + "; ".join(offenders)


def test_canonical_output_dirs_present_for_fe_ii():
    out_root = REPO_ROOT / "outputs"
    for b in (120, 160, 200):
        assert (out_root / f"fe_ion2_{b}").is_dir(), f"missing outputs/fe_ion2_{b}"
