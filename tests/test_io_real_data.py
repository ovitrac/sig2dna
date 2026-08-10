"""
Integration tests against real acquisitions — OPTIONAL and strictly local.

These tests run only when the environment variable ``SIG2DNA_DATA`` points to
a directory containing real chromatograms (searched recursively). No real file
name, sample identifier, or measured value is stored in this repository; the
tests discover files by extension only.

Example (private workstation):
    SIG2DNA_DATA=/path/to/private/data pytest -q tests/test_io_real_data.py
"""

from __future__ import annotations

import glob
import os

import numpy as np
import pytest

from sig2dna_core.io import bin_uniform, read, read_mat

DATA = os.environ.get("SIG2DNA_DATA")

pytestmark = pytest.mark.skipif(
    not DATA, reason="SIG2DNA_DATA not set (private integration tests)"
)


def _find(pattern, limit=3):
    hits = sorted(glob.glob(os.path.join(DATA, "**", pattern), recursive=True))
    return hits[:limit]


def _sibling_mat(path):
    stem = os.path.splitext(path)[0]
    mat = stem + ".mat"
    return mat if os.path.isfile(mat) else None


@pytest.mark.parametrize("pattern", ["*.cdf", "*.CDF", "*.xml"])
def test_raw_reader_smoke(pattern):
    files = _find(pattern)
    if not files:
        pytest.skip(f"no {pattern} under SIG2DNA_DATA")
    for f in files:
        try:
            run = read(f)
        except ValueError as e:
            if pattern == "*.xml" and "not mzXML" in str(e):
                continue  # unrelated xml file
            raise
        assert run.n_scans > 100
        tic = run.get_tic()
        assert np.all(np.isfinite(tic)) and tic.max() > 0


@pytest.mark.parametrize("pattern", ["*.cdf", "*.CDF", "*.xml"])
def test_raw_vs_matlab_oracle(pattern):
    """Binned raw TIC must correlate with Julien's .mat sidecar (> 0.999)."""
    checked = 0
    for f in _find(pattern, limit=10):
        mat = _sibling_mat(f)
        if mat is None:
            continue
        try:
            raw = bin_uniform(read(f))
            oracle = read_mat(mat)
        except ValueError:
            continue
        tic_raw = raw.y.sum(axis=0)
        tic_oracle = oracle.y.sum(axis=0)
        # compare on the overlapping time window, resampled to the oracle grid
        t0 = max(raw.rt[0], oracle.rt[0])
        t1 = min(raw.rt[-1], oracle.rt[-1])
        sel = (oracle.rt >= t0) & (oracle.rt <= t1)
        resampled = np.interp(oracle.rt[sel], raw.rt, tic_raw)
        r = np.corrcoef(resampled, tic_oracle[sel])[0, 1]
        assert r > 0.999, f"TIC mismatch vs oracle (r={r:.5f}) for a {pattern} file"
        checked += 1
    if checked == 0:
        pytest.skip(f"no {pattern}+mat sibling pairs under SIG2DNA_DATA")


def test_agilent_d_smoke():
    folders = sorted(
        glob.glob(os.path.join(DATA, "**", "*.D"), recursive=True)
    )[:3]
    if not folders:
        pytest.skip("no .D folders under SIG2DNA_DATA")
    for d in folders:
        if not os.path.isdir(d):
            continue
        run = read(d)
        assert run.n_scans > 10
        assert run.get_tic().max() > 0
