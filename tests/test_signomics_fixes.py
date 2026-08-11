"""
Regression tests for the five defects fixed in signomics v0.51
(review F1–F5, 2026-08-10). Each test fails on the pre-fix code.

Run: pytest -q tests/test_signomics_fixes.py
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)  # before signomics pulls pyplot

import numpy as np
import pytest

from sig2dna_core.signomics import DNAsignal, peaks


def make_signal():
    p = peaks()
    p.add(x=20, w=2, h=10, name="P1")
    p.add(x=60, w=3, h=5, name="P2")
    return p.to_signal(n=512), p


# F1 — encode_dna must encode ALL requested scales, not just the first
def test_f1_encode_dna_all_scales():
    sig, _ = make_signal()
    d = DNAsignal(sig, encode=True, scales=[1, 2, 4])
    assert {1, 2, 4}.issubset(set(d.codes.keys()))
    # direct call over already-transformed scales must encode every one of them
    d.codes.clear()
    codes = d.encode_dna(scales=[2, 4])
    assert {2, 4}.issubset(set(codes.keys()))


# F2 — constructor with plot=True must not raise on plot_codes
def test_f2_constructor_plot_smoke():
    import matplotlib.pyplot as plt

    sig, _ = make_signal()
    DNAsignal(sig, encode=True, plot=True, scales=[1, 2])
    plt.close("all")


# F3 — sparsify_cwt(inplace=False) returns a copy, original untouched
def test_f3_sparsify_cwt_copy():
    sig, _ = make_signal()
    d = DNAsignal(sig, encode=True, scales=[2])
    before = {s: c.copy() for s, c in d.cwt_coeffs.items()}
    out = d.sparsify_cwt(scale=2, inplace=False)
    assert out is not None and out is not d
    np.testing.assert_array_equal(d.cwt_coeffs[2], before[2])
    n_zero_out = int(np.sum(out.cwt_coeffs[2] == 0.0))
    n_zero_in = int(np.sum(d.cwt_coeffs[2] == 0.0))
    assert n_zero_out >= n_zero_in


# F4 — signal._toDNA must honor its scales argument
def test_f4_todna_forwards_scales():
    sig, _ = make_signal()
    d = sig._toDNA(encode=True, scales=[8])
    assert set(d.codes.keys()) == {8}
    d_default = sig._toDNA(encode=True)
    assert {1, 2, 4, 8, 16, 32}.issubset(set(d_default.codes.keys()))


# F5 — peaks lookup by name returns the record; missing name fails closed
def test_f5_peaks_getitem_str():
    _, p = make_signal()
    hit = p["P1"]
    assert hit is not None and hit["name"] == "P1"
    with pytest.raises(KeyError):
        p["does-not-exist"]
