"""
Tests for sig2dna_core.tools.baseline (ndf / removepeaks port), including
Matlab-oracle parity on the synthetic chromatogram fixture.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

import json
import os

import numpy as np
import pytest

from sig2dna_core.tools.baseline import ndf, removepeaks

HERE = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def oracle():
    with open(os.path.join(HERE, "fixtures", "tools_oracle.json")) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def chromatogram():
    with open(os.path.join(HERE, "fixtures", "tools_oracle_input.json")) as f:
        inp = json.load(f)
    return (
        np.asarray(inp["chromatogram"]["x"]),
        np.asarray(inp["chromatogram"]["y"]),
    )


# ---------------------------------------------------------------------------
# ndf
# ---------------------------------------------------------------------------
class TestNdf:
    def test_first_derivative_polynomial_exact(self):
        # the O(h^6) stencil differentiates low-order polynomials exactly;
        # indices 1..9 may be altered by the reference's boundary
        # glitch-repair heuristic (faithful to ndf.m), so the strict check
        # is on the untouched range
        t = np.linspace(0, 1, 100)
        y = 3 * t**4 - 2 * t**3 + t
        dy = 12 * t**3 - 6 * t**2 + 1
        d = ndf(t, y, 1)
        np.testing.assert_allclose(d[10:], dy[10:], rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(d[0], dy[0], rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(d[:10], dy[:10], rtol=2e-2)

    def test_second_derivative_polynomial_exact(self):
        t = np.linspace(0, 1, 100)
        y = t**4
        np.testing.assert_allclose(ndf(t, y, 2), 12 * t**2, rtol=1e-6, atol=1e-6)

    def test_sine_accuracy(self):
        # same boundary caveat as above: cos reverses direction of its
        # increments near t=0, triggering the reference's glitch repair
        t = np.linspace(0, 2 * np.pi, 1000)
        d = ndf(t, np.sin(t), 1)
        np.testing.assert_allclose(d[10:], np.cos(t)[10:], atol=1e-10)
        np.testing.assert_allclose(d[:10], np.cos(t)[:10], atol=2e-4)

    def test_columnwise(self):
        t = np.linspace(0, 1, 50)
        y = np.column_stack([t**2, t**3])
        d = ndf(t, y, 1)
        assert d.shape == y.shape
        np.testing.assert_allclose(d[:, 0], 2 * t, atol=1e-8)

    def test_too_short(self):
        with pytest.raises(ValueError):
            ndf(np.arange(5.0), np.arange(5.0), 1)

    def test_bad_order(self):
        with pytest.raises(ValueError):
            ndf(np.arange(10.0), np.arange(10.0), 3)


# ---------------------------------------------------------------------------
# removepeaks
# ---------------------------------------------------------------------------
class TestRemovepeaks:
    def test_gaussian_on_linear_baseline(self):
        x = np.linspace(0, 10, 2000)
        baseline = 0.5 + 0.1 * x
        y = baseline + np.exp(-(((x - 5) / 0.3) ** 2))
        ynp = removepeaks(x, y, [[4.0, 6.0]], smooth=5)
        # the bridge must recover the linear baseline under the peak
        inpeak = (x >= 4.2) & (x <= 5.8)
        np.testing.assert_allclose(ynp[inpeak], baseline[inpeak], atol=0.02)

    def test_continuity_at_boundaries(self):
        x = np.linspace(0, 10, 2000)
        y = 0.2 * x + np.exp(-(((x - 5) / 0.3) ** 2))
        ynp = removepeaks(x, y, [[4.0, 6.0]], smooth=5)
        # no jump at the bridge ends
        d = np.abs(np.diff(ynp))
        assert d.max() < 5 * np.median(d[d > 0]) + 1e-3

    def test_multiple_segments(self):
        x = np.linspace(0, 20, 4000)
        y = 1.0 + np.exp(-(((x - 5) / 0.3) ** 2)) + np.exp(-(((x - 15) / 0.4) ** 2))
        ynp = removepeaks(x, y, [[4, 6], [14, 16]], smooth=5)
        np.testing.assert_allclose(ynp, np.ones_like(x), atol=0.02)

    def test_models_output(self):
        x = np.linspace(0, 10, 1000)
        y = 0.1 * x + np.exp(-(((x - 5) / 0.3) ** 2))
        ynp, models = removepeaks(x, y, [[4, 6]], return_models=True)
        assert len(models) == 1 and len(models[0]) == 1
        bad = (x >= 4) & (x <= 6)
        np.testing.assert_allclose(models[0][0](x[bad]), ynp[bad], rtol=1e-12)

    def test_anchor_outside_range_raises(self):
        x = np.linspace(0, 10, 1000)
        y = np.exp(-(((x - 5) / 0.3) ** 2))
        with pytest.raises(ValueError):
            removepeaks(x, y, [[-1.0, 5.0]])


# ---------------------------------------------------------------------------
# Matlab-oracle parity
# ---------------------------------------------------------------------------
class TestRemovepeaksParity:
    def test_chromatogram(self, oracle, chromatogram):
        x, y = chromatogram
        ref = oracle["chromatogram"]["removepeaks"]
        segs = np.asarray(ref["segments"], dtype=float)
        if segs.shape[0] == 2 and segs.shape[1] != 2:
            segs = segs.T
        ynp = removepeaks(x, y, segs, f=0.0, smooth=5)
        np.testing.assert_allclose(
            ynp, np.asarray(ref["ynp"], dtype=float), rtol=1e-7, atol=1e-9
        )
