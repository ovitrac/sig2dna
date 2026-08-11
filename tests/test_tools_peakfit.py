"""
Tests for sig2dna_core.tools.peakfit (monotonepeakfit port): functional
recovery of known mixtures and tolerance-based Matlab-oracle parity
(Nelder-Mead trajectories are close but not bit-identical across
implementations).

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

import json
import os

import numpy as np
import pytest

from sig2dna_core.tools.peakfit import monotonepeakfit
from sig2dna_core.tools.peaks import monotonepeak

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


def synthetic_doublet():
    """Two overlapped 0.6006-Gaussians on a small linear baseline."""
    x = np.linspace(8.0, 17.0, 900)
    true = [(12.0, 0.30, 2.0), (12.6, 0.25, 1.2)]
    y = 0.01 * x + 0.05
    for c, w, a in true:
        y = y + a * np.exp(-(((x - c) / (0.6006 * w)) ** 2))
    return x, y, true


# ---------------------------------------------------------------------------
# functional tests
# ---------------------------------------------------------------------------
class TestMonotonepeakfit:
    def test_doublet_recovery(self):
        x, y, true = synthetic_doublet()
        pk = monotonepeak(y, x, mfilt=10, sort="descend", maxpeak=2, keeporder=True)
        assert len(pk) == 2
        fit = monotonepeakfit(pk, y, x, baseline=True, sort=True, keeporder=True)
        # strategy 2 (free centers) must recover centers, widths, weights
        for i, (c, w, a) in enumerate(true):
            assert abs(fit.position[1, i] - c) < 0.02, f"center {i}"
            assert abs(fit.width[1, i] - w) / w < 0.10, f"width {i}"
            assert abs(fit.weight[i, 1] - a) / a < 0.10, f"weight {i}"

    def test_model_matches_signal(self):
        x, y, _ = synthetic_doublet()
        pk = monotonepeak(y, x, mfilt=10, sort="descend", maxpeak=2, keeporder=True)
        fit = monotonepeakfit(pk, y, x, baseline=True, sort=True)
        model = fit.expansion(x, strategy=2)
        resid = np.linalg.norm(y - model) / np.linalg.norm(y - np.mean(y))
        assert resid < 0.05

    def test_two_strategies_returned(self):
        x, y, _ = synthetic_doublet()
        pk = monotonepeak(y, x, mfilt=10, sort="descend", maxpeak=2, keeporder=True)
        fit = monotonepeakfit(pk, y, x)
        assert fit.position.shape == (2, 2)
        # strategy 1 keeps prescribed centers
        np.testing.assert_allclose(fit.position[0], pk.center)

    def test_sigma_area_conventions(self):
        x, y, _ = synthetic_doublet()
        pk = monotonepeak(y, x, mfilt=10, sort="descend", maxpeak=2, keeporder=True)
        fit = monotonepeakfit(pk, y, x)
        np.testing.assert_allclose(fit.sigma, 0.6006 * fit.width / np.sqrt(2))
        np.testing.assert_allclose(fit.area, 0.6006 * np.sqrt(np.pi) * fit.width)

    def test_lorentzian_smoke(self):
        x, y, _ = synthetic_doublet()
        pk = monotonepeak(y, x, mfilt=10, sort="descend", maxpeak=2, keeporder=True)
        fit = monotonepeakfit(pk, y, x, lorentzian=True)
        assert fit.lorentzian
        assert np.all(np.isfinite(fit.position))

    def test_single_peak(self):
        x = np.linspace(0, 10, 1001)
        y = 3.0 * np.exp(-(((x - 5.0) / (0.6006 * 0.4)) ** 2))
        pk = monotonepeak(y, x, mfilt=5)
        fit = monotonepeakfit(pk, y, x, sort=True)
        assert abs(fit.position[1, 0] - 5.0) < 0.01
        assert abs(fit.weight[0, 1] - 3.0) / 3.0 < 0.05

    def test_empty_peaks_raises(self):
        from sig2dna_core.tools.peaks import PeakTable

        with pytest.raises(ValueError):
            monotonepeakfit(PeakTable.empty(), np.ones(100))


# ---------------------------------------------------------------------------
# Matlab-oracle parity (tolerance-based: same algorithm, same start simplex,
# different fminsearch implementations)
# ---------------------------------------------------------------------------
class TestFitParity:
    @pytest.fixture(scope="class")
    def fitted(self, oracle, chromatogram):
        x, y = chromatogram
        w = oracle["chromatogram"]["fit_window"]
        win = (x >= w["xmin"]) & (x <= w["xmax"])
        xw, yw = x[win], y[win]
        pk = monotonepeak(yw, xw, mfilt=15, sort="descend", maxpeak=2, keeporder=True)
        assert len(pk) == w["npeaks"]
        fit = monotonepeakfit(pk, yw, xw, baseline=True, sort=True)
        return fit, oracle["chromatogram"]["fit"]

    @pytest.mark.parametrize("k,key", [(0, "strategy1"), (1, "strategy2")])
    def test_positions_widths_weights(self, fitted, k, key):
        fit, ref = fitted
        np.testing.assert_allclose(
            fit.position[k], np.atleast_1d(ref[key]["position"]), atol=0.05
        )
        np.testing.assert_allclose(
            fit.width[k], np.atleast_1d(ref[key]["width"]), rtol=0.10
        )
        np.testing.assert_allclose(
            fit.weight[:, k], np.atleast_1d(ref[key]["weight"]), rtol=0.10
        )

    @pytest.mark.parametrize("k,key", [(0, "strategy1"), (1, "strategy2")])
    def test_r2(self, fitted, k, key):
        fit, ref = fitted
        # final expansion r2 must match the oracle's within 1e-3
        np.testing.assert_allclose(
            fit.r2[-1, k], np.atleast_1d(ref[key]["r2"])[-1], atol=1e-3
        )

    def test_rank_order(self, fitted):
        fit, ref = fitted
        np.testing.assert_array_equal(
            fit.rank[:, 1], np.atleast_1d(ref["strategy2"]["rank"]).astype(int)
        )


# ---------------------------------------------------------------------------
# control plot (smoke)
# ---------------------------------------------------------------------------
def test_plot_fit():
    import matplotlib

    matplotlib.use("Agg")
    from sig2dna_core.tools.peakfit import plot_fit

    x, y, _ = synthetic_doublet()
    pk = monotonepeak(y, x, mfilt=10, sort="descend", maxpeak=2, keeporder=True)
    fit = monotonepeakfit(pk, y, x, baseline=True, sort=True)
    ax = plot_fit(y, x, fit)
    assert ax is not None
