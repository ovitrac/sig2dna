"""
Tests for sig2dna_core.tools.cwtarea — peak-area reconstruction from the
Ricker CWT and revalidation of the thesis calibration constants
(Property 5.1 / Figs 3.4 and 5.48 of tel-04194172).

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

import json
import os

import numpy as np
import pytest
from scipy.signal import fftconvolve

from sig2dna_core.icfilter import ricker_kernel
from sig2dna_core.tools.cwtarea import (AMPLIFICATION_FACTOR, LOPT_FACTOR,
                                        calibration_curve,
                                        interpolate_optimal_scale,
                                        optimal_scale, reconstruct_areas,
                                        ricker_gaussian_apex, scale_ridge)
from sig2dna_core.tools.peaks import monotonepeak

HERE = os.path.dirname(os.path.abspath(__file__))
DYADIC = [2.0**k for k in range(0, 9)]


def cwt_stack(y, scales=DYADIC):
    """Multi-scale stack via the exact engine (dict, signomics-compatible)."""
    return {s: fftconvolve(y, ricker_kernel(s), mode="same") for s in scales}


def gaussian_signal(n, peaks):
    """peaks = [(center, sigma, height)] in sample units; returns y, areas."""
    t = np.arange(n, dtype=float)
    y = np.zeros(n)
    for c, s, h in peaks:
        y += h * np.exp(-((t - c) ** 2) / (2.0 * s**2))
    areas = [h * s * np.sqrt(2 * np.pi) for _, s, h in peaks]
    return y, areas


# ---------------------------------------------------------------------------
# closed form vs both engines — the revalidation record
# ---------------------------------------------------------------------------
class TestConstants:
    def test_closed_form_apex_matches_exact_engine(self):
        n, sigma = 8000, 20.0
        y, _ = gaussian_signal(n, [(n / 2, sigma, 1.0)])
        for scale in (10.0, 25.0, 44.72, 80.0):
            w = fftconvolve(y, ricker_kernel(scale), mode="same")
            assert w[n // 2] == pytest.approx(
                ricker_gaussian_apex(sigma, scale), rel=5e-3
            )

    def test_optimal_scale_is_sqrt5_sigma(self):
        n, sigma = 8000, 20.0
        y, _ = gaussian_signal(n, [(n / 2, sigma, 1.0)])
        ls = np.linspace(1.5 * sigma, 3.5 * sigma, 300)
        amps = [fftconvolve(y, ricker_kernel(sc), mode="same")[n // 2] for sc in ls]
        lopt = ls[int(np.argmax(amps))]
        assert lopt / sigma == pytest.approx(LOPT_FACTOR, rel=5e-3)
        assert max(amps) / np.sqrt(sigma) == pytest.approx(
            AMPLIFICATION_FACTOR, rel=5e-3
        )

    def test_pywt_engine_agrees(self):
        pywt = pytest.importorskip("pywt")
        n, sigma = 8000, 20.0
        y, _ = gaussian_signal(n, [(n / 2, sigma, 1.0)])
        for scale in (16.0, 32.0, 64.0):  # dyadic (pywt-safe) scales
            coef, _ = pywt.cwt(y, [scale], "mexh")
            assert coef[0][n // 2] == pytest.approx(
                ricker_gaussian_apex(sigma, scale), rel=5e-3
            )

    def test_thesis_constants_infirmed(self):
        """Revalidation verdict: Property 5.1 constants (2.58, 2.77, 4.79,
        1.588) and Fig 3.4's 2.48 do NOT describe the current engines; the
        closed form (sqrt(5)=2.2361, sqrt(6)=2.4495, sqrt(18)=4.2426,
        1.10597) does — asserted here as the durable record."""
        assert LOPT_FACTOR == pytest.approx(2.2361, abs=1e-4)
        assert AMPLIFICATION_FACTOR == pytest.approx(1.10597, abs=1e-4)
        # the thesis values differ far beyond numerical tolerance
        assert abs(LOPT_FACTOR - 2.58) > 0.3
        assert abs(LOPT_FACTOR - 2.48) > 0.2
        assert abs(AMPLIFICATION_FACTOR - 1.588) > 0.4


# ---------------------------------------------------------------------------
# ridge and interpolation machinery
# ---------------------------------------------------------------------------
class TestRidge:
    def test_ridge_tracks_apex(self):
        n, sigma = 4000, 15.0
        y, _ = gaussian_signal(n, [(n / 2, sigma, 2.0)])
        ridge = scale_ridge(cwt_stack(y), n // 2)
        assert np.all(np.abs(ridge["ipos"] - n // 2) <= 3)

    def test_interpolated_lopt(self):
        n, sigma = 4000, 15.0
        y, _ = gaussian_signal(n, [(n / 2, sigma, 2.0)])
        ridge = scale_ridge(cwt_stack(y), n // 2)
        lopt, wmax, censored = interpolate_optimal_scale(
            ridge["scale"].to_numpy(), ridge["amplitude"].to_numpy()
        )
        assert not censored
        # the log-parabola is a quick estimator with O(10 %) bias on the
        # dyadic grid (asymmetric response); reconstruct_areas uses the
        # exact-shape fit instead, which is bias-free
        assert lopt == pytest.approx(optimal_scale(sigma), rel=0.12)

    def test_censored_when_range_too_narrow(self):
        n, sigma = 4000, 60.0  # l_opt = 134 > max tested scale
        y, _ = gaussian_signal(n, [(n / 2, sigma, 1.0)])
        stack = cwt_stack(y, scales=[1.0, 2.0, 4.0, 8.0])
        df = reconstruct_areas(stack, [n // 2])
        assert bool(df["censored"][0])


# ---------------------------------------------------------------------------
# area reconstruction
# ---------------------------------------------------------------------------
class TestReconstruction:
    def test_clean_multi_gaussian(self):
        n = 20000
        peaks = [(3000.0, 5.0, 1.0), (8000.0, 12.0, 5.0), (15000.0, 30.0, 0.5)]
        y, areas = gaussian_signal(n, peaks)
        pk = monotonepeak(y, mfilt=0, zero=0.0)
        assert len(pk) == 3
        df = reconstruct_areas(cwt_stack(y), pk.icenter)
        assert not df["censored"].any()
        # the closed-form response fit is exact for pure Gaussians
        for true_a, got_a, (_, s, _) in zip(areas, df["area"], peaks):
            assert got_a == pytest.approx(true_a, rel=5e-3), f"sigma={s}"

    def test_noisy_multi_gaussian(self):
        rng = np.random.default_rng(11)
        n = 20000
        peaks = [(4000.0, 8.0, 2.0), (12000.0, 20.0, 1.0)]
        y, areas = gaussian_signal(n, peaks)
        y = y + 0.02 * rng.standard_normal(n)
        pk = monotonepeak(y, mfilt=15, sort="descend", maxpeak=2, keeporder=True)
        df = reconstruct_areas(cwt_stack(y), pk.icenter)
        for true_a, got_a in zip(areas, df["area"]):
            assert got_a == pytest.approx(true_a, rel=0.10)

    def test_chromatogram_fixture(self):
        """End-to-end on the synthetic-chromatogram oracle input: noisy,
        drifting baseline, 6 peaks incl. an overlapped doublet."""
        with open(os.path.join(HERE, "fixtures", "tools_oracle_input.json")) as f:
            inp = json.load(f)
        x = np.asarray(inp["chromatogram"]["x"])
        y = np.asarray(inp["chromatogram"]["y"])
        dx = float(np.median(np.diff(x)))
        true = inp["chromatogram"]["true_peaks"]
        pk = monotonepeak(y, x, mfilt=15, sort="descend", maxpeak=6, keeporder=True)
        assert len(pk) == 6
        df = reconstruct_areas(
            cwt_stack(y, scales=[2.0**k for k in range(0, 8)]),
            pk.icenter,
            dx=dx,
        )
        true_areas = np.array(
            [
                tp["amplitude"]
                * 0.6006
                * tp["width"]
                / np.sqrt(2.0)
                * np.sqrt(2 * np.pi)
                for tp in true
            ]
        )
        is_doublet = np.array([abs(tp["center"] - 12.3) < 1.0 for tp in true])
        got = df["area_x"].to_numpy()
        # isolated peaks: within 2 % despite noise and drifting baseline
        for ta, ga, tp in zip(
            true_areas[~is_doublet],
            got[~is_doublet],
            np.asarray(true, dtype=object)[~is_doublet],
        ):
            assert ga == pytest.approx(ta, rel=0.02), f"peak at {tp['center']}"
        # coeluted doublet (thesis limitation #1): individual areas biased,
        # but mass is conserved — the sum stays within 10 %
        assert got[is_doublet].sum() == pytest.approx(
            true_areas[is_doublet].sum(), rel=0.10
        )
        # the weaker, ridge-truncated member is flagged as censored
        assert bool(df["censored"].to_numpy()[is_doublet][1])

    def test_x_units(self):
        n, sigma = 4000, 10.0
        y, areas = gaussian_signal(n, [(n / 2, sigma, 1.0)])
        x = np.linspace(0.0, 40.0, n)
        dx = float(np.median(np.diff(x)))
        df = reconstruct_areas(cwt_stack(y), [n // 2], x=x)
        assert df["area_x"][0] == pytest.approx(df["area"][0] * dx, rel=1e-12)


# ---------------------------------------------------------------------------
# calibration curve (Fig 5.48)
# ---------------------------------------------------------------------------
class TestCalibrationCurve:
    def test_optimal_locus_bijective(self):
        df = calibration_curve()
        loc = df[df["optimal"]].sort_values("sigma")
        # strictly monotone in both coordinates -> bijective
        assert np.all(np.diff(loc["Ag_over_h"]) > 0)
        assert np.all(np.diff(loc["Acwt_over_h"]) > 0)

    def test_locus_concave_in_log(self):
        # thesis: "la concavité de la courbe de calibration" — verify on the
        # locus in log-log coordinates
        df = calibration_curve(sigmas=np.geomspace(2, 64, 50))
        loc = df[df["optimal"]].sort_values("sigma")
        lx = np.log(loc["Ag_over_h"].to_numpy())
        ly = np.log(loc["Acwt_over_h"].to_numpy())
        slope = np.diff(ly) / np.diff(lx)
        assert np.all(np.diff(slope) <= 1e-9)  # non-increasing slope

    def test_lobe_area_is_scale_pure(self):
        # A_cwt / W(0) = 2 s / sqrt(e): measure on the exact engine
        n, sigma, scale = 8000, 20.0, 40.0
        y, _ = gaussian_signal(n, [(n / 2, sigma, 1.0)])
        w = fftconvolve(y, ricker_kernel(scale), mode="same")
        c = n // 2
        right = np.nonzero(w[c:] < 0)[0][0]
        _trapz = getattr(np, "trapezoid", np.trapz)
        lobe = _trapz(w[c - right : c + right + 1])
        s = np.hypot(sigma, scale)
        assert lobe / w[c] == pytest.approx(2.0 * s / np.sqrt(np.e), rel=0.02)


def test_plot_calibration_smoke():
    import matplotlib

    matplotlib.use("Agg")
    from sig2dna_core.tools.cwtarea import plot_calibration

    ax = plot_calibration()
    assert ax is not None
