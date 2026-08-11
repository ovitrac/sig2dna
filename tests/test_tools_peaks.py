"""
Tests for sig2dna_core.tools.peaks (monotonepeak port), including
Matlab-oracle parity on the synthetic chromatogram fixture.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

import json
import os

import numpy as np
import pytest

from sig2dna_core.tools.peaks import PeakTable, monotonepeak

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
        inp["chromatogram"]["true_peaks"],
    )


def _arr(v, dtype=float):
    return np.atleast_1d(np.asarray(v, dtype=dtype))


# ---------------------------------------------------------------------------
# functional tests
# ---------------------------------------------------------------------------
class TestMonotonepeak:
    def test_single_gaussian(self):
        # odd grid: a sample sits exactly on the apex (an even symmetric
        # grid yields y[n/2-1] == y[n/2], no rise/fall junction — the
        # Matlab reference behaves identically)
        x = np.linspace(0, 10, 1001)
        y = np.exp(-(((x - 5) / 0.5) ** 2))
        pk = monotonepeak(y, x, mfilt=5)
        assert len(pk) == 1
        assert abs(pk.center[0] - 5.0) < 0.05
        assert pk.absheight[0] == pytest.approx(1.0, abs=1e-3)

    def test_finds_all_true_peaks(self, chromatogram):
        x, y, true_peaks = chromatogram
        pk = monotonepeak(y, x, mfilt=15)
        for tp in true_peaks:
            d = np.min(np.abs(pk.center - tp["center"]))
            assert d < 0.2, f"missed true peak at {tp['center']}"

    def test_sort_and_maxpeak(self, chromatogram):
        x, y, _ = chromatogram
        pk = monotonepeak(y, x, mfilt=15, sort="descend", maxpeak=3)
        assert len(pk) == 3
        assert np.all(np.diff(pk.absheight) <= 0)

    def test_keeporder(self, chromatogram):
        x, y, _ = chromatogram
        pk = monotonepeak(y, x, mfilt=15, sort="descend", maxpeak=4, keeporder=True)
        assert np.all(np.diff(pk.center) > 0)

    def test_multi_mfilt_levels(self, chromatogram):
        x, y, _ = chromatogram
        pk1 = monotonepeak(y, x, mfilt=60)
        pk = monotonepeak(y, x, mfilt=[60, 15])
        assert len(pk) >= len(pk1)
        assert set(np.unique(pk.level)) <= {1, 2}

    def test_flat_signal_empty(self):
        pk = monotonepeak(np.zeros(500), mfilt=5)
        assert len(pk) == 0
        assert isinstance(pk, PeakTable)

    def test_apex_is_rise_fall_junction(self):
        x = np.linspace(0, 10, 500)
        y = np.exp(-(((x - 3) / 0.4) ** 2)) + 0.7 * np.exp(-(((x - 7) / 0.6) ** 2))
        pk = monotonepeak(y, x, mfilt=0)  # no smoothing
        assert len(pk) == 2
        for i in range(2):
            ic = pk.icenter[i]
            assert y[ic] >= y[ic - 1] and y[ic] >= y[ic + 1]

    def test_to_frame(self, chromatogram):
        x, y, _ = chromatogram
        df = monotonepeak(y, x, mfilt=15).to_frame()
        assert {"icenter", "center", "width", "tail", "wall"} <= set(df.columns)


# ---------------------------------------------------------------------------
# Matlab-oracle parity
# ---------------------------------------------------------------------------
def assert_peaks_match(pk: PeakTable, ref: dict) -> None:
    np.testing.assert_array_equal(pk.level, _arr(ref["level"], int))
    np.testing.assert_array_equal(pk.icenter, _arr(ref["icenter"], int) - 1)
    np.testing.assert_array_equal(pk.istart, _arr(ref["istart"], int) - 1)
    np.testing.assert_array_equal(pk.istop, _arr(ref["istop"], int) - 1)
    np.testing.assert_array_equal(pk.ibase, _arr(ref["ibase"], int) - 1)
    for f in (
        "center",
        "start",
        "stop",
        "height",
        "absheight",
        "width",
        "ratioheight",
        "ratiowidth",
    ):
        np.testing.assert_allclose(
            getattr(pk, f), _arr(ref[f]), rtol=1e-9, atol=1e-12, err_msg=f
        )
    assert list(pk.tail) == list(ref["tail"])
    assert list(pk.wall) == list(ref["wall"])


class TestPeakParity:
    def test_single_mfilt(self, oracle, chromatogram):
        x, y, _ = chromatogram
        pk = monotonepeak(y, x, mfilt=15)
        assert_peaks_match(pk, oracle["chromatogram"]["peaks_m15"])

    def test_multi_mfilt(self, oracle, chromatogram):
        x, y, _ = chromatogram
        pk = monotonepeak(y, x, mfilt=[60, 15], sort="descend")
        assert_peaks_match(pk, oracle["chromatogram"]["peaks_multi"])


# ---------------------------------------------------------------------------
# control plot (smoke)
# ---------------------------------------------------------------------------
def test_plot_peaks(chromatogram):
    import matplotlib

    matplotlib.use("Agg")
    from sig2dna_core.tools.peaks import plot_peaks

    x, y, _ = chromatogram
    ax = plot_peaks(y, x, mfilt=15)
    assert ax is not None
