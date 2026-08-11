"""
Tests for sig2dna_core.tools.segments (monotone / filtzero port).

Covers the documented Matlab examples and every edge case of the
monotone.m issue history (revs. 2021-03-17, 2021-03-19, 2021-04-16).

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

import numpy as np
import pytest

from sig2dna_core.tools.segments import (MonotoneFull, filtzero, monotone,
                                         monotone_full, to_letters)

EPS = np.finfo(float).eps


def reference_signal():
    """The BASIC EXAMPLE of monotone.m."""
    t = np.linspace(0, 14 * np.pi, 10_000)
    y = np.minimum(0.4, np.maximum(-0.5, np.sin(t) * np.cos(2 * t) * np.sin(3 * t)))
    return t, y


# ---------------------------------------------------------------------------
# monotone — single type
# ---------------------------------------------------------------------------
class TestMonotoneSingle:
    def test_simple_ramp(self):
        seg = monotone(np.array([0.0, 1.0, 2.0, 3.0]), "+")
        assert len(seg) == 1
        assert seg.start[0] == 0
        assert seg.width[0] == 4
        assert seg.height[0] == pytest.approx(3.0)
        assert seg.stop[0] == 3

    def test_v_shape(self):
        x = np.array([2.0, 1.0, 0.0, 1.0, 2.0])
        up = monotone(x, "+")
        dn = monotone(x, "-")
        assert (up.start[0], up.stop[0]) == (2, 4)
        assert (dn.start[0], dn.stop[0]) == (0, 2)
        assert up.height[0] == pytest.approx(2.0)
        assert dn.height[0] == pytest.approx(-2.0)

    def test_constant_segments(self):
        x = np.array([0.0, 0.0, 0.0, 1.0, 1.0])
        seg = monotone(x, "0", zero=1e-9)
        # two constant runs: samples 0-2 and 3-4
        assert len(seg) == 2
        assert list(seg.start) == [0, 3]
        assert list(seg.stop) == [2, 4]

    def test_dead_zone(self):
        # oscillation below the dead zone is constant
        x = np.array([0.0, 1e-12, 0.0, 1e-12, 1.0])
        seg = monotone(x, "+", zero=1e-9)
        assert len(seg) == 1
        assert seg.start[0] == 3

    def test_no_segment(self):
        seg = monotone(np.array([1.0, 1.0, 1.0]), "+")
        assert len(seg) == 0

    def test_bad_kind(self):
        with pytest.raises(ValueError):
            monotone(np.array([0.0, 1.0]), "x")

    def test_signs(self):
        x = np.array([0.0, 1.0, 0.0, 1.0])
        both = monotone(x, "+-")
        assert list(both.sign) == [1, -1, 1]


# ---------------------------------------------------------------------------
# monotone — '+-' merge and keywords
# ---------------------------------------------------------------------------
class TestMonotoneMerge:
    def test_sorted_by_start(self):
        _, y = reference_signal()
        seg = monotone(y, "+-")
        assert np.all(np.diff(seg.start) >= 0)

    def test_leftpriority_shortens_by_one(self):
        x = np.array([0.0, 1.0, 2.0, 1.0, 0.0])
        plain = monotone(x, "+-")
        left = monotone(x, "+-", leftpriority=True)
        assert list(plain.width) == [3, 3]
        assert list(left.width) == [2, 2]
        # heights are those of the full runs, unchanged
        assert list(left.height) == list(plain.height)

    def test_leftpriority_no_zero_width(self):
        # workaround case of rev. 2021-04-16
        for x in (
            [-1.0, 0.0, 1.0, 0.0, 1.0, 2.0],
            [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
        ):
            seg = monotone(np.array(x), "+-", leftpriority=True)
            assert np.all(seg.width >= 1), x

    def test_exclusive_keywords(self):
        with pytest.raises(ValueError):
            monotone(np.array([0.0, 1.0]), "+-", leftpriority=True, nooverlap=True)


# ---------------------------------------------------------------------------
# monotone_full — partition and classes
# ---------------------------------------------------------------------------
class TestMonotoneFull:
    def test_coverage_reference_signal(self):
        _, y = reference_signal()
        full = monotone_full(y)
        assert full.fullclass.size == y.size
        assert np.sum(full.width) == y.size
        assert 1 <= full.nclass <= 16

    @pytest.mark.parametrize(
        "x",
        [
            [0.0, 0.0, 1.0, 1.0, 0.0, 0.0],  # issue 2021-03-17
            [0.0, -2 * EPS, 1.0],  # issue 2021-03-19
            [0.0, -0.1, -1.0, 1e6],  # issue 2021-03-19 (2)
            [-1.0, 0.0, 1.0, 0.0, 1.0, 2.0],  # issue 2021-04-16
            [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
        ],
    )
    def test_matlab_issue_history(self, x):
        full = monotone_full(np.array(x))
        assert isinstance(full, MonotoneFull)
        assert full.fullclass.size == len(x)
        assert np.sum(full.width) == len(x)

    def test_partition_is_contiguous(self):
        rng = np.random.default_rng(0)
        for _ in range(20):
            x = rng.standard_normal(200).cumsum()
            full = monotone_full(x)
            assert np.sum(full.width) == x.size
            # segments sorted, non-overlapping starts
            assert np.all(np.diff(full.start) >= 0)

    def test_attributes_consistency(self):
        _, y = reference_signal()
        full = monotone_full(y)
        attr = full.attributes.to_numpy()
        # every class index maps to a row
        assert full.seg_class.max() < attr.shape[0]
        # monotone segments have nonzero height
        mono = attr[full.seg_class, 0]
        assert np.all(np.abs(full.height[mono]) > 0)

    def test_area_definition(self):
        _, y = reference_signal()
        full = monotone_full(y)
        np.testing.assert_allclose(full.area, full.width * full.height / 2.0)

    def test_too_short(self):
        with pytest.raises(ValueError):
            monotone_full(np.array([0.0, 1.0]))


# ---------------------------------------------------------------------------
# letters (documented example of monotone.m)
# ---------------------------------------------------------------------------
class TestLetters:
    def test_reference_signal_letters(self):
        _, y = reference_signal()
        full = monotone_full(y)
        s, fs = to_letters(full)
        assert len(s) == len(full)
        assert len(fs) == y.size
        assert set(s) <= set("ABCXYZ_")

    def test_lowercase_rule(self):
        _, y = reference_signal()
        full = monotone_full(y)
        s, fs = to_letters(full, lowercase_below_median=True)
        assert any(c.islower() for c in s)
        assert len(fs) == y.size


# ---------------------------------------------------------------------------
# filtzero
# ---------------------------------------------------------------------------
class TestFiltzero:
    def test_identity_below_2(self):
        x = np.arange(10.0)
        np.testing.assert_array_equal(filtzero(x, 1), x)

    def test_zero_phase_preserves_peak_position(self):
        x = np.exp(-((np.arange(1000) - 500.0) ** 2) / (2 * 30.0**2))
        xf = filtzero(x, 21)
        assert abs(int(np.argmax(xf)) - 500) <= 1

    def test_smooths_noise(self):
        rng = np.random.default_rng(1)
        x = np.sin(np.linspace(0, 2 * np.pi, 2000)) + 0.1 * rng.standard_normal(2000)
        xf = filtzero(x, 50)
        clean = np.sin(np.linspace(0, 2 * np.pi, 2000))
        assert np.std(xf - clean) < np.std(x - clean) / 2

    def test_columnwise(self):
        x = np.column_stack([np.arange(300.0), np.arange(300.0)[::-1]])
        xf = filtzero(x, 5)
        assert xf.shape == x.shape

    def test_too_short_raises(self):
        with pytest.raises(ValueError):
            filtzero(np.arange(5.0), 10)

    def test_default_m(self):
        x = np.sin(np.linspace(0, 6.28, 5000))
        xf = filtzero(x)  # m = ceil(5000/100) = 50
        assert xf.shape == x.shape


# ---------------------------------------------------------------------------
# control plots (smoke tests, Agg backend)
# ---------------------------------------------------------------------------
class TestControlPlots:
    def test_plot_segments(self):
        import matplotlib

        matplotlib.use("Agg")
        from sig2dna_core.tools.segments import plot_segments

        _, y = reference_signal()
        ax = plot_segments(y[:500], kind="+")
        assert ax is not None

    def test_plot_classes(self):
        import matplotlib

        matplotlib.use("Agg")
        from sig2dna_core.tools.segments import plot_classes

        _, y = reference_signal()
        ax = plot_classes(y[:500])
        assert ax is not None
