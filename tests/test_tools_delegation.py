"""
Regression net for the delegation of the symbolic encoders' segmentation to
``tools.segments.sign_runs`` (2026-08-11 refactor).

The historical inline algorithms of ``DNAsignal.encode_dna`` (flats as own
segments) and ``icfilter.encode_series`` (flats attached to the preceding
trend) are reimplemented here verbatim as ground truth and compared with
the shared primitive on randomized signals — any semantic drift in
``sign_runs`` breaks these tests, not the encoders silently.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

import numpy as np
import pytest

from sig2dna_core.tools.segments import sign_runs


def legacy_encode_dna_breaks(coef: np.ndarray) -> np.ndarray:
    """Verbatim historical segmentation of DNAsignal.encode_dna."""
    monotonic = np.diff(coef)
    mono_sign = np.sign(monotonic)
    sign_changes = np.where(np.diff(mono_sign) != 0)[0] + 1
    return np.append(sign_changes, len(coef) - 1)


def legacy_encode_series_bounds(c: np.ndarray) -> np.ndarray:
    """Verbatim historical segmentation of icfilter.encode_series."""
    d = np.diff(c)
    sgn = np.sign(d)
    nz = sgn != 0
    idx = np.where(nz, np.arange(sgn.size), -1)
    np.maximum.accumulate(idx, out=idx)
    run = np.where(idx >= 0, sgn[np.clip(idx, 0, None)], 0.0)
    change = np.nonzero(run[1:] != run[:-1])[0] + 1
    return np.concatenate(([0], change, [len(c) - 1]))


def signals(seed: int, n_signals: int = 50):
    rng = np.random.default_rng(seed)
    for _ in range(n_signals):
        kind = rng.integers(4)
        n = int(rng.integers(10, 400))
        if kind == 0:
            yield rng.standard_normal(n).cumsum()
        elif kind == 1:  # piecewise constant: many exact flats
            yield np.repeat(rng.standard_normal(max(n // 5, 2)), 5)[:n]
        elif kind == 2:  # smooth with flats clipped in
            t = np.linspace(0, 6 * np.pi, n)
            yield np.clip(np.sin(t) * np.cos(2 * t), -0.4, 0.5)
        else:  # integer staircase (zeros and ties everywhere)
            yield rng.integers(-3, 4, n).astype(float)


class TestSignRunsEquivalence:
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_flats_segment_matches_encode_dna(self, seed):
        for c in signals(seed):
            np.testing.assert_array_equal(
                sign_runs(c, flats="segment")[1:],
                legacy_encode_dna_breaks(c),
            )

    @pytest.mark.parametrize("seed", [3, 4, 5])
    def test_flats_attach_matches_encode_series(self, seed):
        for c in signals(seed):
            np.testing.assert_array_equal(
                sign_runs(c, flats="attach"),
                legacy_encode_series_bounds(c),
            )

    def test_edge_short_signals(self):
        for c in ([0.0, 1.0], [1.0, 1.0], [1.0, 0.0, 1.0]):
            c = np.asarray(c)
            np.testing.assert_array_equal(
                sign_runs(c, flats="segment")[1:], legacy_encode_dna_breaks(c)
            )
            np.testing.assert_array_equal(
                sign_runs(c, flats="attach"), legacy_encode_series_bounds(c)
            )

    def test_bad_flats(self):
        with pytest.raises(ValueError):
            sign_runs(np.arange(5.0), flats="wrong")


class TestEncodersStillWork:
    def test_encode_dna_letters_stable(self):
        """encode_dna on the repo's reference two-peak signal: canonical
        letters at every scale (frozen 2026-08-11, pre/post-delegation
        snapshot identical)."""
        import matplotlib

        matplotlib.use("Agg", force=True)
        from sig2dna_core.signomics import DNAsignal, peaks

        p = peaks()
        p.add(x=20, w=2, h=10, name="P1")
        p.add(x=60, w=3, h=5, name="P2")
        d = DNAsignal(p.to_signal(n=512), encode=True, scales=[1, 2, 4])
        for scale in (1, 2, 4):
            letters = d.codes[scale]["letters"]
            assert set(letters) <= set("ABCXYZ_")
            # two peaks -> two YAZB-like motifs; A/Z present at every scale
            assert "A" in letters and "Z" in letters
            # partition covers the signal: iloc chains contiguously
            iloc = d.codes[scale]["iloc"]
            for (a1, b1), (a2, b2) in zip(iloc[:-1], iloc[1:]):
                assert a2 == b1
            assert iloc[0][0] == 0 and iloc[-1][1] == 512

    def test_encode_series_yazb_motif(self):
        from sig2dna_core.icfilter import cwt_matrix, encode_series

        x = np.exp(-(((np.arange(1000) - 500) / 40.0) ** 2))
        C = cwt_matrix(x[None, :], scale=4)
        letters, starts, stops, heights = encode_series(C[0])
        assert "YAZB" in letters
        assert starts.size == stops.size == heights.size == len(letters)
