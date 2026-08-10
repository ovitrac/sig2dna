"""
Tests for sig2dna_core.icfilter — per-ion encoding with Poisson rejection
(thesis Prop. 5.1 port). Synthetic fixtures only.
"""

from __future__ import annotations

import numpy as np

from sig2dna_core.icfilter import (
    encode_ic_matrix,
    encode_series,
    noise_threshold,
    sen_filter,
)


def gaussian(x, x0, w, h):
    return h * np.exp(-(((x - x0) / (0.6006 * w)) ** 2))


def test_encode_series_letters_on_ricker_of_gaussian():
    # Ricker transform of a Gaussian ~ negative-lobe / peak / negative-lobe:
    # its monotone letters must contain the YAZB motif
    import pywt

    x = np.arange(1024, dtype=float)
    y = gaussian(x, 512, 20, 1000.0)
    c, _ = pywt.cwt(y, [8], "mexh")
    letters, start, stop, height = encode_series(c[0])
    assert "YAZB" in letters
    assert len(letters) == len(start) == len(stop) == len(height)


def test_sen_filter_rules():
    # realistic context: peak motif enclosed in blanks (as after segmentation
    # of a real ion channel, where baseline plateaus produce '_' runs)
    letters = "_YAZB_"
    start = np.array([0, 5, 10, 20, 30, 40])
    stop = np.array([5, 10, 20, 30, 40, 45])
    height = np.array([0.1, -5.0, 10.0, -10.0, 5.0, 0.1])
    # low threshold: the whole motif survives
    thr = np.full(46, 1.0)
    assert sen_filter(letters, start, stop, height, thr) == "_YAZB_"
    # high threshold: A/Z rejected, enclosed B/Y runs blanked -> full blank-out
    thr = np.full(46, 100.0)
    assert set(sen_filter(letters, start, stop, height, thr)) == {"_"}


def test_sen_filter_restores_flanks():
    # blanks adjacent to surviving A/Z runs are restored (YAZB reconstruction)
    letters = "_YAZB_"
    start = np.array([0, 5, 10, 20, 30, 40])
    stop = np.array([5, 10, 20, 30, 40, 45])
    height = np.array([0.1, -5.0, 50.0, -50.0, 5.0, 0.1])
    thr = np.full(46, 10.0)  # A/Z heights 50 >> 10 survive
    out = sen_filter(letters, start, stop, height, thr)
    assert "AZ" in out and out[1] == "Y" and out[4] == "B"


def test_poisson_rejection_separates_peak_from_noise_ion():
    rng = np.random.default_rng(7)
    n = 4096
    noise_ion = rng.poisson(50, n).astype(float)
    peak_ion = noise_ion.copy() + gaussian(np.arange(n), 2000, 30, 2000.0)
    y = np.vstack([noise_ion, peak_ion])
    enc = encode_ic_matrix(y, mz=np.array([100.0, 149.0]), scale=8)
    surv = enc.survivors()
    kept_noise = enc.letters[0].strip("_")
    kept_peak = enc.letters[1].strip("_")
    assert kept_peak, "true peak must survive the filter"
    assert 149.0 in surv
    # the noise-only ion must be (almost) fully blanked
    az_noise = enc.letters[0].count("A") + enc.letters[0].count("Z")
    az_noise_raw = enc.letters_raw[0].count("A") + enc.letters_raw[0].count("Z")
    assert az_noise <= 0.05 * az_noise_raw
    assert "A" in kept_peak and "Z" in kept_peak


def test_threshold_scales_with_noise_level():
    rng = np.random.default_rng(1)
    low = rng.normal(0, 1.0, (1, 4096))
    high = rng.normal(0, 10.0, (1, 4096))
    import pywt

    t_low = noise_threshold(pywt.cwt(low, [1], "mexh", axis=-1)[0][0])
    t_high = noise_threshold(pywt.cwt(high, [1], "mexh", axis=-1)[0][0])
    assert np.median(t_high) > 5 * np.median(t_low)


def test_text_entropy_matches_reference_examples():
    # textentropy.m: 'abcdefghijklmnop' -> 4.00 (m=1), 3.91 (m=2)
    from sig2dna_core.icfilter import text_entropy

    assert abs(text_entropy("abcdefghijklmnop", 1) - 4.0) < 1e-9
    assert abs(text_entropy("abcdefghijklmnop", 2) - np.log2(15)) < 1e-9
    assert text_entropy("aaaa") == 0.0
    assert text_entropy("", 1) == 0.0


def test_exclusive_entropy_distance():
    from sig2dna_core.icfilter import align_invariants, exclusive_entropy_distance

    # identical sequences: everything invariant, distance 0
    s = "_Y_AZ_B_YAZB__"
    assert exclusive_entropy_distance(s, s) < 1e-12
    # shared motif + exclusive extra peak: distance strictly positive,
    # and smaller than for fully unrelated content
    a = "____YAZB____"
    b = "____YAZB__YAZB__"
    c = "CCCCCCCCCCCC"
    d_shared = exclusive_entropy_distance(a, b)
    d_unrelated = exclusive_entropy_distance(a, c)
    assert 0 <= d_shared < d_unrelated
    inv = align_invariants(a, b)
    assert "YAZB" in inv  # the common peak aligns and cancels
