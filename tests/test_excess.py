"""Tests for the excess mass spectrum (sig2dna_core.excess)."""

from __future__ import annotations

import numpy as np
import pytest

from sig2dna_core.excess import (ExcessSpectrum, excess_spectrum,
                                 negative_envelope, pair_excess)


def _make_encoding(rng, n_ions=60, base_motifs=3, noise=2):
    """Synthetic per-ion letter strings: shared peak motifs + random noise
    letters, mimicking filtered encodings."""
    enc = []
    for _ in range(n_ions):
        s = ""
        for _ in range(base_motifs):
            s += "__YAZB__"
        for _ in range(rng.integers(0, noise + 1)):
            k = rng.integers(0, len(s))
            s = s[:k] + rng.choice(list("BYCX")) + s[k:]
        enc.append(s)
    return enc


def test_pair_excess_zero_for_identical():
    rng = np.random.default_rng(0)
    a = _make_encoding(rng)
    assert np.allclose(pair_excess(a, a), 0.0)


def test_excess_spectrum_flags_spiked_ions_only():
    rng = np.random.default_rng(1)
    refs = [_make_encoding(rng) for _ in range(4)]
    sample = [s for s in refs[0]]
    spiked = [10, 25, 40]
    for i in spiked:
        # a foreign compound: information-rich extra content on selected
        # ions (a repeated clean motif alone adds letters but little
        # entropy — the metric is deliberately insensitive to it)
        sample[i] = sample[i] + "".join(rng.choice(list("YAZBCX_"))
                                        for _ in range(150))
    spec = excess_spectrum(sample, refs)
    sig = set(np.flatnonzero(spec.significant))
    assert set(spiked) <= sig
    # the spiked ions dominate the ranking
    top = [int(m) for m, _, _ in spec.top(3)]
    assert set(top) == set(spiked)


def test_excess_spectrum_clean_sample_mostly_silent():
    rng = np.random.default_rng(2)
    refs = [_make_encoding(rng) for _ in range(4)]
    sample = _make_encoding(rng)   # same population, no spike
    spec = excess_spectrum(sample, refs)
    # no massive excess anywhere: nothing beyond 3x the envelope scale
    env = np.maximum(spec.envelope, 1e-12)
    assert (spec.excess / env).max() < 3.0


def test_negative_envelope_reflects_noise_scale():
    rng = np.random.default_rng(3)
    x = rng.normal(0.0, 1.0, 600)          # pure symmetric noise
    env = negative_envelope(x, window=101, q=0.95)
    # ~95% of positive noise should sit below the reflected envelope
    pos = x[x > 0]
    covered = (pos < env[x > 0]).mean()
    assert covered > 0.85
    # and the envelope tracks the local scale within a factor ~2 of the
    # true 95th percentile of |N(0,1)| (= 1.96)
    assert 1.0 < np.median(env) < 3.0


def test_excess_spectrum_fails_closed():
    rng = np.random.default_rng(4)
    enc = _make_encoding(rng)
    with pytest.raises(ValueError):
        excess_spectrum(enc, [enc])          # one reference: no null
    with pytest.raises(ValueError):
        excess_spectrum(enc[:10], [enc, enc])  # ion-set mismatch
