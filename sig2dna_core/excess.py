"""
Excess mass spectrum — the per-ion readout of foreign chemistry.

Given a sample's per-ion symbolic encoding and a set of reference encodings
(e.g., virgin polymers), the *excess spectrum* is the per-ion aligned
exclusive entropy to the nearest reference, minus the references' own
leave-one-out nearest-neighbour null. It reads like a mass spectrum of the
chemistry that should not be there:

- the shared (expected) chemistry cancels in the alignment (Eq. 4.4 of the
  reference thesis: d = H(A) + H(B) - 2 H(A ∩ B));
- the reference null subtraction removes the residual noise floor, with the
  same statistic on both sides (nearest reference for the sample, strict
  leave-one-out nearest reference for the null);
- the integral of the positive part is the screening score; the *pattern*
  is the etiology (which ions carry the foreign information), directly
  submittable to EI spectral libraries for identification.

Local uncertainty envelope from the negative excess
---------------------------------------------------
Genuine chemistry can only ADD unshared information, so the negative side
of the excess spectrum is signal-free by construction: negative excursions
sample the estimator's noise (including column-bleed misalignment), and
reflecting them yields a local, assumption-light significance envelope —
the same logic as target-decoy false-discovery control, with the negative
half-axis as the decoy. Two documented limits: (i) sparse ions have a
skewed null, so the envelope is estimated on a sliding m/z window with a
robust quantile rather than per ion; (ii) the envelope bounds *random*
error only — one-sided systematic biases (e.g., unequal letter budgets
between acquisition campaigns) leave no trace on the negative side and
must be handled by intensive normalization upstream.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from .icfilter import align_invariants, text_entropy

__all__ = ["ExcessSpectrum", "pair_excess", "negative_envelope",
           "excess_spectrum"]


def pair_excess(a: Sequence[str], b: Sequence[str]) -> np.ndarray:
    """Per-ion aligned exclusive entropy between two encodings (lists of
    letter strings, one per ion). Returns a vector, one value per ion."""
    if len(a) != len(b):
        raise ValueError("encodings must cover the same ion set")
    out = np.empty(len(a))
    for i, (sa, sb) in enumerate(zip(a, b)):
        out[i] = (text_entropy(sa) + text_entropy(sb)
                  - 2.0 * text_entropy(align_invariants(sa, sb)))
    return out


def negative_envelope(excess: np.ndarray, window: int = 51,
                      q: float = 0.95) -> np.ndarray:
    """Local significance envelope from the reflected negative excess.

    For each ion, take the negative excess values within a centered window
    of ``window`` ions, reflect them, and return their ``q``-quantile.
    Windows without negatives fall back to the global envelope; an excess
    with no negative values at all returns zeros (nothing to calibrate on,
    every positive value is then reported as significant — callers should
    treat that case with suspicion)."""
    x = np.asarray(excess, dtype=np.float64)
    neg = -x[x < 0]
    glob = float(np.quantile(neg, q)) if neg.size else 0.0
    half = max(1, window // 2)
    env = np.empty(len(x))
    for i in range(len(x)):
        w = x[max(0, i - half):i + half + 1]
        wneg = -w[w < 0]
        env[i] = float(np.quantile(wneg, q)) if wneg.size >= 5 else glob
    return env


@dataclass
class ExcessSpectrum:
    """Result of :func:`excess_spectrum`."""

    mz: np.ndarray                 # ion axis
    excess: np.ndarray             # per-ion excess (bits), signed
    null_mu: np.ndarray            # per-ion LOO nearest-reference null mean
    null_sd: np.ndarray            # per-ion null SD (Bessel)
    envelope: np.ndarray           # local reflected-negative envelope
    nearest: int = 0               # index of the nearest reference
    total: float = 0.0             # summed distance to nearest reference

    @property
    def z(self) -> np.ndarray:
        """Per-ion z vs the reference null (guarded for zero-SD ions)."""
        return self.excess / np.where(self.null_sd > 0, self.null_sd, np.inf)

    @property
    def significant(self) -> np.ndarray:
        """Ions whose excess exceeds the local envelope."""
        return self.excess > self.envelope

    def top(self, n: int = 10) -> list:
        """The n most excessive significant ions as (mz, excess, z)."""
        order = np.argsort(self.excess)[::-1]
        rows = [(float(self.mz[i]), float(self.excess[i]), float(self.z[i]))
                for i in order if self.significant[i]]
        return rows[:n]


def excess_spectrum(
    sample: Sequence[str],
    references: Sequence[Sequence[str]],
    mz: Optional[np.ndarray] = None,
    envelope_window: int = 51,
    envelope_q: float = 0.95,
) -> ExcessSpectrum:
    """Excess mass spectrum of ``sample`` against ``references``.

    Parameters
    ----------
    sample : per-ion letter strings of the sample.
    references : at least two per-ion encodings of reference materials.
    mz : optional ion axis (defaults to 0..n_ions-1).
    envelope_window, envelope_q : see :func:`negative_envelope`.

    The nearest reference is selected globally (minimum summed distance) so
    the excess is measured against one coherent material; the null uses the
    same statistic (strict leave-one-out nearest reference) per ion.
    """
    refs = list(references)
    if len(refs) < 2:
        raise ValueError("at least two references are required for the null")
    n = len(sample)
    if any(len(r) != n for r in refs):
        raise ValueError("all encodings must cover the same ion set")
    if mz is None:
        mz = np.arange(n, dtype=float)

    d_sample = np.array([pair_excess(sample, r) for r in refs])
    totals = d_sample.sum(axis=1)
    nn = int(np.argmin(totals))

    d_rr = {}
    for i in range(len(refs)):
        for j in range(i + 1, len(refs)):
            d_rr[(i, j)] = d_rr[(j, i)] = pair_excess(refs[i], refs[j])
    null_vectors = []
    for i in range(len(refs)):
        others = [j for j in range(len(refs)) if j != i]
        jn = min(others, key=lambda j: d_rr[(i, j)].sum())
        null_vectors.append(d_rr[(i, jn)])
    null = np.array(null_vectors)
    null_mu = null.mean(axis=0)
    null_sd = null.std(axis=0, ddof=1)

    excess = d_sample[nn] - null_mu
    env = negative_envelope(excess, window=envelope_window, q=envelope_q)
    return ExcessSpectrum(mz=np.asarray(mz, dtype=float), excess=excess,
                          null_mu=null_mu, null_sd=null_sd, envelope=env,
                          nearest=nn, total=float(totals[nn]))
