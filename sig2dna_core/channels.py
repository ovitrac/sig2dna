"""
The channel algebra: retention-time-marginal statistics of event
populations.

Chromatography provides separation; the mass channels carry the chemical
information. Marginalizing the event population over the time axis
(the push-forward onto the channel axis) yields the **channel state**

    ``(B, p)``  with  ``B = sum_m N_m``  and  ``p_m = N_m / B``,

where ``N_m`` counts the resolved events on channel ``m``. ``B`` is the
extensive event budget (kept and reported, never normalized away);
``p`` is the channel composition, robust to overall acquisition
sensitivity. Multiplicities are exactly invariant to displacement of
peaks along the run as long as events remain separately detectable —
the limiting failure is peak merging/splitting, not positional
misalignment. Composition removes a global budget factor; the
MAD-of-log-ratios dispersion removes it *provably* (a per-pair budget
ratio is a constant shift of the log-ratios and MAD is shift-invariant).

These invariances hold within compatible acquisition domains; a change
of observation operator (extraction technique, instrument family) can
replace the observable population entirely — measure a bridge before
comparing (``sig2dna_core.bridges``).
"""
from __future__ import annotations

import numpy as np


def channel_counts(events, n_channels: int) -> np.ndarray:
    """Event multiplicity per channel, ``N_m`` (the RT-marginal)."""
    n = np.zeros(n_channels)
    for e in events:
        ion = e.ion if hasattr(e, "ion") else e[0]
        n[ion] += 1
    return n


def channel_state(events, n_channels: int):
    """The frozen channel state ``(B, p)``: extensive budget + composition."""
    n = channel_counts(events, n_channels)
    b = float(n.sum())
    return b, (n / b if b > 0 else n)


def d_comp(n1: np.ndarray, n2: np.ndarray) -> float:
    """Total-variation distance between channel compositions (0..1).
    0 = identical composition; 1 = disjoint channel usage. Invariant to
    the overall budgets of both runs."""
    p1 = n1 / n1.sum() if n1.sum() else n1
    p2 = n2 / n2.sum() if n2.sum() else n2
    return 0.5 * float(np.abs(p1 - p2).sum())


def mad_logratio(n1: np.ndarray, n2: np.ndarray, kmin: int = 3) -> float:
    """Robust dispersion of per-channel log2 multiplicity ratios over
    channels populated (>= kmin) on both sides. Exactly invariant to a
    global multiplicative budget factor (shift invariance of the MAD);
    a genuinely uniform expansion of chemistry is therefore invisible
    here by design — read it in ``B``."""
    m = (n1 >= kmin) & (n2 >= kmin)
    if m.sum() < 20:
        return float("nan")
    r = np.log2(n1[m] / n2[m])
    return float(1.4826 * np.median(np.abs(r - np.median(r))))
