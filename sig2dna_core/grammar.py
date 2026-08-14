"""
Organization/grammar algebra: the layer connecting the symbolic CWT route
(letters -> texts) with the event/channel algebra (events -> counts).

The idea, taken seriously as an information chain rule::

    I(E, T) = I(E) + I(T | E)

the event layer measures WHAT is present (inventory, extensive); the
symbolic text layer, *conditioned on the event layer*, measures only the
residual organization -- HOW the signal is arranged beyond what the
resolved inventory explains. Practically: a reference grammar is trained
on reference texts, each channel text is scored by its code length under
that grammar, the expected code length given the channel's resolved
event count is subtracted, and the residual is standardized per channel
against reference replicates. Positive residuals (``I_O+``) mean more
surprising organization than the reference; negative (``I_O-``) mean
more predictable/regular organization. Both directions are informative:
perturbations may add texture or erase it.

Validity limits (measured on real campaigns; calibrate on your own data):

- **Conditioning must use the event layer**, never statistics of the
  text itself (self-conditioning absorbs the signal).
- **Censor before summing**: per-channel gates from reference-replicate
  residuals; an ungated sum over hundreds of channels accumulates noise
  faster than any real signal.
- **Levels are acquisition-session-sensitive and profiles only somewhat
  less so.** Never compare organization levels across sessions without a
  same-session reference; floors must come from true replicates of each
  state (replicate-mean floors when contrasting replicate means).
- **The estimator requires a sparse-event regime.** Use
  :func:`rho_e` (median resolved events per channel) as the regime
  coordinate. Two saturation modes exist: scale mismatch (too-fine CWT
  scale -- curable by choosing the scale appropriate to the peak width)
  and event-density saturation in chemically rich matrices (NOT cured by
  coarsening; prefer channel counts there).
- Whether the *sign* of an organization change transports across
  sessions while its magnitude does not is an open question; treat sign
  claims as hypotheses until floored in your own design.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np

#: default gap-run punctuation bounds, in symbolic-segment units; the
#: reference implementation derives them as terciles of the gap-run
#: length distribution of the reference corpus -- recalibrate on yours.
GAP_SHORT = 10
GAP_LONG = 32

VOCAB = ("Y", "A", "Z", "B", "C", "X", ".", ":", "|")


def tokenize(text: str, gap_short: int = GAP_SHORT,
             gap_long: int = GAP_LONG) -> list:
    """Letters + gap punctuation from one channel text.

    Runs of the non-coding letter ``_`` collapse into one punctuation
    token: ``.`` (run <= gap_short), ``:`` (<= gap_long), ``|`` (longer).
    Retention time thereby becomes punctuation, not vocabulary: relative
    spacing is kept, absolute position is not."""
    out, run = [], 0
    for c in text:
        if c == "_":
            run += 1
            continue
        if run:
            out.append("." if run <= gap_short
                       else (":" if run <= gap_long else "|"))
            run = 0
        out.append(c)
    if run:
        out.append("." if run <= gap_short
                   else (":" if run <= gap_long else "|"))
    return out


class MarkovGrammar:
    """Order-2 Markov grammar with Laplace smoothing over the token
    vocabulary. ``codelen`` returns the exact Shannon code length (bits)
    of a token sequence under the trained conditional distributions --
    the reference-relative information of the text."""

    def __init__(self, vocab=VOCAB):
        self.vocab_n = len(vocab)
        self.ctx: dict = defaultdict(Counter)

    def train(self, sequences) -> "MarkovGrammar":
        for seq in sequences:
            s = ["^", "^"] + list(seq) + ["$"]
            for i in range(2, len(s)):
                self.ctx[(s[i - 2], s[i - 1])][s[i]] += 1
        return self

    def codelen(self, seq) -> float:
        s = ["^", "^"] + list(seq) + ["$"]
        bits = 0.0
        for i in range(2, len(s)):
            c = self.ctx.get((s[i - 2], s[i - 1]))
            n = c[s[i]] if c else 0
            tot = sum(c.values()) if c else 0
            bits += -math.log2((n + 1.0) / (tot + self.vocab_n + 1.0))
        return bits


def fit_conditioning(codelens: np.ndarray, counts: np.ndarray):
    """Least-squares fit of expected code length given resolved event
    count: L_hat(N) = b0 + b1*N. Fit ONLY on reference texts scored
    leave-one-out; returns (beta, residuals). Conditioning on the event
    layer is what prevents the organization measure from re-counting
    inventory."""
    A = np.column_stack([np.ones_like(counts, dtype=float), counts])
    beta, *_ = np.linalg.lstsq(A, codelens, rcond=None)
    return beta, codelens - A @ beta


def channel_calibration(residuals: np.ndarray):
    """Per-channel gate statistics from reference residuals shaped
    (n_reference_runs, n_channels): mean and (floored) standard
    deviation per channel."""
    mu = residuals.mean(axis=0)
    sd = np.maximum(residuals.std(axis=0, ddof=1), 1.0)
    return mu, sd


def organization_scores(codelens: np.ndarray, counts: np.ndarray,
                        beta: np.ndarray, mu: np.ndarray, sd: np.ndarray,
                        gate: float = 2.0):
    """Signed, censored organization information of one run.

    Returns ``(z, i_o_pos, i_o_neg)``: per-channel standardized residuals
    and the gated sums of positive (excess organization) and negative
    (over-regularity) deviations, in bits. Channels only contribute
    beyond ``gate`` reference standard deviations -- the analog of local
    event thresholds."""
    z = (codelens - (beta[0] + beta[1] * counts) - mu) / sd
    d = z * sd
    return (z, float(d[z > gate].sum()), float(-d[z < -gate].sum()))


def rho_e(counts: np.ndarray) -> float:
    """Density-regime coordinate: median resolved events per channel.

    Low values -> sparse regime, where conditional-code-length
    organization has been shown measurable; high values -> saturated
    regime, where per-channel grammar residuals are noise-dominated and
    channel counts are the appropriate detector. Bracket the boundary on
    your own corpora; do not assume a universal threshold."""
    return float(np.median(counts))


def profile_distance(z1: np.ndarray, z2: np.ndarray) -> float:
    """Correlation distance between two organization z-profiles.
    Level-free: insensitive to a uniform shift of either profile (the
    component measured to be most session-sensitive)."""
    if z1.std() == 0 or z2.std() == 0:
        return 1.0
    return 1.0 - float(np.corrcoef(z1, z2)[0, 1])


def sign_delta(io_target, io_reference) -> int:
    """Sign of the net organization change target - reference, from two
    ``(i_o_pos, i_o_neg)`` pairs: +1 (net enrichment), -1 (net
    regularization), 0. The sign is the weakest -- and possibly the most
    transportable -- organization statement; treat its cross-session
    stability as a hypothesis to test, not a property to assume."""
    d = (io_target[0] - io_target[1]) - (io_reference[0] - io_reference[1])
    return (d > 0) - (d < 0)


class DensityCalibration:
    """Heteroscedastic fallback calibration: mu(N) and sigma(N) fitted by
    event-count bins on pooled reference residuals.

    Calibration hierarchy (in EVIDENTIAL quality, not numerical size)::

        R3  per-channel true-replicate calibration
        R2  reference-population per-channel (robust) + shrinkage to sigma(N)
        R1  this class: density-conditioned sigma(N) proxy
        R0  one global sigma -- didactic visualization only

    A density-conditioned z must never masquerade as replicate-grade
    significance: it is a proxy surprisal score. Label every score with
    its calibration mode; lack of metrology must make claims harder, not
    easier. The representation can be simple; the calibration must not
    be."""

    def __init__(self, n_bins: int = 10):
        self.n_bins = n_bins

    def fit(self, residuals: np.ndarray,
            counts: np.ndarray) -> "DensityCalibration":
        q = np.quantile(counts, np.linspace(0, 1, self.n_bins + 1))
        q[-1] += 1
        b = np.clip(np.searchsorted(q, counts, side="right") - 1,
                    0, self.n_bins - 1)
        self.edges = q
        self.mu = np.array([residuals[b == k].mean()
                            if (b == k).any() else 0.0
                            for k in range(self.n_bins)])
        self.sd = np.array([max(residuals[b == k].std(ddof=1), 1.0)
                            if (b == k).sum() > 1 else 1.0
                            for k in range(self.n_bins)])
        return self

    def z(self, residuals: np.ndarray, counts: np.ndarray) -> np.ndarray:
        b = np.clip(np.searchsorted(self.edges, counts, side="right") - 1,
                    0, self.n_bins - 1)
        return (residuals - self.mu[b]) / self.sd[b]
