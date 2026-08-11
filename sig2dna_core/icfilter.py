"""
Per-ion (IC-level) symbolic encoding with Poisson noise rejection —
faithful port of the thesis reference implementation.

Provenance of the method (evidence chain):

- Equations: PhD thesis J. Kermorvant (https://theses.hal.science/tel-04194172),
  §5.1.3 (Poisson detector model) and Proposition 5.1 (p. 123): a zero-crossing
  segment is retained iff its height satisfies  h >= 2*sqrt(10 * lambda * dt),
  where lambda*dt is the local Poisson noise level of the ion channel.
- Reference implementation: ``Circular_Daleth`` (GPL v3, J. Kermorvant),
  functions ``samples_compare.m`` / ``sen_filter.m`` / ``monotonetranslate.m``:
  the noise level is estimated as the **moving variance (window 501) of the
  scale-1 Ricker CWT coefficients** of the same ion channel (at scale 1 the
  transform contains essentially noise), i.e. ``filt = 2*sqrt(10*movvar(cfs_1))``,
  and the rejection/restoration rules on the letter string are:
    1. blank A/Z whose |height| < max(filt[start], filt[stop]);
    2. blank any [BCXY]+ run enclosed between blanks;
    3. blank isolated A/Z enclosed between blanks;
    4. restore the original letter for blanks immediately preceding A/B/C or
       immediately following X/Y/Z (rebuilds the YAZB flanks of true peaks).
  This module re-derives those rules in Python (no GPL code copied; the rules
  are method facts, published in the thesis).
- Letters follow the same 6+1 alphabet as ``signomics.DNAsignal._get_letter``:
  A (up, crosses zero), Z (down, crosses zero), B/Y (up/down, negative side),
  C/X (up/down, positive side), '_' blank.

The per-ion architecture is mandatory: on a TIC the summation of hundreds of
ion channels destroys the counting statistics the rejection rule relies on
(verified empirically on D3 data, 2026-08-10).

Ion chromatograms are concatenated in time ("IC1 + IC2 + ..." with '+' the
time-concatenation operator) using '$' as an ion-separator character, matching
the thesis convention (p. 94).

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pywt
from scipy.ndimage import uniform_filter1d
from scipy.signal import fftconvolve

from .tools.segments import sign_runs

__all__ = [
    "ICEncoding",
    "ricker_kernel",
    "cwt_matrix",
    "moving_variance",
    "noise_threshold",
    "encode_series",
    "sen_filter",
    "encode_ic_matrix",
    "text_entropy",
    "align_invariants",
    "exclusive_entropy_distance",
]

SEPARATOR = "$"  # end-of-ion character (thesis p. 94)


# --------------------------------------------------------------------- CWT
def ricker_kernel(scale: float) -> np.ndarray:
    """Analytically sampled Ricker (mexh) wavelet at ``scale`` (in samples),
    L2-normalized, support truncated at +/- 8*scale."""
    n = int(np.ceil(8 * float(scale)))
    t = np.arange(-n, n + 1, dtype=np.float64)
    a = float(scale)
    return (
        (2.0 / (np.sqrt(3.0 * a) * np.pi**0.25))
        * (1 - (t / a) ** 2)
        * np.exp(-(t**2) / (2 * a**2))
    )


def cwt_matrix(
    y: np.ndarray, scale: float, wavelet: str = "mexh", engine: str = "exact"
) -> np.ndarray:
    """Ricker CWT of each row of ``y`` (n_ions x n_time) at one scale.

    Engines
    -------
    ``"exact"`` (default): FFT convolution with the analytically sampled
    Ricker kernel — valid at every scale, and the faster path (cost is
    independent of scale). Only ``wavelet="mexh"`` is implemented; any other
    wavelet raises ``NotImplementedError`` (fail closed).

    ``"pywt"`` (deprecated, kept only to reproduce pre-0.53 results):
    delegates to ``pywt.cwt``, which samples the integrated wavelet on a
    2**10 grid with floor indexing. That sampling is exact only for scales
    with 64/scale integer (1, 2, 4, 8, 16, 32, 64); at any other scale the
    quantization staircase injects high-frequency ripple that silently
    corrupts the monotone-segment structure downstream. Such scales
    therefore raise ``ValueError`` on this engine. Scheduled for removal.

    Coefficients from the two engines are structurally equivalent at the
    grid-exact scales (identical extremum structure) but not bit-identical;
    never mix engines inside one encoding chain.
    """
    y = np.asarray(y, dtype=np.float64)
    if engine == "exact":
        if wavelet != "mexh":
            raise NotImplementedError(
                f"engine='exact' implements the Ricker (mexh) wavelet only, "
                f"got {wavelet!r}"
            )
        return fftconvolve(y, ricker_kernel(scale)[None, :], mode="same", axes=-1)
    if engine == "pywt":
        if not float(64.0 / scale).is_integer():
            raise ValueError(
                f"engine='pywt' is only valid at grid-exact scales "
                f"(64/scale integer); scale {scale} would be silently "
                f"corrupted by the wavelet-grid quantization of pywt.cwt "
                f"— use engine='exact'"
            )
        coefs, _ = pywt.cwt(y, [scale], wavelet, axis=-1)
        return coefs[0]
    raise ValueError(f"unknown engine {engine!r} (use 'exact' or 'pywt')")


def moving_variance(x: np.ndarray, window: int = 501) -> np.ndarray:
    """Centered moving variance along the last axis (Matlab ``movvar`` analog).

    Edge windows are handled by nearest-padding (Matlab shrinks the window;
    the difference affects only the first/last ~window/2 points and is
    documented as an accepted deviation)."""
    m = uniform_filter1d(x, size=window, axis=-1, mode="nearest")
    m2 = uniform_filter1d(x * x, size=window, axis=-1, mode="nearest")
    return np.maximum(m2 - m * m, 0.0)


def noise_threshold(
    cfs_scale1: np.ndarray, window: int = 501, k: float = 2.0, factor: float = 10.0
) -> np.ndarray:
    """Pointwise rejection threshold  k*sqrt(factor * movvar(scale-1 CWT)).

    Thesis Prop. 5.1 with lambda*dt estimated from the scale-1 coefficients
    (reference: samples_compare.m: ``2*sqrt(10*movvar(cfs1f, [250,250]))``)."""
    return k * np.sqrt(factor * moving_variance(cfs_scale1, window))


# ----------------------------------------------------------- segmentation
def encode_series(c: np.ndarray, tol: float = 1e-12):
    """Monotone-segment letters of one transformed series (signomics rules).

    Returns (letters, start_idx, stop_idx, height) — one entry per segment.
    Letter rules identical to ``signomics.DNAsignal._get_letter``."""
    c = np.asarray(c, dtype=np.float64)
    # segment boundaries where the derivative sign changes; plateaus attach
    # to the preceding trend — delegated to the shared primitive
    bounds = sign_runs(c, flats="attach")
    letters, starts, stops, heights = [], [], [], []
    for a, b in zip(bounds[:-1], bounds[1:]):
        if a == b:
            continue
        s, e = c[a], c[b]
        s0 = 0.0 if abs(s) < tol else s
        e0 = 0.0 if abs(e) < tol else e
        if s0 < 0 and e0 > 0:
            L = "A"
        elif s0 > 0 and e0 < 0:
            L = "Z"
        elif s0 <= 0 and e0 <= 0:
            L = "B" if e0 > s0 else ("Y" if e0 < s0 else "_")
        elif s0 >= 0 and e0 >= 0:
            L = "C" if e0 > s0 else ("X" if e0 < s0 else "_")
        else:
            L = "_"
        letters.append(L)
        starts.append(a)
        stops.append(b)
        heights.append(e - s)
    return (
        "".join(letters),
        np.asarray(starts, dtype=np.int64),
        np.asarray(stops, dtype=np.int64),
        np.asarray(heights, dtype=np.float64),
    )


# ------------------------------------------------------------- Prop. 5.1
_RUN_BETWEEN_BLANKS = re.compile(r"(?<=_)[BCXY]+(?=_)")
_LONE_AZ = re.compile(r"(?<=_)[AZ](?=_)")
_RESTORE_BEFORE_UP = re.compile(r"_(?=[ABC])")
_RESTORE_AFTER_DOWN = re.compile(r"(?<=[XYZ])_")


def sen_filter(
    letters: str,
    start: np.ndarray,
    stop: np.ndarray,
    height: np.ndarray,
    threshold: np.ndarray,
) -> str:
    """Literal port of the thesis rejection/restoration rules (sen_filter.m)."""
    sen = np.array(list(letters))
    senf = sen.copy()
    is_az = (senf == "A") | (senf == "Z")
    if is_az.any():
        thr = np.maximum(threshold[start], threshold[stop])
        senf[is_az & (np.abs(height) < thr)] = "_"
    s = "".join(senf)
    # rule 2: [BCXY]+ runs enclosed in blanks (string ends count as enclosed
    # only when actually flanked by '_', as in the reference implementation)
    s = _RUN_BETWEEN_BLANKS.sub(lambda m: "_" * len(m.group(0)), s)
    # rule 3: isolated A/Z between blanks
    s = _LONE_AZ.sub("_", s)
    # rule 4: restore boundary letters (YAZB flank reconstruction)
    out = list(s)
    for m in _RESTORE_BEFORE_UP.finditer(s):
        out[m.start()] = letters[m.start()]
    for m in _RESTORE_AFTER_DOWN.finditer(s):
        out[m.start()] = letters[m.start()]
    return "".join(out)


# -------------------------------------------- aligned mutual-information
def text_entropy(s: str, m: int = 1) -> float:
    """Shannon entropy (bits) of the overlapping m-gram tokens of ``s``.

    Port of ``textentropy.m`` (Circular_Daleth): for m>1 the reference's
    strided base-128 products enumerate exactly the overlapping m-grams
    (its own comment: 'AB|CD|EF|BC|DE'), so the direct m-gram count is the
    faithful equivalent. Empty strings have zero entropy."""
    n = len(s) - m + 1
    if n <= 0:
        return 0.0
    counts: dict = {}
    for k in range(n):
        t = s[k : k + m]
        counts[t] = counts.get(t, 0) + 1
    p = np.array(list(counts.values()), dtype=np.float64) / n
    return float(-(p * np.log2(p)).sum())


def align_invariants(
    a: str, b: str, open_gap: float = -1.0, extend_gap: float = -1.0
) -> str:
    """Global (Needleman-Wunsch) alignment of two letter strings and return
    the concatenated **invariant** characters (exact matches).

    Reference: ``fingerprints_compare.m`` — ``nwalign`` with identity scoring
    matrix (match=1, mismatch=0), gap open/extend penalty 1; invariants are
    the '|' positions of the match line."""
    if not a or not b:
        return ""
    from Bio import Align  # deferred heavy import

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    aln = aligner.align(a, b)[0]
    out = []
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        sa, sb = a[a0:a1], b[b0:b1]
        out.extend(ca for ca, cb in zip(sa, sb) if ca == cb)
    return "".join(out)


def exclusive_entropy_distance(a: str, b: str, m: int = 1) -> float:
    """Thesis Eq. 4.4 (p. 92): d = H(A) + H(B) - 2*H(A∩B), with H(A∩B) the
    m-gram entropy of the aligned invariant characters.

    Shared content (PET is PET: common peaks align and match) cancels out;
    only the mutually exclusive information contributes."""
    h1, h2 = text_entropy(a, m), text_entropy(b, m)
    h12 = text_entropy(align_invariants(a, b), m)
    return h1 + h2 - 2.0 * h12


# ---------------------------------------------------------------- driver
@dataclass
class ICEncoding:
    """Filtered per-ion symbolic encoding of one acquisition at one scale."""

    scale: float
    mz: np.ndarray
    letters: list = field(default_factory=list)  # filtered, one str per ion
    letters_raw: list = field(default_factory=list)  # unfiltered
    n_rejected: int = 0
    n_segments: int = 0

    @property
    def concatenated(self) -> str:
        """IC1 + IC2 + ... with '$' the ion separator (time concatenation)."""
        return SEPARATOR.join(self.letters)

    def survivors(self) -> dict:
        """{mz: filtered letters} for ions with at least one retained letter."""
        return {float(m): s for m, s in zip(self.mz, self.letters) if s.strip("_")}


def encode_ic_matrix(
    y: np.ndarray,
    mz: np.ndarray,
    scale: float,
    window: int = 501,
    k: float = 2.0,
    factor: float = 10.0,
    wavelet: str = "mexh",
    engine: str = "exact",
) -> ICEncoding:
    """Encode every ion chromatogram of ``y`` (n_ions x n_time) at ``scale``
    with Poisson rejection (thesis Prop. 5.1).

    The noise threshold of each ion derives from its own scale-1 CWT
    (scale 1 carries essentially noise); letters are computed on the CWT at
    the analysis scale, then filtered by :func:`sen_filter`. The same CWT
    ``engine`` (see :func:`cwt_matrix`) is used for the threshold and the
    analysis scale — the two stages must never mix engines."""
    y = np.asarray(y, dtype=np.float64)
    cfs1 = cwt_matrix(y, 1.0, wavelet, engine=engine)
    thr = noise_threshold(cfs1, window=window, k=k, factor=factor)
    cfs = cfs1 if scale == 1 else cwt_matrix(y, float(scale), wavelet, engine=engine)

    out = ICEncoding(scale=float(scale), mz=np.asarray(mz))
    for i in range(y.shape[0]):
        letters, start, stop, height = encode_series(cfs[i])
        filtered = sen_filter(letters, start, stop, height, thr[i]) if letters else ""
        out.letters_raw.append(letters)
        out.letters.append(filtered)
        out.n_segments += len(letters)
        out.n_rejected += sum(
            1 for a, b in zip(letters, filtered) if a != b and b == "_"
        )
    return out
