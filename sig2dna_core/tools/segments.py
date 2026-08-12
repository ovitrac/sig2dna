"""
Monotone-segment analysis and zero-phase smoothing.

Python port of the reference Matlab pair ``monotone.m`` / ``filtzero.m``
("MS — Molecular Studio" toolbox, O. Vitrac). ``monotone`` dates back to
2001 (Thermique 1.0, CNRS, WoodOx 1.22 lineage) and was forked 2021 for
the thesis of J. Kermorvant (https://theses.hal.science/tel-04194172),
where the chain ``monotone`` → ``filtzero`` → ``monotonepeak`` is the
signal-analysis foundation; the functions were used by the thesis, not
authored by it.

Semantics follow the post-fork head of ``monotone.m`` (rev. 2021-04-16)
verbatim, with two documented deviations:

- all indices are **0-based** (Matlab is 1-based); a segment spans samples
  ``start .. stop`` inclusive with ``stop = start + width - 1``;
- the thesis fork (Daleth_Circular, frozen 2021-04-16) computes
  ``istop = istart + width`` in the ``'full'`` analysis (pre-fix); the head
  convention ``istop = istart + width - 1`` is used here. The difference is
  one sample on segment stops and is statistically irrelevant on
  chromatographic grids (>10^4 points).

Algorithm (vectorized, O(n)): the dead-zone sign of ``diff(x)`` is computed
in one pass; runs of identical sign are extracted from index gaps, so the
cost is a handful of numpy array operations regardless of the number of
segments — this mirrors the ``find``-based Matlab implementation and its
speed.

The ``'full'`` analysis classifies every segment by four boolean attributes
(monotonic, positive start, increasing, zero-crossing) into at most 16
classes and covers the whole signal (monotone segments are taken with the
``leftpriority`` convention; the remainder forms non-monotonic segments),
so that ``sum(width) == len(x)`` — this invariant is enforced.

Control plots (the Matlab no-output behavior) are provided as explicit
functions :func:`plot_segments` and :func:`plot_classes`.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.signal import filtfilt

__all__ = [
    "filtzero",
    "monotone",
    "monotone_full",
    "sign_runs",
    "MonotoneSegments",
    "MonotoneFull",
    "DEFAULT_LETTER_RULES",
    "to_letters",
    "plot_segments",
    "plot_classes",
]

_EPS = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# filtzero — zero-phase moving average (filtzero.m)
# ---------------------------------------------------------------------------
def filtzero(x: np.ndarray, m: Optional[int] = None) -> np.ndarray:
    """Zero-phase moving-average filtering (port of ``filtzero.m``).

    A boxcar of width ``m`` is applied forward and backward (``filtfilt``),
    yielding exactly zero phase distortion; edge transients are minimized by
    odd reflection over ``3*(m-1)`` samples, as in the Matlab reference.

    Parameters
    ----------
    x : array (n,) or (n, k)
        Signal(s), filtered column-wise.
    m : int, optional
        Filter width in samples; default ``ceil(n/100)``. ``m < 2`` returns
        a copy of ``x`` unchanged.
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    if m is None:
        m = int(np.ceil(n / 100))
    m = int(m)
    if m < 2:
        return x.copy()
    padlen = 3 * (m - 1)
    if n <= padlen:
        raise ValueError(f"filtzero: signal length {n} must exceed 3*(m-1) = {padlen}")
    b = np.ones(m) / m
    return filtfilt(b, [1.0], x, axis=0, padtype="odd", padlen=padlen)


# ---------------------------------------------------------------------------
# monotone — monotonic segments (monotone.m)
# ---------------------------------------------------------------------------
@dataclass
class MonotoneSegments:
    """Monotonic segments of a signal (0-based sample indices).

    A segment spans samples ``start .. stop`` inclusive,
    ``stop = start + width - 1``; ``height = x[stop] - x[start]`` (signed).
    ``sign`` records the segment direction (+1 rising, -1 falling,
    0 constant) — an extension over the Matlab output, useful for plotting.
    """

    start: np.ndarray
    width: np.ndarray
    height: np.ndarray
    sign: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))

    def __post_init__(self) -> None:
        if self.sign.size == 0 and self.start.size:
            self.sign = np.sign(self.height).astype(int)

    @property
    def stop(self) -> np.ndarray:
        return self.start + self.width - 1

    def __len__(self) -> int:
        return int(self.start.size)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "start": self.start,
                "stop": self.stop,
                "width": self.width,
                "height": self.height,
                "sign": self.sign,
            }
        )


def _default_zero(x: np.ndarray) -> float:
    # Matlab default: max(10*eps, (max(X)-min(X))/1e6)
    return max(10 * _EPS, (float(np.max(x)) - float(np.min(x))) / 1e6)


def _runs(idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Start/end sample indices of consecutive runs in a sorted index array
    of ``diff`` positions: ``il`` = run starts, ``ir`` = run ends + 1
    (both as 0-based sample indices into ``x``)."""
    if idx.size == 0:
        e = np.empty(0, dtype=int)
        return e, e
    gap = np.diff(idx) > 1
    il = idx[np.r_[True, gap]]
    ir = idx[np.r_[gap, True]] + 1
    return il, ir


def sign_runs(c: np.ndarray, zero: float = 0.0, flats: str = "segment") -> np.ndarray:
    """Boundaries of constant-trend runs of a series — the shared
    segmentation primitive of the symbolic encoders
    (``signomics.DNAsignal.encode_dna``, ``icfilter.encode_series``).

    The dead-zone sign of ``diff(c)`` is computed (``|diff| <= zero`` counts
    as flat) and the 0-based sample indices bounding each run of constant
    sign are returned as ``bounds``; segment ``k`` spans samples
    ``bounds[k] .. bounds[k+1]`` inclusive, so consecutive segments share
    their boundary sample (unlike :func:`monotone`, whose ``leftpriority``
    convention avoids sharing).

    Parameters
    ----------
    c : array (n,)
        Series (n >= 2 for a meaningful result).
    zero : float
        Dead zone; the encoders use ``0.0`` (exact sign).
    flats : {'segment', 'attach'}
        ``'segment'``: flat runs form segments of their own (``encode_dna``
        behavior). ``'attach'``: flat runs are forward-filled with the
        preceding trend so plateaus attach to it (``encode_series``
        behavior; a leading plateau keeps sign 0).

    Returns
    -------
    bounds : int array
        ``[0, change_1, ..., change_m, n-1]``.
    """
    c = np.asarray(c, dtype=float).ravel()
    n = c.size
    d = np.diff(c)
    s = np.where(d > zero, 1.0, np.where(d < -zero, -1.0, 0.0))
    if flats == "attach":
        nz = s != 0.0
        idx = np.where(nz, np.arange(s.size), -1)
        np.maximum.accumulate(idx, out=idx)
        s = np.where(idx >= 0, s[np.clip(idx, 0, None)], 0.0)
    elif flats != "segment":
        raise ValueError(f"flats must be 'segment' or 'attach', got {flats!r}")
    changes = np.nonzero(s[1:] != s[:-1])[0] + 1
    return np.concatenate(([0], changes, [max(n - 1, 0)]))


def monotone(
    x: np.ndarray,
    kind: str = "+",
    zero: Optional[float] = None,
    leftpriority: bool = False,
    nooverlap: bool = False,
) -> MonotoneSegments:
    """Find monotonic segments in a 1-D signal (port of ``monotone.m``).

    Parameters
    ----------
    x : array (n,)
        Signal.
    kind : {'+', '-', '0', '+-', '-+'}
        Rising, falling, constant, or both rising and falling segments
        (merged and sorted by start; common points belong to '+' by
        construction).
    zero : float, optional
        Dead zone: ``|diff(x)| <= zero`` counts as constant. Default
        ``max(10*eps, (max(x)-min(x))/1e6)`` as in Matlab.
    leftpriority : bool
        Only with ``kind='+-'``: shorten every segment by one sample at its
        end (Matlab keyword ``'leftpriority'``, rev. 2021-04-16), so
        segments never share their boundary point; used by
        :func:`monotone_full` to preserve total coverage.
    nooverlap : bool
        Only with ``kind='+-'``: trim falling segments that share a
        boundary point with a rising one (Matlab keyword ``'nooverlap'``).
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if zero is None:
        zero = _default_zero(x)

    if kind in ("+-", "-+"):
        segp = monotone(x, "+", zero)
        segm = monotone(x, "-", zero)
        pp, lp, ap = segp.start.copy(), segp.width.copy(), segp.height.copy()
        pm, lm, am = segm.start.copy(), segm.width.copy(), segm.height.copy()
        if nooverlap and leftpriority:
            raise ValueError("leftpriority and nooverlap are exclusive")
        if nooverlap:
            # ends of falling segments that are starts of rising ones
            _, _, jover = np.intersect1d(pp, pm + lm - 1, return_indices=True)
            lm[jover] -= 1
            # starts of falling segments that are ends of rising ones
            _, _, jover = np.intersect1d(pp + lp - 1, pm, return_indices=True)
            shift = np.minimum(n - 1, pm[jover] + 1) - pm[jover]
            lm[jover] -= shift
            pm[jover] += shift
            if np.any(lm == 0):
                warnings.warn("monotone: segment of width 0 removed")
                keep = lm != 0
                pm, lm, am = pm[keep], lm[keep], am[keep]
        elif leftpriority:
            lp = lp - 1
            lm = lm - 1
        pos = np.concatenate([pp, pm])
        wid = np.concatenate([lp, lm])
        hei = np.concatenate([ap, am])
        sgn = np.concatenate(
            [np.ones(pp.size, dtype=int), -np.ones(pm.size, dtype=int)]
        )
        isort = np.argsort(pos, kind="stable")
        return MonotoneSegments(pos[isort], wid[isort], hei[isort], sgn[isort])

    d = np.diff(x)
    s = np.where(d > zero, 1, np.where(d < -zero, -1, 0))
    if kind == "+":
        idx = np.nonzero(s > 0)[0]
    elif kind == "-":
        idx = np.nonzero(s < 0)[0]
    elif kind == "0":
        idx = np.nonzero(s == 0)[0]
    else:
        raise ValueError(f"unknown kind={kind!r}")
    il, ir = _runs(idx)
    if il.size == 0:
        e = np.empty(0, dtype=int)
        return MonotoneSegments(e, e, np.empty(0), e)
    sign_val = {"+": 1, "-": -1, "0": 0}[kind]
    return MonotoneSegments(
        il,
        ir - il + 1,
        x[ir] - x[il],
        np.full(il.size, sign_val, dtype=int),
    )


# ---------------------------------------------------------------------------
# monotone 'full' — complete segmentation with class attributes
# ---------------------------------------------------------------------------
@dataclass
class MonotoneFull:
    """Complete segmentation of a signal (``monotone(x, 'full')``).

    Every sample belongs to exactly one segment (``sum(width) == n``);
    segments are classified by four boolean attributes into at most 16
    classes (``attributes`` has one row per class observed, in order of
    first appearance, as Matlab ``unique(..., 'stable')``).
    """

    nclass: int
    fullclass: np.ndarray  # (n,) class index per sample, 0-based
    seg_class: np.ndarray  # (nseg,) class index per segment, 0-based
    start: np.ndarray
    stop: np.ndarray
    width: np.ndarray
    height: np.ndarray
    area: np.ndarray  # signed, base*height/2
    attributes: pd.DataFrame  # ismonotonic, ispositive, isincreasing, iscrossingzero

    def __len__(self) -> int:
        return int(self.start.size)


def monotone_full(x: np.ndarray, zero: Optional[float] = None) -> MonotoneFull:
    """Advanced segmentation (port of ``monotone(X, 'full')``).

    Monotone segments are detected with the ``leftpriority`` convention;
    the uncovered remainder forms non-monotonic segments, so the
    segmentation is a partition of the signal.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 3:
        raise ValueError("monotone_full: provide at least 3 points")
    if zero is None:
        zero = _default_zero(x)

    # truly monotonic segments ('+' or '-'), leftpriority convention
    seg1 = monotone(x, "+-", zero, leftpriority=True)
    istart1, width1, height1 = seg1.start, seg1.width, seg1.height
    istop1 = istart1 + width1 - 1

    # uncovered remainder: cancel covered samples in the ramp t = 1..n
    t = np.arange(1, n + 1, dtype=float)
    for a, b in zip(istart1, istop1):
        t[a : b + 1] = 0.0
    seg2 = monotone(t, "+", 0.0)
    istart2, width2 = seg2.start.copy(), seg2.width.copy()
    iz = t[istart2] == 0.0
    if n >= 2 and t[0] == 1.0 and t[1] == 0.0:  # non-coding left boundary
        istart2 = np.r_[-1, istart2]
        width2 = np.r_[2, width2]
        iz = np.r_[True, iz]
    istart2[iz] = np.minimum(n - 1, istart2[iz] + 1)
    width2[iz] -= 1
    istop2 = istart2 + width2 - 1
    height2 = np.array(
        [x[a : b + 1].max() - x[a : b + 1].min() for a, b in zip(istart2, istop2)]
    )

    # combine and sort
    istart = np.concatenate([istart1, istart2])
    width = np.concatenate([width1, width2])
    height = np.concatenate([height1, height2])
    istop = np.concatenate([istop1, istop2])
    ismono = np.r_[
        np.ones(istart1.size, dtype=bool), np.zeros(istart2.size, dtype=bool)
    ]
    isort = np.argsort(istart, kind="stable")
    istart, width, height, istop, ismono = (
        istart[isort],
        width[isort],
        height[isort],
        istop[isort],
        ismono[isort],
    )

    # attributes and stable-unique classes
    rows = np.column_stack(
        [
            ismono,
            x[istart] > 0,
            x[istop] > x[istart],
            x[istart] * x[istop] < 0,
        ]
    )
    seen: Dict[Tuple[bool, ...], int] = {}
    iclass = np.empty(rows.shape[0], dtype=int)
    uniq = []
    for i, row in enumerate(map(tuple, rows)):
        if row not in seen:
            seen[row] = len(seen)
            uniq.append(row)
        iclass[i] = seen[row]
    attributes = pd.DataFrame(
        np.array(uniq, dtype=bool),
        columns=["ismonotonic", "ispositive", "isincreasing", "iscrossingzero"],
    )

    fullclass = np.repeat(iclass, width)
    if fullclass.size != n:
        raise RuntimeError(
            f"monotone_full: coverage broken (sum(width)={fullclass.size} != n={n})"
        )
    return MonotoneFull(
        nclass=len(uniq),
        fullclass=fullclass,
        seg_class=iclass,
        start=istart,
        stop=istop,
        width=width,
        height=height,
        area=(width * height) / 2.0,
        attributes=attributes,
    )


# ---------------------------------------------------------------------------
# letter coding (documented example of monotone.m)
# ---------------------------------------------------------------------------
#: (ismonotonic, ispositive, isincreasing, iscrossingzero) -> letter
DEFAULT_LETTER_RULES: Dict[Tuple[bool, bool, bool, bool], str] = {
    (True, False, True, True): "A",
    (True, True, True, False): "B",
    (True, False, True, False): "C",
    (True, True, False, True): "Z",
    (True, False, False, False): "Y",
    (True, True, False, False): "X",
}


def to_letters(
    full: MonotoneFull,
    rules: Optional[Dict[Tuple[bool, bool, bool, bool], str]] = None,
    lowercase_below_median: bool = False,
) -> Tuple[str, str]:
    """Translate a full segmentation into letters (monotone.m example).

    Non-monotonic segments map to ``'_'``; monotonic segments map through
    ``rules`` (:data:`DEFAULT_LETTER_RULES` by default). With
    ``lowercase_below_median=True``, monotonic segments whose unsigned area
    is below the median of monotonic areas are lowercased.

    Returns
    -------
    (str_per_segment, str_per_sample)
    """
    if rules is None:
        rules = DEFAULT_LETTER_RULES
    class_letter = []
    for row in map(tuple, full.attributes.to_numpy()):
        if not row[0]:
            class_letter.append("_")
        elif row in rules:
            class_letter.append(rules[row])
        else:
            raise ValueError(f"no letter rule for attributes {row}")
    letters = np.array([class_letter[c] for c in full.seg_class])
    if lowercase_below_median:
        area = np.abs(full.area)
        mono_cls = set(i for i, row in enumerate(full.attributes.to_numpy()) if row[0])
        sel = np.array([c in mono_cls for c in full.seg_class])
        if np.any(sel):
            ref = np.median(area[sel])
            low = sel & (area < ref)
            letters = np.where(low, np.char.lower(letters.astype(str)), letters)
    per_segment = "".join(letters)
    per_sample = "".join(np.repeat(letters, full.width))
    return per_segment, per_sample


# ---------------------------------------------------------------------------
# control plots (Matlab no-output behavior)
# ---------------------------------------------------------------------------
def plot_segments(
    x: np.ndarray,
    segments: Optional[MonotoneSegments] = None,
    kind: str = "+",
    zero: Optional[float] = None,
    ax=None,
):
    """Control plot of :func:`monotone` (Matlab ``monotone(X, ...)`` with no
    output): signal in blue, segment starts as red circles, stops as magenta
    squares, stem of heights as green diamonds."""
    import matplotlib.pyplot as plt

    x = np.asarray(x, dtype=float).ravel()
    if segments is None:
        segments = monotone(x, kind, zero)
    if ax is None:
        _, ax = plt.subplots()
    il, ir = segments.start, segments.stop
    ax.plot(x, "b-", linewidth=0.8)
    ax.plot(il, x[il], "ro", label="start")
    ax.plot(ir, x[ir], "ms", label="stop")
    if il.size:
        ax.stem(il, x[ir] - x[il], linefmt="g-", markerfmt="gd", basefmt=" ")
    ax.set_xlabel("sample index")
    ax.legend(loc="best", fontsize=8)
    return ax


def plot_classes(
    x: np.ndarray,
    full: Optional[MonotoneFull] = None,
    zero: Optional[float] = None,
    ax=None,
    linewidth: float = 2.0,
):
    """Control plot of :func:`monotone_full` (Matlab ``monotone(X, 'full')``
    with no output): each segment class drawn in its own color."""
    import matplotlib.pyplot as plt

    x = np.asarray(x, dtype=float).ravel()
    if full is None:
        full = monotone_full(x, zero)
    if ax is None:
        _, ax = plt.subplots()
    t = np.arange(x.size)
    cmap = plt.get_cmap("tab20" if full.nclass > 10 else "tab10")
    for i in range(full.nclass):
        y = x.copy()
        y[full.fullclass != i] = np.nan
        ax.plot(t, y, linewidth=linewidth, color=cmap(i % cmap.N), label=f"class {i}")
    ax.set_xlabel("sample index")
    ax.legend(loc="best", fontsize=8, ncol=2)
    return ax
