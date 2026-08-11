"""
Peak detection from monotone segments (port of ``monotonepeak.m``).

A peak apex is defined structurally, without derivative thresholds: it is a
sample that terminates a rising monotone segment **and** starts a falling
one on the ``filtzero``-smoothed signal. Detection at several smoothing
bandwidths (``mfilt``) is a discrete scale-space analysis: peaks found at
the coarsest bandwidth are kept, and finer bandwidths contribute only peaks
whose apex is farther than ``itol`` samples from every peak already found.

Reference: ``monotonepeak.m`` ("MS — Molecular Studio" toolbox, O. Vitrac,
rev. 18/01/2018), used by the thesis of J. Kermorvant
(https://theses.hal.science/tel-04194172). All indices are 0-based
(Matlab is 1-based); ``center``/``start``/``stop`` are x-values and
``width`` is in x units, as in the reference.

The control plot (Matlab no-output behavior) is :func:`plot_peaks`:
smoothed segments overlaid in color, numbered apexes, and a stem display of
the sorting criterion.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from .segments import filtzero, monotone

__all__ = ["PeakTable", "monotonepeak", "plot_peaks"]

_MFILT_RATIO_DEFAULT = 1.0 / 100.0


@dataclass
class PeakTable:
    """Detected peaks (one entry per peak, arrays aligned).

    Index fields (0-based samples): ``icenter``, ``istart``, ``istop``,
    ``ibase`` (the lower of the two peak feet, used for baseline anchoring).
    X-valued fields: ``center``, ``start``, ``stop``, ``width``
    (= left-rise + right-fall extent in x units). ``height`` is the larger
    monotone-segment amplitude, ``absheight`` the raw signal value at the
    apex. ``tail``/``wall`` give the asymmetry side by width/height;
    ``ratiowidth``/``ratioheight`` quantify it. ``level`` is the 1-based
    rank of the ``mfilt`` bandwidth that found the peak.
    """

    level: np.ndarray
    icenter: np.ndarray
    istart: np.ndarray
    istop: np.ndarray
    center: np.ndarray
    start: np.ndarray
    stop: np.ndarray
    height: np.ndarray
    absheight: np.ndarray
    width: np.ndarray
    ratioheight: np.ndarray
    ratiowidth: np.ndarray
    ibase: np.ndarray
    tail: np.ndarray
    wall: np.ndarray

    def __len__(self) -> int:
        return int(self.icenter.size)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame({f.name: getattr(self, f.name) for f in fields(self)})

    def take(self, idx: np.ndarray) -> "PeakTable":
        return PeakTable(**{f.name: getattr(self, f.name)[idx] for f in fields(self)})

    @classmethod
    def empty(cls) -> "PeakTable":
        e = np.empty(0, dtype=int)
        f_ = np.empty(0, dtype=float)
        s = np.empty(0, dtype=object)
        return cls(
            e.copy(),
            e.copy(),
            e.copy(),
            e.copy(),
            f_.copy(),
            f_.copy(),
            f_.copy(),
            f_.copy(),
            f_.copy(),
            f_.copy(),
            f_.copy(),
            f_.copy(),
            e.copy(),
            s.copy(),
            s.copy(),
        )

    @classmethod
    def concatenate(cls, tables: Sequence["PeakTable"]) -> "PeakTable":
        return cls(
            **{
                f.name: np.concatenate([getattr(t, f.name) for t in tables])
                for f in fields(cls)
            }
        )


def _detect_single(
    y: np.ndarray,
    x: np.ndarray,
    mfilt: int,
    zero: Optional[float],
    sort: Optional[str],
    sortby: str,
) -> PeakTable:
    """Single-bandwidth detection (the nfilt==1 branch of monotonepeak.m)."""
    yf = filtzero(y, mfilt) if mfilt >= 2 else y
    segp = monotone(yf, "+", zero)
    segm = monotone(yf, "-", zero)
    pp, lp, ap = segp.start, segp.width, segp.height
    pm, lm, am = segm.start, segm.width, -segm.height  # am > 0
    # apex = end of a rising segment that starts a falling one
    p, ip, im = np.intersect1d(pp + lp - 1, pm, return_indices=True)
    if p.size == 0:
        return PeakTable.empty()
    lpx = x[pp + lp - 1] - x[pp]  # rise extents (x units)
    lmx = x[pm + lm - 1] - x[pm]  # fall extents (x units)
    sides = np.array(["left", "right"], dtype=object)
    istop = pm[im] + lm[im] - 1
    peaks = PeakTable(
        level=np.ones(p.size, dtype=int),
        icenter=p,
        istart=pp[ip],
        istop=istop,
        center=x[p],
        start=x[pp[ip]],
        stop=x[istop],
        height=np.maximum(ap[ip], am[im]),
        absheight=y[p],
        width=lpx[ip] + lmx[im],
        ratioheight=np.maximum(ap[ip], am[im]) / np.minimum(ap[ip], am[im]),
        ratiowidth=np.maximum(lpx[ip], lmx[im]) / np.minimum(lpx[ip], lmx[im]),
        ibase=np.where(y[istop] < y[pp[ip]], istop, pp[ip]),
        tail=sides[(lpx[ip] < lmx[im]).astype(int)],
        wall=sides[(ap[ip] < am[im]).astype(int)],
    )
    return _sorted(peaks, sort, sortby)


def _sorted(peaks: PeakTable, sort: Optional[str], sortby: str) -> PeakTable:
    if not sort or len(peaks) == 0:
        return peaks
    key = np.asarray(getattr(peaks, sortby), dtype=float)
    if sort == "descend":
        idx = np.argsort(-key, kind="stable")
    elif sort == "ascend":
        idx = np.argsort(key, kind="stable")
    else:
        raise ValueError(f"sort must be 'ascend' or 'descend', got {sort!r}")
    return peaks.take(idx)


def monotonepeak(
    y: np.ndarray,
    x: Optional[np.ndarray] = None,
    mfilt: Union[None, int, Sequence[int]] = None,
    zero: Optional[float] = None,
    xtol: Optional[float] = None,
    itol: int = 10,
    sort: Optional[str] = None,
    sortby: str = "absheight",
    mfiltsort: str = "descend",
    maxpeak: Optional[int] = None,
    keeporder: bool = False,
) -> PeakTable:
    """Find peaks by monotone-segment analysis (port of ``monotonepeak.m``).

    Parameters
    ----------
    y : array (n,)
        Signal.
    x : array (n,), optional
        Abscissae (default: sample indices ``0..n-1``).
    mfilt : int or sequence of int, optional
        ``filtzero`` bandwidth(s); default ``round(n/100)``. With several
        values (sorted per ``mfiltsort``, coarsest first by default), finer
        bandwidths only contribute peaks whose apex differs by more than
        ``itol`` samples from every peak already found; ``level`` records
        the 1-based bandwidth rank.
    zero : float, optional
        Dead zone forwarded to :func:`~.segments.monotone`.
    xtol : float, optional
        Tolerance in x units for peak identity; overrides ``itol``
        (``itol = max(round(xtol/mean(diff(x))), 1)``).
    itol : int
        Tolerance in samples for peak identity across bandwidths.
    sort : {'ascend', 'descend', None}
        Sort peaks by ``sortby`` (default: natural apex order).
    sortby : str
        Sorting field, default ``'absheight'``.
    mfiltsort : {'descend', 'ascend'}
        Processing order of the ``mfilt`` values.
    maxpeak : int, optional
        Keep at most this number of peaks (after sorting).
    keeporder : bool
        Re-sort the final list by increasing apex position.
    """
    y = np.asarray(y, dtype=float).ravel()
    n = y.size
    if n == 0:
        raise ValueError("y is empty")
    x = np.arange(n, dtype=float) if x is None else np.asarray(x, dtype=float).ravel()
    if x.size != n:
        raise ValueError("x and y must be of the same size")
    if mfilt is None:
        mfilt = _MFILT_RATIO_DEFAULT * n
    mfilt_arr = np.atleast_1d(np.asarray(mfilt, dtype=float))
    mfilt_arr = np.maximum(0, np.round(mfilt_arr)).astype(int)
    if xtol is not None:
        itol = max(int(round(xtol / float(np.mean(np.diff(x))))), 1)

    if mfilt_arr.size == 1:
        peaks = _detect_single(y, x, int(mfilt_arr[0]), zero, sort, sortby)
    else:
        if mfiltsort == "descend":
            mfilt_arr = np.sort(mfilt_arr)[::-1]
        else:
            mfilt_arr = np.sort(mfilt_arr)
        # per-bandwidth candidates are merged in natural (ascending-apex)
        # order — as the Matlab reference does — so that the first of two
        # nearby candidates wins the itol dedup; the requested sort is
        # applied to the final list only
        peaks = _detect_single(y, x, int(mfilt_arr[0]), zero, None, sortby)
        for lev in range(2, mfilt_arr.size + 1):
            ptmp = _detect_single(y, x, int(mfilt_arr[lev - 1]), zero, None, sortby)
            for j in range(len(ptmp)):
                if (
                    len(peaks) == 0
                    or np.min(np.abs(ptmp.icenter[j] - peaks.icenter)) > itol
                ):
                    new = ptmp.take(np.array([j]))
                    new.level[:] = lev
                    peaks = PeakTable.concatenate([peaks, new])
        peaks = _sorted(peaks, sort, sortby)

    if maxpeak is not None and len(peaks) > maxpeak:
        peaks = peaks.take(np.arange(maxpeak))
    if keeporder and len(peaks):
        peaks = peaks.take(np.argsort(peaks.center, kind="stable"))
    return peaks


def plot_peaks(
    y: np.ndarray,
    x: Optional[np.ndarray] = None,
    peaks: Optional[PeakTable] = None,
    sortby: str = "absheight",
    ax=None,
    **kwargs,
):
    """Control plot of :func:`monotonepeak` (Matlab no-output behavior):
    signal as a thin black line, each peak's span overlaid in color with a
    numbered apex, and a stem display of the sorting criterion from the
    bottom of the axes."""
    import matplotlib.pyplot as plt

    y = np.asarray(y, dtype=float).ravel()
    x = (
        np.arange(y.size, dtype=float)
        if x is None
        else np.asarray(x, dtype=float).ravel()
    )
    if peaks is None:
        peaks = monotonepeak(y, x, sort="descend", sortby=sortby, **kwargs)
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y, "k-", linewidth=0.5)
    npk = len(peaks)
    if npk == 0:
        return ax
    cmap = plt.get_cmap("Set2")
    y0 = float(np.min(y))
    crit = np.asarray(getattr(peaks, sortby), dtype=float)
    for i in range(npk):
        col = cmap(i % cmap.N)
        sl = slice(peaks.istart[i], peaks.istop[i] + 1)
        ax.plot(x[sl], y[sl], "-", linewidth=2, color=col)
        ax.annotate(
            f"{i + 1}",
            (peaks.center[i], y[peaks.icenter[i]]),
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )
        ax.plot([peaks.center[i]] * 2, [y0, y0 + crit[i]], "-", color=col, linewidth=2)
        ax.plot(peaks.center[i], y0 + crit[i], "o", markersize=5, color=col)
        ax.plot(
            [peaks.start[i], peaks.center[i], peaks.stop[i]],
            [y0, y0 + crit[i], y0],
            "-",
            color=col,
            linewidth=1,
        )
    ax.set_xlabel("x")
    return ax
