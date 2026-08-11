"""
Baseline reconstruction by peak removal (port of ``removepeaks.m`` /
``ndf.m``).

``removepeaks`` replaces each peak region by a closed-form regularized
interpolant — a cubic (optionally quartic, coefficient ``f``) polynomial
matching the smoothed signal value **and** first derivative at both region
boundaries — thereby rebuilding the baseline under the peaks. When several
left-shifted trial windows are given (``xtrials``), the one minimizing the
bending energy $\\int (d^2y/dx^2)^2$ of the bridge is retained.

``ndf`` is the companion high-order numerical differentiator: 7-point
central stencils of order $O(h^6)$ in the interior with matching 6-point
one-sided stencils at the boundaries (first and second derivative), plus a
sign-consistency repair of early-boundary glitches.

Reference: ``removepeaks.m`` (MS 2.1, 2012, O. Vitrac & X. Fang) and
``ndf.m`` (MS-MATLAB 1.0, 2004, O. Vitrac) — "MS — Molecular Studio"
toolbox; used by the thesis of J. Kermorvant
(https://theses.hal.science/tel-04194172).

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

from typing import Callable, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.interpolate import PchipInterpolator

from .segments import filtzero

__all__ = ["ndf", "removepeaks"]


# ---------------------------------------------------------------------------
# ndf — high-order numerical differentiation (ndf.m)
# ---------------------------------------------------------------------------
def _ndf_1d(
    t: np.ndarray,
    y: np.ndarray,
    order: int,
    dydt0: Optional[float],
    method: str,
) -> np.ndarray:
    n = y.size
    if n < 8:
        raise ValueError(f"ndf: at least 8 points are required, got {n}")
    dt = t[1] - t[0]  # uniform sampling assumed, as in the reference
    if order == 1:
        D = np.array([-1.0, 9.0, -45.0, 0.0, 45.0, -9.0, 1.0]) / (60.0 * dt)
        D0 = np.array([-137.0, 300.0, -300.0, 200.0, -75.0, 12.0]) / (60.0 * dt)
        D1 = -D0[::-1]
    elif order == 2:
        D = np.array([2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0]) / (180.0 * dt**2)
        D0 = np.array([45.0, -154.0, 214.0, -156.0, 61.0, -10.0]) / (12.0 * dt**2)
        D1 = D0[::-1]
    else:
        raise ValueError("ndf: order must be 1 or 2")

    # boundary rows (6-point one-sided) + interior rows (7-point centered)
    head = np.array([np.dot(y[k : k + 6], D0) for k in range(3)])
    tail = np.array([np.dot(y[n - 8 + k : n - 2 + k], D1) for k in range(3)])
    interior = np.convolve(y, D[::-1], mode="valid")  # length n - 6
    dydt = np.concatenate([head, interior, tail])

    if dydt0 is not None:
        dydt[0] = dydt0

    # repair sign-inconsistent glitches near the left boundary
    # (Matlab: i = intersect(find(dy*sign(mean(sign(dy)))<0), 2:10))
    dy = np.diff(dydt)
    dominant = np.sign(np.mean(np.sign(dy)))
    bad = np.nonzero(dy * dominant < 0)[0]
    bad = bad[(bad >= 1) & (bad <= 9)]
    if bad.size:
        keep = np.setdiff1d(np.arange(n), bad)
        if method != "pchip":
            raise NotImplementedError("ndf: only method='pchip' is ported")
        dydt[bad] = PchipInterpolator(t[keep], dydt[keep], extrapolate=True)(t[bad])
    return dydt


def ndf(
    t: np.ndarray,
    y: np.ndarray,
    order: int = 1,
    dydt0: Optional[float] = None,
    method: str = "pchip",
    makeuniform: bool = False,
    resolution: int = 10_000,
) -> np.ndarray:
    """High-order numerical derivative (port of ``ndf.m``).

    First (``order=1``, $O(h^6)$ interior accuracy) or second
    (``order=2``) derivative of ``y`` sampled on the uniform grid ``t``
    (uniformity is assumed, not checked — as in the reference; use
    ``makeuniform=True`` for non-uniform grids to interpolate on a uniform
    grid of ``resolution`` points first). ``y`` may be 1-D or 2-D
    (column-wise derivation).
    """
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float)
    if y.ndim == 2:
        return np.column_stack(
            [
                ndf(t, y[:, j], order, dydt0, method, makeuniform, resolution)
                for j in range(y.shape[1])
            ]
        )
    y = y.ravel()
    if y.size != t.size:
        raise ValueError("ndf: t and y must have the same length")
    if makeuniform:
        t0 = t
        if np.unique(t0).size == 1:
            return np.full_like(y, np.inf)
        tu = np.linspace(t0.min(), t0.max(), int(resolution))
        yu = PchipInterpolator(t0, y)(tu)
        dydt = _ndf_1d(tu, yu, order, dydt0, method)
        return PchipInterpolator(tu, dydt)(t0)
    return _ndf_1d(t, y, order, dydt0, method)


# ---------------------------------------------------------------------------
# removepeaks — closed-form Hermite bridge across peak regions (removepeaks.m)
# ---------------------------------------------------------------------------
def _interpolant(f, x, x1, x2, y1, y2, yy1, yy2):
    """Closed-form cubic/quartic polynomial with prescribed values (y1, y2)
    and first derivatives (yy1, yy2) at x1 and x2 (Matlab symbolic
    derivation reproduced verbatim in removepeaks.m)."""
    return (
        y2
        + x * yy2
        - x2 * yy2
        + f * x**4
        - 2.0 * f * x**3 * x1
        - 2.0 * f * x**3 * x2
        - 2.0 * (x - x2) ** 3 / (x1 - x2) ** 3 * (y1 - y2)
        - (x - x2) ** 2 * (yy1 + 2.0 * yy2) / (x1 - x2)
        + f * x**2 * x1**2
        + f * x**2 * x2**2
        + f * x1**2 * x2**2
        + (x - x2) ** 2
        / (x1 - x2) ** 2
        * (3.0 * y1 - 3.0 * y2 + x * yy1 + x * yy2 - x2 * yy1 - x2 * yy2)
        - 2.0 * f * x * x1 * x2**2
        - 2.0 * f * x * x1**2 * x2
        + 4.0 * f * x**2 * x1 * x2
    )


def _interp_strict(xg: np.ndarray, yg: np.ndarray, xq: float) -> float:
    """Linear interpolation, failing closed outside the data range
    (Matlab interp1 would return NaN there)."""
    if xq < xg[0] or xq > xg[-1]:
        raise ValueError(
            f"removepeaks: anchor x={xq} outside the peak-free data range "
            f"[{xg[0]}, {xg[-1]}]"
        )
    return float(np.interp(xq, xg, yg))


def removepeaks(
    x: np.ndarray,
    y: np.ndarray,
    segments: Union[np.ndarray, Sequence[Sequence[float]]],
    f: float = 0.0,
    smooth: int = 5,
    xtrials: Union[float, Sequence[float]] = 0.0,
    return_models: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, List[List[Callable]]]]:
    """Remove peak regions and rebuild the baseline (port of
    ``removepeaks.m``).

    Parameters
    ----------
    x : array (m,)
        Abscissae.
    y : array (m,) or (m, k)
        Signal(s), processed column-wise.
    segments : array (nseg, 2)
        Peak regions ``[start, stop]`` in x units (e.g. the
        ``start``/``stop`` fields of :func:`~.peaks.monotonepeak`).
    f : float
        Quartic regularization coefficient (0 → cubic bridge).
    smooth : int
        ``filtzero`` bandwidth applied before anchoring (the returned
        signal is the smoothed one outside the bridged regions).
    xtrials : float or sequence
        Left-shift trials of the region start; the bridge minimizing the
        bending energy $\\sum (d^2y/dx^2)^2$ is retained.
    return_models : bool
        Also return the per-segment/per-column bridge polynomials as
        callables.
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float)
    squeeze = y.ndim == 1
    if squeeze:
        y = y[:, None]
    m, ncol = y.shape
    if x.size != m:
        raise ValueError("x and y must have the same number of rows")
    segments = np.atleast_2d(np.asarray(segments, dtype=float))
    if segments.shape == (2, 1):
        segments = segments.T
    if segments.shape[1] != 2:
        raise ValueError("segments must be an (nseg, 2) array")
    xtrials_arr = np.atleast_1d(np.asarray(xtrials, dtype=float))
    ntrials = xtrials_arr.size

    yf = filtzero(y, smooth)
    dyfdx = ndf(x, yf, 1)

    ynp = yf.copy()
    good = np.ones((m, ntrials), dtype=bool)
    models: List[List[Callable]] = []
    for x1, x2 in segments:
        bad0 = (x >= x1) & (x <= x2)
        bad = np.zeros((m, ntrials), dtype=bool)
        for k in range(ntrials):
            bad[:, k] = (x >= x1 + xtrials_arr[k]) & (x <= x2)
            good[:, k] &= ~bad[:, k]
        seg_models: List[Callable] = []
        for j in range(ncol):
            trials = []
            crit = np.zeros(ntrials)
            anchors = []
            for k in range(ntrials):
                g = good[:, k]
                y1 = _interp_strict(x[g], yf[g, j], x1)
                y2 = _interp_strict(x[g], yf[g, j], x2)
                yy1 = _interp_strict(x[g], dyfdx[g, j], x1)
                yy2 = _interp_strict(x[g], dyfdx[g, j], x2)
                bridge = _interpolant(f, x[bad0], x1, x2, y1, y2, yy1, yy2)
                trials.append(bridge)
                anchors.append((y1, y2, yy1, yy2))
                crit[k] = float(np.sum(ndf(x[bad0], bridge, 2) ** 2))
            k = int(np.argmin(crit))
            ynp[bad0, j] = trials[k]
            if return_models:
                y1, y2, yy1, yy2 = anchors[k]
                seg_models.append(
                    lambda xq, x1=x1, x2=x2, y1=y1, y2=y2, yy1=yy1, yy2=yy2: (
                        _interpolant(
                            f, np.asarray(xq, dtype=float), x1, x2, y1, y2, yy1, yy2
                        )
                    )
                )
        if return_models:
            models.append(seg_models)
    out = ynp[:, 0] if squeeze else ynp
    if return_models:
        return out, models
    return out
