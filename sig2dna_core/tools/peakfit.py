"""
Deconvolution of detected peaks as a Gaussian/Lorentzian mixture (port of
``monotonepeakfit.m``).

Two fits are always performed and both returned:

- **strategy 1** — peak centers prescribed at the detected apexes
  (:func:`~.peaks.monotonepeak`), only widths are optimized;
- **strategy 2** — centers and widths free, softly constrained by
  tanh-penalty terms (shift/width buffers and Lagrange-like scales); after
  the fit, the peak↔guess assignment is resolved by scoring all
  permutations (n ≤ 7) against ordering/shift/width violations, the
  identity permutation first as tie-breaker (Matlab ``fliplr(perms(...))``).

Kernels (the reference's simplified parameterization)::

    gaussian(x; c, w)   = exp(-((x-c)/(0.6006 w))^2)
        area = 0.6006*sqrt(pi)*w,  sigma = 0.6006/sqrt(2)*w
    lorentzian(x; c, w) = 1/(1+((x-c)/(0.5 w))^2)   (same half-height value)

Extension over the Matlab reference (``emg=True``): an exponentially
modified Gaussian kernel for chromatographic tailing — same 0.6006 Gaussian
core, convolved with an exponential of **signed** time constant ``tau``
(positive tails right, negative mirrors to fronting), fitted per peak in
both strategies and initialized from the detected asymmetry side
(``PeakTable.tail``). EMG kernels are normalized to unit maximum so weights
remain peak amplitudes; their area is reported as ``1/norm`` (the EMG
density has unit area).

Peak amplitudes (weights) are resolved at every objective evaluation by
linear least squares ``|G\\y|`` (or NNLS when ``nonneg`` and G is
rank-deficient). ``scipy``'s Nelder–Mead uses the same initial-simplex
construction and coefficients as Matlab ``fminsearch``, so both strategies
track the reference closely.

Ported verbatim, including two documented reference quirks: (i) ``r2max``
is computed from the objective value, penalties included; (ii) in the
sorted-output block, the window-center denominator sums the weights up to
the *first strategy's* significant-peak count for both strategies (Matlab
uses a vector as a colon bound there and silently takes its first element,
with a warning since R2024a).

Reference: ``monotonepeakfit.m`` (MS 2.1, 2013, O. Vitrac,
rev. 07/07/2018), "MS — Molecular Studio" toolbox; used by the thesis of
J. Kermorvant (https://theses.hal.science/tel-04194172).

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy.optimize import minimize, nnls
from scipy.special import erfcx
from scipy.stats import norm as _norm

from .peaks import PeakTable

__all__ = ["FitResult", "monotonepeakfit", "plot_fit"]

_MAXPEAKS_FOR_CONSTRAINTS = 7
_ALPHA = 0.05


def _gaussian_kernel(x, position, width):
    return np.exp(-(((x - position) / (0.6006 * width)) ** 2))


def _lorentzian_kernel(x, position, width):
    return 1.0 / (1.0 + ((x - position) / (0.5 * width)) ** 2)


def _emg_pdf(x, position, width, tau):
    """Exponentially modified Gaussian density (unnormalized amplitude).

    Gaussian core follows the 0.6006-width convention
    (``sigma = 0.6006*width/sqrt(2)``); ``tau`` is the exponential tail
    time constant, **signed**: ``tau > 0`` tails to the right (classical GC
    tailing), ``tau < 0`` mirrors the shape to tail left (fronting).

    Evaluated in the numerically stable scaled form
    ``pdf = exp(-xi^2/(2 sigma^2)) * erfcx(v) / (2|tau|)`` with
    ``v = (sigma/|tau| - xi/sigma)/sqrt(2)``, which degrades gracefully to
    the Gaussian shape as ``tau -> 0`` (a relative floor on ``|tau|``
    avoids the singular limit).
    """
    sigma = 0.6006 * width / np.sqrt(2.0)
    t = np.maximum(np.abs(tau), 1e-8 * width)
    xi = np.where(tau >= 0, x - position, position - x)
    v = (sigma / t - xi / sigma) / np.sqrt(2.0)
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        core = np.exp(-(xi**2) / (2.0 * sigma**2)) * erfcx(v)
        # far decay side (v << 0): erfc(v) ~= 2 exactly; the exponent below
        # is negative wherever this branch is selected, so no overflow
        far = 2.0 * np.exp(sigma**2 / (2.0 * t**2) - xi / t)
        return np.where(v > -25.0, core, far) / (2.0 * t)


def _argpad(v, n: int) -> np.ndarray:
    """Pad/truncate to n elements, repeating the last one (Matlab argpad)."""
    v = np.atleast_1d(np.asarray(v, dtype=float)).ravel()
    if v.size >= n:
        return v[:n].copy()
    return np.r_[v, np.full(n - v.size, v[-1])]


def _nearestpointunique(val: np.ndarray, lst: np.ndarray) -> np.ndarray:
    """Unique nearest-point assignment, closest pairs served first
    (port of nearestpointunique.m, mode 'nearest')."""
    val = np.asarray(val, dtype=float)
    lst = np.asarray(lst, dtype=float)
    if val.size > lst.size:
        raise ValueError("val must not have more elements than list")
    crit = np.abs(val - lst[np.argmin(np.abs(val[:, None] - lst[None, :]), axis=1)])
    order = np.argsort(crit, kind="stable")
    ind = np.zeros(val.size, dtype=int)
    ok = np.ones(lst.size, dtype=bool)
    for i in order:
        jvalid = np.nonzero(ok)[0]
        ind[i] = jvalid[np.argmin(np.abs(lst[jvalid] - val[i]))]
        ok[ind[i]] = False
    return ind


def _matlab_perms(n: int):
    """Rows of Matlab fliplr(perms(1:n)): identity first."""
    return [
        tuple(reversed(p)) for p in reversed(list(itertools.permutations(range(n))))
    ]


@dataclass
class FitResult:
    """Result of :func:`monotonepeakfit` (both strategies).

    Row 0 of the (2, n) arrays is strategy 1 (prescribed centers), row 1 is
    strategy 2 (free centers). ``weight``/``relativeweight``/``rank``/``r2``
    are (n, 2), peak-wise per strategy, as in the reference. When the
    result is sorted (``sort=True``), peaks are ordered by decreasing
    weight and ``relativeweight`` holds the cumulative relative weight.
    """

    position: np.ndarray
    width: np.ndarray
    sigma: np.ndarray
    area: np.ndarray
    weight: np.ndarray
    relativeweight: np.ndarray
    rank: np.ndarray
    r2: np.ndarray
    r2max: np.ndarray
    numrank: np.ndarray
    baseline: Tuple[float, float]
    isbaseline: Optional[np.ndarray]
    window: dict
    nsignificantpeaks: np.ndarray
    lorentzian: bool = False
    tau: Optional[np.ndarray] = None  # (2, n), EMG only (signed tail constant)
    norm: Optional[np.ndarray] = None  # (2, n), EMG unit-max normalization
    _x: np.ndarray = field(default_factory=lambda: np.empty(0), repr=False)

    @property
    def npeaks(self) -> int:
        return int(self.position.shape[1])

    def kernel(self, i: int, strategy: int = 2) -> Callable:
        """Unit-amplitude kernel of peak ``i`` (0-based) for a strategy (1|2)."""
        k = strategy - 1
        c, w = self.position[k, i], self.width[k, i]
        if self.tau is not None:
            t, nrm = self.tau[k, i], self.norm[k, i]
            return lambda x: _emg_pdf(np.asarray(x, dtype=float), c, w, t) / nrm
        ker = _lorentzian_kernel if self.lorentzian else _gaussian_kernel
        return lambda x: ker(np.asarray(x, dtype=float), c, w)

    def expansion(
        self,
        x,
        npeaks: Optional[int] = None,
        strategy: int = 2,
        with_baseline: bool = True,
    ) -> np.ndarray:
        """Mixture model truncated to the first ``npeaks`` peaks
        (Matlab ``outsum``/``outsumbaseline``)."""
        x = np.asarray(x, dtype=float)
        k = strategy - 1
        n = self.npeaks if npeaks is None else int(npeaks)
        if self.tau is not None:
            G = (
                _emg_pdf(
                    x[:, None],
                    self.position[k, :n][None, :],
                    self.width[k, :n][None, :],
                    self.tau[k, :n][None, :],
                )
                / self.norm[k, :n][None, :]
            )
        else:
            ker = _lorentzian_kernel if self.lorentzian else _gaussian_kernel
            G = ker(
                x[:, None], self.position[k, :n][None, :], self.width[k, :n][None, :]
            )
        y = G @ self.weight[:n, k]
        if with_baseline:
            y = y + np.polyval(self.baseline, x)
        return y

    def to_frame(self, strategy: int = 2) -> pd.DataFrame:
        k = strategy - 1
        cols = {
            "position": self.position[k],
            "width": self.width[k],
            "sigma": self.sigma[k],
            "area": self.area[k],
            "weight": self.weight[:, k],
            "relativeweight": self.relativeweight[:, k],
            "rank": self.rank[:, k],
            "r2": self.r2[:, k],
        }
        if self.tau is not None:
            cols["tau"] = self.tau[k]
        return pd.DataFrame(cols)


def monotonepeakfit(
    p: PeakTable,
    y: np.ndarray,
    x: Optional[np.ndarray] = None,
    significant: float = 0.95,
    minpointsinbaseline: int = 5,
    shiftmax: Union[None, float, Sequence[float]] = None,
    shiftbuffer: Union[None, float, Sequence[float]] = None,
    shiftpenaltyscale: Union[None, float, Sequence[float]] = None,
    widthmax: Union[float, Sequence[float]] = np.inf,
    widthmin: Union[float, Sequence[float]] = -np.inf,
    widthbuffer: Union[None, float, Sequence[float]] = None,
    widthpenaltyscale: Union[None, float, Sequence[float]] = None,
    preject: float = 5.0,
    nonneg: bool = False,
    baseline: bool = False,
    sort: bool = False,
    lorentzian: bool = False,
    emg: bool = False,
    endforced: bool = False,
    keeporder: bool = False,
    keepinitialorder: bool = False,
    keepfittedorder: bool = False,
    independent: bool = False,
    maxiter: int = 1000,
    tol: float = 1e-6,
) -> FitResult:
    """Fit ``(x, y)`` as a sum of Gaussian/Lorentzian peaks whose apexes
    were identified by :func:`~.peaks.monotonepeak` (port of
    ``monotonepeakfit.m``). See the module docstring for the two-strategy
    scheme and the constraint machinery.
    """
    n = len(p)
    if n == 0:
        raise ValueError("p contains no peak")
    y = np.asarray(y, dtype=float).ravel()
    m = y.size
    if m == 0:
        raise ValueError("y is empty")
    x = np.arange(m, dtype=float) if x is None else np.asarray(x, dtype=float).ravel()
    if x.size != m:
        raise ValueError("x and y must be of the same size")

    sce = float(np.sum((y - np.mean(y)) ** 2))  # var(y)*(m-1)
    penaltyscale = np.sqrt(sce) / n
    xwidth = float(np.max(x) - np.min(x))
    dx = abs(float(np.median(np.diff(x))) / 4.0)

    def H(v):  # smoothed heaviside forcing positive widths
        return 0.5 * (1.0 + np.tanh(v / dx))

    if emg and lorentzian:
        raise ValueError("emg and lorentzian are exclusive")
    kernel = _lorentzian_kernel if lorentzian else _gaussian_kernel

    # constraint defaults (strategy 2 only)
    shiftmax_v = np.abs(_argpad(np.min(p.width) if shiftmax is None else shiftmax, n))
    shiftbuffer_v = _argpad(dx if shiftbuffer is None else shiftbuffer, n)
    shiftpenalty_v = _argpad(
        penaltyscale if shiftpenaltyscale is None else shiftpenaltyscale, n
    )
    widthmax_v = _argpad(widthmax, n)
    widthmin_v = _argpad(widthmin, n)
    widthbuffer_v = _argpad(shiftbuffer_v if widthbuffer is None else widthbuffer, n)
    widthpenalty_v = _argpad(
        shiftpenalty_v if widthpenaltyscale is None else widthpenaltyscale, n
    )

    def Hshift(v):
        return 0.5 * (1.0 + np.tanh(v / shiftbuffer_v))

    def Hwidth(v):
        return 0.5 * (1.0 + np.tanh(v / widthbuffer_v))

    # baseline removal
    isbaseline: Optional[np.ndarray] = None
    if baseline:
        isbaseline = np.ones(m, dtype=bool)
        for i in range(n):
            isbaseline[(x >= p.start[i]) & (x <= p.stop[i])] = False
        if np.count_nonzero(isbaseline) < minpointsinbaseline:
            isbaseline[[0, -1]] = True
            isbaseline[p.ibase] = True
        if endforced:
            rej = np.percentile(y, preject)
            above = np.nonzero(y > rej)[0]
            if above.size:
                isbaseline[: above[0] + 1] = True
                isbaseline[above[-1] :] = True
        b = np.polyfit(x[isbaseline], y[isbaseline], 1)
        y = y - np.polyval(b, x)
    else:
        b = np.array([0.0, 0.0])

    widthguess = np.asarray(p.width, dtype=float) / 4.0
    positionguess = np.asarray(p.center, dtype=float)
    # EMG tail guess: sign from the detected asymmetry side, |tau| ~ width/4
    tauguess = (
        np.where(np.asarray(p.tail, dtype=object) == "right", 1.0, -1.0) * widthguess
        if emg
        else None
    )

    def mixture(c, s, t=None):
        """Mixture evaluation: kernels + least-squares amplitudes.
        Returns (yfit, weights, G, norms); norms is the per-column unit-max
        normalization (ones unless EMG)."""
        c = np.atleast_1d(c)[None, :]
        s = np.atleast_1d(s)[None, :]
        if emg:
            G = _emg_pdf(x[:, None], c, s, np.atleast_1d(t)[None, :])
            norms = G.max(axis=0)
            norms[norms <= 0] = 1.0
            G = G / norms[None, :]
        else:
            G = kernel(x[:, None], c, s)
            norms = np.ones(n)
        if nonneg and np.linalg.matrix_rank(G) < n:
            W, _ = nnls(G, y)
        else:
            W = np.abs(np.linalg.lstsq(G, y, rcond=None)[0])
        return G @ W, W, G, norms

    def fit_fixedpos(z):
        if emg:
            s, t = z.reshape(2, n, order="F")
        else:
            s, t = z, None
        yfit, _, _, _ = mixture(positionguess, s, t)
        return np.linalg.norm(yfit - y) + np.linalg.norm((1.0 - H(s)) * penaltyscale)

    def fit_free(z):
        ps = z.reshape(3 if emg else 2, n, order="F")
        t = ps[2] if emg else None
        WC = float(
            np.sum(
                (Hwidth(ps[1] - widthmax_v) + Hwidth(widthmin_v - ps[1]))
                * widthpenalty_v
            )
        )
        SC = float(
            np.sum(Hshift(np.abs(ps[0] - positionguess) - shiftmax_v) * shiftpenalty_v)
        )
        yfit, _, _, _ = mixture(ps[0], ps[1], t)
        return (
            np.linalg.norm(yfit - y)
            + np.linalg.norm((1.0 - H(ps[1])) * penaltyscale)
            + WC
            + SC
        )

    nvar1 = 2 * n if emg else n
    nm_opts = dict(xatol=tol, fatol=tol, maxiter=maxiter, maxfev=200 * max(nvar1, 1))
    position = np.full((2, n), np.nan)
    width = np.full((2, n), np.nan)
    tau = np.full((2, n), np.nan) if emg else None
    critfit = np.zeros(2)

    # STRATEGY 1: widths (and EMG taus) only, prescribed centers
    z1 = np.vstack([widthguess, tauguess]).reshape(-1, order="F") if emg else widthguess
    res1 = minimize(fit_fixedpos, z1, method="Nelder-Mead", options=nm_opts)
    if emg:
        width[0], tau[0] = res1.x.reshape(2, n, order="F")
    else:
        width[0] = res1.x
    position[0] = positionguess
    critfit[0] = res1.fun

    # STRATEGY 2: centers and widths (and EMG taus) free
    rows2 = [positionguess, widthguess if independent else width[0]]
    if emg:
        rows2.append(tauguess if independent else tau[0])
    z0 = np.vstack(rows2).reshape(-1, order="F")
    nm_opts2 = dict(nm_opts, maxfev=200 * len(z0))
    res2 = minimize(fit_free, z0, method="Nelder-Mead", options=nm_opts2)
    ps2 = res2.x.reshape(3 if emg else 2, n, order="F")
    position[1], width[1] = ps2[0], ps2[1]
    if emg:
        tau[1] = ps2[2]
    critfit[1] = res2.fun

    # permutation matching of strategy 2 against the guesses
    if n <= _MAXPEAKS_FOR_CONSTRAINTS:
        perms = _matlab_perms(n)
        nperm = len(perms)
        violations = np.zeros((nperm, 7))
        p0order = np.argsort(positionguess, kind="stable")
        porder = np.argsort(position[1], kind="stable")
        widthmax_fitted = float(np.max(width[1]))
        for i, pi in enumerate(perms):
            pi = np.asarray(pi)
            shift = np.abs(position[1, pi] - positionguess)
            violations[i, 0] = np.sum(np.abs(porder[pi] - p0order))
            violations[i, 1] = np.sum(shift)
            violations[i, 2] = np.sum(np.abs(width[0] - width[1, pi]))
            violations[i, 3] = np.sum(shift > shiftmax_v) / n
            violations[i, 4] = np.sum(width[1, pi] > widthmax_v) / n
            violations[i, 5] = np.sum(width[1, pi] > 3.0 * 4.0 * widthguess) / n
            violations[i, 6] = np.sum(width[1, pi] < widthguess) / n
        with np.errstate(invalid="ignore", divide="ignore"):
            maxv = violations.max(axis=0)
            violations[:, :3] = violations[:, :3] / maxv[:3]
        if widthmax_fitted > xwidth / 2:
            vw = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0])
        elif widthmax_fitted > xwidth / 4:
            vw = np.array([0.5, 0.5, 1.0, 0.5, 1.0, 1.0, 1.0])
        else:
            vw = np.array([1.5, 1.0, 1.0, 1.0, 0.3, 0.3, 0.3])
        score = violations @ vw
        if np.all(np.isnan(score)):
            ibest = 0
        else:
            ibest = int(np.nanargmin(score))
        best = np.asarray(perms[ibest])
        position[1] = position[1, best]
        width[1] = width[1, best]
        if emg:
            tau[1] = tau[1, best]

    # weights, ranks, fitted models
    weight = np.full((n, 2), np.nan)
    numrank = np.zeros(2, dtype=int)
    order = np.zeros((n, 2), dtype=int)
    rank_ = np.zeros((n, 2), dtype=int)
    norm_ = np.ones((2, n)) if emg else None
    for k in range(2):
        _, W, G, norms = mixture(position[k], width[k], tau[k] if emg else None)
        weight[:, k] = W
        numrank[k] = np.linalg.matrix_rank(G)
        order[:, k] = np.argsort(-weight[:, k], kind="stable")
        rank_[order[:, k], k] = np.arange(1, n + 1)
        if emg:
            norm_[k] = norms

    sigma = (0.6006 * width) / np.sqrt(2.0)
    # unit-max kernel area: EMG density has unit area, so area = 1/norm
    area = 1.0 / norm_ if emg else 0.6006 * np.sqrt(np.pi) * width
    relativeweight = weight / np.sum(weight, axis=0, keepdims=True)
    r2max = 1.0 - critfit**2 / sce
    windowwidth = float(x[-1] - x[0])
    windowcenter = np.array(
        [float(position[k] @ weight[:, k] / np.sum(weight[:, k])) for k in range(2)]
    )
    r2 = np.full((n, 2), np.nan)
    nsignificantpeaks = np.zeros(2, dtype=int)

    window = {
        "center": windowcenter.copy(),
        "width": windowwidth,
        "alpha": _ALPHA,
        "widthalpha": np.array(
            [
                float(
                    np.diff(
                        _norm.ppf(
                            [_ALPHA / 2, 1 - _ALPHA / 2],
                            loc=windowcenter[k],
                            scale=(0.6006 * windowwidth) / np.sqrt(2.0),
                        )
                    )[0]
                )
                for k in range(2)
            ]
        ),
        "nsignificantpeaks": nsignificantpeaks,
        "significantproba": significant,
    }

    if sort:
        initialorder = np.zeros((2, n), dtype=int)
        for k in range(2):
            initialorder[k, np.argsort(position[k], kind="stable")] = np.arange(n)
        cumweight = np.full((n, 2), np.nan)

        def _permute(k, idx):
            position[k] = position[k, idx]
            width[k] = width[k, idx]
            sigma[k] = sigma[k, idx]
            area[k] = area[k, idx]
            weight[:, k] = weight[idx, k]
            rank_[:, k] = rank_[idx, k]
            if emg:
                tau[k] = tau[k, idx]
                norm_[k] = norm_[k, idx]

        for k in range(2):
            idx = order[:, k]
            _permute(k, idx)
            cumweight[:, k] = np.cumsum(relativeweight[idx, k])
            if not np.all(np.isnan(cumweight[:, k])):
                nsignificantpeaks[k] = (
                    int(np.nonzero(cumweight[:, k] >= significant)[0][0]) + 1
                )
            # reference quirk: the denominator uses the FIRST strategy's
            # significant count for both strategies (vector colon bound)
            ns_num = nsignificantpeaks[k]
            ns_den = nsignificantpeaks[0] if nsignificantpeaks[0] else ns_num
            windowcenter[k] = float(
                position[k, :ns_num] @ weight[:ns_num, k] / np.sum(weight[:ns_den, k])
            )
        relativeweight = cumweight

        if keepinitialorder:
            for k in range(2):
                idx = _nearestpointunique(positionguess, position[k])
                _permute(k, idx)
                relativeweight[:, k] = relativeweight[idx, k]
        elif keepfittedorder:
            for k in range(2):
                srt = np.argsort(position[k], kind="stable")
                idx = srt[initialorder[k]]
                _permute(k, idx)
                relativeweight[:, k] = relativeweight[idx, k]
        elif keeporder:
            idx = np.argsort(position[0], kind="stable")
            for k in range(2):
                _permute(k, idx)
            relativeweight = relativeweight[idx, :]

        # refresh r2 with expansion order i (verbatim, baseline included)
        for k in range(2):
            if emg:
                Gk = (
                    _emg_pdf(
                        x[:, None],
                        position[k][None, :],
                        width[k][None, :],
                        tau[k][None, :],
                    )
                    / norm_[k][None, :]
                )
            else:
                Gk = kernel(x[:, None], position[k][None, :], width[k][None, :])
            for i in range(n):
                model = Gk[:, : i + 1] @ weight[: i + 1, k] + np.polyval(b, x)
                r2[i, k] = 1.0 - float(np.linalg.norm(y - model) ** 2) / sce
            window["widthalpha"][k] = float(
                np.diff(
                    _norm.ppf(
                        [_ALPHA / 2, 1 - _ALPHA / 2],
                        loc=windowcenter[k],
                        scale=(0.6006 * windowwidth) / np.sqrt(2.0),
                    )
                )[0]
            )
        window["center"] = windowcenter.copy()

    return FitResult(
        position=position,
        width=width,
        sigma=sigma,
        area=area,
        weight=weight,
        relativeweight=relativeweight,
        rank=rank_,
        r2=r2,
        r2max=r2max,
        numrank=numrank,
        baseline=(float(b[0]), float(b[1])),
        isbaseline=isbaseline,
        window=window,
        nsignificantpeaks=nsignificantpeaks,
        lorentzian=lorentzian,
        tau=tau,
        norm=norm_,
        _x=x,
    )


def plot_fit(
    y: np.ndarray,
    x: np.ndarray,
    fit: FitResult,
    strategy: int = 2,
    ax=None,
):
    """Control plot: measured signal, mixture model, and individual peak
    contributions (the manual plots of the monotonepeakfit.m examples)."""
    import matplotlib.pyplot as plt

    y = np.asarray(y, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y, "k-", linewidth=0.8, label="measured")
    ax.plot(
        x,
        fit.expansion(x, strategy=strategy),
        "-",
        linewidth=1.5,
        label=f"model (strategy {strategy})",
    )
    k = strategy - 1
    base = np.polyval(fit.baseline, x)
    for i in range(fit.npeaks):
        ax.plot(
            x,
            fit.weight[i, k] * fit.kernel(i, strategy)(x) + base,
            "--",
            linewidth=0.8,
        )
    ax.legend(loc="best", fontsize=8)
    ax.set_xlabel("x")
    return ax
