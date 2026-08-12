"""
Peak-area reconstruction from the Ricker CWT (thesis §5.5.3.3, Figs
5.46–5.48 of https://theses.hal.science/tel-04194172).

Closed form. For a Gaussian peak ``g(t) = h exp(-t^2/(2 sigma^2))`` and the
L2-normalized Ricker wavelet at scale ``l`` (the exact engine of
:mod:`..icfilter`; pywt's ``mexh`` used by ``signomics.compute_cwt`` matches
it to <0.3 %), the transform is again Ricker-shaped with effective scale
``s = sqrt(sigma^2 + l^2)``::

    W(t; h, sigma, l) = K h sigma l^{5/2} / s^3 * (1 - t^2/s^2) exp(-t^2/(2 s^2))
    K = 2 sqrt(2 pi) / (sqrt(3) pi^{1/4}) = 2.17418...

Consequences (all validated numerically against both engines):

- optimal scale (max apex response):  ``l_opt = sqrt(5) sigma = 2.2361 sigma``
- at ``l_opt``: zeros at ``+/- s = sqrt(6) sigma = 2.4495 sigma``,
  minima at ``+/- sqrt(3) s = 4.2426 sigma``
- amplification ``W(0; l_opt)/h = K 5^{5/4}/6^{3/2} sqrt(sigma)
  = 1.10597 sqrt(sigma)`` (sigma in samples)
- unsigned central-lobe area over apex amplitude is scale-pure:
  ``A_cwt / W(0) = 2 s / sqrt(e)``

**Revalidation note (2026-08-12).** Thesis Property 5.1 states
``l_opt = 2.58 sigma``, zeros ``+/-2.77 sigma``, minima ``+/-4.79 sigma``,
amplification ``1.588 sqrt(sigma)``; Fig 3.4 reads ``l_opt = 2.48 x width``.
Direct measurement with both current engines reproduces the closed form
above and **not** those constants (which are mutually consistent but match
neither engine; they presumably derive from a different normalization
convention of the thesis-era tooling). The *structure* of the method —
existence of ``l_opt`` proportional to ``sigma`` and the bijective
calibration of Fig 5.48 with ``l_opt`` acting as the Lagrange multiplier —
is confirmed; the numerical constants are superseded for the current
implementation.

Reconstruction (the Fig 5.48 bijection): the apex of each detected peak is
followed across scales (simplified ridge filiation after Du et al. 2006);
``l_opt`` is interpolated log-parabolically from the apex amplitudes; then::

    sigma = l_opt / sqrt(5)
    h     = W_max / (1.10597 sqrt(sigma))
    A_g   = h sigma sqrt(2 pi)

Compatible inputs: ``signomics.DNAsignal.cwt_coeffs`` (dict scale -> 1-D
coefficient array, pywt path) and stacks built from
``icfilter.cwt_matrix`` (exact path). Scales are in samples; pass ``dx``
(or ``x``) to obtain areas in physical x units.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

__all__ = [
    "K_RICKER",
    "LOPT_FACTOR",
    "AMPLIFICATION_FACTOR",
    "ricker_gaussian_apex",
    "optimal_scale",
    "sigma_from_scale",
    "scale_ridge",
    "interpolate_optimal_scale",
    "reconstruct_areas",
    "calibration_curve",
    "plot_calibration",
]

#: closed-form constant K = 2 sqrt(2 pi)/(sqrt(3) pi^{1/4})
K_RICKER = 2.0 * np.sqrt(2.0 * np.pi) / (np.sqrt(3.0) * np.pi**0.25)
#: l_opt / sigma (supersedes thesis 2.58 / Fig 3.4's 2.48 — see module note)
LOPT_FACTOR = float(np.sqrt(5.0))
#: W(0; l_opt) / (h sqrt(sigma)) (supersedes thesis 1.588 — see module note)
AMPLIFICATION_FACTOR = float(K_RICKER * 5.0**1.25 / 6.0**1.5)
#: A_cwt / W(0) = LOBE_FACTOR * s with s = sqrt(sigma^2 + l^2)
LOBE_FACTOR = 2.0 / np.sqrt(np.e)


def ricker_gaussian_apex(sigma: float, scale: float, height: float = 1.0) -> float:
    """Apex amplitude ``W(0)`` of the transform of a Gaussian (closed form)."""
    s2 = sigma**2 + scale**2
    return float(K_RICKER * height * sigma * scale**2.5 / s2**1.5)


def optimal_scale(sigma: float) -> float:
    """``l_opt = sqrt(5) sigma`` (sigma in samples)."""
    return LOPT_FACTOR * sigma


def sigma_from_scale(lopt: float) -> float:
    """Inverse of :func:`optimal_scale`."""
    return lopt / LOPT_FACTOR


def _as_stack(
    coeffs: Union[Mapping[float, np.ndarray], Tuple[Sequence[float], np.ndarray]],
) -> Tuple[np.ndarray, np.ndarray]:
    """Normalize input to (scales ascending, matrix (nscales, n))."""
    if isinstance(coeffs, Mapping):
        scales = np.array(sorted(coeffs), dtype=float)
        mat = np.vstack([np.asarray(coeffs[s], dtype=float).ravel() for s in scales])
    else:
        scales, mat = coeffs
        scales = np.asarray(scales, dtype=float)
        mat = np.asarray(mat, dtype=float)
        idx = np.argsort(scales)
        scales, mat = scales[idx], mat[idx]
    if mat.ndim != 2 or mat.shape[0] != scales.size:
        raise ValueError("coeffs must map each scale to one coefficient row")
    return scales, mat


def scale_ridge(
    coeffs: Union[Mapping[float, np.ndarray], Tuple[Sequence[float], np.ndarray]],
    icenter: int,
    search: float = 1.0,
) -> pd.DataFrame:
    """Follow one peak apex across scales (simplified filiation after
    Du et al. 2006: at each scale the local maximum nearest to the position
    tracked at the previous scale, within ``+/- max(3, search*scale)``
    samples).

    Returns a DataFrame with columns ``scale``, ``ipos``, ``amplitude``.
    """
    scales, mat = _as_stack(coeffs)
    n = mat.shape[1]
    pos = int(icenter)
    rows = []
    for k, sc in enumerate(scales):
        w = max(3, int(round(search * sc)))
        lo, hi = max(0, pos - w), min(n, pos + w + 1)
        j = lo + int(np.argmax(mat[k, lo:hi]))
        pos = j
        rows.append((sc, j, mat[k, j]))
    return pd.DataFrame(rows, columns=["scale", "ipos", "amplitude"])


def interpolate_optimal_scale(
    scales: np.ndarray, amplitudes: np.ndarray
) -> Tuple[float, float, bool]:
    """Log-parabolic interpolation of the apex-amplitude maximum over
    scales (the thesis interpolates between the dyadic transforms).

    Returns ``(l_opt, w_max, censored)``; ``censored`` is True when the
    maximum sits on the scale-range boundary (no interpolation possible —
    widen the scale range for a reliable estimate).
    """
    scales = np.asarray(scales, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    k = int(np.argmax(amplitudes))
    if k == 0 or k == scales.size - 1:
        return float(scales[k]), float(amplitudes[k]), True
    u = np.log(scales[k - 1 : k + 2])
    a = amplitudes[k - 1 : k + 2]
    d1, d2 = u[1] - u[0], u[2] - u[1]
    s1, s2 = (a[1] - a[0]) / d1, (a[2] - a[1]) / d2
    ca = (s2 - s1) / (d1 + d2)  # second divided difference (curvature)
    if ca >= 0:  # degenerate curvature — keep the grid maximum
        return float(scales[k]), float(a[1]), True
    # Newton form f(u) = a0 + s1 (u-u0) + ca (u-u0)(u-u1); f'(u*) = 0
    ustar = 0.5 * (u[0] + u[1]) - s1 / (2.0 * ca)
    wmax = a[0] + s1 * (ustar - u[0]) + ca * (ustar - u[0]) * (ustar - u[1])
    return float(np.exp(ustar)), float(wmax), False


def _fit_ridge_model(
    scales: np.ndarray, amps: np.ndarray, valid: np.ndarray
) -> Tuple[float, float, float, bool]:
    """Fit the closed-form scale response ``W(l) = C l^{5/2}/(sigma^2+l^2)^{3/2}``
    to the valid ridge amplitudes (amplitude given ``sigma`` is linear in
    ``C``, so ``C`` is projected out and only ``sigma`` is optimized).

    Returns ``(sigma, C, wmax_model, censored)``. ``censored`` is True when
    the response maximum is not bracketed by the valid scales (fit
    extrapolates — widen the scale range or check for interference).
    """
    sc_all = np.asarray(scales, dtype=float)
    am_all = np.asarray(amps, dtype=float)
    sc, am = sc_all[valid], am_all[valid]
    if sc.size < 3:
        sc, am = sc_all, am_all  # degraded fallback: use everything
    k = int(np.argmax(am))
    # contiguous block around the max with amplitude >= 35 % of it
    strong = am >= 0.35 * am[k]
    lo = k
    while lo > 0 and strong[lo - 1]:
        lo -= 1
    hi = k
    while hi < am.size - 1 and strong[hi + 1]:
        hi += 1
    scf, amf = sc[lo : hi + 1], am[lo : hi + 1]
    censored = k == 0 or k == am.size - 1 or scf.size < 3

    def sse(sigma: float) -> float:
        g = scf**2.5 / (sigma**2 + scf**2) ** 1.5
        c = float(amf @ g) / float(g @ g)
        return float(np.sum((amf - c * g) ** 2))

    lo_b = sc[max(k - 1, 0)] / (2.0 * LOPT_FACTOR)
    hi_b = sc[min(k + 1, sc.size - 1)] * 2.0 / LOPT_FACTOR
    res = minimize_scalar(sse, bounds=(lo_b, hi_b), method="bounded")
    sigma = float(res.x)
    g = scf**2.5 / (sigma**2 + scf**2) ** 1.5
    C = float(amf @ g) / float(g @ g)
    lopt = optimal_scale(sigma)
    wmax = C * lopt**2.5 / (sigma**2 + lopt**2) ** 1.5
    return sigma, C, wmax, censored


def reconstruct_areas(
    coeffs: Union[Mapping[float, np.ndarray], Tuple[Sequence[float], np.ndarray]],
    icenter: Union[Sequence[int], np.ndarray],
    dx: Optional[float] = None,
    x: Optional[np.ndarray] = None,
    search: float = 1.0,
) -> pd.DataFrame:
    """Reconstruct Gaussian peak areas from a multi-scale Ricker CWT
    (the Fig 5.48 bijection, ``l_opt`` as Lagrange multiplier).

    Parameters
    ----------
    coeffs : dict {scale: 1-D array} or (scales, matrix)
        Multi-scale transform of one signal — ``DNAsignal.cwt_coeffs`` or a
        stack of ``icfilter.cwt_matrix`` rows. Scales in samples, spanning
        the expected ``l_opt`` range (``~2.24 sigma``) on both sides.
    icenter : sequence of int
        Apex sample indices of the peaks to quantify (e.g.
        ``monotonepeak(...).icenter``).
    dx : float, optional
        Sample step; taken from ``x`` when given. When available, the
        physical-unit columns ``sigma_x`` and ``area_x`` are added.
    search : float
        Ridge search-window factor (see :func:`scale_ridge`).

    Method
    ------
    For each peak the apex ridge is extracted, then rows whose apex drifted
    more than ``max(3, 0.25*scale)`` samples from the detected position are
    discarded (the Ricker transform preserves the position of an isolated
    Gaussian, so drift marks interference with a neighbor — coelution — or
    baseline takeover at coarse scales), and the closed-form scale response
    is fitted; the area follows directly from the fitted prefactor,
    ``A_g = C sqrt(2 pi)/K``.

    Returns
    -------
    DataFrame with per-peak columns: ``icenter``, ``lopt`` (samples),
    ``wmax`` (modeled transform apex amplitude), ``sigma`` (samples),
    ``height``, ``area`` (sample units), ``censored`` (True when the
    response maximum was not bracketed by usable scales — widen the scale
    range, or the peak is coeluted), and ``sigma_x``/``area_x`` when ``dx``
    is known.
    """
    if x is not None and dx is None:
        x = np.asarray(x, dtype=float).ravel()
        dx = float(np.median(np.diff(x)))
    rows = []
    for ic in np.atleast_1d(np.asarray(icenter, dtype=int)):
        ridge = scale_ridge(coeffs, int(ic), search=search)
        sc = ridge["scale"].to_numpy()
        am = ridge["amplitude"].to_numpy()
        drift = np.abs(ridge["ipos"].to_numpy() - int(ic))
        valid = drift <= np.maximum(3.0, 0.25 * sc)
        if not valid.all():  # truncate at first contaminated scale
            first_bad = int(np.argmin(valid))
            valid[first_bad:] = False
        sigma, C, wmax, censored = _fit_ridge_model(sc, am, valid)
        height = C / (K_RICKER * sigma)
        area = C * np.sqrt(2.0 * np.pi) / K_RICKER
        rows.append(
            (int(ic), optimal_scale(sigma), wmax, sigma, height, area, censored)
        )
    out = pd.DataFrame(
        rows,
        columns=["icenter", "lopt", "wmax", "sigma", "height", "area", "censored"],
    )
    if dx is not None:
        out["sigma_x"] = out["sigma"] * dx
        out["area_x"] = out["area"] * dx
    return out


def calibration_curve(
    sigmas: Optional[np.ndarray] = None,
    scales: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Regenerate the Fig 5.48 model analytically.

    For each ``(sigma, l)``: the normalized unsigned central-lobe area of
    the transform ``A_cwt/W(0) = 2 s/sqrt(e)`` versus the normalized
    Gaussian area ``A_g/W(0) = sqrt(2 pi) s^3/(K l^{5/2})`` (both per unit
    transform amplitude, as in the thesis figure). Rows with
    ``optimal=True`` form the bijective locus ``l = sqrt(5) sigma``.

    Returns a DataFrame with columns ``sigma``, ``scale``, ``s``,
    ``Ag_over_h``, ``Acwt_over_h``, ``optimal``.
    """
    if sigmas is None:
        sigmas = np.geomspace(2.0, 64.0, 25)
    if scales is None:
        scales = 2.0 ** np.arange(0, 9)
    rows = []
    for sig in np.asarray(sigmas, dtype=float):
        for sc in np.asarray(scales, dtype=float):
            s = np.hypot(sig, sc)
            w0 = ricker_gaussian_apex(sig, sc)
            rows.append(
                (
                    sig,
                    sc,
                    s,
                    sig * np.sqrt(2.0 * np.pi) / w0,
                    LOBE_FACTOR * s,
                    False,
                )
            )
        lo = optimal_scale(sig)
        s = np.hypot(sig, lo)
        w0 = ricker_gaussian_apex(sig, lo)
        rows.append(
            (sig, lo, s, sig * np.sqrt(2.0 * np.pi) / w0, LOBE_FACTOR * s, True)
        )
    return pd.DataFrame(
        rows, columns=["sigma", "scale", "s", "Ag_over_h", "Acwt_over_h", "optimal"]
    )


def plot_calibration(df: Optional[pd.DataFrame] = None, ax=None):
    """Control plot reproducing Fig 5.48: per-scale isocontours (dashed)
    and the bijective optimal-scale locus (solid)."""
    import matplotlib.pyplot as plt

    if df is None:
        df = calibration_curve()
    if ax is None:
        _, ax = plt.subplots()
    for sc, grp in df[~df["optimal"]].groupby("scale"):
        g = grp.sort_values("Ag_over_h")
        ax.plot(
            g["Ag_over_h"],
            g["Acwt_over_h"],
            "--",
            linewidth=0.8,
            label=f"$\\ell$={sc:g}",
        )
    loc = df[df["optimal"]].sort_values("Ag_over_h")
    ax.plot(
        loc["Ag_over_h"],
        loc["Acwt_over_h"],
        "r-",
        linewidth=2,
        label="optimal locus $\\ell_{opt}=\\sqrt{5}\\,\\sigma$",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$A_g / h_{cwt}$")
    ax.set_ylabel(r"$A_{cwt} / h_{cwt}$")
    ax.legend(loc="best", fontsize=7, ncol=2)
    return ax
