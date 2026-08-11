"""
Chemical-space analysis on symbolic distances — generic tools extracted from
the D3 validation campaign (2026-08): PCoA (weighted Gower), L-curve
dimension selection, truncated-space distances, nearest-reference scoring
with leave-one-out nulls, and z-threshold rules.

These functions operate on any symmetric distance matrix (e.g., summed
per-IC Levenshtein or exclusive-entropy distances from
``sig2dna_core.icfilter``) and carry no dataset-specific logic.

Method notes
------------
- ``pcoa`` is Principal Coordinate Analysis (Torgerson/Gower), NOT PCA: it
  consumes distances (no feature matrix), double-centers the squared
  distances, B = -1/2 * J D^2 J, and eigendecomposes B exactly. Optional
  weighted centering places the origin at a weighted centroid (e.g., a set
  of reference samples) - Gower (1966) generalized form.
- ``lcurve_elbow`` selects the number of significant dimensions as the point
  of maximum perpendicular distance to the chord of the eigenvalue curve.
  Distances computed in the truncated space shed the extensive noise offset
  that accumulates over many ion channels and residual dimensions.
- ``nearest_reference_scores`` implements the nearest-known-reference
  criterion with strict leave-one-out for reference samples, and
  ``z_threshold`` the k-sigma rule on the references' own nearest-neighbour
  null (k = 2..4 by the inverse normal law).

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

__all__ = [
    "pcoa",
    "lcurve_elbow",
    "fraction_of_distance",
    "truncated_distances",
    "nearest_reference_scores",
    "reference_null",
    "z_threshold",
]


def pcoa(D: np.ndarray, weights: Optional[np.ndarray] = None):
    """Principal Coordinate Analysis of a symmetric distance matrix.

    Parameters
    ----------
    D : (n, n) symmetric distance matrix, zero diagonal.
    weights : optional (n,) nonnegative centering weights (sum normalized
        internally). Uniform weights give classical Gower centering; putting
        the mass on a reference subset anchors the origin on it.

    Returns
    -------
    X : (n, p) coordinates for the p positive eigenvalues (columns ordered
        by decreasing eigenvalue).
    eigenvalues : (p,) positive eigenvalues.
    negative_inertia : |sum of negative eigenvalues| / sum of positive ones
        (0 for Euclidean-embeddable distances).
    """
    D = np.asarray(D, dtype=np.float64)
    n = len(D)
    if D.shape != (n, n):
        raise ValueError("D must be square")
    if not np.allclose(D, D.T, atol=1e-8):
        raise ValueError("D must be symmetric")
    if weights is None:
        w = np.full(n, 1.0 / n)
    else:
        w = np.asarray(weights, dtype=np.float64)
        if w.min() < 0 or w.sum() <= 0:
            raise ValueError("weights must be nonnegative with positive sum")
        w = w / w.sum()
    Jl = np.eye(n) - np.outer(np.ones(n), w)
    B = -0.5 * Jl @ (D ** 2) @ Jl.T
    lam, U = np.linalg.eigh((B + B.T) / 2)
    order = np.argsort(lam)[::-1]
    lam, U = lam[order], U[:, order]
    pos = lam > 1e-9 * max(lam.max(), 1.0)
    neg = float(-lam[lam < 0].sum() / lam[pos].sum()) if pos.any() else 0.0
    return U[:, pos] * np.sqrt(lam[pos]), lam[pos], neg


def lcurve_elbow(eigenvalues: Sequence[float]) -> int:
    """Number of significant dimensions by the L-curve criterion: the index
    of maximum perpendicular distance to the chord joining the first and
    last points of the (normalized) eigenvalue curve."""
    lam = np.asarray(eigenvalues, dtype=np.float64)
    if len(lam) < 3:
        return len(lam)
    x = np.arange(1, len(lam) + 1, dtype=np.float64)
    y = lam / lam[0]
    num = np.abs((y[-1] - y[0]) * x - (x[-1] - x[0]) * y
                 + x[-1] * y[0] - y[-1] * x[0])
    den = float(np.hypot(y[-1] - y[0], x[-1] - x[0]))
    return int(np.argmax(num / den) + 1)


def fraction_of_distance(D: np.ndarray, X: np.ndarray, axis: int) -> float:
    """Share of the total pairwise distance carried by one coordinate axis:
    sum of pairwise |coordinate gaps| on that axis / sum of pairwise D."""
    iu = np.triu_indices(len(D), 1)
    dk = np.abs(X[:, axis][:, None] - X[:, axis][None, :])[iu]
    return float(dk.sum() / np.asarray(D)[iu].sum())


def truncated_distances(X: np.ndarray, m: int) -> np.ndarray:
    """Euclidean pairwise distances using the first ``m`` coordinates only
    (the denoised chemical space)."""
    Xm = np.asarray(X)[:, :m]
    return np.sqrt(((Xm[:, None, :] - Xm[None, :, :]) ** 2).sum(-1))


def nearest_reference_scores(
    dist: np.ndarray, ref_idx: Sequence[int]
) -> np.ndarray:
    """Distance of every sample to its nearest reference, with strict
    leave-one-out: a sample that is itself a reference is scored against
    the other references only."""
    dist = np.asarray(dist)
    ref = list(ref_idx)
    out = np.empty(len(dist))
    for i in range(len(dist)):
        refs = [r for r in ref if r != i]
        if not refs:
            raise ValueError("at least two references are required")
        out[i] = dist[i, refs].min()
    return out


def reference_null(dist: np.ndarray, ref_idx: Sequence[int],
                   exclude: Optional[int] = None):
    """Mean and SD (Bessel) of the references' leave-one-out
    nearest-neighbour distances; ``exclude`` removes one reference from the
    null entirely (strict LOO when evaluating that reference)."""
    dist = np.asarray(dist)
    ref = [r for r in ref_idx if r != exclude]
    vals = []
    for r in ref:
        others = [u for u in ref if u != r]
        vals.append(dist[r, others].min())
    return float(np.mean(vals)), float(np.std(vals, ddof=1))


def z_threshold(mu: float, sd: float, k: float = 3.0) -> float:
    """k-sigma decision threshold on the reference null (k = 2..4 per the
    inverse normal law: ~98 %, ~99.9 %, ~99.997 % one-sided)."""
    return mu + k * sd
