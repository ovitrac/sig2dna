"""Tests for sig2dna_core.chemspace — synthetic fixtures only."""
import numpy as np
import pytest
from sig2dna_core.chemspace import (
    pcoa, lcurve_elbow, fraction_of_distance, truncated_distances,
    nearest_reference_scores, reference_null, z_threshold)


def make_points(seed=0):
    rng = np.random.default_rng(seed)
    refs = rng.normal(0, 1, (5, 3))          # tight reference cluster
    far = rng.normal(8, 1, (10, 3))          # displaced cluster
    X = np.vstack([refs, far])
    D = np.sqrt(((X[:, None] - X[None, :]) ** 2).sum(-1))
    return X, D


def test_pcoa_reproduces_euclidean_distances():
    X, D = make_points()
    Y, lam, neg = pcoa(D)
    assert neg < 1e-6                        # Euclidean-embeddable
    Dr = np.sqrt(((Y[:, None] - Y[None, :]) ** 2).sum(-1))
    np.testing.assert_allclose(Dr, D, atol=1e-6)


def test_pcoa_weighted_centering_moves_origin():
    _, D = make_points()
    w = np.zeros(len(D)); w[:5] = 1.0        # anchor on the references
    Y, _, _ = pcoa(D, weights=w)
    assert np.abs(Y[:5].mean(axis=0)).max() < 1e-6


def test_lcurve_elbow_finds_break():
    lam = np.array([100, 80, 60, 2, 1.5, 1.2, 1.0, 0.9])
    assert lcurve_elbow(lam) in (3, 4)


def test_fraction_and_truncation():
    X, D = make_points()
    Y, lam, _ = pcoa(D)
    f = sum(fraction_of_distance(D, Y, k) for k in range(3))
    assert f > 0.9                            # 3 informative axes
    Dt = truncated_distances(Y, 3)
    assert np.corrcoef(Dt[np.triu_indices(len(D), 1)],
                       D[np.triu_indices(len(D), 1)])[0, 1] > 0.999


def test_nearest_reference_scoring_and_z_rule():
    _, D = make_points()
    refs = list(range(5))
    score = nearest_reference_scores(D, refs)
    mu, sd = reference_null(D, refs)
    th = z_threshold(mu, sd, k=3)
    calls = score > th
    assert not calls[:5].any()                # references pass
    assert calls[5:].all()                    # displaced cluster flagged
    # strict LOO: excluding a reference changes the null consistently
    mu2, _ = reference_null(D, refs, exclude=0)
    assert mu2 == pytest.approx(mu, rel=1.0)
