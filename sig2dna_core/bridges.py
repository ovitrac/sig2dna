"""
Bridge metrology: when may two acquisition domains be compared?

A chemical fingerprint need not be invariant across analytical
operators; what must be *measurable* is the information lost, preserved
or transformed between them. A **bridge** is the same set of physical
materials measured under two acquisition conditions. Its edge record
carries two coordinates:

* ``R`` — transport cost: median same-material cross-condition distance
  divided by a replicate **floor** (the distance between repeat
  measurements of one material within a domain). Symmetric bridges use
  the conservative ``max`` of both floors; if one side has no replicates
  the ratio is explicitly *one-sided*.
* ``Q`` — structure preserved, in levels: ``Q0`` none, ``Q1``
  group/stream structure, ``Q2`` material neighborhood (significant rank
  enrichment), ``Q3`` material identity resolved.

Cost alone cannot characterize an edge (a moderate ``R`` can carry no
structure at all); both coordinates are always reported. Comparisons
between domains are legitimate only across measured edges.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np


def replicate_floor(replicate_pairs, dist) -> float:
    """Median distance over ``(x, y)`` same-material repeat pairs."""
    return float(np.median([dist(x, y) for x, y in replicate_pairs]))


def measure_bridge(side_a: dict, side_b: dict, dist, floor: float,
                   groups: dict | None = None, one_sided: bool = False,
                   n_perm: int = 10000, seed: int = 0):
    """Edge record for a same-material bridge.

    ``side_a``/``side_b`` map material -> representation (same keys);
    ``dist`` is the domain-valid distance; ``floor`` per the module
    doctrine; ``groups`` optionally maps material -> group label for the
    stream endpoint. Returns the two-coordinate record with permutation
    p-values (never a bare pass/fail).
    """
    mats = sorted(side_a)
    n = len(mats)
    D = np.array([[dist(side_a[i], side_b[j]) for j in mats] for i in mats])
    diag = np.diag(D)
    ranks = [int(1 + np.sum(D[i] < D[i, i])) for i in range(n)]
    top1 = float(np.mean([r == 1 for r in ranks]))
    mrr = float(np.mean([1.0 / r for r in ranks]))
    rng = np.random.default_rng(seed)
    null = [np.mean([1.0 / (1 + np.sum(D[i] < D[i, p[i]]))
                     for i in range(n)])
            for p in (rng.permutation(n) for _ in range(n_perm))]
    p_mat = float(np.mean(np.array(null) >= mrr))
    rec = dict(materials=n, floor=float(floor), one_sided=bool(one_sided),
               median_true_pair=float(np.median(diag)),
               R=float(np.median(diag) / floor), ranks=ranks, top1=top1,
               mrr=mrr, p_material=p_mat)
    if groups:
        g = [groups[m] for m in mats]
        within = [D[i, j] for i in range(n) for j in range(n)
                  if g[i] == g[j]]
        cross = [D[i, j] for i in range(n) for j in range(n)
                 if g[i] != g[j]]
        T = float(np.mean(cross) - np.mean(within))
        labels = sorted(set(g))
        if len(labels) == 2:
            k = g.count(labels[0])
            perms = []
            for c in combinations(range(n), k):
                gg = [labels[0] if i in c else labels[1] for i in range(n)]
                w = [D[i, j] for i in range(n) for j in range(n)
                     if gg[i] == gg[j]]
                x = [D[i, j] for i in range(n) for j in range(n)
                     if gg[i] != gg[j]]
                perms.append(np.mean(x) - np.mean(w))
            p_stream = float(np.mean(np.array(perms) >= T))
        else:
            p_stream = float("nan")
        rec.update(T_stream=T, p_stream=p_stream)
    rec["Q"] = ("Q3" if p_mat < 0.05 and top1 >= 0.5 else
                "Q2" if p_mat < 0.05 else
                "Q1" if groups and rec.get("p_stream", 1.0) < 0.05 else
                "Q0")
    rec["edge_class"] = {"Q3": "material identity resolved",
                         "Q2": "partial material / neighborhood transport",
                         "Q1": "stream-level transport",
                         "Q0": "operator-dominates"}[rec["Q"]]
    return rec
