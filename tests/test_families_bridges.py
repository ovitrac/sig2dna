"""Synthetic tests for compound families and bridge metrology."""
import numpy as np
import pytest

from sig2dna_core.families import (cosine, dedup, group_families,
                                   make_decoys, recognize)
from sig2dna_core.bridges import measure_bridge, replicate_floor
from sig2dna_core.channels import d_comp


def _events(triples):
    return [(ion, apex, 0, 0, a, 1.0) for ion, apex, a in triples]


def test_grouping_and_dedup():
    ev = _events([(10, 1000, 5.0), (20, 1004, 3.0), (30, 1008, 1.0),
                  (10, 2000, 4.0), (20, 2003, 2.0),   # only 2 ions
                  (40, 3000, 9.0), (41, 3001, 8.0), (42, 3002, 7.0)])
    fams = group_families(ev, gap=12, min_ions=3)
    assert len(fams) == 2                     # the 2-ion group is dropped
    twice = dedup(fams + [dict(f) for f in fams])
    assert len(twice) == len(fams)            # identity rule merges clones


def test_recognition_beats_decoys():
    rng = np.random.default_rng(7)
    ev = _events([(m, 1000 + m % 5, 10.0 - 0.01 * m) for m in range(50, 60)])
    fams = group_families(ev, gap=12, min_ions=3)
    assert recognize(fams[0], ev)             # self-recognition
    decoys = make_decoys(fams, n_time=10000, rng=rng)
    assert decoys and all(not recognize(d, ev) for d in decoys)
    assert cosine(fams[0]["spec"], fams[0]["spec"]) == pytest.approx(1.0)


def test_bridge_perfect_and_shuffled():
    rng = np.random.default_rng(2)
    base = {m: rng.poisson(30.0, 200).astype(float) + 1 for m in range(8)}
    noisy = {m: v + rng.poisson(2.0, 200) for m, v in base.items()}
    reps = [(base[m] + rng.poisson(2.0, 200), base[m] + rng.poisson(2.0, 200))
            for m in base]
    floor = replicate_floor(reps, d_comp)
    rec = measure_bridge(base, noisy, d_comp, floor,
                         groups={m: m % 2 for m in base}, n_perm=2000)
    assert rec["top1"] == 1.0 and rec["p_material"] < 0.05
    assert rec["Q"] == "Q3" and rec["R"] < 3.0
    # shuffled sides -> no material structure
    keys = list(base)
    perm = dict(zip(keys, [noisy[k] for k in np.roll(keys, 3)]))
    rec2 = measure_bridge(base, perm, d_comp, floor, n_perm=2000)
    assert rec2["top1"] < 0.5
    assert rec2["Q"] in ("Q0", "Q1", "Q2")
    assert "edge_class" in rec2 and "one_sided" in rec2
