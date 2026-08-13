"""Exact laws of the event/channel information algebra (synthetic only)."""
import numpy as np
import pytest

from sig2dna_core.events import (cluster_classes, d_event, escore,
                                 extract_events, is_recycled_pet,
                                 validate_step0)
from sig2dna_core.channels import (channel_counts, channel_state, d_comp,
                                   mad_logratio)


def _gauss_run(apexes, amp=2000.0, n=8000, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    y = rng.poisson(30.0, n).astype(float)
    for a in apexes:
        y += amp * np.exp(-((t - a) / 10.0) ** 2)
    ev, _, _ = extract_events(y[None, :])
    return ev


def test_step0_one_gaussian_one_event():
    assert validate_step0(verbose=False)


def test_laplace_gain_exact_value():
    # one class present only in the target, absent from n references
    tgt = [(0, 4000, 3900, 4100, 100.0, 1.0)]
    pops = {"X": tgt}
    n_ref = 7
    for r in range(n_ref):
        pops[f"R{r}"] = [(0, 6000, 5900, 6100, 100.0, 1.0)]
    occ = cluster_classes(pops, tau=36)
    eg, el = escore(occ, "X", [f"R{r}" for r in range(n_ref)])
    # gained class: k=0 -> -log2(1/(n+2)); lost class: k=n -> -log2(1-q)
    assert eg == pytest.approx(-np.log2(1.0 / (n_ref + 2)), abs=1e-12)
    assert el == pytest.approx(-np.log2(1.0 - (n_ref + 1) / (n_ref + 2)),
                               abs=1e-12)


def test_cancellation_of_shared_population():
    # adding one disjoint shared class to target AND references changes nothing
    base = {"X": [(0, 4000, 0, 0, 9.0, 1.0)],
            "R0": [(0, 6000, 0, 0, 9.0, 1.0)],
            "R1": [(0, 6000, 0, 0, 9.0, 1.0)]}
    occ0 = cluster_classes(base, tau=36)
    e0 = escore(occ0, "X", ["R0", "R1"])
    shared = (5, 2000, 0, 0, 9.0, 1.0)
    aug = {k: v + [shared] for k, v in base.items()}
    occ1 = cluster_classes(aug, tau=36)
    e1 = escore(occ1, "X", ["R0", "R1"])
    # the shared class is ubiquitous (q=1-eps) -> contributes -log2 q > 0,
    # identical on both sides of any target-vs-target comparison; the
    # exact cancellation law holds for DIFFERENCES between targets:
    aug2 = {"Y": base["X"] + [shared], "X": base["X"] + [shared],
            "R0": aug["R0"], "R1": aug["R1"]}
    occ2 = cluster_classes(aug2, tau=36)
    ex = escore(occ2, "X", ["R0", "R1"])
    ey = escore(occ2, "Y", ["R0", "R1"])
    assert ex == pytest.approx(ey, abs=1e-12)
    assert d_event(occ2, "X", "Y") == 0
    assert e1[0] - e0[0] == pytest.approx(-np.log2(3.0 / 4.0), abs=1e-12)


def test_channel_counts_displacement_invariant():
    ev1 = _gauss_run([2000, 4000, 6000])
    ev2 = _gauss_run([2300, 4300, 6300])   # peaks moved, still resolved
    n1 = channel_counts(ev1, 1)
    n2 = channel_counts(ev2, 1)
    assert np.array_equal(n1, n2)


def test_composition_and_mad_budget_invariance():
    rng = np.random.default_rng(1)
    n1 = rng.poisson(20.0, 300).astype(float) + 3
    n2 = 2.7 * n1                            # pure budget change
    assert d_comp(n1, n2) == pytest.approx(0.0, abs=1e-12)
    assert mad_logratio(n1, n2) == pytest.approx(0.0, abs=1e-9)
    b1, p1 = channel_state([], 4)
    assert b1 == 0.0


def test_is_recycled_pet_structured():
    tgt = _gauss_run([2000, 4000, 6000], seed=3)
    refs = {f"V{k}": _gauss_run([2000, 4000], seed=10 + k) for k in range(4)}
    out = is_recycled_pet(tgt, refs, history_score=5.0)
    assert set(out) == {"recycled_history", "authenticity_departure",
                        "direction", "domain"}
    assert out["direction"]["E_gain_bits"] > 0
    assert not isinstance(out, float)        # never a scalar, by design
