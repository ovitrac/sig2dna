"""Exact laws of the organization/grammar algebra and the sig2text
round-trip (synthetic only)."""
import numpy as np
import pytest

from sig2dna_core.grammar import (MarkovGrammar, channel_calibration,
                                  fit_conditioning, organization_scores,
                                  profile_distance, rho_e, sign_delta,
                                  tokenize)
from sig2dna_core.sig2text import Claim, Document, parse


def test_tokenize_gap_classes():
    toks = tokenize("YA" + "_" * 5 + "Z" + "_" * 20 + "B" + "_" * 100 + "Y")
    assert toks == ["Y", "A", ".", "Z", ":", "B", "|", "Y"]


def test_codelen_exact_laplace():
    # one training sequence 'YA'; scoring the same sequence:
    # steps: ^^->Y (ctx count 1/ (1+10)), ^Y->A (1/11), YA->$ (1/11)
    g = MarkovGrammar().train([["Y", "A"]])
    expect = 3 * -np.log2(2.0 / 11.0)
    assert g.codelen(["Y", "A"]) == pytest.approx(expect, abs=1e-12)
    # first step: known context (^,^) but unseen symbol Z -> 1/11;
    # second step: unseen context (^,Z) -> uniform Laplace floor 1/10
    assert g.codelen(["Z"]) == pytest.approx(
        -np.log2(1.0 / 11.0) - np.log2(1.0 / 10.0), abs=1e-12)


def test_conditioning_removes_multiplicity():
    rng = np.random.default_rng(0)
    n = rng.integers(1, 40, 400).astype(float)
    L = 7.0 + 3.5 * n + rng.normal(0, 0.1, 400)
    beta, res = fit_conditioning(L, n)
    assert beta[1] == pytest.approx(3.5, abs=0.01)
    assert abs(np.corrcoef(res, n)[0, 1]) < 0.05


def _score_setup(n_ref=8, n_ch=50, seed=1):
    rng = np.random.default_rng(seed)
    base = [["Y", "A", "Z", "B"] * int(k) for k in rng.integers(2, 6, n_ch)]
    g = MarkovGrammar().train(base)
    N = np.array([len(s) / 4 for s in base])
    refs = []
    for _ in range(n_ref):
        L = np.array([g.codelen(s) for s in base]) \
            + rng.normal(0, 0.3, n_ch)
        refs.append(L)
    beta, _ = fit_conditioning(np.concatenate(refs),
                               np.tile(N, n_ref))
    res = np.stack([r - (beta[0] + beta[1] * N) for r in refs])
    mu, sd = channel_calibration(res)
    return g, N, beta, mu, sd, base


def test_signed_directions():
    g, N, beta, mu, sd, base = _score_setup()
    # regular target: identical to training -> no gated deviation
    L0 = np.array([g.codelen(s) for s in base])
    z0, p0, n0 = organization_scores(L0, N, beta, mu, sd)
    # enriched target: inject surprising letters -> I_O+ dominates
    rich = [s[:2] + ["C", "X", "|"] * 4 + s[2:] for s in base]
    L1 = np.array([g.codelen(s) for s in rich])
    z1, p1, n1 = organization_scores(L1, N, beta, mu, sd)
    assert p1 > p0 and p1 > 10 * max(n1, 1e-9)
    assert sign_delta((p1, n1), (p0, n0)) == +1
    # regularized target: texture erased -- one motif left while the
    # resolved inventory N is unchanged (the erasure signature) -> I_O-
    reg = [["Y", "A", "Z", "B"] for _ in base]
    L2 = np.array([g.codelen(s) for s in reg])
    z2, p2, n2 = organization_scores(L2, N, beta, mu, sd)
    assert n2 > p2
    assert sign_delta((p2, n2), (p0, n0)) == -1


def test_profile_distance_level_free():
    rng = np.random.default_rng(3)
    z = rng.normal(0, 1, 200)
    assert profile_distance(z, z + 5.0) == pytest.approx(0.0, abs=1e-12)
    assert rho_e(np.array([1, 2, 100])) == 2.0


def test_gain_invariance_of_texts():
    # per-channel amplitude gains leave the symbolic text unchanged:
    # CWT is linear and the rejection threshold scales with the channel's
    # own noise -- the structural non-double-counting law.
    from sig2dna_core.icfilter import encode_ic_matrix
    rng = np.random.default_rng(4)
    t = np.arange(4000)
    y = rng.poisson(20.0, (6, 4000)).astype(float)
    for i, a in enumerate((900, 1400, 2600, 3200, 1900, 700)):
        y[i] += 3000.0 * np.exp(-((t - a) / 12.0) ** 2)
    gains = np.exp(rng.uniform(np.log(0.5), np.log(3.0), 6))
    e0 = encode_ic_matrix(y, np.arange(6), 16)
    e1 = encode_ic_matrix(y * gains[:, None], np.arange(6), 16)
    assert e0.letters == e1.letters


def test_sig2text_roundtrip_and_masks():
    doc = Document(sample_id="SYN-001", domain="synthetic-demo")
    doc.claims = [
        Claim("phi", "H", state="+", value=3.7, scope="same-domain"),
        Claim("phi", "R", state="?", value=None, scope="same-domain",
              status="candidate", mask_reason="below reference floor"),
        Claim("mu", "ORIGIN", state="", value=None, scope="cohort",
              status="suppressed",
              mask_reason="negative finding: not identifiable"),
    ]
    back = parse(doc.serialize())
    assert back == doc                      # the round-trip law
    with pytest.raises(ValueError):
        Claim("phi", "X", value=None)       # missing without reason
    with pytest.raises(ValueError):
        Claim("bad", "X", value=1.0)        # unknown namespace
