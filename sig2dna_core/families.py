"""
Compound families: co-elution grouping, pseudo-EI spectra, and
decoy-calibrated recognition.

Events that co-elute (apex gap <= ``gap`` samples) are grouped into a
**compound family** summarized by a pseudo-EI spectrum (per-channel
maximum event amplitude). Families are deduplicated, and recognized in
other runs, with one frozen identity rule: apex distance within ``tau``
AND spectral cosine >= ``cos_min``.

Recognition is only ever claimed against a **decoy null**: real spectra
displaced to wrong positions. A detection must beat what chance
coincidence achieves in equally crowded data. Family identity is a
*within-acquisition-domain* notion: across sessions it degrades (large
chained families first) and across observation operators it does not
transfer — dictionaries interpret measurements, they never define them.
"""
from __future__ import annotations

import numpy as np


def cosine(s1: dict, s2: dict) -> float:
    """Cosine similarity of two sparse spectra ``{channel: intensity}``."""
    ions = set(s1) | set(s2)
    v1 = np.array([s1.get(i, 0.0) for i in ions])
    v2 = np.array([s2.get(i, 0.0) for i in ions])
    d = np.linalg.norm(v1) * np.linalg.norm(v2)
    return float(v1 @ v2 / d) if d > 0 else 0.0


def group_families(events, gap: int = 12, min_ions: int = 3):
    """Co-elution grouping of one run's events into compound families.

    Returns ``[{"apex": int, "spec": {channel: max amplitude}, "n_ions"}]``
    keeping only families with >= ``min_ions`` distinct channels."""
    evs = sorted(((e.ion, e.apex_idx, e.a_stat) if hasattr(e, "ion")
                  else (e[0], e[1], e[4]) for e in events),
                 key=lambda t: t[1])
    groups, cur, last = [], [], None
    for e in evs:
        if cur and e[1] - last > gap:
            groups.append(cur)
            cur = []
        cur.append(e)
        last = e[1]
    if cur:
        groups.append(cur)
    out = []
    for g in groups:
        spec = {}
        for ion, apex, a in g:
            spec[ion] = max(spec.get(ion, 0.0), float(a))
        if len(spec) >= min_ions:
            out.append(dict(apex=int(np.median([e[1] for e in g])),
                            spec=spec, n_ions=len(spec)))
    return out


def dedup(fams, tau: float = 36.0, cos_min: float = 0.7):
    """Merge families that satisfy the identity rule (apex within tau AND
    cosine >= cos_min), keeping per-channel maxima."""
    fams = sorted(fams, key=lambda f: f["apex"])
    out = []
    for f in fams:
        for g in out:
            if (abs(f["apex"] - g["apex"]) <= tau
                    and cosine(f["spec"], g["spec"]) >= cos_min):
                for i, a in f["spec"].items():
                    g["spec"][i] = max(g["spec"].get(i, 0.0), a)
                g["n_ions"] = len(g["spec"])
                break
        else:
            out.append(dict(f))
    return out


def make_decoys(fams, n_time: int, rng, n_sets: int = 3,
                shift_lo: int = 100, shift_hi: int = 400, tau: float = 36.0):
    """Decoy families: real spectra displaced by random apex offsets
    (|shift| in [shift_lo, shift_hi]), rejected if any of the family's
    top-5 channels lands within 2*tau of a true same-channel family.
    The decoy recognition rate is the chance floor any claimed detection
    must beat."""
    occupied = {}
    for f in fams:
        for m in sorted(f["spec"], key=lambda k: -f["spec"][k])[:5]:
            occupied.setdefault(m, []).append(f["apex"])
    decoys = []
    for _ in range(n_sets):
        for f in fams:
            top5 = sorted(f["spec"], key=lambda k: -f["spec"][k])[:5]
            for _try in range(300):
                d = int(rng.integers(shift_lo, shift_hi + 1))
                if rng.random() < 0.5:
                    d = -d
                apex = f["apex"] + d
                if not (0 < apex < n_time):
                    continue
                if any(abs(apex - x) < 2 * tau for m in top5
                       for x in occupied.get(m, [])):
                    continue
                decoys.append(dict(apex=apex, spec=f["spec"],
                                   n_ions=f["n_ions"]))
                break
    return decoys


def recognize(family, run_events, tau: float = 36.0, m_min: int = 3):
    """Frozen matcher: ``family`` is RECOVERED in a run iff at least
    ``m_min`` of its top-5 channels have a same-channel event within
    ``tau`` of the family apex. Use identically on true families and
    decoys; report both rates."""
    by_ion = {}
    for e in run_events:
        ion, apex = (e.ion, e.apex_idx) if hasattr(e, "ion") else (e[0], e[1])
        by_ion.setdefault(ion, []).append(apex)
    top5 = sorted(family["spec"], key=lambda k: -family["spec"][k])[:5]
    hits = sum(1 for m in top5
               if any(abs(x - family["apex"]) <= tau
                      for x in by_ion.get(m, [])))
    return hits >= m_min
