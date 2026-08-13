"""
Peak-event extraction and the event information algebra.

An **event** is a chromatographic peak detected on one ion channel: a
contiguous run of informative symbolic letters whose apex is the stable
adjacent "A|Z" zero-crossing pair of the CWT encoding. Detection is
against a *local* noise threshold (moving window), never a global cutoff,
and record-edge artifacts are excluded.

The algebra layer clusters events of several runs into per-ion event
classes (1-D clustering with a tolerance ``tau``), estimates each class's
reference occurrence probability with Laplace smoothing, and scores a
target run against a reference population:

* ``E_gain``  — summed surprisal ``-log2(q)`` of classes present in the
  target (optionally amplitude-censored against detectability),
* ``E_loss``  — summed surprisal ``-log2(1-q)`` of reference classes the
  target lacks although it could have detected them,
* ``d_event`` — a symmetric censored disagreement count usable as a
  transfer statistic.

Exact laws (used as unit tests): one isolated Gaussian above threshold
yields exactly one event across amplitudes; adding the same disjoint
population to target and references cancels exactly; a class absent from
``n`` references contributes exactly ``-log2(1/(n+2))`` bits.

Validity limits: single-column-family chromatography, EI single-quadrupole
style channel data; all scores are conditional on the chosen reference
population and on a common acquisition domain (see channels/bridges for
what survives acquisition changes).
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from sig2dna_core.icfilter import (cwt_matrix, encode_series, noise_threshold,
                                   sen_filter)


@dataclass
class PeakEvent:
    ion: int              # row index in the y matrix
    apex_idx: int         # sample index of the A|Z zero crossing
    lo: int               # event extent (samples, inclusive)
    hi: int
    a_stat: float         # detection statistic: max |height| of apex pair
    thr_at_apex: float    # local Prop. 5.1 threshold at the apex


def _split_run(letters, starts, stops, heights, idx):
    """Split one contiguous informative run (list of segment indices) at
    inter-apex boundaries; yield one (apex position, member indices,
    a_stat) triple per adjacent 'AZ' pair."""
    apexes = [k for k in range(len(idx) - 1)
              if letters[idx[k]] == "A" and letters[idx[k + 1]] == "Z"
              and stops[idx[k]] == starts[idx[k + 1]]]
    if not apexes:
        return
    # boundaries between consecutive apexes: split at the segment gap
    cuts = [0]
    for k1, k2 in zip(apexes[:-1], apexes[1:]):
        cuts.append((k1 + 1 + k2) // 2 + 1)
    cuts.append(len(idx))
    for j, k in enumerate(apexes):
        members = idx[cuts[j]:cuts[j + 1]]
        a = max(abs(heights[idx[k]]), abs(heights[idx[k + 1]]))
        yield stops[idx[k]], members, a


def extract_events(y: np.ndarray, scale: float = 16.0, window: int = 501,
                   edge: int | None = None):
    """PeakEvents of every ion of ``y`` (n_ions x n_time), full record.

    Returns (events, thr_ds, ds) where thr_ds is the threshold field
    downsampled by ds along time (for cross-run censoring lookups)."""
    y = np.asarray(y, dtype=np.float64)
    n_t = y.shape[1]
    if edge is None:
        edge = int(8 * scale)
    cfs1 = cwt_matrix(y, 1.0)
    thr = noise_threshold(cfs1, window=window)
    cfs = cfs1 if scale == 1 else cwt_matrix(y, float(scale))
    ds = 10
    thr_ds = thr[:, ::ds].astype(np.float32)
    events = []
    for i in range(y.shape[0]):
        letters, starts, stops, heights = encode_series(cfs[i])
        if not letters:
            continue
        filt = sen_filter(letters, starts, stops, heights, thr[i])
        # contiguous informative runs
        run = []
        for k, L in enumerate(filt):
            if L != "_":
                if run and starts[k] != stops[run[-1]]:
                    _flush(run, filt, starts, stops, heights, thr[i], i,
                           events, edge, n_t)
                    run = []
                run.append(k)
            elif run:
                _flush(run, filt, starts, stops, heights, thr[i], i,
                       events, edge, n_t)
                run = []
        if run:
            _flush(run, filt, starts, stops, heights, thr[i], i,
                   events, edge, n_t)
    return events, thr_ds, ds


def _flush(run, letters, starts, stops, heights, thr_i, ion, events,
           edge, n_t):
    lo, hi = starts[run[0]], stops[run[-1]]
    if lo < edge or hi > n_t - edge:      # record-boundary artifact
        return
    for apex, members, a in _split_run(letters, starts, stops, heights, run):
        events.append(PeakEvent(
            ion=ion, apex_idx=int(apex),
            lo=int(starts[members[0]]), hi=int(stops[members[-1]]),
            a_stat=float(a), thr_at_apex=float(thr_i[apex])))


# ---------------------------------------------------------------- step 0
def validate_step0(verbose: bool = True) -> bool:
    """The law: one isolated Gaussian above threshold => exactly ONE event
    with a stable apex, across >2 orders of magnitude of amplitude, for
    several noise seeds. Also: no peak => no event (edge exclusion works),
    and a resolved doublet => exactly TWO events."""
    ok = True
    for seed in (5, 11, 23):
        rng = np.random.default_rng(seed)
        n = 8000
        t = np.arange(n)
        base = rng.poisson(30.0, n).astype(float)
        apexes = {}
        for amp in (50, 100, 400, 1600, 6400, 25600):
            y = base + amp * np.exp(-((t - 4000) / 10.0) ** 2)
            ev, _, _ = extract_events(y[None, :])
            n_ev = len(ev)
            if n_ev != 1:
                ok = False
            if ev:
                apexes[amp] = ev[0].apex_idx
            if verbose:
                print(f"  seed {seed} amp {amp:6d}: events={n_ev} "
                      f"apex={[e.apex_idx for e in ev]}")
        drift = (max(apexes.values()) - min(apexes.values())) if apexes else 99
        if drift > 3:
            ok = False
        # no peak -> no events
        ev0, _, _ = extract_events(base[None, :])
        if len(ev0) != 0:
            ok = False
        if verbose:
            print(f"  seed {seed}: apex drift {drift} samples; "
                  f"no-peak events={len(ev0)}")
        # resolved doublet -> two events
        y2 = base + 2000 * (np.exp(-((t - 3800) / 10.0) ** 2)
                            + np.exp(-((t - 4200) / 10.0) ** 2))
        ev2, _, _ = extract_events(y2[None, :])
        if len(ev2) != 2:
            ok = False
        if verbose:
            print(f"  seed {seed}: doublet events={len(ev2)} "
                  f"apexes={[e.apex_idx for e in ev2]}")
    return ok


if __name__ == "__main__":
    print("STEP-0 VALIDATION (single-Gaussian one-event law)")
    print("PASS" if validate_step0() else "FAIL")


# ---------------------------------------------------------------- algebra
def cluster_classes(populations: dict, tau: float):
    """Cluster the events of several runs into per-ion event classes.

    ``populations`` maps run name -> iterable of PeakEvent (or tuples
    ``(ion, apex_idx, lo, hi, a_stat, thr_at_apex)``). Returns a dict
    ``(ion, class_id) -> {run: (a_stat, thr_at_apex, apex_idx)}`` where a
    new class starts whenever the apex gap exceeds ``tau`` samples.
    """
    occ = {}
    pts_by_ion = {}
    for name, evs in populations.items():
        for e in evs:
            ion, apex, a, th = ((e.ion, e.apex_idx, e.a_stat, e.thr_at_apex)
                                if isinstance(e, PeakEvent)
                                else (e[0], e[1], e[4], e[5]))
            pts_by_ion.setdefault(ion, []).append((apex, name, a, th))
    for ion, pts in pts_by_ion.items():
        pts.sort()
        cid, prev = 0, None
        for apex, name, a, th in pts:
            if prev is not None and apex - prev > tau:
                cid += 1
            prev = apex
            d = occ.setdefault((ion, cid), {})
            if name not in d or a > d[name][0]:
                d[name] = (a, th, apex)
    return occ


def escore(occ, target: str, references: list):
    """(E_gain, E_loss) of ``target`` against ``references`` (bits).

    Laplace-smoothed occurrence ``q = (k+1)/(n_ref+2)``; a class present
    in the target adds ``-log2 q``; a reference class the target lacks
    adds ``-log2(1-q)``. Censoring against detectability thresholds is
    the caller's concern (drop non-detectable classes from ``occ`` first);
    scores are conditional on the reference population.
    """
    n = len(references)
    eg = el = 0.0
    for d in occ.values():
        k = sum(1 for r in references if r in d)
        q = (k + 1.0) / (n + 2.0)
        if target in d:
            eg += -np.log2(q)
        elif k > 0:
            el += -np.log2(1.0 - q)
    return eg, el


def d_event(occ, a: str, b: str):
    """Symmetric event disagreement between two runs: the number of
    classes present in exactly one of them (censoring, if any, applied
    upstream). A transfer-stable statistic within one acquisition domain."""
    return sum(1 for d in occ.values() if (a in d) != (b in d))


def is_recycled_pet(target_events, reference_events: dict,
                    history_score: float | None = None, tau: float = 36.0):
    """Structured verdict for one sample against virgin references.

    Never returns a single scalar: recycling *history* and *present
    chemistry* are independent memories (a strongly recycled PET can be
    chemically stripped below the virgin floor, and a virgin resin can
    carry an atypical formulation).

    ``history_score`` is the caller's dedicated history coordinate (e.g.
    a degradation-entropy z-score from the symbolic route); without it
    the history field is reported as ``"not evaluated"``.
    """
    pops = {"__target__": target_events}
    pops.update(reference_events)
    refs = list(reference_events)
    occ = cluster_classes(pops, tau)
    eg, el = escore(occ, "__target__", refs)
    dts = [d_event(occ, "__target__", r) for r in refs]
    null = [d_event(occ, a, b) for a, b in combinations(refs, 2)]
    mu = float(np.mean(null)) if null else float("nan")
    sd = float(np.std(null, ddof=1)) if len(null) > 1 else float("nan")
    z = (float(np.median(dts)) - mu) / sd if sd and sd == sd else float("nan")
    return {
        "recycled_history": ({"score": history_score,
                              "verdict": "recycled-history evidence"
                              if history_score is not None and history_score > 3
                              else "virgin-compatible history"}
                             if history_score is not None else "not evaluated"),
        "authenticity_departure": {"d_event_median": float(np.median(dts)),
                                   "z_vs_reference_null": z},
        "direction": {"E_gain_bits": eg, "E_loss_bits": el},
        "domain": ("valid only within one acquisition domain and against "
                   "this reference population; cross-domain comparisons "
                   "require measured bridges (sig2dna_core.bridges)"),
    }
