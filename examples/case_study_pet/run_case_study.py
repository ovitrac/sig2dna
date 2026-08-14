#!/usr/bin/env python3
"""
Case study — one virgin and one recycled PET, four readings of the same
two signals: CWT letters (shape/organization), event/channel algebra
(inventory/amount), the organization/grammar algebra connecting them,
and pseudo-EI spectra reconstructed from co-eluting differential events
(recognizable chemistry: the two faces of recycling).

Route 4 is blind-first: the event algebra nominates the event classes
gained by the recycled run, chain-clusters their apexes into co-eluting
families (|t_m - t*| <= tau), and reconstructs each family's spectrum
from the raw local-baseline-corrected ion amplitudes (base peak = 100;
event surprisal bits are never used as intensities). Only afterwards is
each spectrum compared with the frozen candidate-marker dictionaries
(degradation/history vs contamination/cleanliness); exemplars are
chosen at the poles, gated on genuine amplitude excess
(X = sum S+ / sum S_R >= 0.25) so apex/class shifts are named as such,
and the strongest excess families outside both dictionaries are shown
as their own category instead of being forced onto the nearest marker.

Fully replayable: everything in report.html is recomputed from the
shipped anonymized data (case_vr.npz: two full-scan GC-MS runs,
n_channels x n_time integer intensities; 0.14 s/scan; m/z 40-650;
all identifiers, timestamps and provenance removed). The report ships
with light and dark themes (system default, manual toggle, print=light).

CALIBRATION DEMONSTRATION: the dataset ships, alongside the two runs,
an anonymized reference-population calibration (per-channel location and
spread of grammar residuals and event counts from nine anonymized virgin
references, leave-one-out; mode = reference_population, R2). The grammar
section scores the SAME residuals under three calibrations — R0 global
proxy, R1 density-conditioned sigma(N) fitted on the single shipped
reference, R2 shipped reference-population — and shows how the apparent
low-m/z concentration of outliers under R0 is a heteroscedasticity
artifact that calibration removes. No replicate floors exist in this
example; even R2 remains reference-population grade, not
technical-replicate grade (R3).

Run:  python run_case_study.py     (writes report.html next to itself)
"""
import base64
import io
import os
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sig2dna_core.events import cluster_classes, escore, extract_events
from sig2dna_core.channels import channel_counts, d_comp
from sig2dna_core.grammar import (DensityCalibration, MarkovGrammar,
                                  fit_conditioning, organization_scores,
                                  profile_distance, rho_e, sign_delta,
                                  tokenize)
from sig2dna_core.icfilter import encode_ic_matrix, text_entropy

HERE = os.path.dirname(os.path.abspath(__file__))
TAU = 36.0
SCALE = 16.0
GATE = 2.0
# route 4 (events -> chemistry): pseudo-EI reconstruction
MIN_FRAG = 3       # a co-eluting family needs >= this many distinct ions
TOP_FAM = 12       # reconstruct at most this many families (by support)
BASE_W = 350       # samples each side of t* for the local baseline median
N_EXEMPLAR = 3     # pole exemplars per axis, at most
X_MIN = 0.25       # exemplars need this excess-support fraction
# Candidate marker dictionaries for the two faces of recycling
# (family-level EI fragments; candidate assignments pending library
# confirmation -- labels deliberately stay generic, never molecule names).
DICT_D = {45: "dioxolane/acetaldehyde-EG-like", 87: "dioxolane-like",
          96: "furanic-like", 109: "methylfuranic-like",
          148: "anhydride-like", 572: "cyclic-oligomer-like"}
DICT_C = {55: "alkene-like", 56: "alkene-like", 59: "oxygenate-like",
          69: "terpene/alkene-like", 72: "oxygenate-like",
          74: "fatty-ester-like", 86: "oxygenate-like",
          88: "oxygenate-like", 91: "alkylbenzene-like",
          121: "terpene-like", 165: "PAH-like"}

THEMES = {
    "light": dict(bg="#ffffff", ink="#1c2321", mut="#5c6663",
                  grid="#e3e0d8", v="#2e7d54", r="#b3542e",
                  acc="#22536b", band="#00000012"),
    "dark": dict(bg="#1d2326", ink="#e6e2d8", mut="#9aa4a0",
                 grid="#333a3d", v="#5fae86", r="#d08159",
                 acc="#7fb4cf", band="#ffffff14"),
}


def figbytes(fig):
    """Render once as SVG; the same bytes feed the HTML (base64-inline)
    and the Markdown report (committed files under figures/)."""
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def themed_axes(theme, figsize=(9, 2.8)):
    t = THEMES[theme]
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor(t["bg"])
    ax.set_facecolor(t["bg"])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(t["mut"])
    ax.tick_params(colors=t["mut"], labelcolor=t["ink"])
    ax.xaxis.label.set_color(t["ink"])
    ax.yaxis.label.set_color(t["ink"])
    return fig, ax, t


def legend(ax, t):
    leg = ax.legend(frameon=False)
    for txt in leg.get_texts():
        txt.set_color(t["ink"])


def local_amp(y, ion, apex, t_star):
    """Raw amplitude at ``apex`` minus the local baseline (median of a
    window centred on the family apex t*). Floored at zero. The spectrum
    intensities come from the SIGNAL, never from event surprisal bits."""
    lo, hi = max(0, t_star - BASE_W), min(y.shape[1], t_star + BASE_W)
    return max(float(y[ion, apex] - np.median(y[ion, lo:hi])), 0.0)


def ref_amp(y, ion, t_star, tau):
    """Best baseline-corrected response of ``ion`` in the reference run
    within +-tau of t* (tolerates small retention shifts)."""
    lo = max(0, int(t_star - tau))
    hi = min(y.shape[1], int(t_star + tau) + 1)
    wlo, whi = max(0, t_star - BASE_W), min(y.shape[1], t_star + BASE_W)
    return max(float(y[ion, lo:hi].max()
                     - np.median(y[ion, wlo:whi])), 0.0)


def dict_frac(spec_mz, keys, top=10):
    """Fraction of the spectrum's SIGNATURE intensity (its ``top`` most
    intense fragments, standard EI main-peak practice) carried by
    dictionary channels — a full-spectrum fraction would bury rich
    spectra under their minor fragments. Matching is EXACT unit-mass: a
    +-1 tolerance conflates chemically opposite markers (e.g.
    phthalic-anhydride 148, a degradation marker, with the
    dialkyl-phthalate 149 of plasticizer contamination)."""
    sig = dict(sorted(spec_mz.items(), key=lambda kv: kv[1])[-top:])
    tot = sum(sig.values())
    if tot <= 0:
        return 0.0
    return sum(v for m, v in sig.items() if m in keys) / tot


def main():
    t0 = time.time()
    z = np.load(os.path.join(HERE, "case_vr.npz"))
    yV = z["y_virgin"].astype(float)
    yR = z["y_recycled"].astype(float)
    dt = float(z["dt_s"])
    mz0 = int(z["mz_first"])
    n_ch, n_t = yV.shape
    mz = np.arange(n_ch) + mz0
    rt_min = np.arange(n_t) * dt / 60.0
    print(f"loaded: {n_ch} channels x {n_t} scans, dt={dt:.2f}s")

    stats = {}

    # ---------------- route 1: CWT letters -------------------------------
    print("route 1: encoding letters...", flush=True)
    encV = encode_ic_matrix(yV, mz, SCALE)
    encR = encode_ic_matrix(yR, mz, SCALE)
    entV = np.array([text_entropy(s) for s in encV.letters])
    entR = np.array([text_entropy(s) for s in encR.letters])

    def lcount(letters):
        c = {}
        for s in letters:
            for ch in s:
                c[ch] = c.get(ch, 0) + 1
        return c

    stats["lettersV"], stats["lettersR"] = lcount(encV.letters), lcount(
        encR.letters)
    demo_ch = int(np.argmax([len(s.replace("_", ""))
                             for s in encV.letters]))
    stats["demo_mz"] = int(mz[demo_ch])
    stats["textV"] = encV.letters[demo_ch]
    stats["textR"] = encR.letters[demo_ch]

    # ---------------- route 2: events ------------------------------------
    print("route 2: extracting events...", flush=True)
    evV, _, _ = extract_events(yV, scale=SCALE, window=501)
    evR, _, _ = extract_events(yR, scale=SCALE, window=501)
    nV = channel_counts(evV, n_ch)
    nR = channel_counts(evR, n_ch)
    stats["BV"], stats["BR"] = int(nV.sum()), int(nR.sum())
    stats["rhoV"], stats["rhoR"] = rho_e(nV), rho_e(nR)
    stats["dcomp"] = d_comp(nV, nR)

    def tup(evs):
        return [(e.ion, e.apex_idx, e.lo, e.hi, e.a_stat, e.thr_at_apex)
                for e in evs]

    occ = cluster_classes({"R": tup(evR), "V": tup(evV)}, tau=TAU)
    eg, el = escore(occ, "R", ["V"])
    stats["Egain"], stats["Eloss"] = eg, el
    gain_ch = np.maximum(nR - nV, 0)
    top_idx = np.argsort(-gain_ch)[:12]
    stats["topgain"] = [(int(mz[i]), int(gain_ch[i])) for i in top_idx
                        if gain_ch[i] > 0]

    # ---------------- route 4: events -> chemistry -----------------------
    # Blind-first: the event algebra nominates gained co-eluting families;
    # spectra are reconstructed from RAW amplitudes; only afterwards are
    # they compared with the candidate marker dictionaries.
    print("route 4: reconstructing pseudo-EI spectra...", flush=True)
    ticV, ticR = yV.sum(0), yR.sum(0)
    ok = (ticV > 0) & (ticR > 0)
    alpha = float(np.median(ticR[ok] / ticV[ok]))  # response calibration
    stats["alpha"] = alpha
    gained = sorted((d["R"][2], ion, d["R"][0], d["R"][1])
                    for (ion, cid), d in occ.items()
                    if "R" in d and "V" not in d)
    stats["n_gained_cls"] = len(gained)
    fams = []                       # chain clustering of gained apexes
    for rec in gained:
        if fams and rec[0] - fams[-1][-1][0] <= TAU:
            fams[-1].append(rec)
        else:
            fams.append([rec])
    fams = [f for f in fams if len({r[1] for r in f}) >= MIN_FRAG]
    fams.sort(key=lambda f: (-len({r[1] for r in f}),
                             -sum(r[2] for r in f)))
    stats["n_fam_all"] = len(fams)
    fams = fams[:TOP_FAM]
    families = []
    for f in fams:
        t_star = max(f, key=lambda r: r[2])[0]   # strongest gained apex
        gained_ions = {r[1] for r in f}
        sup = {}      # ALL co-eluting R events, strongest per ion
        for e in evR:
            if abs(e.apex_idx - t_star) <= TAU:
                if e.ion not in sup or e.a_stat > sup[e.ion].a_stat:
                    sup[e.ion] = e
        SR = {}
        for ion, e in sup.items():
            a = local_amp(yR, ion, e.apex_idx, t_star)
            if a > 0:
                SR[ion] = a
        if len(SR) < MIN_FRAG:
            continue
        SV = {ion: ref_amp(yV, ion, t_star, TAU) for ion in SR}
        SP = {ion: max(SR[ion] - alpha * SV[ion], 0.0) for ion in SR}
        smz = {int(mz[i]): v for i, v in SR.items()}
        D_c = dict_frac(smz, DICT_D)
        C_c = dict_frac(smz, DICT_C)
        base_mz = max(smz, key=smz.get)
        # generic family label from the most intense matched marker of
        # the signature (top-10 fragments), consistent with dict_frac
        sig = dict(sorted(smz.items(), key=lambda kv: kv[1])[-10:])
        lbl, best_v = "", 0.0
        for dic in (DICT_D, DICT_C):
            for m, v in sig.items():
                if m in dic and v > best_v:
                    lbl, best_v = dic[m], v
        if D_c >= 2 * C_c and D_c >= 0.15:
            reading = "degradation-like"
        elif C_c >= 2 * D_c and C_c >= 0.15:
            reading = "contamination-like"
        else:
            reading = "mixed / unresolved"
        # excess support: fraction of the R spectrum's intensity that is
        # genuine excess after response matching. Low X means the family
        # was nominated by apex/class shifts (inventory change without
        # amplitude excess), not by new chemistry.
        X = sum(SP.values()) / sum(SR.values())
        families.append(dict(
            t_star=int(t_star), rt=t_star * dt / 60.0, SR=SR, SV=SV, SP=SP,
            gained=gained_ions, D=D_c, C=C_c, X=X, base_mz=int(base_mz),
            n_frag=len(SR), n_gained=len(gained_ions), label=lbl,
            reading=reading))
    stats["families"] = families
    deg = sorted((f for f in families if f["reading"] == "degradation-like"
                  and f["X"] >= X_MIN),
                 key=lambda f: -(f["D"] - f["C"]))[:N_EXEMPLAR]
    con = sorted((f for f in families
                  if f["reading"] == "contamination-like"
                  and f["X"] >= X_MIN),
                 key=lambda f: -(f["C"] - f["D"]))[:N_EXEMPLAR]
    # strongest genuine-excess families that neither dictionary spans:
    # shown as their own category rather than forced onto a pole
    unm = sorted((f for f in families if f["reading"] == "mixed / unresolved"
                  and f["X"] >= 0.5), key=lambda f: -f["X"])[:N_EXEMPLAR]
    stats["exemplars"] = deg + con + unm
    stats["n_deg"], stats["n_con"] = len(deg), len(con)
    stats["n_unm"] = len(unm)
    stats["has_himass"] = any(f["base_mz"] > 500 or
                              any(int(mz[i]) > 500 for i in f["SR"])
                              for f in families)
    print(f"  {len(gained)} gained classes -> {stats['n_fam_all']} families "
          f"(>= {MIN_FRAG} ions); reconstructed {len(families)}; "
          f"exemplars {len(deg)} degradation-like / {len(con)} "
          f"contamination-like; alpha={alpha:.3f}", flush=True)

    # ---------------- route 3: grammar, three calibrations ---------------
    print("route 3: grammar...", flush=True)
    toksV = [tokenize(s) for s in encV.letters]
    toksR = [tokenize(s) for s in encR.letters]
    # shipped reference grammar (trained on the same anonymized virgin
    # references as the shipped calibration -- aggregate transition
    # counts only, part of the calibration bundle)
    import gzip as _gzip
    import json as _json
    G = MarkovGrammar()
    with _gzip.open(os.path.join(HERE, "case_grammar.json.gz"), "rt") as f:
        for a, b, cnt in _json.load(f):
            for tok, n in cnt.items():
                G.ctx[(a, b)][tok] += n
    LV = np.array([G.codelen(s) for s in toksV])
    LR = np.array([G.codelen(s) for s in toksR])
    beta = z["beta_ref"].astype(float)      # shipped reference conditioning
    stats["beta"] = beta
    rV = LV - (beta[0] + beta[1] * nV)
    rR = LR - (beta[0] + beta[1] * nR)
    stats["R2"] = 1 - rV.var() / LV.var()
    mu_L, sd_L = z["mu_L"].astype(float), z["sd_L"].astype(float)
    mu_N, sd_N = z["mu_N"].astype(float), z["sd_N"].astype(float)
    stats["n_ref"] = int(z["n_ref"])
    # R0: one global sigma (didactic visualization only)
    z0 = (rR - rV.mean()) / max(rV.std(ddof=1), 1.0)
    # R1: density-conditioned sigma(N) fitted on the single reference
    cal1 = DensityCalibration(n_bins=8).fit(rV, nV)
    z1 = cal1.z(rR, nR)
    # R2: shipped reference-population per-channel calibration
    z2 = (rR - mu_L) / sd_L
    z2V = (rV - mu_L) / sd_L
    low = mz < 200
    def ratio(zz):
        hit = np.abs(zz) > GATE
        return hit[low].mean(), hit[~low].mean()
    stats["ladder"] = {name: ratio(zz) for name, zz in
                       (("R0 global_proxy", z0),
                        ("R1 density_conditioned", z1),
                        ("R2 reference_population", z2))}
    d2 = z2 * sd_L
    d2V = z2V * sd_L
    stats["ioR"] = (float(d2[z2 > GATE].sum()),
                    float(-d2[z2 < -GATE].sum()))
    stats["ioV"] = (float(d2V[z2V > GATE].sum()),
                    float(-d2V[z2V < -GATE].sum()))
    stats["sign"] = sign_delta(stats["ioR"], stats["ioV"])
    stats["pdist"] = profile_distance(z2V, z2)
    # channel-state quadrants (ratified display): inventory vs organization
    zN = (nR - mu_N) / sd_N
    inv = np.abs(zN) > GATE
    org = np.abs(z2) > GATE
    stats["quad"] = {}
    for lbl, m_ in (("stable", ~inv & ~org), ("inventory-only", inv & ~org),
                    ("organization-only", ~inv & org), ("both", inv & org)):
        stats["quad"][lbl] = (float(m_[low].mean()), float(m_[~low].mean()))
    stats["zN"], stats["z0"], stats["z2"] = zN, z0, z2

    # ---------------- figures, both themes -------------------------------
    print("rendering figures (light + dark)...", flush=True)
    figs = {}
    for theme in THEMES:
        fig, ax, t = themed_axes(theme)
        ax.plot(rt_min, yV.sum(0), color=t["v"], lw=0.7,
                label="virgin (V)")
        ax.plot(rt_min, yR.sum(0), color=t["r"], lw=0.7, alpha=0.85,
                label="recycled (R)")
        ax.set_xlabel("retention time (min)")
        ax.set_ylabel("total ion current")
        legend(ax, t)
        figs[("tic", theme)] = figbytes(fig)

        fig, ax, t = themed_axes(theme)
        ax.plot(mz, entV, color=t["v"], lw=0.8, label="virgin")
        ax.plot(mz, entR, color=t["r"], lw=0.8, alpha=0.85,
                label="recycled")
        ax.set_xlabel("m/z channel")
        ax.set_ylabel("text entropy (bits)")
        legend(ax, t)
        figs[("entropy", theme)] = figbytes(fig)

        fig, ax, t = themed_axes(theme)
        ax.plot(mz, nV, color=t["v"], lw=0.8, label="virgin")
        ax.plot(mz, nR, color=t["r"], lw=0.8, alpha=0.85,
                label="recycled")
        ax.set_xlabel("m/z channel")
        ax.set_ylabel("resolved events per channel")
        legend(ax, t)
        figs[("counts", theme)] = figbytes(fig)

        fig, axs = plt.subplots(2, 1, figsize=(9, 4.6), sharex=True)
        tt = THEMES[theme]
        fig.patch.set_facecolor(tt["bg"])
        for ax, zz, lbl in ((axs[0], stats["z0"],
                             "R0: global proxy (didactic only)"),
                            (axs[1], stats["z2"],
                             "R2: shipped reference-population "
                             "calibration")):
            ax.set_facecolor(tt["bg"])
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color(tt["mut"])
            ax.tick_params(colors=tt["mut"], labelcolor=tt["ink"])
            ax.axhspan(-GATE, GATE, color=tt["band"], lw=0)
            ax.plot(mz, zz, color=tt["r"], lw=0.8)
            ax.axhline(0, color=tt["mut"], lw=0.6)
            ax.set_ylabel("z", color=tt["ink"])
            ax.set_title(lbl, fontsize=9, color=tt["ink"], loc="left")
        axs[1].set_xlabel("m/z channel", color=tt["ink"])
        figs[("grammar", theme)] = figbytes(fig)

        fig, ax, t = themed_axes(theme, figsize=(5.6, 4.2))
        m_low = mz < 200
        ax.scatter(stats["zN"][~m_low], stats["z2"][~m_low], s=8,
                   c=t["mut"], alpha=.55, label="m/z \u2265 200")
        ax.scatter(stats["zN"][m_low], stats["z2"][m_low], s=10,
                   c=t["r"], alpha=.8, label="m/z < 200")
        for g in (-GATE, GATE):
            ax.axhline(g, color=t["grid"], lw=0.8)
            ax.axvline(g, color=t["grid"], lw=0.8)
        ax.set_xlabel("inventory z  (event count vs reference)")
        ax.set_ylabel("organization z  (R2-calibrated)")
        legend(ax, t)
        figs[("quad", theme)] = figbytes(fig)

        # route 4 figures: D-C family map + three-spectra exemplars
        fig, ax, t = themed_axes(theme, figsize=(5.6, 4.2))
        col = {"degradation-like": t["acc"], "contamination-like": t["r"],
               "mixed / unresolved": t["mut"]}
        seen = set()
        for f in stats["families"]:
            lab = f["reading"] if f["reading"] not in seen else None
            seen.add(f["reading"])
            ax.scatter(f["C"], f["D"], s=28 + 3 * f["n_frag"],
                       c=col[f["reading"]], alpha=.85, label=lab)
        for f in stats["exemplars"]:
            ax.annotate(f"{f['rt']:.1f} min", (f["C"], f["D"]),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=7, color=t["ink"])
        lim = max([0.3] + [f["D"] for f in stats["families"]]
                  + [f["C"] for f in stats["families"]]) * 1.15
        ax.plot([0, lim], [0, lim], color=t["grid"], lw=0.8)
        ax.set_xlim(-0.01, lim)
        ax.set_ylim(-0.01, lim)
        ax.set_xlabel("C  (contamination-dictionary intensity fraction)")
        ax.set_ylabel("D  (degradation-dictionary intensity fraction)")
        legend(ax, t)
        figs[("fammap", theme)] = figbytes(fig)

        for k, f in enumerate(stats["exemplars"]):
            smax = max(f["SR"].values())
            fam_mz = [int(mz[i]) for i in f["SR"]]
            xlo = max(min(fam_mz) - 8, int(mz[0]) - 2)
            xhi = max(fam_mz) + 14
            fig, axs = plt.subplots(3, 1, figsize=(9, 5.2), sharex=True)
            tt = THEMES[theme]
            fig.patch.set_facecolor(tt["bg"])
            rows = ((axs[0], f["SV"], tt["v"],
                     "virgin reference, response-matched "
                     "(\u00d7\u03b1, same scale)", alpha),
                    (axs[1], f["SR"], tt["r"],
                     "recycled raw pseudo-EI \u2014 the chemical object",
                     1.0),
                    (axs[2], f["SP"], tt["acc"],
                     "excess  [S_R \u2212 \u03b1\u00b7S_V]\u208a \u2014 "
                     "explanatory overlay", 1.0))
            for ax, S, c, lab, mult in rows:
                ax.set_facecolor(tt["bg"])
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
                for sp in ("left", "bottom"):
                    ax.spines[sp].set_color(tt["mut"])
                ax.tick_params(colors=tt["mut"], labelcolor=tt["ink"])
                vals = {int(mz[i]): 100.0 * mult * v / smax
                        for i, v in S.items()}
                if vals:
                    ax.vlines(list(vals), 0, list(vals.values()),
                              color=c, lw=1.6)
                top5 = sorted(vals, key=vals.get)[-5:]
                for m in top5:
                    if vals[m] > 4:
                        ax.annotate(str(m), (m, vals[m]),
                                    textcoords="offset points",
                                    xytext=(0, 2), ha="center",
                                    fontsize=7, color=tt["ink"])
                ax.set_ylim(0, 112)
                ax.set_ylabel("% base", color=tt["ink"], fontsize=8)
                ax.set_title(lab, fontsize=8.5, color=tt["ink"],
                             loc="left")
            axs[2].set_xlabel("m/z", color=tt["ink"])
            axs[0].set_xlim(xlo, xhi)
            beyond = f["reading"] == "mixed / unresolved"
            cat = "beyond the frozen dictionaries" if beyond else f["reading"]
            fig.suptitle(
                cat + (f" \u00b7 {f['label']}" if f["label"] and not beyond
                       else "")
                + f" \u00b7 RT {f['rt']:.1f} min \u00b7 "
                  f"{f['n_frag']} fragments ({f['n_gained']} gained)",
                fontsize=9.5, color=tt["ink"], x=0.12, ha="left")
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            figs[(f"spec{k}", theme)] = figbytes(fig)

    write_report(figs, stats, n_ch, n_t, dt, time.time() - t0)
    write_markdown(figs, stats, n_ch, n_t, dt)


def img_pair(figs, name, alt):
    b64l = base64.b64encode(figs[(name, "light")]).decode()
    b64d = base64.b64encode(figs[(name, "dark")]).decode()
    return (f'<img class="fl" src="data:image/svg+xml;base64,'
            f'{b64l}" alt="{alt}">'
            f'<img class="fd" src="data:image/svg+xml;base64,'
            f'{b64d}" alt="{alt}">')


def md_pic(name, alt):
    """GitHub-rendered dual-theme figure: the <picture> element with a
    prefers-color-scheme source is supported in GitHub Markdown."""
    return (f'<picture>\n'
            f'<source media="(prefers-color-scheme: dark)" '
            f'srcset="figures/{name}_dark.svg">\n'
            f'<img src="figures/{name}_light.svg" alt="{alt}">\n'
            f'</picture>')


def write_report(figs, s, n_ch, n_t, dt, elapsed):
    lv = s["lettersV"]
    lr = s["lettersR"]
    keys = [k for k in "YAZBCX_" if k in lv or k in lr]
    lrow_v = "".join(f"<td>{lv.get(k, 0):,}</td>" for k in keys)
    lrow_r = "".join(f"<td>{lr.get(k, 0):,}</td>" for k in keys)
    lhead = "".join(f"<th><code>{k}</code></th>" for k in keys)
    topg = " ".join(f"<span class='chip'>m/z {m} (+{g})</span>"
                    for m, g in s["topgain"])
    sign_word = {1: "net enrichment (I<sub>O</sub><sup>+</sup> dominated)",
                 -1: "net regularization (I<sub>O</sub><sup>&minus;</sup> "
                     "dominated)", 0: "balanced"}[s["sign"]]

    def _ratio(lo, hi):
        if hi == 0:
            return ("&infin; (no channel above the gate at m/z &ge; 200)"
                    if lo > 0 else "&mdash;")
        return f"{lo / hi:.2f}"
    ladder_rows = "".join(
        f"<tr><td><code>{name}</code></td><td>{lo:.3f}</td>"
        f"<td>{hi:.3f}</td><td><b>{_ratio(lo, hi)}</b></td></tr>"
        for name, (lo, hi) in s["ladder"].items())
    quad_rows = "".join(
        f"<tr><td>{k}</td><td>{a:.3f}</td><td>{b:.3f}</td></tr>"
        for k, (a, b) in s["quad"].items())
    img_grammar = img_pair(figs, "grammar", "organization residuals, R0 vs R2")
    img_quad = img_pair(figs, "quad", "channel-state quadrants")
    fam_rows = "".join(
        f"<tr><td>{f['rt']:.1f}</td><td>{f['n_frag']}</td>"
        f"<td>{f['n_gained']}</td><td>{f['base_mz']}</td>"
        f"<td>{f['D']:.2f}</td><td>{f['C']:.2f}</td>"
        f"<td>{f['X']:.2f}</td>"
        f"<td>{f['reading']}"
        f"{' &mdash; ' + f['label'] if f['label'] else ''}"
        f"{'' if f['X'] >= 0.25 else ' <span class=chip>class shift</span>'}"
        f"</td></tr>"
        for f in s["families"])
    cat_head = {"degradation-like": "Degradation face",
                "contamination-like": "Contamination face",
                "mixed / unresolved": "Beyond the frozen dictionaries"}
    spec_html = "".join(
        f"<h3>{cat_head[f['reading']]} &mdash; "
        + (f["label"] if f["label"]
           and f["reading"] != "mixed / unresolved"
           else "signature outside both dictionaries")
        + f" (RT {f['rt']:.1f}&thinsp;min, base peak m/z {f['base_mz']})</h3>"
        + img_pair(figs, f"spec{k}", "three-spectra pseudo-EI reconstruction")
        for k, f in enumerate(s["exemplars"]))
    himass_note = ("" if s["has_himass"] else
                   " No credible high-mass family (m/z &gt; 500, e.g. a "
                   "cyclic-oligomer-type marker) was nominated by the "
                   "events of this particular pair; the panel reports what "
                   "the data support and does not manufacture known cohort "
                   "markers.")
    tau_s = TAU * dt
    n_fam_rec = len(s["families"])

    def trim(t, n=90):
        t = t.strip("_")
        return (t[:n] + "&hellip;") if len(t) > n else (t or "(empty)")

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>sig2dna case study — one virgin and one recycled PET, four readings</title>
<style>
  :root {{ --ink:#1c2321; --mut:#5c6663; --bg:#fbfaf7; --card:#ffffff;
          --line:#e3e0d8; --v:#2e7d54; --r:#b3542e; --acc:#22536b;
          --notebg:#f3efe4; --noteln:#c8b37a; --chip:#eee9dc; }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{ --ink:#e6e2d8; --mut:#9aa4a0;
      --bg:#14181a; --card:#1d2326; --line:#333a3d; --v:#5fae86;
      --r:#d08159; --acc:#7fb4cf; --notebg:#242a21; --noteln:#8a7c4f;
      --chip:#2a3134; }}
  }}
  :root[data-theme="dark"] {{ --ink:#e6e2d8; --mut:#9aa4a0;
      --bg:#14181a; --card:#1d2326; --line:#333a3d; --v:#5fae86;
      --r:#d08159; --acc:#7fb4cf; --notebg:#242a21; --noteln:#8a7c4f;
      --chip:#2a3134; }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--bg); color:var(--ink);
         font:16px/1.55 Georgia,'Times New Roman',serif; }}
  main {{ max-width:60rem; margin:0 auto; padding:2.5rem 1.25rem 4rem; }}
  h1 {{ font-size:1.7rem; line-height:1.25; margin:0 0 .3rem;
       text-wrap:balance; }}
  h2 {{ font-size:1.15rem; margin:2.4rem 0 .6rem; color:var(--acc);
       border-bottom:1px solid var(--line); padding-bottom:.3rem; }}
  p, li {{ max-width:46rem; }}
  .sub {{ color:var(--mut); margin:0 0 1.4rem; }}
  .note {{ background:var(--notebg); border-left:3px solid var(--noteln);
          padding:.7rem 1rem; font-size:.92rem; max-width:46rem; }}
  img {{ max-width:100%; border:1px solid var(--line); border-radius:4px; }}
  .fd {{ display:none; }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) .fl {{ display:none; }}
    :root:not([data-theme="light"]) .fd {{ display:block; }} }}
  :root[data-theme="dark"] .fl {{ display:none; }}
  :root[data-theme="dark"] .fd {{ display:block; }}
  table {{ border-collapse:collapse; font-size:.92rem;
          font-variant-numeric:tabular-nums; }}
  th, td {{ border:1px solid var(--line); padding:.35rem .7rem;
           text-align:right; }}
  th:first-child, td:first-child {{ text-align:left; }}
  code, .mono {{ font:.85rem/1.5 ui-monospace,Consolas,monospace; }}
  .text-sample {{ background:var(--card); border:1px solid var(--line);
    border-radius:4px; padding:.6rem .8rem; overflow-x:auto;
    white-space:nowrap; font:.8rem/1.7 ui-monospace,monospace; }}
  .v {{ color:var(--v); font-weight:bold; }} .r {{ color:var(--r);
    font-weight:bold; }}
  .chip {{ display:inline-block; background:var(--chip); border-radius:9px;
          padding:.06rem .55rem; font-size:.82rem; margin:.1rem; }}
  .verdict {{ background:var(--card); border:1px solid var(--line);
    border-left:4px solid var(--acc); padding:.9rem 1.1rem;
    max-width:46rem; }}
  footer {{ margin-top:3rem; color:var(--mut); font-size:.85rem; }}
  #themebtn {{ position:fixed; top:.8rem; right:.9rem; cursor:pointer;
    background:var(--card); color:var(--ink); border:1px solid var(--line);
    border-radius:50%; width:2.2rem; height:2.2rem; font-size:1.05rem;
    line-height:1; }}
  @media print {{ #themebtn {{ display:none; }}
    :root, :root[data-theme="dark"] {{ --ink:#1c2321; --mut:#5c6663;
      --bg:#ffffff; --card:#ffffff; --line:#e3e0d8; --v:#2e7d54;
      --r:#b3542e; --acc:#22536b; --notebg:#f3efe4; --noteln:#c8b37a;
      --chip:#eee9dc; }}
    .fd {{ display:none !important; }} .fl {{ display:block !important; }} }}
</style></head><body>
<button id="themebtn" title="toggle light/dark" aria-label="toggle theme"
 onclick="(function(){{var h=document.documentElement,
  c=h.getAttribute('data-theme'),
  s=window.matchMedia('(prefers-color-scheme: dark)').matches,
  n=c?(c==='dark'?'light':'dark'):(s?'light':'dark');
  h.setAttribute('data-theme',n);
  try{{localStorage.setItem('csTheme',n);}}catch(e){{}}}})()">&#9681;</button>
<script>try{{var t=localStorage.getItem('csTheme');
if(t)document.documentElement.setAttribute('data-theme',t);}}catch(e){{}}
</script>
<main>

<h1>One virgin and one recycled PET &mdash; four readings of the same two signals</h1>
<p class="sub">sig2dna case study &middot; regenerated by
<code>run_case_study.py</code> from the shipped anonymized data
({n_ch} channels &times; {n_t:,} scans, {dt:.2f}&thinsp;s/scan,
m/z 40&ndash;650) &middot; compute time {elapsed:.0f}&thinsp;s</p>

<p class="note"><b>Scope of this example.</b> The two chromatograms are
real, anonymized measurements (identifiers, timestamps and provenance
removed). The dataset also ships an anonymized
<b>reference-population calibration</b> (per-channel location and spread
of grammar residuals and event counts, from {s['n_ref']} anonymized
virgin references, leave-one-out; mode = <code>reference_population</code>).
No technical-replicate floors exist here: even the calibrated scores are
reference-population grade, not replicate grade. One further lesson is
shipped deliberately: the calibration bundle was built by the same
software stack that generates this report &mdash; marginal letters can
flicker between library versions, so per-channel grammar calibration is
<b>pipeline-local</b>; recomputing on a different stack may shift
individual z values. The numbers illustrate the three representations;
they are not a validated classification.</p>

<h2>0 &middot; The two signals</h2>
{img_pair(figs, 'tic', 'TIC overlay')}
<p>Summed over all channels the two runs look broadly similar &mdash;
the differences that matter live per channel, and each route below reads
them differently.</p>

<h2>1 &middot; Route 1: CWT letters &mdash; shape &amp; organization</h2>
<p>Each channel is wavelet-encoded into letters
(<code>Y</code>&thinsp;rise, <code>A</code>&thinsp;apex,
<code>Z</code>&thinsp;fall, <code>B</code>&thinsp;return,
<code>C</code>/<code>X</code>&thinsp;merged forms,
<code>_</code>&thinsp;non-coding), after Poisson rejection of noise
segments. Amplitude is deliberately forgotten &mdash; the letters read
the <em>morphology</em>.</p>
<p class="mono">m/z {s['demo_mz']} &mdash; most-coded channel of the
virgin run:</p>
<div class="text-sample"><span class="v">V</span>&nbsp; {trim(s['textV'])}</div>
<div class="text-sample"><span class="r">R</span>&nbsp; {trim(s['textR'])}</div>
<table>
<tr><th>letter counts (all channels)</th>{lhead}</tr>
<tr><td class="v">virgin</td>{lrow_v}</tr>
<tr><td class="r">recycled</td>{lrow_r}</tr>
</table>
<p></p>
{img_pair(figs, 'entropy', 'per-channel text entropy')}
<p>Per-channel text entropy is the intensive organization statistic:
where the recycled trace departs, its channel texts are organized
differently at similar coding density.</p>

<h2>2 &middot; Route 2: events &mdash; inventory &amp; amount</h2>
<p>Peaks above each channel's <em>local</em> noise threshold become
events; counting them per channel gives the extensive inventory.
Population comparison uses Laplace surprisal: an event class absent from
the reference weighs &minus;log&#8322;(1/(n+2)) bits.</p>
{img_pair(figs, 'counts', 'per-channel event counts')}
<table>
<tr><th></th><th>event budget B</th><th>&rho;<sub>E</sub> (median
events/channel)</th></tr>
<tr><td class="v">virgin</td><td>{s['BV']:,}</td><td>{s['rhoV']:.0f}</td></tr>
<tr><td class="r">recycled</td><td>{s['BR']:,}</td><td>{s['rhoR']:.0f}</td></tr>
</table>
<p>Composition distance d<sub>comp</sub>(V,&thinsp;R) =
<b>{s['dcomp']:.4f}</b>; directed information of R against the
single-reference V: E<sup>+</sup> = <b>{s['Egain']:.0f} bits</b>
(gained classes), E<sup>&minus;</sup> = <b>{s['Eloss']:.0f} bits</b>
(lost classes). Channels gaining most events: {topg}</p>

<h2>3 &middot; Route 3: grammar &mdash; organization conditioned on
inventory, and why calibration matters</h2>
<p>The chain rule I(E,&thinsp;T) = I(E) + I(T&thinsp;|&thinsp;E) in
practice: an order-2 Markov grammar is trained on the virgin channel
texts (gap runs become punctuation &mdash; retention time is punctuation,
not vocabulary); each channel text gets a code length in bits; the
expectation given that channel's <em>resolved event count</em> is
subtracted (shipped reference fit:
L&#770; = {s['beta'][0]:.1f} + {s['beta'][1]:.2f}&thinsp;N); residuals
are standardized and censored at |z| &gt; 2. <b>How they are
standardized decides what you see</b> &mdash; the same residuals under
three calibrations of increasing evidential quality:</p>
<table>
<tr><th>calibration of z</th><th>P(|z|&gt;2), m/z&lt;200</th>
<th>P(|z|&gt;2), m/z&ge;200</th><th>low/high ratio</th></tr>
{ladder_rows}
</table>
<p>Under the global proxy (R0) the flagged channels are
<b>exclusively low-mass</b> &mdash; a <b>heteroscedasticity
artifact</b>: grammar-residual variance grows with channel event
density, and one global &sigma; preferentially pushes dense low-mass
channels past the gate while starving the sparse high-mass channels.
Density-conditioned &sigma;(N) (R1, fitted on the single shipped
reference) flattens the mass dependence almost completely. The shipped
per-channel reference-population calibration (R2) is the properly
scaled score: moderate exceedance on both sides with a residual
low-mass lean for this particular pair &mdash; a tendency, not a
partition of chemical information by mass. <em>The representation can
be simple; the calibration must not be.</em></p>
{img_grammar}
<table>
<tr><th>R2-calibrated</th><th>I<sub>O</sub><sup>+</sup> (bits)</th>
<th>I<sub>O</sub><sup>&minus;</sup> (bits)</th><th>net</th></tr>
<tr><td class="v">virgin (vs shipped reference population)</td>
<td>{s['ioV'][0]:.0f}</td><td>{s['ioV'][1]:.0f}</td>
<td>{s['ioV'][0] - s['ioV'][1]:+.0f}</td></tr>
<tr><td class="r">recycled</td><td>{s['ioR'][0]:.0f}</td>
<td>{s['ioR'][1]:.0f}</td><td>{s['ioR'][0] - s['ioR'][1]:+.0f}</td></tr>
</table>
<p>Direction of the organization change R vs V: <b>{sign_word}</b>.
Profile distance (level-free) = {s['pdist']:.3f}.
Calibration mode of every number in this section:
<code>reference_population</code> &mdash; never to be read as
replicate-grade significance.</p>

<h3>Channel states: what kind of difference does each channel carry?</h3>
<p>Crossing the inventory deviation (event count vs the shipped
reference statistics) with the R2-calibrated organization deviation
classifies every channel into four readable states:</p>
{img_quad}
<table>
<tr><th>fraction of channels</th><th>m/z &lt; 200</th>
<th>m/z &ge; 200</th></tr>
{quad_rows}
</table>
<p>Inventory change and organization change are largely carried by
<em>different</em> channels &mdash; the two routes localize different
aspects of the same perturbation, which is exactly what the chain-rule
construction intends. Organization-only channels appear at both ends of
the mass axis: a tendency, not a partition of chemical information by
mass.</p>

<h2>4 &middot; From differential events to chemistry &mdash; the two
faces of recycling</h2>
<p>The event ontology extends one step further, from statistics back to
<em>recognizable chemistry</em>:</p>
<p class="mono" style="text-align:center">letters &rarr; fragment-peak
events &rarr; co-eluting fragments &rarr; pseudo-EI spectrum &rarr;
candidate chemistry</p>
<p>The selection is <b>blind-first</b>: the event algebra nominates the
event classes <em>gained</em> by the recycled run (present in R, absent
in V), chain-clusters their apexes into co-eluting families
(|t<sub>m</sub>&nbsp;&minus;&nbsp;t<sup>*</sup>| &le; &tau; =
{tau_s:.1f}&thinsp;s), and reconstructs each family's spectrum from the
<b>raw, local-baseline-corrected ion amplitudes</b> of <em>all</em>
co-eluting events, normalized to base peak = 100. The event machinery
supplies the candidate peak; the original signal supplies the spectral
intensities &mdash; event surprisal bits are never used as intensities.
Only <em>after</em> reconstruction is each spectrum compared with two
frozen candidate-marker dictionaries: <b>&#119967;</b>
(degradation/history: furanic-like 96/109, anhydride-like 148,
dioxolane/acetaldehyde-EG-like 87/45, cyclic-oligomer-like 572) and
<b>&#119966;</b> (contamination/cleanliness: fatty-ester-like 74,
alkylbenzene-like 91, terpene/PAH/alkene-like 121/165/55/56/69,
oxygenate-like 59/72/86/88). D and C are the fractions of a spectrum's
<em>signature</em> intensity (its 10 most intense fragments, standard
EI main-peak practice) carried by each dictionary, with exact unit-mass
matching &mdash; a &plusmn;1 tolerance would conflate chemically
opposite markers such as anhydride 148 (degradation) and
dialkyl-phthalate 149 (plasticizer contamination);
exemplars are chosen only near the poles (D&nbsp;&Gt;&nbsp;C or
C&nbsp;&Gt;&nbsp;D) <em>and</em> with genuine amplitude excess:
X&nbsp;=&nbsp;&Sigma;S<sup>+</sup>/&Sigma;S<sub>R</sub>&nbsp;&ge;&nbsp;0.25,
where a low X unmasks families the event layer flagged as gained but
whose fragments the response-matched virgin also carries &mdash;
apex/class <em>shifts</em> (inventory change without amplitude excess),
not new chemistry: <b>event-class gain &ne; chemical excess</b>. The
route is thus two successive questions &mdash; <b>4a</b>: is there
actual chemical excess? (X answers) &mdash; <b>4b</b>: what kind of
chemistry is that excess compatible with? (the dictionaries answer,
and only for families that passed 4a). In full: gained event &rarr;
co-eluting family &rarr; pseudo-EI spectrum &rarr; excess gate &rarr;
chemical interpretation. The three-spectra display below makes the 4a
distinction visible instead of hiding it.</p>
<p>{s['n_gained_cls']} gained event classes cluster into
{s['n_fam_all']} co-eluting families with &ge;&thinsp;3 fragments;
the {n_fam_rec} strongest were reconstructed. At the poles with genuine
excess (X&nbsp;&ge;&nbsp;0.25): {s['n_deg']} degradation-like and
{s['n_con']} contamination-like famil{'ies' if s['n_con'] != 1 else 'y'};
in addition, the {s['n_unm']} strongest excess families
(X&nbsp;&ge;&nbsp;0.5) whose signatures <em>neither</em> dictionary
spans are shown as their own category &mdash; the dictionaries were
frozen before looking, and the panel does not stretch them after the
fact. Response calibration between the two runs:
&alpha; = {s['alpha']:.3f} (median per-scan TIC
ratio).{himass_note}</p>
{img_pair(figs, 'fammap', 'degradation vs contamination family map')}
<table>
<tr><th>RT (min)</th><th>fragments</th><th>gained</th>
<th>base peak m/z</th><th>D</th><th>C</th><th>X</th>
<th>reading</th></tr>
{fam_rows}
</table>
<p>For each pole exemplar, three spectra: the <b>virgin reference</b>
around the same retention region (response-matched by &alpha;, drawn on
the same scale), the <b>recycled raw pseudo-EI spectrum</b> &mdash;
this is the chemical object &mdash; and the <b>excess spectrum</b>
[S<sub>R</sub>&nbsp;&minus;&nbsp;&alpha;S<sub>V</sub>]<sub>+</sub>, the
explanatory overlay highlighting the fragments responsible for the
difference.</p>
{spec_html}
<p>What this particular pair says, read blind: the one clean
degradation-pole exemplar carries the
dioxolane/acetaldehyde&ndash;ethylene-glycol signature (base peak
m/z&thinsp;45) &mdash; a recycling-<em>history</em> marker;
{('no contamination-pole family survives the excess gate, the '
  'contamination-flagged families being apex/class shifts rather than '
  'new chemistry') if s['n_con'] == 0 else
 (str(s['n_con']) + ' contamination-pole '
  + ('family survives' if s['n_con'] == 1 else 'families survive')
  + ' the excess gate')}. The largest genuine excess lies
<em>beyond</em> the frozen dictionaries: a
<b>dimethyl-benzenedicarboxylate-like aromatic ester</b>
(terephthalate-type expected in PET context; positional isomer
unresolved; 163/194/135/103) &mdash; a recycling/<em>process</em>-history
reading &mdash; <span class="chip">E2 &middot; MSP+RI-library supported
candidate</span>, and a <b>heavy dialkyl-phthalate-like plasticizer
signature</b> (149/167 over an alkyl series)
<span class="chip">E2 &middot; MSP+RI-library supported candidate</span>.
Both were promoted from fragment-pattern candidates by post-hoc library
screening (spectral match <em>plus</em> an independent library
retention region); neither is <em>confirmed</em> &mdash; see the
evidentiality ladder below. Blind-first reporting keeps such surprises
visible instead of forcing them onto the nearest known marker.</p>
<div class="verdict">
<p style="margin-bottom:0"><b>Recycling leaves at least two chemically
distinct memories.</b> One is polymer/process history, illustrated
here by a genuine dioxolane/acetaldehyde&ndash;EG-type excess. The
other is foreign contamination chemistry, which the frozen dictionary
can detect but which is not expressed as a genuine excess in this
decontaminated example. The largest additional spectral information in
this sample lies instead in previously unencoded aromatic-ester and
ester-plasticizer-type families, illustrating why the dictionary must
remain open to unknown chemistry.</p>
</div>
<p class="note"><b>Naming caution &mdash; the evidentiality ladder.</b>
Every chemical label in this panel carries an explicit evidence level:
<b>E0</b>&nbsp;fragment-pattern candidate &middot;
<b>E1</b>&nbsp;external spectral-library supported &middot;
<b>E2</b>&nbsp;spectral + independent-library retention-region
supported &middot; <b>E3</b>&nbsp;measured retention index and/or
authentic-standard confirmation. The dictionary labels above
(&laquo;&thinsp;furanic-like&thinsp;&raquo;,
&laquo;&thinsp;fatty-ester-like&thinsp;&raquo;&hellip;) are E0
family-level candidate EI assignments; the two promoted families are
E2 &mdash; <em>supported, not confirmed</em>: their retention evidence
is a library retention <em>region</em>, not a retention index measured
for these peaks, and E3 requires exactly that (or an authentic
standard). A standing warning from this screen: a very high unit-mass
EI match score is <em>not</em> molecular identity &mdash; structurally
unrelated molecules can share the same nominal mass and the same
nominal neutral loss. A contamination spectrum reads
roughly as <em>foreign chemistry was added</em>; a degradation spectrum
as <em>the polymer/process generated or transformed chemistry</em>
&mdash; and the event algebra reminds us that a transformation
A&thinsp;&rarr;&thinsp;B writes both signs at once
(E<sup>+</sup>&thinsp;&gt;&thinsp;0 and
E<sup>&minus;</sup>&thinsp;&gt;&thinsp;0), which is why the
two-dimensional (E<sup>+</sup>,&thinsp;E<sup>&minus;</sup>) reading is
mechanistically richer than any single scalar.</p>

<h2>5 &middot; Four readings, one verdict structure</h2>
<div class="verdict">
<p><b>Letters</b> say how the chemistry is <em>arranged</em> (shape,
history). <b>Events</b> say <em>what and how much</em> is there
(inventory, in additive bits). <b>Grammar</b> says what remains unusual
about the arrangement <em>once the inventory is known</em> &mdash;
signed: perturbations may add texture (I<sub>O</sub><sup>+</sup>) or
erase it (I<sub>O</sub><sup>&minus;</sup>). <b>Reconstructed spectra</b>
say whether the gained chemistry looks like degradation or
contamination &mdash; the two faces of recycling.</p>
<p style="margin-bottom:0">No single scalar decides; a defensible
verdict is structured &mdash; history &times; inventory departure
&times; organization &times; chemistry &times; validity domain &mdash;
and every claim carries the scope it was measured in.</p>
</div>

<footer>Generated deterministically by <code>run_case_study.py</code>
(sig2dna, MIT). Data: anonymized full-scan GC-MS, shipped as
<code>case_vr.npz</code>. No language model was involved in producing
any number in this report.</footer>
</main></body></html>"""
    out = os.path.join(HERE, "report.html")
    with open(out, "w") as f:
        f.write(html)
    print(f"written {out} ({os.path.getsize(out) / 1e6:.1f} MB)")


def write_markdown(figs, s, n_ch, n_t, dt):
    """GitHub-rendered companion of report.html: same generator, same
    numbers; figures are committed dual-theme SVGs under figures/."""
    fdir = os.path.join(HERE, "figures")
    os.makedirs(fdir, exist_ok=True)
    for (name, theme), data in figs.items():
        with open(os.path.join(fdir, f"{name}_{theme}.svg"), "wb") as fh:
            fh.write(data)

    lv, lr = s["lettersV"], s["lettersR"]
    keys = [k for k in "YAZBCX_" if k in lv or k in lr]
    lhead = " | ".join(f"`{k}`" for k in keys)
    lsep = " | ".join("---:" for _ in keys)
    lrow_v = " | ".join(f"{lv.get(k, 0):,}" for k in keys)
    lrow_r = " | ".join(f"{lr.get(k, 0):,}" for k in keys)
    topg = " · ".join(f"m/z {m} (+{g})" for m, g in s["topgain"])
    sign_word = {1: "net enrichment (I_O⁺ dominated)",
                 -1: "net regularization (I_O⁻ dominated)",
                 0: "balanced"}[s["sign"]]

    def _ratio(lo, hi):
        if hi == 0:
            return ("∞ (no channel above the gate at m/z ≥ 200)"
                    if lo > 0 else "—")
        return f"{lo / hi:.2f}"
    ladder_rows = "\n".join(
        f"| `{name}` | {lo:.3f} | {hi:.3f} | **{_ratio(lo, hi)}** |"
        for name, (lo, hi) in s["ladder"].items())
    quad_rows = "\n".join(f"| {k} | {a:.3f} | {b:.3f} |"
                          for k, (a, b) in s["quad"].items())
    fam_rows = "\n".join(
        f"| {f['rt']:.1f} | {f['n_frag']} | {f['n_gained']} "
        f"| {f['base_mz']} | {f['D']:.2f} | {f['C']:.2f} | {f['X']:.2f} "
        f"| {f['reading']}"
        + (f" — {f['label']}" if f["label"] else "")
        + ("" if f["X"] >= 0.25 else " *(class shift)*") + " |"
        for f in s["families"])
    cat_head = {"degradation-like": "Degradation face",
                "contamination-like": "Contamination face",
                "mixed / unresolved": "Beyond the frozen dictionaries"}
    spec_md = "\n\n".join(
        f"#### {cat_head[f['reading']]} — "
        + (f["label"] if f["label"]
           and f["reading"] != "mixed / unresolved"
           else "signature outside both dictionaries")
        + f" (RT {f['rt']:.1f} min, base peak m/z {f['base_mz']})\n\n"
        + md_pic(f"spec{k}", "three-spectra pseudo-EI reconstruction")
        for k, f in enumerate(s["exemplars"]))
    himass = ("" if s["has_himass"] else
              " No credible high-mass family (m/z > 500, e.g. a "
              "cyclic-oligomer-type marker) was nominated by the events "
              "of this particular pair; the panel reports what the data "
              "support and does not manufacture known cohort markers.")
    con_clause = (
        "no contamination-pole family survives the excess gate, the "
        "contamination-flagged families being apex/class shifts rather "
        "than new chemistry" if s["n_con"] == 0 else
        f"{s['n_con']} contamination-pole "
        + ("family survives" if s["n_con"] == 1 else "families survive")
        + " the excess gate")

    def trim(t, n=90):
        t = t.strip("_")
        return (t[:n] + "…") if len(t) > n else (t or "(empty)")

    md = f"""# One virgin and one recycled PET — four readings of the same two signals

*sig2dna case study · GitHub-rendered companion of
[`report.html`](report.html) — same generator
([`run_case_study.py`](run_case_study.py)), same numbers, regenerated
from the shipped anonymized data ({n_ch} channels × {n_t:,} scans,
{dt:.2f} s/scan, m/z 40–650). Figures are committed dual-theme SVGs.*

> **Scope of this example.** The two chromatograms are real, anonymized
> measurements (identifiers, timestamps and provenance removed). The
> dataset also ships an anonymized **reference-population calibration**
> (per-channel location and spread of grammar residuals and event
> counts, from {s['n_ref']} anonymized virgin references,
> leave-one-out; mode = `reference_population`). No technical-replicate
> floors exist here. The calibration bundle was built by the same
> software stack that generates this report — per-channel grammar
> calibration is **pipeline-local**; recomputing on a different stack
> may shift individual z values. The numbers illustrate the four
> representations; they are not a validated classification.

## 0 · The two signals

{md_pic('tic', 'TIC overlay')}

Summed over all channels the two runs look broadly similar — the
differences that matter live per channel, and each route below reads
them differently.

## 1 · Route 1: CWT letters — shape & organization

Each channel is wavelet-encoded into letters (`Y` rise, `A` apex,
`Z` fall, `B` return, `C`/`X` merged forms, `_` non-coding) after
Poisson rejection of noise segments; amplitude is deliberately
forgotten — the letters read the *morphology*. m/z {s['demo_mz']},
most-coded channel of the virgin run:

```text
V  {trim(s['textV'])}
R  {trim(s['textR'])}
```

| letter counts (all channels) | {lhead} |
|---|{lsep}|
| virgin | {lrow_v} |
| recycled | {lrow_r} |

{md_pic('entropy', 'per-channel text entropy')}

Per-channel text entropy is the intensive organization statistic:
where the recycled trace departs, its channel texts are organized
differently at similar coding density.

## 2 · Route 2: events — inventory & amount

Peaks above each channel's *local* noise threshold become events;
counting them per channel gives the extensive inventory. Population
comparison uses Laplace surprisal: an event class absent from the
reference weighs −log₂(1/(n+2)) bits.

{md_pic('counts', 'per-channel event counts')}

| | event budget B | ρ_E (median events/channel) |
|---|---:|---:|
| virgin | {s['BV']:,} | {s['rhoV']:.0f} |
| recycled | {s['BR']:,} | {s['rhoR']:.0f} |

Composition distance d_comp(V, R) = **{s['dcomp']:.4f}**; directed
information of R against the single-reference V:
E⁺ = **{s['Egain']:.0f} bits** (gained classes),
E⁻ = **{s['Eloss']:.0f} bits** (lost classes).
Channels gaining most events: {topg}.

## 3 · Route 3: grammar — organization conditioned on inventory, and why calibration matters

The chain rule I(E, T) = I(E) + I(T | E) in practice: an order-2
Markov grammar trained on the virgin channel texts (gap runs become
punctuation — retention time is punctuation, not vocabulary) assigns
each channel text a code length in bits; the expectation given that
channel's *resolved event count* is subtracted (shipped reference fit:
L̂ = {s['beta'][0]:.1f} + {s['beta'][1]:.2f} N); residuals are
standardized and censored at |z| > 2. **How they are standardized
decides what you see** — the same residuals under three calibrations
of increasing evidential quality:

| calibration of z | P(\\|z\\|>2), m/z<200 | P(\\|z\\|>2), m/z≥200 | low/high ratio |
|---|---:|---:|---:|
{ladder_rows}

Under the global proxy (R0) the flagged channels are **exclusively
low-mass** — a **heteroscedasticity artifact**: grammar-residual
variance grows with channel event density, and one global σ
preferentially pushes dense low-mass channels past the gate while
starving the sparse high-mass channels. Density-conditioned σ(N) (R1)
flattens the mass dependence almost completely; the shipped
per-channel reference-population calibration (R2) is the properly
scaled score. *The representation can be simple; the calibration must
not be.*

{md_pic('grammar', 'organization residuals, R0 vs R2')}

| R2-calibrated | I_O⁺ (bits) | I_O⁻ (bits) | net |
|---|---:|---:|---:|
| virgin (vs shipped reference population) | {s['ioV'][0]:.0f} \
| {s['ioV'][1]:.0f} | {s['ioV'][0] - s['ioV'][1]:+.0f} |
| recycled | {s['ioR'][0]:.0f} | {s['ioR'][1]:.0f} | {s['ioR'][0] - s['ioR'][1]:+.0f} |

Direction of the organization change R vs V: **{sign_word}**. Profile
distance (level-free) = {s['pdist']:.3f}. Calibration mode of every
number in this section: `reference_population` — never to be read as
replicate-grade significance.

### Channel states: what kind of difference does each channel carry?

{md_pic('quad', 'channel-state quadrants')}

| fraction of channels | m/z < 200 | m/z ≥ 200 |
|---|---:|---:|
{quad_rows}

Inventory change and organization change are largely carried by
*different* channels — the two routes localize different aspects of
the same perturbation, which is exactly what the chain-rule
construction intends.

## 4 · From differential events to chemistry — the two faces of recycling

letters → fragment-peak events → co-eluting fragments → pseudo-EI
spectrum → candidate chemistry

The selection is **blind-first**: the event algebra nominates the
event classes *gained* by the recycled run, chain-clusters their
apexes into co-eluting families (\\|t_m − t\\*\\| ≤ τ =
{TAU * dt:.1f} s), and reconstructs each family's spectrum from the
**raw, local-baseline-corrected ion amplitudes** of *all* co-eluting
events, normalized to base peak = 100 — event surprisal bits are never
used as intensities. Only *after* reconstruction is each spectrum
compared with two frozen candidate-marker dictionaries: **𝒟**
(degradation/history: furanic-like 96/109, anhydride-like 148,
dioxolane/acetaldehyde-EG-like 87/45, cyclic-oligomer-like 572) and
**𝒞** (contamination/cleanliness: fatty-ester-like 74,
alkylbenzene-like 91, terpene/PAH/alkene-like 121/165/55/56/69,
oxygenate-like 59/72/86/88). D and C are the fractions of a spectrum's
*signature* intensity (10 most intense fragments, exact unit-mass
matching — a ±1 tolerance would conflate anhydride 148 with
dialkyl-phthalate 149). Exemplars need a pole (D ≫ C or C ≫ D) *and*
genuine amplitude excess X = ΣS⁺/ΣS_R ≥ 0.25:
**event-class gain ≠ chemical excess** — the route asks **4a** is
there actual chemical excess? then **4b** what chemistry is it
compatible with?

{s['n_gained_cls']} gained event classes cluster into
{s['n_fam_all']} co-eluting families with ≥ 3 fragments; the
{len(s['families'])} strongest were reconstructed
({s['n_deg']} degradation-like, {s['n_con']} contamination-like at the
poles with X ≥ 0.25; the {s['n_unm']} strongest excess families whose
signatures *neither* dictionary spans are shown as their own
category). Response calibration between the two runs:
α = {s['alpha']:.3f} (median per-scan TIC ratio).{himass}

{md_pic('fammap', 'degradation vs contamination family map')}

| RT (min) | fragments | gained | base peak m/z | D | C | X | reading |
|---:|---:|---:|---:|---:|---:|---:|---|
{fam_rows}

For each exemplar, three spectra: the **virgin reference** around the
same retention region (response-matched by α, same scale), the
**recycled raw pseudo-EI spectrum** — the chemical object — and the
**excess spectrum** [S_R − α·S_V]₊, the explanatory overlay.

{spec_md}

What this particular pair says, read blind: the one clean
degradation-pole exemplar carries the dioxolane/acetaldehyde–
ethylene-glycol signature (base peak m/z 45) — a recycling-*history*
marker; {con_clause}. The largest genuine excess lies *beyond* the
frozen dictionaries: a **dimethyl-benzenedicarboxylate-like aromatic
ester** (terephthalate-type expected in PET context; positional isomer
unresolved; 163/194/135/103) — a recycling/*process*-history reading —
`E2 · MSP+RI-library supported candidate`, and a **heavy
dialkyl-phthalate-like plasticizer signature** (149/167 over an alkyl
series) `E2 · MSP+RI-library supported candidate`. Neither is
*confirmed* — see the evidentiality ladder below.

> **Recycling leaves at least two chemically distinct memories.** One
> is polymer/process history, illustrated here by a genuine
> dioxolane/acetaldehyde–EG-type excess. The other is foreign
> contamination chemistry, which the frozen dictionary can detect but
> which is not expressed as a genuine excess in this decontaminated
> example. The largest additional spectral information in this sample
> lies instead in previously unencoded aromatic-ester and
> ester-plasticizer-type families, illustrating why the dictionary
> must remain open to unknown chemistry.

> **Naming caution — the evidentiality ladder.** **E0** fragment-
> pattern candidate · **E1** external spectral-library supported ·
> **E2** spectral + independent-library retention-region supported ·
> **E3** measured retention index and/or authentic-standard
> confirmation. Dictionary labels above are E0; the two promoted
> families are E2 — *supported, not confirmed*: their retention
> evidence is a library retention *region*, not a retention index
> measured for these peaks. A standing warning: a very high unit-mass
> EI match score is *not* molecular identity — structurally unrelated
> molecules can share the same nominal mass and the same nominal
> neutral loss.

## 5 · Four readings, one verdict structure

**Letters** say how the chemistry is *arranged* (shape, history).
**Events** say *what and how much* is there (inventory, in additive
bits). **Grammar** says what remains unusual about the arrangement
*once the inventory is known* — signed: perturbations may add texture
(I_O⁺) or erase it (I_O⁻). **Reconstructed spectra** say whether the
gained chemistry looks like degradation or contamination — the two
faces of recycling. No single scalar decides; a defensible verdict is
structured — history × inventory departure × organization × chemistry
× validity domain — and every claim carries the scope it was measured
in.

---

*Generated deterministically by `run_case_study.py` (sig2dna, MIT).
Data: anonymized full-scan GC-MS, shipped as `case_vr.npz`. No
language model was involved in producing any number in this report.*
"""
    out = os.path.join(HERE, "report.md")
    with open(out, "w") as f:
        f.write(md)
    n_svg = len(figs)
    tot = sum(os.path.getsize(os.path.join(fdir, p))
              for p in os.listdir(fdir))
    print(f"written {out} ({os.path.getsize(out) / 1e3:.0f} kB) + "
          f"{n_svg} SVGs ({tot / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
