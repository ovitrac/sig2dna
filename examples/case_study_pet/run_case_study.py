#!/usr/bin/env python3
"""
Case study — one virgin and one recycled PET, three readings of the same
two signals: CWT letters (shape/organization), event/channel algebra
(inventory/amount), and the organization/grammar algebra connecting them.

Fully replayable: everything in report.html is recomputed from the
shipped anonymized data (case_vr.npz: two full-scan GC-MS runs,
n_channels x n_time integer intensities; 0.14 s/scan; m/z 40-650;
all identifiers, timestamps and provenance removed). The report ships
with light and dark themes (system default, manual toggle, print=light).

DEMO-GRADE CALIBRATION NOTICE: with a single reference run there are no
replicate floors; the grammar gates here use a global residual-spread
proxy and the virgin run is scored in-sample as a visual baseline. A
real deployment calibrates per channel on true replicates and states
scopes explicitly (see module docstrings).

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
from sig2dna_core.grammar import (MarkovGrammar, fit_conditioning,
                                  organization_scores, profile_distance,
                                  rho_e, sign_delta, tokenize)
from sig2dna_core.icfilter import encode_ic_matrix, text_entropy

HERE = os.path.dirname(os.path.abspath(__file__))
TAU = 36.0
SCALE = 16.0
GATE = 2.0

THEMES = {
    "light": dict(bg="#ffffff", ink="#1c2321", mut="#5c6663",
                  grid="#e3e0d8", v="#2e7d54", r="#b3542e",
                  band="#00000012"),
    "dark": dict(bg="#1d2326", ink="#e6e2d8", mut="#9aa4a0",
                 grid="#333a3d", v="#5fae86", r="#d08159",
                 band="#ffffff14"),
}


def fig64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


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

    # ---------------- route 3: grammar -----------------------------------
    print("route 3: grammar...", flush=True)
    toksV = [tokenize(s) for s in encV.letters]
    toksR = [tokenize(s) for s in encR.letters]
    G = MarkovGrammar().train(toksV)
    LV = np.array([G.codelen(s) for s in toksV])
    LR = np.array([G.codelen(s) for s in toksR])
    beta, res = fit_conditioning(LV, nV)
    stats["beta"] = beta
    stats["R2"] = 1 - res.var() / LV.var()
    # demo-grade calibration: global residual spread (see notice above)
    mu = np.zeros(n_ch)
    sd = np.full(n_ch, max(np.std(res, ddof=1), 1.0))
    zV, pV, nnV = organization_scores(LV, nV, beta, mu, sd, gate=GATE)
    zR, pR, nnR = organization_scores(LR, nR, beta, mu, sd, gate=GATE)
    stats["ioV"], stats["ioR"] = (pV, nnV), (pR, nnR)
    stats["sign"] = sign_delta((pR, nnR), (pV, nnV))
    stats["pdist"] = profile_distance(zV, zR)

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
        figs[("tic", theme)] = fig64(fig)

        fig, ax, t = themed_axes(theme)
        ax.plot(mz, entV, color=t["v"], lw=0.8, label="virgin")
        ax.plot(mz, entR, color=t["r"], lw=0.8, alpha=0.85,
                label="recycled")
        ax.set_xlabel("m/z channel")
        ax.set_ylabel("text entropy (bits)")
        legend(ax, t)
        figs[("entropy", theme)] = fig64(fig)

        fig, ax, t = themed_axes(theme)
        ax.plot(mz, nV, color=t["v"], lw=0.8, label="virgin")
        ax.plot(mz, nR, color=t["r"], lw=0.8, alpha=0.85,
                label="recycled")
        ax.set_xlabel("m/z channel")
        ax.set_ylabel("resolved events per channel")
        legend(ax, t)
        figs[("counts", theme)] = fig64(fig)

        fig, ax, t = themed_axes(theme)
        ax.axhspan(-GATE, GATE, color=t["band"], lw=0)
        ax.plot(mz, zR, color=t["r"], lw=0.8,
                label="recycled, z vs virgin grammar")
        ax.axhline(0, color=t["mut"], lw=0.6)
        ax.set_xlabel("m/z channel")
        ax.set_ylabel("organization residual z")
        legend(ax, t)
        figs[("grammar", theme)] = fig64(fig)

    write_report(figs, stats, n_ch, n_t, dt, time.time() - t0)


def img_pair(figs, name, alt):
    return (f'<img class="fl" src="data:image/png;base64,'
            f'{figs[(name, "light")]}" alt="{alt}">'
            f'<img class="fd" src="data:image/png;base64,'
            f'{figs[(name, "dark")]}" alt="{alt}">')


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

    def trim(t, n=90):
        t = t.strip("_")
        return (t[:n] + "&hellip;") if len(t) > n else (t or "(empty)")

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>sig2dna case study — one virgin and one recycled PET, three readings</title>
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

<h1>One virgin and one recycled PET &mdash; three readings of the same two signals</h1>
<p class="sub">sig2dna case study &middot; regenerated by
<code>run_case_study.py</code> from the shipped anonymized data
({n_ch} channels &times; {n_t:,} scans, {dt:.2f}&thinsp;s/scan,
m/z 40&ndash;650) &middot; compute time {elapsed:.0f}&thinsp;s</p>

<p class="note"><b>Scope of this example.</b> The two chromatograms are
real, anonymized measurements (identifiers, timestamps and provenance
removed). With a single reference run there are no replicate floors:
grammar gates use a global residual-spread proxy and the virgin run is
scored in-sample as a visual baseline. The numbers illustrate the three
representations; they are not a validated classification. Real
deployments calibrate per channel on true replicates and state validity
scopes explicitly (see the module docstrings).</p>

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
inventory</h2>
<p>The chain rule I(E,&thinsp;T) = I(E) + I(T&thinsp;|&thinsp;E) in
practice: an order-2 Markov grammar is trained on the virgin channel
texts (gap runs become punctuation &mdash; retention time is punctuation,
not vocabulary); each recycled channel text gets a code length in bits;
the expectation given that channel's <em>resolved event count</em>
(fit: L&#770; = {s['beta'][0]:.1f} + {s['beta'][1]:.2f}&thinsp;N,
R&sup2; = {s['R2']:.2f}) is subtracted; residuals are standardized and
censored at |z| &gt; 2.</p>
{img_pair(figs, 'grammar', 'organization residuals')}
<table>
<tr><th></th><th>I<sub>O</sub><sup>+</sup> (bits)</th>
<th>I<sub>O</sub><sup>&minus;</sup> (bits)</th><th>net</th></tr>
<tr><td class="v">virgin (in-sample baseline)</td>
<td>{s['ioV'][0]:.0f}</td><td>{s['ioV'][1]:.0f}</td>
<td>{s['ioV'][0] - s['ioV'][1]:+.0f}</td></tr>
<tr><td class="r">recycled</td><td>{s['ioR'][0]:.0f}</td>
<td>{s['ioR'][1]:.0f}</td><td>{s['ioR'][0] - s['ioR'][1]:+.0f}</td></tr>
</table>
<p>Direction of the organization change R vs V: <b>{sign_word}</b>.
Profile distance (level-free) = {s['pdist']:.3f}.</p>

<h2>4 &middot; Three readings, one verdict structure</h2>
<div class="verdict">
<p><b>Letters</b> say how the chemistry is <em>arranged</em> (shape,
history). <b>Events</b> say <em>what and how much</em> is there
(inventory, in additive bits). <b>Grammar</b> says what remains unusual
about the arrangement <em>once the inventory is known</em> &mdash;
signed: perturbations may add texture (I<sub>O</sub><sup>+</sup>) or
erase it (I<sub>O</sub><sup>&minus;</sup>).</p>
<p style="margin-bottom:0">No single scalar decides; a defensible
verdict is structured &mdash; history &times; inventory departure
&times; organization &times; validity domain &mdash; and every claim
carries the scope it was measured in.</p>
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


if __name__ == "__main__":
    main()
