"""
plausibility_scaling.py
=======================
Final biological-plausibility stress test of the surviving titration mechanism.

Replaces the arbitrary universal S_site=0.00634 nM with explicit copy-number
accounting: for a diploid locus the concentration is derived per nuclear volume
from Avogadro's number, and active Cas9 is specified as a molecule count (also
converted to nM per volume). We then ask: does the transient on-target kinetic
delay survive when competitor burden AND Cas9 abundance are constrained to
physically plausible molecule counts?

Kinetics are unchanged (same epsilon / forward rates). DEFENSIBLE composition
only. No panel tuning.

Physics note: total site conc = M*S_site = 2M/(NA*V); Cas9 conc = Nc/(NA*V).
Their RATIO = 2M/Nc is volume-INDEPENDENT; volume only rescales absolute time.
So the phase boundary should track the molecule-count ratio 2M/Nc.
"""

import os, json, csv, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import go_no_go as g
import audit_panels as ap   # for build_panel_det (deterministic DEFENSIBLE)

try:
    import cupy as cp
    _HAVE_CUPY = True
except Exception:
    cp = None; _HAVE_CUPY = False

NA = 6.02214076e23
OUT = os.path.join(g.BASE_DIR, "artifacts", "plausibility")
os.makedirs(OUT, exist_ok=True)
EPS, BASE_FWD = g.load_parameters()

VOLUMES_FL = [250, 500, 1000]
NC_LIST = [10, 30, 100, 300, 1000, 3000, 10000, 30000]     # active Cas9 molecules
M_LIST = [10, 30, 100, 300, 1000, 3000, 10000, 30000]      # occ-competent loci


def s_site(V_fL):
    return 2.0 / (NA * V_fL * 1e-15) * 1e9         # nM, diploid locus


def cas_conc(Nc, V_fL):
    return Nc / (NA * V_fL * 1e-15) * 1e9          # nM


# ---------------------------------------------------------------------------
# Fast competitive solver: warm-started Illinois root find for C_free(t).
# g(cf) = cf - (C_total - bound(cf)) is monotone increasing -> unique root.
# ---------------------------------------------------------------------------
def _bound_of(xp, P, cf, fwd, bck, h, theta):
    Pt = g._cn_step(xp, P, cf, fwd, bck, h, theta)
    return Pt, float(Pt[:, 1:].sum())


def _solve_cfree(xp, P, C_total, fwd, bck, h, cf_warm, theta=0.5,
                 tol=1e-7, maxit=20):
    def gf(cf):
        Pt, b = _bound_of(xp, P, cf, fwd, bck, h, theta)
        return cf - (C_total - b), Pt
    a = 0.0; fa, _ = gf(a)
    if fa >= 0.0:
        _, Pt = gf(0.0); return Pt, 0.0
    b = C_total; fb, Ptb = gf(b)
    if fb <= 0.0:
        return Ptb, C_total
    # warm-start bracket split using previous cf
    if 0.0 < cf_warm < C_total:
        fw, Ptw = gf(cf_warm)
        if fw == 0.0:
            return Ptw, cf_warm
        if fw > 0.0:
            b, fb, Ptb = cf_warm, fw, Ptw
        else:
            a, fa = cf_warm, fw
    Ptc = Ptb; c = b
    for _ in range(maxit):
        c = (a * fb - b * fa) / (fb - fa)          # regula falsi
        fc, Ptc = gf(c)
        if abs(fc) <= tol * C_total or abs(b - a) <= tol * C_total:
            return Ptc, c
        if fc * fb < 0.0:
            a, fa = b, fb
            b, fb = c, fc
        else:
            fa *= 0.5                               # Illinois down-weight
            b, fb = c, fc
    return Ptc, c


def simulate_fast(C_total, mms, t, S_site, use_gpu=True, theta=0.5):
    xp = cp if (use_gpu and _HAVE_CUPY) else np
    M = len(mms)
    fwd = np.zeros((M, g.NSTATES)); bck = np.zeros((M, g.NSTATES))
    for i, mm in enumerate(mms):
        fwd[i], bck[i] = g.get_rate_pair(EPS, BASE_FWD, mm)
    fwd = xp.asarray(fwd); bck = xp.asarray(bck)
    P = xp.zeros((M, g.NSTATES)); P[:, 0] = S_site
    T = len(t)
    Cf = np.zeros(T); Cf[0] = C_total
    on = np.zeros(T)                    # on-target cleaved fraction (site 0)
    cf = C_total
    on_gpu = (xp is cp)
    for k in range(T - 1):
        h = t[k + 1] - t[k]
        P, cf = _solve_cfree(xp, P, C_total, fwd, bck, h, cf, theta)
        P = xp.maximum(P, 0.0)
        Cf[k + 1] = cf
        s0 = float(P[0].sum())
        on[k + 1] = (S_site - s0) / S_site
    return dict(Cf=Cf, on=on)


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------
def t_reach(t, curve, frac):
    """Absolute time to reach `frac` cleaved fraction (interp). NaN if never."""
    if curve[-1] < frac:
        return np.nan
    idx = np.searchsorted(curve, frac)
    if idx <= 0:
        return t[0]
    t0, t1, c0, c1 = t[idx - 1], t[idx], curve[idx - 1], curve[idx]
    return t1 if c1 == c0 else t0 + (frac - c0) * (t1 - t0) / (c1 - c0)


def duration_above(t, series, thresh):
    """Total time (s) for which series > thresh, via midpoint rule on segments."""
    dur = 0.0
    for k in range(len(t) - 1):
        a, b = series[k] > thresh, series[k + 1] > thresh
        if a and b:
            dur += t[k + 1] - t[k]
        elif a != b:
            dur += 0.5 * (t[k + 1] - t[k])
    return dur


# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
def run_grid(t):
    rows = []
    # cache Nature on-target per unique C_total (independent of M, S_site)
    nat_cache = {}
    panel_cache = {}
    total = len(VOLUMES_FL) * len(NC_LIST) * len(M_LIST)
    n = 0
    t0 = time.time()
    for V in VOLUMES_FL:
        S = s_site(V)
        for Nc in NC_LIST:
            C = cas_conc(Nc, V)
            if C not in nat_cache:
                on_A = g.nature_cleavage(C, [], t, EPS, BASE_FWD)
                nat_cache[C] = on_A
            on_A = nat_cache[C]
            tA = {p: t_reach(t, on_A, p) for p in (0.1, 0.5, 0.9)}
            for M in M_LIST:
                if M not in panel_cache:
                    mms, _ = ap.build_panel_det("DEFENSIBLE", M)
                    panel_cache[M] = mms
                mms = panel_cache[M]
                res = simulate_fast(C, mms, t, S, use_gpu=True)
                on_B, Cf = res["on"], res["Cf"]
                depl_series = 1.0 - Cf / C
                d = (on_B - on_A) * 100.0
                d_max = float(np.max(np.abs(d)))
                d_final = float(d[-1])
                tB = {p: t_reach(t, on_B, p) for p in (0.1, 0.5, 0.9)}
                def dl(p):
                    a, b = tA[p], tB[p]
                    if a and np.isfinite(a) and a > 0 and np.isfinite(b):
                        return float((b - a) / a * 100.0)
                    return np.nan
                depl_max = float(np.max(depl_series) * 100.0)
                dur5 = duration_above(t, depl_series, 0.05)
                ratio = 2.0 * M / Nc
                rows.append(dict(
                    V_fL=V, S_site_nM=S, Nc_molecules=Nc, C_total_nM=C, M_loci=M,
                    site_copies=2 * M, total_site_nM=M * S, ratio_copies_per_cas9=ratio,
                    d_on_max_pp=d_max, d_on_final_pp=d_final,
                    t10_delay_pct=dl(0.1), t50_delay_pct=dl(0.5), t90_delay_pct=dl(0.9),
                    depletion_max_pct=depl_max, dur_gt5pct_depl_s=dur5,
                    t50_A_s=float(tA[0.5]) if np.isfinite(tA[0.5]) else None,
                    t50_B_s=float(tB[0.5]) if np.isfinite(tB[0.5]) else None,
                    on_B_final=float(on_B[-1]), on_A_final=float(on_A[-1]),
                ))
                n += 1
                if n % 8 == 0 or n == total:
                    print(f"  {n:3d}/{total} V={V} Nc={Nc} M={M} "
                          f"dOnMax={d_max:6.2f}pp t50d={dl(0.5) if np.isfinite(dl(0.5)) else float('nan'):7.1f}% "
                          f"depl={depl_max:5.1f}% ratio={ratio:.1f} ({time.time()-t0:.0f}s)",
                          flush=True)
    return rows


# ---------------------------------------------------------------------------
# Phase maps
# ---------------------------------------------------------------------------
def phase_maps(rows):
    def grid_of(V, field):
        Z = np.full((len(NC_LIST), len(M_LIST)), np.nan)
        for r in rows:
            if r["V_fL"] != V:
                continue
            i = NC_LIST.index(r["Nc_molecules"]); j = M_LIST.index(r["M_loci"])
            Z[i, j] = r[field]
        return Z

    # transient dOn phase map (binned)
    bounds = [0, 1, 5, 10, 1e9]
    cmap = ListedColormap(["#2c7fb8", "#7fcdbb", "#fec44f", "#d95f0e"])
    norm = BoundaryNorm(bounds, cmap.N)
    for V in VOLUMES_FL:
        Z = grid_of(V, "d_on_max_pp")
        fig, ax = plt.subplots(figsize=(7, 5.5))
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap=cmap, norm=norm)
        _annot(ax, Z, fmt="{:.1f}")
        ax.set_xticks(range(len(M_LIST))); ax.set_xticklabels(M_LIST, rotation=45)
        ax.set_yticks(range(len(NC_LIST))); ax.set_yticklabels(NC_LIST)
        ax.set_xlabel("effective competitor loci  M"); ax.set_ylabel("active Cas9 molecules")
        ax.set_title(f"Transient on-target |Δ| (pp)  —  V={V} fL\n"
                     "bins: <1 / 1-5 / 5-10 / >10 pp")
        cb = fig.colorbar(im, ax=ax, ticks=[0.5, 3, 7.5, 50])
        cb.ax.set_yticklabels(["<1", "1-5", "5-10", ">10"])
        plt.tight_layout(); plt.savefig(os.path.join(OUT, f"phase_dOn_V{V}.png"), dpi=150)
        plt.close()

    # t50 delay maps (>10%, >25%)
    tb = [-1e9, 10, 25, 1e9]
    tcmap = ListedColormap(["#edf8b1", "#7fcdbb", "#2c7fb8"])
    tnorm = BoundaryNorm(tb, tcmap.N)
    for V in VOLUMES_FL:
        Z = grid_of(V, "t50_delay_pct")
        Zp = np.where(np.isfinite(Z), Z, 1e8)   # censored (never reached) -> huge
        fig, ax = plt.subplots(figsize=(7, 5.5))
        im = ax.imshow(Zp, origin="lower", aspect="auto", cmap=tcmap, norm=tnorm)
        _annot(ax, Z, fmt="{:.0f}")
        ax.set_xticks(range(len(M_LIST))); ax.set_xticklabels(M_LIST, rotation=45)
        ax.set_yticks(range(len(NC_LIST))); ax.set_yticklabels(NC_LIST)
        ax.set_xlabel("effective competitor loci  M"); ax.set_ylabel("active Cas9 molecules")
        ax.set_title(f"t50 delay (%)  —  V={V} fL\nbins: <10 / 10-25 / >25 %")
        cb = fig.colorbar(im, ax=ax, ticks=[5, 17, 60])
        cb.ax.set_yticklabels(["<10", "10-25", ">25"])
        plt.tight_layout(); plt.savefig(os.path.join(OUT, f"phase_t50delay_V{V}.png"), dpi=150)
        plt.close()


def _annot(ax, Z, fmt):
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if np.isfinite(Z[i, j]):
                ax.text(j, i, fmt.format(Z[i, j]), ha="center", va="center",
                        fontsize=6.5, color="black")


# ---------------------------------------------------------------------------
# Minimum competitor thresholds
# ---------------------------------------------------------------------------
def min_competitor_thresholds(rows):
    out = []
    for V in VOLUMES_FL:
        for Nc in NC_LIST:
            sub = sorted([r for r in rows if r["V_fL"] == V and r["Nc_molecules"] == Nc],
                         key=lambda r: r["M_loci"])
            def min_M(pred):
                for r in sub:
                    if pred(r):
                        return r["M_loci"]
                return None
            out.append(dict(
                V_fL=V, Nc_molecules=Nc, C_total_nM=cas_conc(Nc, V),
                min_M_1pp=min_M(lambda r: r["d_on_max_pp"] >= 1.0),
                min_M_5pp=min_M(lambda r: r["d_on_max_pp"] >= 5.0),
                min_M_t50_10pct=min_M(lambda r: np.isfinite(r["t50_delay_pct"]) and r["t50_delay_pct"] >= 10.0),
                min_M_t50_25pct=min_M(lambda r: np.isfinite(r["t50_delay_pct"]) and r["t50_delay_pct"] >= 25.0),
            ))
    return out


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
def decide(rows):
    if any(not np.isfinite(r["d_on_max_pp"]) for r in rows):
        return "NUMERICALLY_INCONCLUSIVE", {}
    # Plausibility envelope (adversarial, deliberately generous ceilings):
    #  - active nuclear Cas9-RNP: plausibly 1e3-1e5 molecules; scarce end ~1e2.
    #  - occupancy-competent off-target loci per guide: plausibly <=~300;
    #    1e3 is high; 1e4+ is not credible.
    PLAUS_M = 1000          # generous ceiling on occupancy-competent loci
    PLAUS_NC_MIN = 100      # below ~100 active Cas9 molecules is scarce/edge
    def plausible(r):
        return r["M_loci"] <= PLAUS_M and r["Nc_molecules"] >= PLAUS_NC_MIN

    strong = [r for r in rows if r["d_on_max_pp"] >= 5.0
              or (np.isfinite(r["t50_delay_pct"]) and r["t50_delay_pct"] >= 25.0)]
    strong_plaus = [r for r in strong if plausible(r)]
    weak = [r for r in rows if r["d_on_max_pp"] >= 1.0
            or (np.isfinite(r["t50_delay_pct"]) and r["t50_delay_pct"] >= 10.0)]
    weak_plaus = [r for r in weak if plausible(r)]

    info = dict(n_strong=len(strong), n_strong_plausible=len(strong_plaus),
                n_weak_plausible=len(weak_plaus),
                plausibility_ceiling_M=PLAUS_M, plausibility_floor_Nc=PLAUS_NC_MIN)
    if len(strong_plaus) >= 6:
        return "CONTINUE_BROAD_PLAUSIBLE_REGIME", info
    if strong_plaus or weak_plaus:
        return "PIVOT_NARROW_PLAUSIBLE_REGIME", info
    return "ABORT_PHYSICALLY_IMPLAUSIBLE", info


# ---------------------------------------------------------------------------
def write_report(rows, thresholds, verdict, vinfo):
    L = ["# Physical-Scaling Plausibility Stress Test\n",
         f"**Verdict: `{verdict}`**\n"]
    L.append("## Copy-number accounting\n")
    L.append("Diploid locus concentration derived per nuclear volume "
             "(S_site = 2 / (N_A · V)); Cas9 given as molecule count.\n")
    L.append("| V (fL) | S_site (nM) | 1 Cas9 molecule (nM) | 100 nM Cas9 = N molecules |")
    L.append("|---|---|---|---|")
    for V in VOLUMES_FL:
        L.append(f"| {V} | {s_site(V):.5f} | {cas_conc(1,V):.5f} | {100/cas_conc(1,V):.0f} |")
    L.append(f"\n_The original universal S_site=0.00634 nM corresponds to a "
             f"~{2/(NA*0.00634e-9)*1e15:.0f} fL nucleus (i.e. the 500 fL case). "
             f"It is NOT reused universally here._\n")
    L.append("**Volume (near-)invariance:** total-site/Cas9 ratio = 2M/Nc is "
             "volume-independent; volume only rescales absolute time. The phase "
             "boundary therefore tracks the molecule-count ratio 2M/Nc.\n")

    # headline extremes
    strong = [r for r in rows if r["d_on_max_pp"] >= 5.0]
    def smallest(rs, key):
        return min(rs, key=key) if rs else None
    s5 = smallest([r for r in rows if r["d_on_max_pp"] >= 5.0],
                  key=lambda r: (r["M_loci"], -r["Nc_molecules"]))
    t25 = smallest([r for r in rows if np.isfinite(r["t50_delay_pct"]) and r["t50_delay_pct"] >= 25.0],
                   key=lambda r: (r["M_loci"], -r["Nc_molecules"]))
    L.append("## Smallest regimes crossing thresholds\n")
    if s5:
        L.append(f"- **>5 pp transient:** smallest at M={s5['M_loci']} loci "
                 f"({s5['site_copies']} copies), Cas9={s5['Nc_molecules']} molecules "
                 f"({s5['C_total_nM']:.3g} nM), V={s5['V_fL']} fL "
                 f"(ratio 2M/Nc={s5['ratio_copies_per_cas9']:.1f}); dOnMax={s5['d_on_max_pp']:.1f} pp.")
    else:
        L.append("- **>5 pp transient:** never reached anywhere on the grid.")
    if t25:
        L.append(f"- **>25% t50 delay:** smallest at M={t25['M_loci']} loci, "
                 f"Cas9={t25['Nc_molecules']} molecules ({t25['C_total_nM']:.3g} nM), "
                 f"V={t25['V_fL']} fL (ratio {t25['ratio_copies_per_cas9']:.1f}).")
    else:
        L.append("- **>25% t50 delay:** never reached anywhere on the grid.")
    L.append("")

    L.append("## Minimum competitor loci to cross each threshold\n")
    L.append("(None = threshold not reached even at M=30000)\n")
    L.append("| V (fL) | Cas9 (molecules / nM) | M for 1pp | M for 5pp | M for 10% t50 | M for 25% t50 |")
    L.append("|---|---|---|---|---|---|")
    for r in thresholds:
        L.append(f"| {r['V_fL']} | {r['Nc_molecules']} / {r['C_total_nM']:.3g} | "
                 f"{r['min_M_1pp']} | {r['min_M_5pp']} | {r['min_M_t50_10pct']} | {r['min_M_t50_25pct']} |")
    L.append("")

    L.append("## Honesty audit: PAM sites vs occupancy-competent competitors\n")
    L.append("- **Total NGG PAM sites** in a human diploid genome are ~4×10^8 "
             "(one every ~8 bp, both strands). This is NOT the competitor count.")
    L.append("- **Occupancy-competent competitors** (M here) are the far smaller "
             "subset of loci that actually bind and hold Cas9 long enough to "
             "sequester it. Empirical off-target assays (GUIDE-seq / CIRCLE-seq) "
             "typically detect tens to low-hundreds of cleavage-active off-targets "
             "per guide; the number that transiently *sequester* Cas9 is unknown "
             "but is bounded by that order, not by the 10^8 PAM count.")
    L.append("- **We do NOT assume 10^4 effective competitors is realistic.** "
             "M=10^4 occupancy-competent loci for a single guide has no empirical "
             "support and is treated here as an upper-bound stress point only.")
    L.append("- **Controlling unknown:** the effect is governed by the ratio "
             "2M/Nc = (occupancy-competent off-target copies) / (active nuclear "
             "Cas9-RNP molecules). Both numerator and denominator are unmeasured "
             "in most experiments. Relevance now hinges on ONE composite unknown: "
             "**whether occupancy-competent off-target loci can rival or exceed the "
             "active nuclear Cas9-RNP molecule count.** Absolute timescale (set by "
             "volume + Cas9 conc) is a secondary honesty flag: at scarce Cas9 the "
             "delays play out over 10^5-10^7 s, longer than RNP lifetime / cell "
             "cycle, so a large % delay may be biologically moot.\n")

    L.append(f"## Verdict rationale\n")
    L.append(f"Plausibility envelope: occupancy-competent loci M ≤ {vinfo.get('plausibility_ceiling_M')}, "
             f"active Cas9 ≥ {vinfo.get('plausibility_floor_Nc')} molecules. "
             f"Grid cells with a strong effect (>5 pp or >25% t50): {vinfo.get('n_strong')}; "
             f"of those inside the plausibility envelope: {vinfo.get('n_strong_plausible')}; "
             f"weak-effect cells inside envelope: {vinfo.get('n_weak_plausible')}.\n")

    L.append("## Full grid\n")
    L.append("| V | Nc | C(nM) | M | copies | ratio2M/Nc | dOnMax(pp) | dOnFin(pp) | t10% | t50% | t90% | depl% | dur>5%(s) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        def f(x): return f"{x:.1f}" if (x is not None and np.isfinite(x)) else "—"
        L.append(f"| {r['V_fL']} | {r['Nc_molecules']} | {r['C_total_nM']:.3g} | {r['M_loci']} | "
                 f"{r['site_copies']} | {r['ratio_copies_per_cas9']:.2g} | {r['d_on_max_pp']:.2f} | "
                 f"{r['d_on_final_pp']:.3f} | {f(r['t10_delay_pct'])} | {f(r['t50_delay_pct'])} | "
                 f"{f(r['t90_delay_pct'])} | {r['depletion_max_pct']:.1f} | {r['dur_gt5pct_depl_s']:.2g} |")
    with open(os.path.join(OUT, "plausibility_report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L))


def main():
    if not g.gpu_banner():
        raise SystemExit("GPU unavailable")
    # time grid must reach on-target completion even at scarce Cas9 (slow binding)
    t = np.insert(np.logspace(-2, 7.5, 120), 0, 0.0)
    print("grid:", len(VOLUMES_FL), "vols x", len(NC_LIST), "Cas9 x", len(M_LIST), "M")
    t0 = time.time()
    rows = run_grid(t)
    print(f"grid done in {time.time()-t0:.0f}s")
    thresholds = min_competitor_thresholds(rows)
    phase_maps(rows)
    verdict, vinfo = decide(rows)

    with open(os.path.join(OUT, "plausibility_grid.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(OUT, "plausibility.json"), "w") as f:
        json.dump(dict(grid=rows, thresholds=thresholds, verdict=verdict,
                       verdict_info=vinfo,
                       accounting=dict(NA=NA, volumes_fL=VOLUMES_FL,
                                       s_site_nM={V: s_site(V) for V in VOLUMES_FL})),
                  f, indent=2, default=float)
    write_report(rows, thresholds, verdict, vinfo)
    print("VERDICT:", verdict, vinfo)


if __name__ == "__main__":
    main()
