"""
audit_panels.py
===============
Forensic audit: is the surviving M=1000-10000 competitive signal biologically
real, or an artifact of how the CONTROLLED_STRESS / DEFENSIBLE competitor panels
were constructed?

No new large campaign, no parameter tuning. Targeted ablations only.

Strategy
--------
1. Document exact panel construction (patterns, multiplicities, positions).
2. Quantify per-class occupancy + Cas9 sequestration (who actually holds the
   scarce Cas9) using time-integrated bound-state occupancy.
3. Ablate site classes at M=1000 and 10000 (C=0.1 nM, the 17.29 pp / 59% corner)
   and test whether the transient on-target gap and timing delay survive.
   Adversarial control: remove a RANDOM 10% as well as the strongest 10%.
4. Audit the isolated ranking inversion at M=1000, C=1 nM.

The panel realization is PINNED with fixed integer seeds (go_no_go.build_panel
seeds off Python's per-process randomized str hash, so the campaign's exact draw
is not reproducible; at M>=1000 the class aggregates are draw-independent by
large numbers, which we verify by reproducing the baseline magnitude).
"""

import os, json, csv, time
import numpy as np
import go_no_go as g

OUT = os.path.join(g.BASE_DIR, "artifacts", "go_no_go")
os.makedirs(OUT, exist_ok=True)
EPS, BASE_FWD = g.load_parameters()

# Fixed, reproducible name->seed (replaces randomized hash(name))
NAME_SEED = {"CONTROLLED_STRESS": 101, "DEFENSIBLE": 202}
GUIDE = g.GUIDE_LEN
DISTAL_MIN = 13     # PAM-distal region start (nt 13-20)
SEED_MAX = 8        # PAM-proximal seed region (nt 1-8)


# ---------------------------------------------------------------------------
# Reproducible replica of go_no_go.build_panel (identical rules, fixed seed)
# ---------------------------------------------------------------------------
def _mm(rng, lo, hi, k):
    hi = min(hi, GUIDE)
    k = min(k, hi - lo + 1)
    return sorted(int(x) for x in rng.choice(range(lo, hi + 1), size=k, replace=False))


def build_panel_det(name, M):
    rng = np.random.default_rng(NAME_SEED[name] + M)
    mms = [[]]; cls = ["on_target"]
    n_off = M - 1
    if n_off <= 0:
        return mms, cls
    if name == "CONTROLLED_STRESS":
        counts = g._counts(dict(near_perfect=0.20, seed=0.20, mid=0.20,
                                distal=0.20, weak=0.20), n_off)
        for _ in range(counts["near_perfect"]):
            mms.append(_mm(rng, 18, 20, 1)); cls.append("near_perfect")
        for _ in range(counts["seed"]):
            mms.append(_mm(rng, 1, 8, rng.integers(1, 3))); cls.append("seed")
        for _ in range(counts["mid"]):
            mms.append(_mm(rng, 9, 12, rng.integers(1, 3))); cls.append("mid")
        for _ in range(counts["distal"]):
            mms.append(_mm(rng, 13, 20, rng.integers(1, 3))); cls.append("distal")
        for _ in range(counts["weak"]):
            mms.append(_mm(rng, 1, 20, rng.integers(4, 6))); cls.append("weak")
    elif name == "DEFENSIBLE":
        counts = g._counts(dict(sink=0.02, relevant=0.18, weak=0.80), n_off)
        for _ in range(counts["sink"]):
            mms.append(_mm(rng, 15, 20, rng.integers(1, 3))); cls.append("sink")
        for _ in range(counts["relevant"]):
            mms.append(_mm(rng, 9, 20, rng.integers(3, 5))); cls.append("relevant")
        for _ in range(counts["weak"]):
            if rng.random() < 0.6:
                sm = _mm(rng, 1, 8, rng.integers(1, 3))
                dm = _mm(rng, 9, 20, 3)
                mms.append(sorted(set(sm + dm)))
            else:
                mms.append(_mm(rng, 1, 20, 5))
            cls.append("weak")
    mms = mms[:M]; cls = cls[:M]
    while len(mms) < M:
        mms.append(_mm(rng, 1, 20, 5)); cls.append("weak")
    return mms, cls


# ---------------------------------------------------------------------------
# Residence time of an off-target (independent-site strength proxy)
# MFPT to leave the bound manifold back toward solution (excludes cleavage sink)
# ---------------------------------------------------------------------------
def residence_time(mm):
    fwd, bck = g.get_rate_pair(EPS, BASE_FWD, mm)
    fwd_d = fwd.copy(); fwd_d[-1] = 0.0
    diag = -(fwd_d[1:] + bck[1:])
    K_sub = np.diag(diag) + np.diag(bck[2:], k=1) + np.diag(fwd_d[1:-1], k=-1)
    try:
        tau = np.linalg.inv(-K_sub)[:, 0].sum()
    except np.linalg.LinAlgError:
        tau = np.inf
    return float(tau)


def mm_count_bin(n):
    return "0" if n == 0 else ("1" if n == 1 else ("2" if n == 2 else
           ("3" if n == 3 else ("4" if n == 4 else "5+"))))


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------
def t_half(t, curve):
    final = curve[-1]
    if final <= 1e-9:
        return np.nan
    target = 0.5 * final
    idx = np.searchsorted(curve, target)
    if idx <= 0 or idx >= len(curve):
        return np.nan
    t0, t1, c0, c1 = t[idx - 1], t[idx], curve[idx - 1], curve[idx]
    return t1 if c1 == c0 else t0 + (target - c0) * (t1 - t0) / (c1 - c0)


def run_condition(C, mms, t, want_occ=False):
    """Return competitive on-target curve, Cf, depletion, and optional occ."""
    resB = g.simulate_competitive(C, mms, t, EPS, BASE_FWD, use_gpu=True,
                                  track_occ=want_occ)
    return resB


# ---------------------------------------------------------------------------
# PART 1+2: documentation & class contributions
# ---------------------------------------------------------------------------
def document_and_contributions(panels, Ms, C_ref, t):
    doc = {}
    on_A = g.nature_cleavage(C_ref, [], t, EPS, BASE_FWD)   # Nature on-target
    th_A = t_half(t, on_A)
    for panel in panels:
        for M in Ms:
            mms, cls = build_panel_det(panel, M)
            cls = np.array(cls)
            counts = [len(m) for m in mms]
            # histogram by mismatch count
            binhist = {}
            for c in counts:
                b = mm_count_bin(c)
                binhist[b] = binhist.get(b, 0) + 1
            # class multiplicities
            clsmult = {c: int((cls == c).sum()) for c in sorted(set(cls))}
            # position histogram (off-targets only)
            poshist = {p: 0 for p in range(1, GUIDE + 1)}
            for m in mms[1:]:
                for p in m:
                    poshist[p] += 1
            near_perfect = int(sum(1 for c in counts[1:] if c <= 1))
            distal_only = int(sum(1 for m in mms[1:]
                                  if len(m) > 0 and all(p >= DISTAL_MIN for p in m)))

            # independent occupancy (capacity) + competitive sequestration (actual)
            resI = g.simulate_competitive(C_ref, mms, t, EPS, BASE_FWD,
                                          use_gpu=True, constant=True, track_occ=True)
            resC = g.simulate_competitive(C_ref, mms, t, EPS, BASE_FWD,
                                          use_gpu=True, constant=False, track_occ=True)
            occI = resI["occ"]; occC = resC["occ"]
            # aggregate by class
            def agg(vec, labels):
                d = {}
                for lab in sorted(set(labels)):
                    d[lab] = float(vec[labels == lab].sum())
                return d
            occI_cls = agg(occI, cls); occC_cls = agg(occC, cls)
            # aggregate by mismatch-count bin
            cbin = np.array([mm_count_bin(c) for c in counts])
            occI_bin = agg(occI, cbin); occC_bin = agg(occC, cbin)
            tot_seq = float(occC[1:].sum()) + 1e-30   # off-target sequestration
            depl = float((1 - resC["Cf"].min() / C_ref) * 100)
            th_B = t_half(t, resC["clv"][:, 0])
            timing = float((th_B - th_A) / th_A * 100) if th_A and np.isfinite(th_B) else np.nan
            don = float(np.max(np.abs(resC["clv"][:, 0] - on_A)) * 100)

            doc[f"{panel}_M{M}"] = dict(
                panel=panel, M=M, C_ref=C_ref,
                mm_count_hist=binhist, class_mult=clsmult,
                n_near_perfect_le1mm=near_perfect, n_distal_only=distal_only,
                position_hist=poshist,
                occ_independent_by_class=occI_cls,
                occ_independent_by_mmbin=occI_bin,
                sequestration_competitive_by_class=occC_cls,
                sequestration_competitive_by_mmbin=occC_bin,
                sequestration_frac_by_class={k: v / tot_seq
                                             for k, v in occC_cls.items()},
                depletion_pct=depl, timing_pct=timing, d_on_max_pp=don,
            )
            print(f"  [doc] {panel} M={M}: depl={depl:.2f}% timing={timing:.1f}% "
                  f"dOn={don:.2f}pp  seq_frac={ {k: round(v/tot_seq,3) for k,v in occC_cls.items()} }",
                  flush=True)
    return doc, dict(on_A=on_A.tolist(), th_A=float(th_A))


# ---------------------------------------------------------------------------
# PART 3+4: ablations
# ---------------------------------------------------------------------------
def ablate(panel, M, C, t, on_A, th_A):
    mms, cls = build_panel_det(panel, M)
    cls = np.array(cls)
    counts = np.array([len(m) for m in mms])
    # independent occupancy for strength ranking
    resI = g.simulate_competitive(C, mms, t, EPS, BASE_FWD, use_gpu=True,
                                  constant=True, track_occ=True)
    occI = resI["occ"].copy()
    occI[0] = -np.inf  # never rank/remove the on-target
    off_order = np.argsort(-occI)   # strongest sink first (on-target excluded via -inf sorts last)

    n_off = M - 1
    rng = np.random.default_rng(7)

    def strongest_frac(frac):
        k = max(1, int(round(frac * n_off)))
        return set(off_order[:k].tolist())

    def keep_mask(drop_idx):
        keep = np.ones(M, dtype=bool)
        keep[list(drop_idx)] = False
        keep[0] = True  # always keep on-target
        return keep

    distal_only = set(i for i in range(1, M)
                      if len(mms[i]) > 0 and all(p >= DISTAL_MIN for p in mms[i]))
    variants = {
        "baseline": set(),
        "drop_near_perfect_le1mm": set(i for i in range(1, M) if counts[i] <= 1),
        "drop_1_2mm": set(i for i in range(1, M) if counts[i] in (1, 2)),
        "drop_pam_distal_only": distal_only,
        "drop_top1pct_sink": strongest_frac(0.01),
        "drop_top5pct_sink": strongest_frac(0.05),
        "drop_top10pct_sink": strongest_frac(0.10),
        "drop_random10pct_CTRL": set(rng.choice(range(1, M),
                                     size=max(1, int(0.10 * n_off)), replace=False).tolist()),
    }
    rows = []
    for name, drop in variants.items():
        keep = keep_mask(drop)
        sub = [mms[i] for i in range(M) if keep[i]]
        resB = g.simulate_competitive(C, sub, t, EPS, BASE_FWD, use_gpu=True)
        onB = resB["clv"][:, 0]
        don_max = float(np.max(np.abs(onB - on_A)) * 100)
        don_final = float((onB[-1] - on_A[-1]) * 100)
        th_B = t_half(t, onB)
        timing = float((th_B - th_A) / th_A * 100) if th_A and np.isfinite(th_B) else np.nan
        depl = float((1 - resB["Cf"].min() / C) * 100)
        rows.append(dict(panel=panel, M=M, C=C, ablation=name,
                         n_removed=int(M - keep.sum()), M_eff=int(keep.sum()),
                         d_on_max_pp=don_max, d_on_final_pp=don_final,
                         timing_pct=timing, depletion_pct=depl))
        print(f"  [ablate] {panel} M={M} C={C} {name:26s} rm={M-keep.sum():5d} "
              f"dOnMax={don_max:6.2f}pp timing={timing:6.1f}% depl={depl:5.2f}%",
              flush=True)
    return rows


# ---------------------------------------------------------------------------
# PART 5: ranking-inversion audit at M=1000, C=1 nM
# ---------------------------------------------------------------------------
def audit_inversion(panel, M, C, t):
    mms, cls = build_panel_det(panel, M)
    cls = np.array(cls)
    off_idx = np.where(cls != "on_target")[0]
    resA = g.simulate_competitive(C, mms, t, EPS, BASE_FWD, use_gpu=True, constant=True)
    resB = g.simulate_competitive(C, mms, t, EPS, BASE_FWD, use_gpu=True, constant=False)
    a = resA["clv"][-1, off_idx]     # independent final cleavage
    b = resB["clv"][-1, off_idx]     # competitive final cleavage
    # spearman
    sp = g.spearman(a, b)
    # find inversions among top competitors with margin
    order = np.argsort(-a)[:80]
    aa, bb = a[order], b[order]
    inv = []
    for i in range(len(aa)):
        for j in range(i + 1, len(aa)):
            if (aa[i] - aa[j]) > 1e-9 and (bb[i] - bb[j]) < -1e-9:
                inv.append((order[i], order[j], float(aa[i] - aa[j]), float(bb[i] - bb[j])))
    details = []
    for (i, j, da, db) in inv[:20]:
        gi, gj = off_idx[i], off_idx[j]
        details.append(dict(
            site_i=int(gi), site_j=int(gj),
            mm_i=mms[gi], mm_j=mms[gj], class_i=str(cls[gi]), class_j=str(cls[gj]),
            A_i=float(a[i]), A_j=float(a[j]), B_i=float(b[i]), B_j=float(b[j]),
            A_gap=da, B_gap=db, gap_vs_tol=abs(da) / 1e-6))
    return dict(panel=panel, M=M, C=C, spearman=float(sp),
                n_inversions_margin1e9=len(inv),
                numerical_tol=1e-6, inversions=details,
                max_abs_final_gap=float(np.max(np.abs(a - b))))


# ---------------------------------------------------------------------------
def main():
    ok = g.gpu_banner()
    if not ok:
        raise SystemExit("GPU unavailable")
    t_abl = np.insert(np.logspace(-2, 6, 110), 0, 0.0)
    C_ref = 0.1
    panels = ["DEFENSIBLE", "CONTROLLED_STRESS"]

    print("\n== PART 1/2: documentation + class contributions ==")
    doc, natref = document_and_contributions(panels, [1000, 10000], C_ref, t_abl)

    # baseline reproduction check vs campaign headline
    repro = {k: dict(d_on_max_pp=doc[k]["d_on_max_pp"], timing_pct=doc[k]["timing_pct"],
                     depletion_pct=doc[k]["depletion_pct"])
             for k in ("DEFENSIBLE_M10000", "DEFENSIBLE_M1000")}

    print("\n== PART 3/4: ablations (C=0.1 nM) ==")
    on_A = np.array(natref["on_A"]); th_A = natref["th_A"]
    ablation_rows = []
    for panel, M in [("DEFENSIBLE", 10000), ("DEFENSIBLE", 1000),
                     ("CONTROLLED_STRESS", 10000)]:
        ablation_rows += ablate(panel, M, C_ref, t_abl, on_A, th_A)

    print("\n== PART 5: ranking-inversion audit (M=1000, C=1 nM) ==")
    t_inv = np.insert(np.logspace(-2, 6, 150), 0, 0.0)
    inv_audit = audit_inversion("DEFENSIBLE", 1000, 1.0, t_inv)
    print(f"  spearman={inv_audit['spearman']:.6f} "
          f"n_inv={inv_audit['n_inversions_margin1e9']} "
          f"max_final_gap={inv_audit['max_abs_final_gap']:.3e}", flush=True)

    # ---- verdict ----
    verdict = decide(doc, ablation_rows, inv_audit)

    # ---- save ----
    with open(os.path.join(OUT, "panel_audit.json"), "w") as f:
        json.dump(dict(construction_doc=doc, baseline_reproduction=repro,
                       nature_ref=dict(th_A=th_A),
                       inversion_audit=inv_audit, verdict=verdict),
                  f, indent=2, default=float)
    with open(os.path.join(OUT, "ablation_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(ablation_rows[0].keys()))
        w.writeheader()
        for r in ablation_rows:
            w.writerow(r)
    write_report(doc, repro, ablation_rows, inv_audit, verdict)
    print("\nVERDICT:", verdict)


def decide(doc, ablation_rows, inv_audit):
    """Adversarial verdict on whether the DEFENSIBLE-panel signal is a
    panel-construction artifact.

    Discriminator hierarchy:
      * RARE_STRONG_SINKS: removing the few injected near-cognate/strong-sink
        sites (near-perfect <=1 mm, 1-2 mm, or the top ~1% by occupancy)
        collapses the signal.
      * ROBUST_TO_PANEL_COMPOSITION: those removals barely change it, random
        subsetting reduces it ~proportionally, and the majority survives even
        removing the strongest 10% -> a broad mass-action burden effect that
        does not depend on composition.
      * PANEL_ARTIFACT: neither -- the signal hinges on one specific,
        non-rare, over-represented synthetic class.
    """
    if any(not np.isfinite(r["d_on_max_pp"]) for r in ablation_rows):
        return "NUMERICALLY_INCONCLUSIVE"
    if inv_audit.get("spearman") is not None and np.isfinite(inv_audit["spearman"]) \
            and inv_audit["spearman"] < 0.9:
        return "NUMERICALLY_INCONCLUSIVE"
    R = {(r["panel"], r["M"], r["ablation"]): r for r in ablation_rows}
    key = ("DEFENSIBLE", 10000)
    base = R[(key[0], key[1], "baseline")]
    if base["d_on_max_pp"] < 1.0:
        return "NUMERICALLY_INCONCLUSIVE"

    def retain(ablation, metric="d_on_max_pp"):
        return R[(key[0], key[1], ablation)][metric] / (base[metric] + 1e-12)

    r_np = retain("drop_near_perfect_le1mm")
    r_12 = retain("drop_1_2mm")
    r_top1 = retain("drop_top1pct_sink")
    r_top10 = retain("drop_top10pct_sink")
    r_rand10 = retain("drop_random10pct_CTRL")

    # Rare strong sinks / near-cognates drive it?
    if r_np < 0.6 or r_12 < 0.5 or r_top1 < 0.6:
        return "SIGNAL_DRIVEN_BY_RARE_STRONG_SINKS"
    # Robust broad burden: injected sinks irrelevant, random ~ proportional,
    # majority survives strongest-10% removal.
    if r_np > 0.9 and r_12 > 0.85 and r_rand10 > 0.85 and r_top10 > 0.55:
        return "SIGNAL_ROBUST_TO_PANEL_COMPOSITION"
    return "SIGNAL_PANEL_ARTIFACT"


def write_report(doc, repro, ablation_rows, inv_audit, verdict):
    L = ["# Panel-Construction Audit of the M=1000-10000 Competitive Signal\n",
         f"**Verdict: `{verdict}`**\n"]
    L.append("## Baseline reproduction (fresh deterministic panel realization)\n")
    L.append("| condition | dOnMax(pp) | timing% | depletion% |")
    L.append("|---|---|---|---|")
    for k, v in repro.items():
        L.append(f"| {k} | {v['d_on_max_pp']:.2f} | {v['timing_pct']:.1f} | {v['depletion_pct']:.2f} |")
    L.append("\n_Campaign headline was 17.29 pp / 59% / 38% at DEFENSIBLE M=10000, C=0.1._\n")

    L.append("## Sequestration fraction by class (who holds the scarce Cas9, C=0.1 nM)\n")
    for k in ("DEFENSIBLE_M10000", "DEFENSIBLE_M1000",
              "CONTROLLED_STRESS_M10000", "CONTROLLED_STRESS_M1000"):
        if k not in doc:
            continue
        d = doc[k]
        fr = d["sequestration_frac_by_class"]
        mult = d["class_mult"]
        L.append(f"### {k}")
        L.append(f"- class multiplicities: {mult}")
        L.append(f"- mismatch-count histogram: {d['mm_count_hist']}")
        L.append(f"- near-perfect (<=1 mm) sites: {d['n_near_perfect_le1mm']}; "
                 f"distal-only near-cognates: {d['n_distal_only']}")
        L.append(f"- **sequestration fraction by class:** "
                 f"{ {kk: round(vv,3) for kk,vv in fr.items()} }")
        L.append("")

    L.append("## Ablation results (C=0.1 nM)\n")
    L.append("| panel | M | ablation | removed | M_eff | dOnMax(pp) | dOnFinal(pp) | timing% | depl% |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for r in ablation_rows:
        L.append(f"| {r['panel']} | {r['M']} | {r['ablation']} | {r['n_removed']} | "
                 f"{r['M_eff']} | {r['d_on_max_pp']:.2f} | {r['d_on_final_pp']:.3f} | "
                 f"{r['timing_pct']:.1f} | {r['depletion_pct']:.2f} |")
    L.append("")

    L.append("## Ranking-inversion audit (DEFENSIBLE, M=1000, C=1 nM)\n")
    L.append(f"- Spearman(off-target final cleavage, A vs B) = {inv_audit['spearman']:.6f}")
    L.append(f"- inversions above margin 1e-9: {inv_audit['n_inversions_margin1e9']}")
    L.append(f"- max |A-B| final off-target cleavage = {inv_audit['max_abs_final_gap']:.3e}")
    L.append(f"- numerical tolerance = {inv_audit['numerical_tol']:.0e}")
    if inv_audit["inversions"]:
        L.append("\n| site_i | site_j | mm_i | mm_j | A_gap | B_gap | gap/tol |")
        L.append("|---|---|---|---|---|---|---|")
        for d in inv_audit["inversions"]:
            L.append(f"| {d['site_i']} | {d['site_j']} | {d['mm_i']} | {d['mm_j']} | "
                     f"{d['A_gap']:.2e} | {d['B_gap']:.2e} | {d['gap_vs_tol']:.2f} |")
    else:
        L.append("\n_No inversions above a 1e-9 margin in this realization._")
    L.append("")

    with open(os.path.join(OUT, "panel_audit_report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"total {time.time()-t0:.1f}s")
