"""
diagnostic_finite_lifetime.py
=============================
Minimal, targeted diagnostic (NOT a grid) for the pivot audit.

Question: if active Cas9 exists only for a finite window [0, T_cas] (transient
delivery / finite RNP lifetime, step model), does off-target burden reduce the
FINAL on-target editing yield in a PLAUSIBLE regime -- as opposed to the killed
steady-state (t->inf) hypothesis where both models reach ~1.0?

Method: the competition dynamics up to T_cas are exactly the constant-total-Cas9
finite-pool dynamics we already trust (fast Illinois solver). Reading the curves
at a hard cutoff T_cas gives the editing yield if the enzyme vanishes at T_cas.
Rigorous bound: deficit(T_cas) = on_A(T_cas) - on_B(T_cas) <= max_t transient dOn.

No parameter tuning. A handful of conditions only.
"""
import json, numpy as np
import go_no_go as g
import audit_panels as ap
import plausibility_scaling as ps

EPS, BASE_FWD = g.load_parameters()
H = 3600.0
CUTOFFS_H = [6, 12, 24, 48, 72]          # realistic RNP active windows (hours)

# (label, M loci, Nc molecules, V fL, plausibility)
CONDS = [
    ("plausible_max  M1000/Cas100/250fL", 1000, 100, 250, "plausible (envelope corner)"),
    ("plausible_mod  M1000/Cas1000/500fL", 1000, 1000, 500, "plausible (moderate Cas9)"),
    ("plausible_low  M300/Cas300/500fL",   300, 300, 500, "plausible (defensible burden)"),
    ("IMPLAUSIBLE    M3000/Cas100/250fL",  3000, 100, 250, "IMPLAUSIBLE burden (contrast)"),
    ("IMPLAUSIBLE    M10000/Cas100/500fL", 10000, 100, 500, "IMPLAUSIBLE burden (contrast)"),
]


def t_reach(t, curve, frac):
    if curve[-1] < frac:
        return np.nan
    i = np.searchsorted(curve, frac)
    if i <= 0:
        return t[0]
    t0, t1, c0, c1 = t[i-1], t[i], curve[i-1], curve[i]
    return t1 if c1 == c0 else t0 + (frac - c0)*(t1-t0)/(c1-c0)


def at(t, curve, T):
    return float(np.interp(T, t, curve))


def main():
    g.gpu_banner()
    t = np.insert(np.logspace(-2, 6.5, 240), 0, 0.0)
    out = []
    for label, M, Nc, V, plaus in CONDS:
        S = ps.s_site(V)
        C = ps.cas_conc(Nc, V)
        mms, _ = ap.build_panel_det("DEFENSIBLE", M)
        on_A = g.nature_cleavage(C, [], t, EPS, BASE_FWD)          # no competition
        res = ps.simulate_fast(C, mms, t, S, use_gpu=True)
        on_B = res["on"]                                            # with burden
        d = on_A - on_B
        kmax = int(np.argmax(d))
        rec = dict(label=label, plausibility=plaus, M=M, Nc=Nc, V_fL=V,
                   C_total_nM=C, S_site_nM=S, ratio_2M_over_Nc=2.0*M/Nc,
                   max_transient_pp=float(d[kmax]*100),
                   time_of_max_transient_s=float(t[kmax]),
                   time_of_max_transient_h=float(t[kmax]/H),
                   on_A_final=float(on_A[-1]), on_B_final=float(on_B[-1]),
                   cutoffs={})
        for Th in CUTOFFS_H:
            T = Th*H
            a, b = at(t, on_A, T), at(t, on_B, T)
            rec["cutoffs"][f"{Th}h"] = dict(
                on_A_yield=round(a, 4), on_B_yield=round(b, 4),
                deficit_pp=round((a-b)*100, 3),
                relative_loss_pct=round((a-b)/a*100 if a > 1e-9 else 0.0, 2))
        out.append(rec)
        cs = rec["cutoffs"]
        print(f"\n{label}  [{plaus}]  C={C:.3g}nM ratio2M/Nc={2*M/Nc:.1f}")
        print(f"  max transient {rec['max_transient_pp']:.2f} pp at "
              f"t={rec['time_of_max_transient_h']:.1f} h; "
              f"final A={rec['on_A_final']:.4f} B={rec['on_B_final']:.4f}")
        for Th in CUTOFFS_H:
            c = cs[f"{Th}h"]
            print(f"  T_cas={Th:3d}h: yield noComp={c['on_A_yield']:.3f} "
                  f"withBurden={c['on_B_yield']:.3f}  deficit={c['deficit_pp']:6.2f} pp "
                  f"({c['relative_loss_pct']:.1f}% rel)")
    json.dump(out, open("artifacts/plausibility/finite_lifetime_diagnostic.json", "w"),
              indent=2, default=float)
    print("\nsaved artifacts/plausibility/finite_lifetime_diagnostic.json")


if __name__ == "__main__":
    main()
