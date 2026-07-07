"""
simulate_competitive_cuda.py
=======================
Multi-site competitive SpCas9 kinetic model, accelerated with CuPy (CUDA).
Extends the Nature Communications (2022) single-site Markov chain to a
coupled non-linear dynamical system where all genomic sites share a finite
Cas9 pool C_free(t).

Solver: Vectorized Crank-Nicolson with Brent root-finding for C_free.
        Calculations for large M are offloaded to GPU.

Author: Priyansh Kuniyal / Antigravity
"""

import argparse
import os, json, time
import numpy as np
import cupy as cp
import scipy.linalg as linalg
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# -- aesthetics --------------------------------------------------------------
sns.set_style("ticks")
plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "figure.titlesize": 14,
    "legend.fontsize": 9, "font.family": "sans-serif",
    "axes.spines.top": False, "axes.spines.right": False
})


def _gpu_name():
    device_count = cp.cuda.runtime.getDeviceCount()
    if device_count < 1:
        raise RuntimeError("No CUDA-capable GPU detected by CuPy")
    props = cp.cuda.runtime.getDeviceProperties(0)
    name = props.get("name", b"unknown")
    if isinstance(name, (bytes, bytearray)):
        name = name.decode("utf-8", errors="ignore")
    return name, device_count

# -- paths --------------------------------------------------------------------
BASE_DIR      = os.path.dirname(os.path.abspath(__file__))
EPS_PATH      = os.path.join(BASE_DIR, "dashboard_repo", "parameters", "SpCas9_epsilon.txt")
RATES_PATH    = os.path.join(BASE_DIR, "dashboard_repo", "parameters", "SpCas9_forward_rates.txt")
ARTIFACTS_DIR = os.path.join(BASE_DIR, "artifacts")
os.makedirs(ARTIFACTS_DIR, exist_ok=True)

# -- physical constant --------------------------------------------------------
S_DIPLOID = 0.00634   # 6.34 pM - diploid concentration of a unique genomic locus (nM)

# ----------------------------------------------------------------------------
# 1. PARAMETER LOADING
# ----------------------------------------------------------------------------
def load_parameters():
    epsilon       = np.loadtxt(EPS_PATH)
    forward_rates = np.loadtxt(RATES_PATH)
    return epsilon, forward_rates

def get_step_energies(epsilon, mismatch_positions, guide_length=20):
    energies = -1.0 * epsilon[0:(guide_length + 1)]
    energies[0] = epsilon[0]
    if len(mismatch_positions) > 0:
        epsI = epsilon[(guide_length + 1):]
        for pos in mismatch_positions:
            if 1 <= pos <= guide_length:
                energies[pos] += epsI[pos - 1]
    return energies

def get_rate_pair(epsilon, base_fwd, mismatch_positions, guide_length=20):
    energies = get_step_energies(epsilon, mismatch_positions, guide_length)
    bck = np.zeros_like(base_fwd)
    bck[1:] = base_fwd[:-1] * np.exp(energies)
    bck[0]  = 0.0
    return base_fwd.copy(), bck

# ----------------------------------------------------------------------------
# 2. VECTORIZED CRANK-NICOLSON TRIDIAGONAL SOLVER (GPU)
# ----------------------------------------------------------------------------
def _tri_solve_gpu(L, D, U, B):
    n = D.shape[1]
    cp_arr = cp.zeros_like(L)
    dp_arr = cp.zeros_like(D)

    cp_arr[:, 0] = U[:, 0] / D[:, 0]
    dp_arr[:, 0] = B[:, 0] / D[:, 0]

    for j in range(1, n - 1):
        denom    = D[:, j] - L[:, j-1] * cp_arr[:, j-1]
        cp_arr[:, j] = U[:, j] / denom
        dp_arr[:, j] = (B[:, j] - L[:, j-1] * dp_arr[:, j-1]) / denom

    denom      = D[:, n-1] - L[:, n-2] * cp_arr[:, n-2]
    dp_arr[:, n-1] = (B[:, n-1] - L[:, n-2] * dp_arr[:, n-2]) / denom

    X = cp.zeros_like(D)
    X[:, n-1] = dp_arr[:, n-1]
    for j in range(n - 2, -1, -1):
        X[:, j] = dp_arr[:, j] - cp_arr[:, j] * X[:, j+1]
    return X

def _rhs_times_P_gpu(P, kf_eff, kb, h, theta):
    alpha  = (1.0 - theta) * h
    diag_K = -(kf_eff + kb)
    R = P + alpha * diag_K * P
    R[:, :-1] += alpha * kb[:, 1:]  * P[:, 1:]
    R[:, 1:]  += alpha * kf_eff[:, :-1] * P[:, :-1]
    return R

def _lhs_matrices_gpu(kf_eff, kb, h, theta):
    diag_K = -(kf_eff + kb)
    D = 1.0 - theta * h * diag_K
    U = -theta * h * kb[:, 1:]
    L = -theta * h * kf_eff[:, :-1]
    return L, D, U

def cn_step_gpu(P, c_free, fwd, bck, h, theta=0.5):
    kf_eff          = fwd.copy()
    kf_eff[:, 0]   *= c_free
    RHS = _rhs_times_P_gpu(P, kf_eff, bck, h, theta)
    L, D, U = _lhs_matrices_gpu(kf_eff, bck, h, theta)
    return _tri_solve_gpu(L, D, U, RHS)

# ----------------------------------------------------------------------------
# 3. COMPETITIVE STEP
# ----------------------------------------------------------------------------
def competitive_step_gpu(P, C_total, fwd, bck, h, theta=0.5, S_site=S_DIPLOID):
    def g(cf):
        P_trial = cn_step_gpu(P, cf, fwd, bck, h, theta)
        C_bound = P_trial[:, 1:].sum()
        return float(cf - (C_total - float(C_bound)))

    try:
        c_free = brentq(g, 0.0, C_total, xtol=1e-14, rtol=1e-10, maxiter=100)
    except ValueError:
        g0, gC = g(0.0), g(C_total)
        c_free = 0.0 if abs(g0) < abs(gC) else C_total

    P_next = cn_step_gpu(P, c_free, fwd, bck, h, theta)
    return P_next, c_free

# ----------------------------------------------------------------------------
# 4. FULL MULTI-SITE SIMULATION
# ----------------------------------------------------------------------------
def simulate_gpu(C_total, site_mismatches, t_points, S_site=S_DIPLOID, theta=0.5):
    epsilon, base_fwd = load_parameters()
    M = len(site_mismatches)

    fwd_all = np.zeros((M, 22))
    bck_all = np.zeros((M, 22))
    for i, mm in enumerate(site_mismatches):
        fwd_all[i], bck_all[i] = get_rate_pair(epsilon, base_fwd, mm)

    # Move to GPU
    fwd_gpu = cp.array(fwd_all)
    bck_gpu = cp.array(bck_all)

    P_gpu = cp.zeros((M, 22))
    P_gpu[:, 0] = S_site

    T = len(t_points)
    P_hist  = np.zeros((T, M, 22))
    Cf_hist = np.zeros(T)
    clv_hist= np.zeros((T, M))

    # CPU side arrays to store history
    P_hist[0]  = P_gpu.get()
    Cf_hist[0] = C_total
    clv_hist[0]= 0.0

    for k in range(T - 1):
        h = t_points[k+1] - t_points[k]
        P_gpu, cf = competitive_step_gpu(P_gpu, C_total, fwd_gpu, bck_gpu, h, theta, S_site)
        P_cpu = P_gpu.get()
        P_hist[k+1]  = P_cpu
        Cf_hist[k+1] = cf
        clv_hist[k+1] = (S_site - P_cpu.sum(axis=1)) / S_site

    return P_hist, Cf_hist, clv_hist

# ----------------------------------------------------------------------------
# 5. ORIGINAL NATURE MODEL  (CPU, since it's just 1 site)
# ----------------------------------------------------------------------------
def nature_model(C_total, mismatch_positions, t_points):
    epsilon, base_fwd = load_parameters()
    fwd, bck = get_rate_pair(epsilon, base_fwd, mismatch_positions)

    diag = -(fwd + bck)
    K = np.diag(diag) + np.diag(bck[1:], k=1) + np.diag(fwd[:-1], k=-1)
    K[0, 0] *= C_total
    K[1, 0] *= C_total

    P0 = np.zeros(22); P0[0] = 1.0
    P_clv = np.zeros(len(t_points))
    for i, t in enumerate(t_points):
        if t > 0:
            Pt = linalg.expm(K * t) @ P0
            P_clv[i] = 1.0 - Pt.sum()
    return P_clv

def nature_off_mean(mismatches, C_total, t_pts):
    off_clvs = []
    for mm in mismatches[1:]:
        off_clvs.append(nature_model(C_total, mm, t_pts))
    return np.mean(off_clvs, axis=0) if off_clvs else np.zeros(len(t_pts))

# ----------------------------------------------------------------------------
# 6. POPULATION GENERATION
# ----------------------------------------------------------------------------
def compute_residence_time(mm_positions, epsilon, base_fwd, guide_length=20):
    fwd, bck = get_rate_pair(epsilon, base_fwd, mm_positions, guide_length)
    fwd_d = fwd.copy(); fwd_d[-1] = 0.0
    diag = -(fwd_d[1:] + bck[1:])
    K_sub = np.diag(diag) + np.diag(bck[2:], k=1) + np.diag(fwd_d[1:-1], k=-1)
    try:
        Minv = np.linalg.inv(-K_sub)
        tau  = Minv[:, 0].sum()
    except np.linalg.LinAlgError:
        tau = 1e10
    return tau

def generate_population(M, rng_seed=42):
    rng = np.random.default_rng(rng_seed)
    epsilon, base_fwd = load_parameters()

    mismatches     = [[]]
    classifications= ["on_target"]
    res_times      = [compute_residence_time([], epsilon, base_fwd)]

    n_sinks    = max(1, int(0.02 * M))
    n_relevant = max(1, int(0.18 * M))
    n_weak     = M - 1 - n_sinks - n_relevant

    for _ in range(n_sinks):
        k   = rng.integers(1, 3)
        mm  = list(rng.choice(range(15, 21), size=k, replace=False))
        mismatches.append(mm); classifications.append("sink")
        res_times.append(compute_residence_time(mm, epsilon, base_fwd))

    for _ in range(n_relevant):
        k   = rng.integers(3, 5)
        mm  = list(rng.choice(range(9, 21), size=k, replace=False))
        mismatches.append(mm); classifications.append("relevant")
        res_times.append(compute_residence_time(mm, epsilon, base_fwd))

    for _ in range(n_weak):
        if rng.random() < 0.6:
            ns  = rng.integers(1, 3)
            sm  = list(rng.choice(range(1, 9),  size=ns, replace=False))
            dm  = list(rng.choice(range(9, 21), size=3, replace=False))
            mm  = list(set(sm + dm))
        else:
            mm  = list(rng.choice(range(1, 21), size=5, replace=False))
        mismatches.append(mm); classifications.append("weak")
        res_times.append(compute_residence_time(mm, epsilon, base_fwd))

    return mismatches, classifications, np.array(res_times)

# ----------------------------------------------------------------------------
# 7. CAMPAIGN RUNNER
# ----------------------------------------------------------------------------
def run_campaign(quick=False):
    print("\n" + "="*60)
    print("PART 4: SIMULATION CAMPAIGN (GPU ACCELERATED)")
    print("="*60)

    if quick:
        M_values = [10, 100]
        C_values = [1.0, 10.0]
        t_pts = np.insert(np.logspace(-2, 2, 60), 0, 0.0)
    else:
        M_values = [10, 100, 1000, 10000]
        C_values = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
        t_pts = np.insert(np.logspace(-2, 4.5, 300), 0, 0.0)
    results = {}

    # Warmup GPU
    print("Warming up GPU...")
    _ = simulate_gpu(1.0, [[]], np.array([0.0, 1.0]))

    for M in M_values:
        print(f"\n  M = {M:>6} sites", end="", flush=True)
        t0 = time.time()

        mismatches, classes, res_times = generate_population(M)
        results[M] = {"classes": classes, "res_times": res_times.tolist(),
                      "C_sweeps": {}}

        for C in C_values:
            print(".", end="", flush=True)
            _, Cf_hist, clv_hist = simulate_gpu(C, mismatches, t_pts)

            on_clv  = clv_hist[:, 0]
            off_clv = clv_hist[:, 1:].mean(axis=1) if M > 1 else np.zeros_like(on_clv)
            spec    = on_clv / (off_clv + 1e-12)

            results[M]["C_sweeps"][C] = {
                "C_free"   : Cf_hist.tolist(),
                "on_clv"   : on_clv.tolist(),
                "off_clv"  : off_clv.tolist(),
                "spec"     : spec.tolist(),
            }

        elapsed = time.time() - t0
        print(f"  ({elapsed:.1f} s)")

    results["t_pts"] = t_pts.tolist()
    return results

# ----------------------------------------------------------------------------
# 8. FIGURES & AUDIT
# ----------------------------------------------------------------------------
def make_figures(results):
    print("\n" + "="*60)
    print("PART 8: GENERATING FIGURES")
    print("="*60)
    # Re-using the same figure logic ...
    # (Just basic figures to save time / ensure completeness)
    t_pts  = np.array(results["t_pts"])
    M_vals = sorted(k for k in results.keys() if isinstance(k, int))
    C_vals = sorted(results[M_vals[0]]["C_sweeps"].keys())
    M_colors = {10:"#74c476", 100:"#41ab5d", 1000:"#238b45", 10000:"#005a32"}
    reference_c = 1.0 if 1.0 in results[M_vals[0]]["C_sweeps"] else C_vals[0]

    # Fig 1
    fig, ax = plt.subplots(figsize=(8, 5))
    nat_on_1nM = nature_model(1.0, [], t_pts)
    ax.plot(t_pts, nat_on_1nM, "k--", lw=2.5, label="Nature model")
    for M in M_vals:
        on_clv = np.array(results[M]["C_sweeps"][reference_c]["on_clv"])
        ax.plot(t_pts, on_clv, lw=2, color=M_colors[M], label=f"M={M}")
    ax.set_xscale("log"); ax.set_xlabel("Time (s)"); ax.set_ylabel("On-Target Cleaved Fraction")
    ax.legend(fontsize=8); ax.grid(True, which="both", ls="--", alpha=0.4)
    plt.tight_layout(); plt.savefig(os.path.join(ARTIFACTS_DIR, "fig1_on_target_cleavage.png"), dpi=300); plt.close()

    # Just creating this script to run the solver, saving the results
    with open(os.path.join(ARTIFACTS_DIR, "results_cuda.json"), "w") as f:
        json.dump(results, f)
    print("Results saved and Figure 1 generated.")

def run_audit(results):
    M_vals = sorted(k for k in results.keys() if isinstance(k, int))
    C_vals = sorted(results[M_vals[0]]["C_sweeps"].keys())
    t_pts  = np.array(results["t_pts"])
    records = []
    for M in M_vals:
        mm_M, _, _ = generate_population(M)
        for C in C_vals:
            nat_on   = nature_model(C, [], t_pts)[-1]
            nat_off  = nature_off_mean(mm_M, C, t_pts)[-1]
            nat_spec = nat_on / (nat_off + 1e-12)
            comp_on   = results[M]["C_sweeps"][C]["on_clv"][-1]
            comp_off  = results[M]["C_sweeps"][C]["off_clv"][-1]
            comp_spec = comp_on / (comp_off + 1e-12)
            d_on   = (comp_on  - nat_on)  * 100.0
            d_spec = np.log2((comp_spec + 1e-15) / (nat_spec + 1e-15))
            records.append(dict(M=M, C=C, d_on_pp=d_on, d_spec_log2=d_spec))
            
    with open(os.path.join(ARTIFACTS_DIR, "audit_report_cuda.json"), "w") as f:
        json.dump(records, f, indent=2)
    print("Audit completed.")

def main():
    parser = argparse.ArgumentParser(description="GPU-accelerated competitive SpCas9 simulation")
    parser.add_argument("--quick", action="store_true", help="Run a smaller GPU smoke test campaign")
    args = parser.parse_args()

    gpu_name, device_count = _gpu_name()
    print(f"Detected {device_count} CUDA device(s); using GPU 0: {gpu_name}")

    t_wall = time.time()
    results = run_campaign(quick=args.quick)
    make_figures(results)
    run_audit(results)
    print(f"\nCOMPLETE - total wall time: {time.time()-t_wall:.1f} s")

if __name__ == "__main__":
    main()
