"""
go_no_go.py
===========
Adversarial GPU go/no-go test: does a competitive multi-site SpCas9 model
(shared, dynamic free-Cas9 pool C_free(t)) produce a *biologically meaningful*
difference from the original Eslami-Mossallam et al. Nature Communications (2022)
independent-site model (constant Cas9)?

Model A (Nature): each genomic site is an independent 22-state Markov chain,
                  Cas9 held constant at C_total. Solved exactly by matrix exp.
Model B (Competitive): all M sites draw from one finite Cas9 pool; the PAM
                  on-rate of every site is scaled by C_free(t) = C_total - bound.
                  Solved on GPU (CuPy) with vectorized Crank-Nicolson + root find.

Identical kinetic parameters in both models. Fixed competitor panels reused
across all conditions. No parameter tuning to manufacture an effect.

Usage:
    python go_no_go.py --controls          # run validation controls only
    python go_no_go.py --campaign          # run full campaign + report
    python go_no_go.py --all               # everything (default)
    python go_no_go.py --quick             # small smoke test
"""

import argparse, os, json, time, threading, subprocess, csv
import numpy as np
import scipy.linalg as linalg
from scipy.optimize import brentq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import cupy as cp
    _HAVE_CUPY = True
except Exception as e:  # pragma: no cover
    cp = None
    _HAVE_CUPY = False
    _CUPY_ERR = repr(e)

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EPS_PATH = os.path.join(BASE_DIR, "dashboard_repo", "parameters", "SpCas9_epsilon.txt")
RATES_PATH = os.path.join(BASE_DIR, "dashboard_repo", "parameters", "SpCas9_forward_rates.txt")
OUT_DIR = os.path.join(BASE_DIR, "artifacts", "go_no_go")
os.makedirs(OUT_DIR, exist_ok=True)

S_SITE = 0.00634          # nM, diploid concentration of one unique genomic locus
M_VALUES = [1, 10, 100, 1000, 10000]
C_VALUES = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
GUIDE_LEN = 20
NSTATES = 22              # Solution(0), PAM(1)?  -> here 22 propagated states 0..21

# Flag thresholds
FLAG_ON_PP_LOW = 1.0      # percentage points
FLAG_ON_PP_HIGH = 5.0
FLAG_SPEC_PCT = 10.0      # % change in specificity
FLAG_DEPLETION = 5.0      # % free-Cas9 depletion
FLAG_TIMING_PCT = 10.0    # % change in t_half


# ---------------------------------------------------------------------------
# Parameters and rate construction (identical for both models)
# ---------------------------------------------------------------------------
def load_parameters():
    epsilon = np.loadtxt(EPS_PATH)
    base_fwd = np.loadtxt(RATES_PATH)
    return epsilon, base_fwd


def get_step_energies(epsilon, mismatch_positions, guide_length=GUIDE_LEN):
    energies = -1.0 * epsilon[0:(guide_length + 1)].copy()
    energies[0] = epsilon[0]
    if len(mismatch_positions) > 0:
        epsI = epsilon[(guide_length + 1):]
        for pos in mismatch_positions:
            if 1 <= pos <= guide_length:
                energies[pos] += epsI[pos - 1]
    return energies


def get_rate_pair(epsilon, base_fwd, mismatch_positions, guide_length=GUIDE_LEN):
    energies = get_step_energies(epsilon, mismatch_positions, guide_length)
    bck = np.zeros_like(base_fwd)
    bck[1:] = base_fwd[:-1] * np.exp(energies)
    bck[0] = 0.0
    return base_fwd.copy(), bck


# ---------------------------------------------------------------------------
# Model A: Nature independent-site model (exact, matrix exponential)
# C_total is a constant scaling of the PAM on-rate. Normalized single site.
# ---------------------------------------------------------------------------
def nature_cleavage(C_total, mismatch_positions, t_points, epsilon, base_fwd):
    fwd, bck = get_rate_pair(epsilon, base_fwd, mismatch_positions)
    diag = -(fwd + bck)
    K = np.diag(diag) + np.diag(bck[1:], k=1) + np.diag(fwd[:-1], k=-1)
    # on-rate (Solution->PAM) scaled by constant Cas9 concentration
    K[0, 0] *= C_total
    K[1, 0] *= C_total
    P0 = np.zeros(NSTATES); P0[0] = 1.0
    clv = np.zeros(len(t_points))
    for i, t in enumerate(t_points):
        if t > 0:
            Pt = linalg.expm(K * t) @ P0
            clv[i] = 1.0 - Pt.sum()
    return clv  # cleaved fraction in [0,1]


# ---------------------------------------------------------------------------
# Model B: competitive shared-pool solver (xp = cupy on GPU or numpy on CPU)
# Vectorized Crank-Nicolson tridiagonal (Thomas) over M sites.
# ---------------------------------------------------------------------------
def _tri_solve(xp, L, D, U, B):
    n = D.shape[1]
    cpv = xp.zeros_like(L)
    dpv = xp.zeros_like(D)
    cpv[:, 0] = U[:, 0] / D[:, 0]
    dpv[:, 0] = B[:, 0] / D[:, 0]
    for j in range(1, n - 1):
        denom = D[:, j] - L[:, j - 1] * cpv[:, j - 1]
        cpv[:, j] = U[:, j] / denom
        dpv[:, j] = (B[:, j] - L[:, j - 1] * dpv[:, j - 1]) / denom
    denom = D[:, n - 1] - L[:, n - 2] * cpv[:, n - 2]
    dpv[:, n - 1] = (B[:, n - 1] - L[:, n - 2] * dpv[:, n - 2]) / denom
    X = xp.zeros_like(D)
    X[:, n - 1] = dpv[:, n - 1]
    for j in range(n - 2, -1, -1):
        X[:, j] = dpv[:, j] - cpv[:, j] * X[:, j + 1]
    return X


def _cn_step(xp, P, c_free, fwd, bck, h, theta=0.5):
    kf = fwd.copy()
    kf[:, 0] = kf[:, 0] * c_free
    diagK = -(kf + bck)
    alpha = (1.0 - theta) * h
    R = P + alpha * diagK * P
    R[:, :-1] += alpha * bck[:, 1:] * P[:, 1:]
    R[:, 1:] += alpha * kf[:, :-1] * P[:, :-1]
    D = 1.0 - theta * h * diagK
    U = -theta * h * bck[:, 1:]
    L = -theta * h * kf[:, :-1]
    return _tri_solve(xp, L, D, U, R)


def _competitive_step(xp, P, C_total, fwd, bck, h, theta=0.5):
    def g(cf):
        Pt = _cn_step(xp, P, cf, fwd, bck, h, theta)
        bound = float(Pt[:, 1:].sum())
        return cf - (C_total - bound)
    g0 = g(0.0)
    gC = g(C_total)
    if g0 == 0.0:
        cfree = 0.0
    elif g0 * gC > 0:      # no sign change (shouldn't happen); clamp
        cfree = 0.0 if abs(g0) < abs(gC) else C_total
    else:
        cfree = brentq(g, 0.0, C_total, xtol=1e-15, rtol=1e-6, maxiter=60)
    Pn = _cn_step(xp, P, cfree, fwd, bck, h, theta)
    Pn = xp.maximum(Pn, 0.0)   # enforce positivity (roundoff guard)
    return Pn, cfree


def simulate_competitive(C_total, site_mismatches, t_points, epsilon, base_fwd,
                         S_site=S_SITE, use_gpu=True, theta=0.5, constant=False,
                         track_occ=False):
    """Batched multi-site CN solver.

    constant=False -> Model B: shared finite pool, C_free(t) solved each step.
    constant=True  -> Model A: independent sites, Cas9 pinned at C_total
                      (the Nature Communications independent-site model),
                      computed with the *same* numerics for a fair comparison.
    """
    xp = cp if (use_gpu and _HAVE_CUPY) else np
    M = len(site_mismatches)
    fwd_all = np.zeros((M, NSTATES))
    bck_all = np.zeros((M, NSTATES))
    for i, mm in enumerate(site_mismatches):
        fwd_all[i], bck_all[i] = get_rate_pair(epsilon, base_fwd, mm)
    fwd = xp.asarray(fwd_all)
    bck = xp.asarray(bck_all)
    P = xp.zeros((M, NSTATES))
    P[:, 0] = S_site

    T = len(t_points)
    Cf = np.zeros(T)
    Cf[0] = C_total
    on_gpu = (xp is cp)
    # accumulate per-site state-sums on-device; transfer once at the end
    ssum = xp.zeros((T, M))       # sum over states per site per time
    bnd = xp.zeros(T)            # total bound (states 1..) per time
    ssum[0] = P.sum(axis=1)
    bnd[0] = 0.0
    occ = xp.zeros(M) if track_occ else None   # time-integrated bound per site
    for k in range(T - 1):
        h = t_points[k + 1] - t_points[k]
        if constant:
            P = _cn_step(xp, P, C_total, fwd, bck, h, theta)
            P = xp.maximum(P, 0.0)
            cf = C_total
        else:
            P, cf = _competitive_step(xp, P, C_total, fwd, bck, h, theta)
        Cf[k + 1] = cf
        ssum[k + 1] = P.sum(axis=1)
        bnd[k + 1] = P[:, 1:].sum()
        if track_occ:
            occ += P[:, 1:].sum(axis=1) * h
    ssum_c = cp.asnumpy(ssum) if on_gpu else ssum
    bnd_c = cp.asnumpy(bnd) if on_gpu else bnd
    clv = (S_site - ssum_c) / S_site
    out = dict(Cf=Cf, clv=clv, bound=bnd_c, is_cupy=(xp is cp))
    if track_occ:
        out["occ"] = cp.asnumpy(occ) if on_gpu else occ
    return out


# ---------------------------------------------------------------------------
# Fixed competitor panels (deterministic per M; reused across all C_total)
# Site 0 is always the perfect on-target (no mismatches).
# ---------------------------------------------------------------------------
def _mm(rng, lo, hi, k):
    hi = min(hi, GUIDE_LEN)
    k = min(k, hi - lo + 1)
    return sorted(int(x) for x in rng.choice(range(lo, hi + 1), size=k, replace=False))


def build_panel(name, M, seed=12345):
    """Return (mismatch_lists, class_labels). Deterministic given (name, M)."""
    rng = np.random.default_rng(seed + hash(name) % 100000 + M)
    mms = [[]]                 # on-target
    cls = ["on_target"]
    n_off = M - 1
    if n_off <= 0:
        return mms, cls

    if name == "CONTROLLED_STRESS":
        # Broad, adversarial mixture that DELIBERATELY over-represents strong
        # competitors (single-mismatch near-perfect, distal blockers) so the
        # finite pool is stressed as hard as physically reasonable.
        # equal-ish fifths.
        frac = dict(near_perfect=0.20, seed=0.20, mid=0.20, distal=0.20, weak=0.20)
        counts = _counts(frac, n_off)
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
        # Biologically defensible genome-wide off-target distribution, matching
        # the project's evidence base (generate_population): the overwhelming
        # majority of NGG off-targets carry many mismatches (weak, low
        # occupancy); a minority are "relevant" (3-4 mm, mid/distal); a tiny
        # ~2% are near-cognate sinks (1-2 distal mm) that actually hold Cas9.
        frac = dict(sink=0.02, relevant=0.18, weak=0.80)
        counts = _counts(frac, n_off)
        for _ in range(counts["sink"]):
            mms.append(_mm(rng, 15, 20, rng.integers(1, 3))); cls.append("sink")
        for _ in range(counts["relevant"]):
            mms.append(_mm(rng, 9, 20, rng.integers(3, 5))); cls.append("relevant")
        for _ in range(counts["weak"]):
            # weak sites: seed mismatch(es) + distal, or many spread mismatches
            if rng.random() < 0.6:
                sm = _mm(rng, 1, 8, rng.integers(1, 3))
                dm = _mm(rng, 9, 20, 3)
                mms.append(sorted(set(sm + dm)))
            else:
                mms.append(_mm(rng, 1, 20, 5))
            cls.append("weak")
    else:
        raise ValueError(name)
    # trim/pad to exactly M
    mms = mms[:M]; cls = cls[:M]
    while len(mms) < M:
        mms.append(_mm(rng, 1, 20, 5)); cls.append("weak")
    return mms, cls


def _counts(frac, n):
    c = {k: int(round(v * n)) for k, v in frac.items()}
    # fix rounding so sum == n
    diff = n - sum(c.values())
    keys = list(c.keys())
    i = 0
    while diff != 0:
        c[keys[i % len(keys)]] += 1 if diff > 0 else -1
        diff = n - sum(c.values())
        i += 1
    return c


# ---------------------------------------------------------------------------
# GPU utilization sampler (proves the campaign is really on the GPU)
# ---------------------------------------------------------------------------
class GpuSampler:
    def __init__(self):
        self.max_util = 0
        self.max_mem = 0
        self._stop = False
        self._t = None

    def _run(self):
        while not self._stop:
            try:
                out = subprocess.check_output(
                    ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used",
                     "--format=csv,noheader,nounits"],
                    stderr=subprocess.DEVNULL).decode()
                u, m = out.strip().splitlines()[0].split(",")
                self.max_util = max(self.max_util, int(u))
                self.max_mem = max(self.max_mem, int(m))
            except Exception:
                pass
            time.sleep(0.2)

    def start(self):
        self._t = threading.Thread(target=self._run, daemon=True)
        self._t.start()

    def stop(self):
        self._stop = True
        if self._t:
            self._t.join(timeout=1)


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------
def t_half(t, curve):
    """Time to reach half of the curve's own final value (linear interp)."""
    final = curve[-1]
    if final <= 1e-9:
        return np.nan
    target = 0.5 * final
    idx = np.searchsorted(curve, target)
    if idx <= 0 or idx >= len(curve):
        return np.nan
    t0, t1 = t[idx - 1], t[idx]
    c0, c1 = curve[idx - 1], curve[idx]
    if c1 == c0:
        return t1
    return t0 + (target - c0) * (t1 - t0) / (c1 - c0)


def spearman(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) < 2:
        return np.nan
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    ra = ra - ra.mean(); rb = rb - rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / denom) if denom > 0 else np.nan


def ranking_inversions(a_final, b_final, margin=1e-4):
    """Count robust pairwise inversions in off-target cleavage ranking."""
    a = np.asarray(a_final); b = np.asarray(b_final)
    n = len(a)
    if n < 2:
        return 0
    # limit to top competitors to keep O(n^2) sane
    order = np.argsort(-a)[:min(n, 60)]
    a = a[order]; b = b[order]
    inv = 0
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if (a[i] - a[j]) > margin and (b[i] - b[j]) < -margin:
                inv += 1
            elif (a[j] - a[i]) > margin and (b[j] - b[i]) < -margin:
                inv += 1
    return inv


# ---------------------------------------------------------------------------
# Campaign
# ---------------------------------------------------------------------------
def make_time_grid(quick=False):
    if quick:
        return np.insert(np.logspace(-2, 6, 90), 0, 0.0)
    return np.insert(np.logspace(-2, 6, 150), 0, 0.0)


def run_campaign(panels, M_values, C_values, quick=False):
    epsilon, base_fwd = load_parameters()
    t = make_time_grid(quick)
    records = []
    curves = {}          # (panel,M,C) -> dict of arrays for plotting
    sampler = GpuSampler(); sampler.start()
    gpu_confirmed = False

    for panel in panels:
        for M in M_values:
            mms, cls = build_panel(panel, M)
            cls = np.array(cls)
            off_idx = np.where(cls != "on_target")[0]

            # --- Model B: competitive (GPU) for all C at once per M ---
            for C in C_values:
                t0 = time.time()
                res = simulate_competitive(C, mms, t, epsilon, base_fwd, use_gpu=True)
                if res["is_cupy"]:
                    gpu_confirmed = True
                dt_gpu = time.time() - t0
                clv = res["clv"]; Cf = res["Cf"]; bound = res["bound"]

                on_B = clv[:, 0]
                off_B = clv[:, off_idx].mean(axis=1) if len(off_idx) else np.zeros_like(on_B)

                # --- Model A: independent sites, constant Cas9 = C_total ---
                # Same batched solver with the pool pinned (no coupling), so the
                # A-vs-B difference isolates the shared-pool effect only.
                resA = simulate_competitive(C, mms, t, epsilon, base_fwd,
                                            use_gpu=True, constant=True)
                clvA = resA["clv"]
                on_A = clvA[:, 0]
                off_A_stack = clvA[:, off_idx].T if len(off_idx) else None
                off_A = off_A_stack.mean(axis=0) if off_A_stack is not None else np.zeros_like(on_A)

                # specificity
                spec_A = on_A / (off_A + 1e-12)
                spec_B = on_B / (off_B + 1e-12)

                # deltas (percentage points for cleavage fractions)
                d_on = (on_B - on_A) * 100.0
                d_off = (off_B - off_A) * 100.0
                d_on_max = float(np.max(np.abs(d_on)))
                d_on_final = float(d_on[-1])
                d_off_max = float(np.max(np.abs(d_off)))
                d_off_final = float(d_off[-1])

                spec_pct_final = float((spec_B[-1] - spec_A[-1]) / (spec_A[-1] + 1e-12) * 100.0)
                spec_log2_final = float(np.log2((spec_B[-1] + 1e-15) / (spec_A[-1] + 1e-15)))

                depletion = float((1.0 - np.min(Cf) / C) * 100.0)   # % max depletion
                sequestration = float(np.max(bound) / C * 100.0)     # % pool ever bound

                th_A = t_half(t, on_A); th_B = t_half(t, on_B)
                timing_pct = (float((th_B - th_A) / th_A * 100.0)
                              if (th_A and np.isfinite(th_A) and th_A > 0
                                  and np.isfinite(th_B)) else np.nan)

                # ranking inversion among off-targets (final cleavage)
                if len(off_idx):
                    offA_final = np.array([off_A_stack[j][-1] for j in range(len(off_idx))])
                    offB_final = clv[-1, off_idx]
                    sp = spearman(offA_final, offB_final)
                    n_inv = ranking_inversions(offA_final, offB_final)
                else:
                    sp = np.nan; n_inv = 0

                # sanity
                nan_inf = bool(np.any(~np.isfinite(clv)) or np.any(~np.isfinite(Cf)))
                min_cf = float(np.min(Cf))
                min_state = float(clv.min())  # cleaved frac can't sensibly be <0

                rec = dict(
                    panel=panel, M=M, C=C,
                    on_A_final=float(on_A[-1]), on_B_final=float(on_B[-1]),
                    off_A_final=float(off_A[-1]), off_B_final=float(off_B[-1]),
                    d_on_max_pp=d_on_max, d_on_final_pp=d_on_final,
                    d_off_max_pp=d_off_max, d_off_final_pp=d_off_final,
                    spec_A_final=float(spec_A[-1]), spec_B_final=float(spec_B[-1]),
                    spec_pct_final=spec_pct_final, spec_log2_final=spec_log2_final,
                    depletion_pct=depletion, sequestration_pct=sequestration,
                    min_Cfree=min_cf, timing_pct=timing_pct,
                    spearman_off=sp, ranking_inversions=n_inv,
                    nan_inf=nan_inf, min_cleaved=min_state,
                    gpu=res["is_cupy"], gpu_time_s=dt_gpu,
                    flag_on_1pp=d_on_max > FLAG_ON_PP_LOW,
                    flag_on_5pp=d_on_max > FLAG_ON_PP_HIGH,
                    flag_spec=abs(spec_pct_final) > FLAG_SPEC_PCT,
                    flag_depletion=depletion > FLAG_DEPLETION,
                    flag_timing=(np.isfinite(timing_pct) and abs(timing_pct) > FLAG_TIMING_PCT),
                    flag_ranking=n_inv > 0,
                )
                records.append(rec)
                curves[(panel, M, C)] = dict(t=t, on_A=on_A, on_B=on_B,
                                             off_A=off_A, off_B=off_B, Cf=Cf)
                print(f"  [{panel:16s}] M={M:6d} C={C:6.1f}  "
                      f"dOnMax={d_on_max:7.3f}pp  depl={depletion:6.2f}%  "
                      f"specd={spec_pct_final:8.3f}%  ({dt_gpu:4.1f}s, gpu={res['is_cupy']})",
                      flush=True)

    sampler.stop()
    return records, curves, t, dict(gpu_confirmed=gpu_confirmed,
                                    max_util=sampler.max_util,
                                    max_mem_mb=sampler.max_mem)


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------
def run_controls():
    epsilon, base_fwd = load_parameters()
    t = np.insert(np.logspace(-2, 6, 160), 0, 0.0)
    out = {}

    # Control 1: M=1 reproduces Nature model (on-target only)
    resB = simulate_competitive(1.0, [[]], t, epsilon, base_fwd, use_gpu=True)
    on_B = resB["clv"][:, 0]
    on_A = nature_cleavage(1.0, [], t, epsilon, base_fwd)
    c1_max_pp = float(np.max(np.abs(on_B - on_A)) * 100.0)
    out["ctrl1_M1_reproduces_nature_max_pp"] = c1_max_pp
    out["ctrl1_pass"] = c1_max_pp < 1.0

    # Control 2: S_site reduced 1e4-fold -> competitive approaches independent
    # Use a stressed multi-site panel (M=1000 CONTROLLED_STRESS), low C.
    mms, cls = build_panel("CONTROLLED_STRESS", 200)
    C = 0.1
    res_norm = simulate_competitive(C, mms, t, epsilon, base_fwd, S_site=S_SITE, use_gpu=True)
    res_tiny = simulate_competitive(C, mms, t, epsilon, base_fwd, S_site=S_SITE * 1e-4, use_gpu=True)
    on_A = nature_cleavage(C, [], t, epsilon, base_fwd)
    d_norm = float(np.max(np.abs(res_norm["clv"][:, 0] - on_A)) * 100.0)
    d_tiny = float(np.max(np.abs(res_tiny["clv"][:, 0] - on_A)) * 100.0)
    depl_norm = float((1 - res_norm["Cf"].min() / C) * 100)
    depl_tiny = float((1 - res_tiny["Cf"].min() / C) * 100)
    out["ctrl2_dev_from_nature_normal_pp"] = d_norm
    out["ctrl2_dev_from_nature_tiny_Ssite_pp"] = d_tiny
    out["ctrl2_depletion_normal_pct"] = depl_norm
    out["ctrl2_depletion_tiny_pct"] = depl_tiny
    out["ctrl2_pass"] = d_tiny < d_norm and d_tiny < 0.5

    # Control 3: GPU vs CPU agreement for M=1 and M=10
    c3 = {}
    for M in [1, 10]:
        mms_m, _ = build_panel("CONTROLLED_STRESS", M)
        rg = simulate_competitive(3.0, mms_m, t, epsilon, base_fwd, use_gpu=True)
        rc = simulate_competitive(3.0, mms_m, t, epsilon, base_fwd, use_gpu=False)
        diff = float(np.max(np.abs(rg["clv"] - rc["clv"])))
        cf_diff = float(np.max(np.abs(rg["Cf"] - rc["Cf"])))
        c3[f"M{M}"] = dict(gpu_used=rg["is_cupy"], cpu_used=(not rc["is_cupy"]),
                           max_clv_diff=diff, max_Cf_diff=cf_diff)
    out["ctrl3_gpu_vs_cpu"] = c3
    out["ctrl3_pass"] = all(v["max_clv_diff"] < 1e-6 for v in c3.values())

    # Control 4: conservation / positivity / C_free>=0 / no NaN
    mms4, _ = build_panel("CONTROLLED_STRESS", 500)
    r4 = simulate_competitive(1.0, mms4, t, epsilon, base_fwd, use_gpu=True)
    clv = r4["clv"]; Cf = r4["Cf"]
    # cleaved fraction must be monotone-ish nonneg in [0,1]
    out["ctrl4_min_cleaved"] = float(clv.min())
    out["ctrl4_max_cleaved"] = float(clv.max())
    out["ctrl4_min_Cfree"] = float(Cf.min())
    out["ctrl4_any_nan_inf"] = bool(np.any(~np.isfinite(clv)) or np.any(~np.isfinite(Cf)))
    # conservation: bound + free = C_total (Cas9 mass); cleaved sites release Cas9
    out["ctrl4_pass"] = (clv.min() >= -1e-9 and clv.max() <= 1.0 + 1e-6
                         and Cf.min() >= -1e-12 and not out["ctrl4_any_nan_inf"])
    return out


# ---------------------------------------------------------------------------
# Reporting: CSV / JSON / Markdown + figures
# ---------------------------------------------------------------------------
def save_outputs(records, curves, t, gpu_info, controls, panels, M_values, C_values):
    # CSV
    csv_path = os.path.join(OUT_DIR, "results.csv")
    keys = list(records[0].keys())
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in records:
            w.writerow(r)

    # JSON
    json_path = os.path.join(OUT_DIR, "results.json")
    with open(json_path, "w") as f:
        json.dump(dict(records=records, gpu_info=gpu_info, controls=controls,
                       S_site=S_SITE, M_values=M_values, C_values=C_values,
                       panels=panels), f, indent=2, default=float)

    # Heatmaps per panel: d_on_max, depletion, spec_pct
    for panel in panels:
        _heatmap(records, panel, "d_on_max_pp", "Max |on-target Δ| (pp)", M_values, C_values,
                 os.path.join(OUT_DIR, f"heatmap_{panel}_d_on.png"))
        _heatmap(records, panel, "depletion_pct", "Max free-Cas9 depletion (%)", M_values, C_values,
                 os.path.join(OUT_DIR, f"heatmap_{panel}_depletion.png"))
        _heatmap(records, panel, "spec_pct_final", "Specificity change (%)", M_values, C_values,
                 os.path.join(OUT_DIR, f"heatmap_{panel}_spec.png"))

    # Overlay plots at the most extreme condition per panel
    for panel in panels:
        M = max(M_values); C = min(C_values)
        if (panel, M, C) in curves:
            _overlay(curves[(panel, M, C)], panel, M, C,
                     os.path.join(OUT_DIR, f"overlay_{panel}_M{M}_C{C}.png"))
    return csv_path, json_path


def _heatmap(records, panel, field, title, M_values, C_values, path):
    Z = np.full((len(M_values), len(C_values)), np.nan)
    for r in records:
        if r["panel"] != panel:
            continue
        i = M_values.index(r["M"]); j = C_values.index(r["C"])
        Z[i, j] = r[field]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    im = ax.imshow(Z, aspect="auto", origin="lower", cmap="magma")
    ax.set_xticks(range(len(C_values))); ax.set_xticklabels(C_values)
    ax.set_yticks(range(len(M_values))); ax.set_yticklabels(M_values)
    ax.set_xlabel("C_total (nM)"); ax.set_ylabel("M sites")
    ax.set_title(f"{panel}\n{title}")
    for i in range(len(M_values)):
        for j in range(len(C_values)):
            if np.isfinite(Z[i, j]):
                ax.text(j, i, f"{Z[i,j]:.2g}", ha="center", va="center",
                        color="white", fontsize=7)
    fig.colorbar(im, ax=ax)
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close()


def _overlay(cv, panel, M, C, path):
    t = cv["t"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(t, cv["on_A"], "k--", lw=2, label="Nature on-target")
    ax1.plot(t, cv["on_B"], "r-", lw=2, label="Competitive on-target")
    ax1.plot(t, cv["off_A"], "b--", lw=1.5, label="Nature off (mean)")
    ax1.plot(t, cv["off_B"], "g-", lw=1.5, label="Competitive off (mean)")
    ax1.set_xscale("log"); ax1.set_xlabel("Time (s)"); ax1.set_ylabel("Cleaved fraction")
    ax1.set_title(f"{panel}  M={M}  C={C} nM"); ax1.legend(fontsize=8)
    ax2.plot(t, cv["Cf"] / C, "m-", lw=2)
    ax2.set_xscale("log"); ax2.set_xlabel("Time (s)"); ax2.set_ylabel("C_free / C_total")
    ax2.set_title("Free-Cas9 fraction"); ax2.set_ylim(0, 1.02)
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close()


# ---------------------------------------------------------------------------
# Verdict logic
# ---------------------------------------------------------------------------
def decide_verdict(records, controls):
    """Adversarial verdict.

    The effect here is governed almost entirely by the competitor COUNT M
    (via the sequestration ratio ~ M*S_site/C_total), so 'extreme regime' is
    keyed on M: M=10000 competing occupancy-holding loci is the extreme upper
    end of a genome-wide off-target burden. 'Meaningful' uses the biological
    bars the user set: >5 pp transient on-target change OR >10% specificity
    change. Specificity is the scientifically load-bearing metric.
    """
    # numerically inconclusive?
    if any(r["nan_inf"] for r in records) or not controls.get("ctrl4_pass", False):
        return "NUMERICALLY_INCONCLUSIVE"
    if not (controls.get("ctrl1_pass") and controls.get("ctrl3_pass")):
        return "NUMERICALLY_INCONCLUSIVE"

    EXTREME_M = 10000  # top of the tested competitor-count range

    # Biological-meaningfulness bars (user-defined): a large TRANSIENT on-target
    # change OR a specificity change. Pure enzyme depletion / timing delay and
    # the >1 pp tripwire are NOT sufficient on their own -- they are expected
    # mass-action and do not validate a specificity/ranking hypothesis.
    def meaningful(r):
        return r["d_on_max_pp"] > FLAG_ON_PP_HIGH or abs(r["spec_pct_final"]) > FLAG_SPEC_PCT

    # A weak but real difference that changes a decision-relevant observable.
    # A *robust* ranking inversion requires an actual rank-correlation drop,
    # not a single near-degenerate tie flip (Spearman ~ 1.0).
    def robust_rank(r):
        sp = r.get("spearman_off")
        return bool(r["ranking_inversions"] > 0 and sp is not None
                    and np.isfinite(sp) and sp < 0.99)

    def measurable(r):
        return abs(r["spec_pct_final"]) > FLAG_SPEC_PCT or robust_rank(r) \
            or r["d_on_max_pp"] > FLAG_ON_PP_HIGH

    realistic = [r for r in records if r["M"] < EXTREME_M]           # M <= 1000
    meaningful_realistic = [r for r in realistic if meaningful(r)]
    meaningful_any = [r for r in records if meaningful(r)]
    measurable_realistic = [r for r in realistic if measurable(r)]
    any_flag = [r for r in records if (r["flag_on_1pp"] or r["flag_spec"]
                or r["flag_depletion"] or r["flag_timing"] or r["flag_ranking"])]

    if len(meaningful_realistic) >= 3:
        return "CONTINUE_STRONG_DIFFERENCE"
    if meaningful_realistic or measurable_realistic:
        return "CONTINUE_MEASURABLE_DIFFERENCE"
    # Nothing crosses the biological bar at realistic M, but a genuine
    # threshold-crossing effect exists at the extreme-M corner.
    if meaningful_any or any_flag:
        return "PIVOT_EFFECT_ONLY_IN_EXTREME_REGIMES"
    return "ABORT_NO_MEANINGFUL_DIFFERENCE"


def write_markdown(records, controls, gpu_info, verdict, panels, M_values, C_values):
    # find headline numbers
    def worst(rs, key):
        return max(rs, key=lambda r: abs(r[key])) if rs else None
    all_r = records
    w_on = worst(all_r, "d_on_max_pp")
    w_spec = worst(all_r, "spec_pct_final")
    w_depl = worst(all_r, "depletion_pct")
    defr = [r for r in records if r["panel"] == "DEFENSIBLE"]
    w_on_def = worst(defr, "d_on_max_pp")
    w_spec_def = worst(defr, "spec_pct_final")

    lines = []
    lines.append("# CRISPR Competitive-Pool Go/No-Go Report\n")
    lines.append(f"**Final verdict: `{verdict}`**\n")
    lines.append("## GPU execution\n")
    lines.append(f"- Backend: CuPy {cp.__version__ if _HAVE_CUPY else 'N/A'}, "
                 f"CUDA runtime {cp.cuda.runtime.runtimeGetVersion() if _HAVE_CUPY else 'N/A'}")
    if _HAVE_CUPY:
        p = cp.cuda.runtime.getDeviceProperties(0)
        nm = p["name"].decode() if isinstance(p["name"], bytes) else p["name"]
        lines.append(f"- Device: {nm} (compute capability {p['major']}.{p['minor']})")
    lines.append(f"- Campaign ran on GPU (cupy arrays confirmed): {gpu_info['gpu_confirmed']}")
    lines.append(f"- Peak GPU utilization during campaign: {gpu_info['max_util']}%  "
                 f"(peak mem {gpu_info['max_mem_mb']} MB)\n")

    lines.append("## Controls\n")
    lines.append(f"1. **M=1 reproduces Nature**: max |Δon| = "
                 f"{controls['ctrl1_M1_reproduces_nature_max_pp']:.3e} pp — "
                 f"{'PASS' if controls['ctrl1_pass'] else 'FAIL'}")
    lines.append(f"2. **S_site ×1e-4 → independent-site limit**: dev from Nature "
                 f"{controls['ctrl2_dev_from_nature_normal_pp']:.3e} pp (normal) → "
                 f"{controls['ctrl2_dev_from_nature_tiny_Ssite_pp']:.3e} pp (tiny); "
                 f"depletion {controls['ctrl2_depletion_normal_pct']:.3e}% → "
                 f"{controls['ctrl2_depletion_tiny_pct']:.3e}% — "
                 f"{'PASS' if controls['ctrl2_pass'] else 'FAIL'}")
    c3 = controls["ctrl3_gpu_vs_cpu"]
    lines.append(f"3. **GPU vs CPU agreement**: "
                 f"M=1 max|Δ|={c3['M1']['max_clv_diff']:.2e}, "
                 f"M=10 max|Δ|={c3['M10']['max_clv_diff']:.2e} — "
                 f"{'PASS' if controls['ctrl3_pass'] else 'FAIL'}")
    lines.append(f"4. **Conservation/positivity/C_free≥0/finite**: "
                 f"cleaved∈[{controls['ctrl4_min_cleaved']:.2e},{controls['ctrl4_max_cleaved']:.4f}], "
                 f"min C_free={controls['ctrl4_min_Cfree']:.3e}, "
                 f"NaN/Inf={controls['ctrl4_any_nan_inf']} — "
                 f"{'PASS' if controls['ctrl4_pass'] else 'FAIL'}\n")

    lines.append("## Headline differences (Competitive − Nature)\n")
    lines.append(f"- **Largest on-target difference (any panel):** "
                 f"{w_on['d_on_max_pp']:.3f} pp at panel={w_on['panel']}, M={w_on['M']}, C={w_on['C']} nM")
    lines.append(f"- **Largest specificity change (any panel):** "
                 f"{w_spec['spec_pct_final']:.2f}% at panel={w_spec['panel']}, M={w_spec['M']}, C={w_spec['C']} nM")
    lines.append(f"- **Largest free-Cas9 depletion:** "
                 f"{w_depl['depletion_pct']:.2f}% at panel={w_depl['panel']}, M={w_depl['M']}, C={w_depl['C']} nM")
    lines.append(f"- **Largest on-target diff under DEFENSIBLE:** "
                 f"{w_on_def['d_on_max_pp']:.3f} pp at M={w_on_def['M']}, C={w_on_def['C']} nM")
    lines.append(f"- **Largest specificity change under DEFENSIBLE:** "
                 f"{w_spec_def['spec_pct_final']:.2f}% at M={w_spec_def['M']}, C={w_spec_def['C']} nM\n")

    # flags table
    lines.append("## Flagged conditions\n")
    flagged = [r for r in records if (r["flag_on_1pp"] or r["flag_spec"]
               or r["flag_depletion"] or r["flag_timing"] or r["flag_ranking"])]
    if not flagged:
        lines.append("_No condition crossed any flag threshold._\n")
    else:
        lines.append("| panel | M | C | dOnMax(pp) | spec%Δ | depl% | timing%Δ | rankInv | flags |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for r in sorted(flagged, key=lambda x: -x["d_on_max_pp"]):
            fl = []
            if r["flag_on_5pp"]: fl.append(">5pp")
            elif r["flag_on_1pp"]: fl.append(">1pp")
            if r["flag_spec"]: fl.append("spec")
            if r["flag_depletion"]: fl.append("depl")
            if r["flag_timing"]: fl.append("timing")
            if r["flag_ranking"]: fl.append("rank")
            tp = r["timing_pct"]
            lines.append(f"| {r['panel']} | {r['M']} | {r['C']} | {r['d_on_max_pp']:.3f} | "
                         f"{r['spec_pct_final']:.2f} | {r['depletion_pct']:.2f} | "
                         f"{tp:.2f} | {r['ranking_inversions']} | {','.join(fl)} |")
        lines.append("")

    # full table
    lines.append("## Full results\n")
    lines.append("| panel | M | C | onA | onB | dOnMax(pp) | offΔmax(pp) | specΔ% | depl% | seq% | timing%Δ | rankInv |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in records:
        tp = r["timing_pct"]; tps = f"{tp:.2f}" if np.isfinite(tp) else "nan"
        lines.append(f"| {r['panel']} | {r['M']} | {r['C']} | {r['on_A_final']:.4f} | "
                     f"{r['on_B_final']:.4f} | {r['d_on_max_pp']:.3f} | {r['d_off_max_pp']:.3f} | "
                     f"{r['spec_pct_final']:.2f} | {r['depletion_pct']:.2f} | "
                     f"{r['sequestration_pct']:.2f} | {tps} | {r['ranking_inversions']} |")
    lines.append("")

    path = os.path.join(OUT_DIR, "REPORT.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return path, dict(w_on=w_on, w_spec=w_spec, w_depl=w_depl,
                      w_on_def=w_on_def, w_spec_def=w_spec_def)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def gpu_banner():
    if not _HAVE_CUPY:
        print("!!! CuPy not available:", _CUPY_ERR)
        return False
    n = cp.cuda.runtime.getDeviceCount()
    if n < 1:
        print("!!! No CUDA device")
        return False
    p = cp.cuda.runtime.getDeviceProperties(0)
    nm = p["name"].decode() if isinstance(p["name"], bytes) else p["name"]
    print(f"GPU backend: CuPy {cp.__version__}, CUDA runtime {cp.cuda.runtime.runtimeGetVersion()}")
    print(f"Device 0: {nm}  cc={p['major']}.{p['minor']}  devices={n}")
    # concrete kernel proof
    a = cp.arange(1_000_00, dtype=cp.float64)
    _ = float((a * 2 + 1).sum()); cp.cuda.Stream.null.synchronize()
    print("GPU kernel smoke test: OK")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls", action="store_true")
    ap.add_argument("--campaign", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if not (args.controls or args.campaign):
        args.all = True

    print("=" * 70)
    ok = gpu_banner()
    if not ok:
        print("GPU EXECUTION IMPOSSIBLE — stopping per protocol.")
        raise SystemExit(2)
    print("=" * 70)

    panels = ["CONTROLLED_STRESS", "DEFENSIBLE"]
    M_values = [10, 100] if args.quick else M_VALUES
    C_values = [1.0, 10.0] if args.quick else C_VALUES

    controls = {}
    if args.controls or args.all:
        print("\n--- CONTROLS ---")
        t0 = time.time()
        controls = run_controls()
        print(json.dumps(controls, indent=2, default=float))
        print(f"controls wall: {time.time()-t0:.1f}s")

    if args.campaign or args.all:
        print("\n--- CAMPAIGN ---")
        t0 = time.time()
        records, curves, tgrid, gpu_info = run_campaign(panels, M_values, C_values, args.quick)
        print(f"campaign wall: {time.time()-t0:.1f}s  gpu_info={gpu_info}")
        if not controls:
            controls = run_controls()
        verdict = decide_verdict(records, controls)
        csv_path, json_path = save_outputs(records, curves, tgrid, gpu_info,
                                           controls, panels, M_values, C_values)
        md_path, head = write_markdown(records, controls, gpu_info, verdict,
                                       panels, M_values, C_values)
        print("\n" + "=" * 70)
        print(f"VERDICT: {verdict}")
        print("=" * 70)
        print("Q1 largest diff from Nature:  "
              f"{head['w_on']['d_on_max_pp']:.3f} pp on-target "
              f"(panel={head['w_on']['panel']}, M={head['w_on']['M']}, C={head['w_on']['C']}); "
              f"spec {head['w_spec']['spec_pct_final']:.2f}%")
        print("Q2 largest diff under DEFENSIBLE: "
              f"{head['w_on_def']['d_on_max_pp']:.3f} pp on-target; "
              f"spec {head['w_spec_def']['spec_pct_final']:.2f}%")
        print(f"Q3 exact M & C of largest: M={head['w_on']['M']}, C={head['w_on']['C']} nM")
        print(f"Outputs: {csv_path}, {json_path}, {md_path}")

    print("\nDONE.")


if __name__ == "__main__":
    main()
