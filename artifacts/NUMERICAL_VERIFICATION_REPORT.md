# NUMERICAL VERIFICATION REPORT
## Competitive Cas9 Multi-Site Kinetic Model

**Date:** June 21, 2026  
**Project:** Biophysical audit of finite Cas9 competition in CRISPR-Cas9 targeting  
**Status:** NUMERICAL VERIFICATION PHASE (Steps 1–5 complete)

---

## EXECUTIVE SUMMARY

The competitive Cas9 solver has been implemented and partially validated. The solver passes the critical M=1 validation test against the original Nature Communications (2022) reference model with acceptable error (max 2.17e-7, requirement < 1e-6). Mass conservation is excellent (error < 1e-10). However, **state-concentration normalization shows systematic error accumulation (~7.6% relative loss)**, exceeding the strict target of < 1e-10. This blocks progression to the full campaign until resolved.

---

## WORK COMPLETED

### 1. Solver Implementation

**File:** `simulate_competitive.py`

**Changes:**
- Replaced unstable vectorized Thomas tridiagonal algorithm with robust dense batched solver (`np.linalg.solve`)
- Removed Unicode arrows that triggered encoding issues
- Solver now constructs full 22×22 coefficient matrices per site and solves via LAPACK
- **Rationale:** Thomas algorithm exhibited catastrophic numerical instability (errors ~4e14) due to ill-conditioned matrices at large time steps; dense solve trades efficiency for correctness

**Code Structure:**
```python
def _tri_solve(L, D, U, B):
    # Builds full tridiagonal matrices and uses np.linalg.solve
    # Input: L (M, 21), D (M, 22), U (M, 21), B (M, 22)
    # Output: X (M, 22)
```

### 2. Numerical Verification Framework

**Implemented five diagnostic steps:**

1. **M=1 Validation:** Compare competitive solver against matrix exponentiation reference
2. **Positivity Test:** Verify all state concentrations remain ≥ 0
3. **Normalization Test:** Check Σ_states per site = S_site throughout simulation
4. **Mass Conservation Test:** Verify C_total = C_free + C_bound throughout
5. **M=10 Stability Test:** Stress test with 1 on-target + 9 random off-targets

**Input Parameters:**
- Site concentration: S_site = 0.00634 nM (corrected diploid locus concentration)
- Cas9 concentration: C_total = 1.0 nM
- Time grid: 300 log-spaced points (validation), 200 points (M=10)
- Time range: 10⁻² to 10³ seconds
- Integrator: Crank-Nicolson (θ=0.5, implicit)

---

## NUMERICAL VERIFICATION RESULTS

### STEP 1: M=1 VALIDATION

**Test:** Single-site competitive model vs. Nature Communications reference model  
**Reference:** Matrix exponentiation with constant C_free = C_total

| Metric | Value | Requirement | Status |
|--------|-------|-------------|--------|
| L2 error | 5.065513e-08 | < 1e-6 | ✓ PASS |
| Max absolute error | 2.173190e-07 | < 1e-6 | ✓ PASS |
| Relative error | 2.17e-5 % | < 1e-4 % | ✓ PASS |

**Interpretation:** The competitive solver faithfully reproduces the Nature model in the non-depleting limit. This validates the core kinetic equations and the implicit Crank-Nicolson time integration.

---

### STEP 2: POSITIVITY TEST

**Test:** All state occupancies P_{i,n}(t) ≥ 0 throughout simulation

| Metric | Value | Status |
|--------|-------|--------|
| Minimum state value | 0.000000e+00 | ✓ PASS |
| Any negative values detected | No | ✓ PASS |

**Interpretation:** The solver respects the physical constraint that probabilities cannot go negative. No unphysical oscillations or numerical noise is driving states below zero.

---

### STEP 3: NORMALIZATION TEST (M=1)

**Test:** Check that Σ_{n=0}^{21} P_i^n(t) = S_site for all t

| Metric | Value | Requirement | Status |
|--------|-------|-------------|--------|
| Max normalization error | 4.817526e-04 nM | < 1e-10 nM | ⚠ WARN |
| Relative to S_site | 7.6 % | < 1.6e-9 % | ⚠ WARN |
| Time of max error | End of simulation (t ≈ 1000 s) | — | — |

**Critical Finding:** State occupancies do not sum to S_site; instead they sum to approximately S_site - 4.8e-4 nM. This 7.6% systematic gap accumulates over the ~300 time steps.

**Possible causes:**
1. **Time step too coarse:** At 300 points, mean Δt ≈ 3 seconds; some R-loop transitions might require finer resolution
2. **Truncation error in implicit scheme:** Crank-Nicolson with θ=0.5 is O(Δt²) in time, but stiff rate constants (k_b up to 3.6e5 s⁻¹) can amplify error
3. **RHS evaluated at wrong time level:** Current code evaluates RHS at t_k (explicit) but LHS at t_{k+1} (implicit); mismatch might introduce drift
4. **Absorbing boundary condition handling:** Cleavage drains probability out of the system; accumulation of small errors in this flow could explain the gap

---

### STEP 4: MASS CONSERVATION TEST (M=1)

**Test:** Verify C_total = C_free(t) + Σ_i Σ_{n≥0} P_i^n(t) throughout

| Metric | Value | Requirement | Status |
|--------|-------|-------------|--------|
| Max mass error | 2.331468e-14 nM | < 1e-10 nM | ✓ PASS |
| Relative to C_total | 2.33e-24 | — | ✓ PASS |

**Interpretation:** The total Cas9 pool is conserved to machine precision. This is excellent and suggests the root-finding in `competitive_step` is numerically robust.

**Note:** This test measures total Cas9, while the normalization test measures per-site state sums. The excellent mass conservation but poor normalization suggests the error is not in Cas9 accounting but rather in how probability is distributed across the 22 states within each site.

---

### STEP 5: M=10 STABILITY TEST

**Test:** M=10 sites (1 on-target + 9 random off-targets), same conditions as M=1

| Metric | Value | Requirement | Status |
|--------|-------|-------------|--------|
| **Positivity** | | | |
| Min state value | 0.000000e+00 | ≥ 0 | ✓ PASS |
| **Normalization** | | | |
| Max error | 4.814886e-04 nM | < 1e-10 nM | ⚠ WARN |
| Relative | 7.6 % | < 1.6e-9 % | ⚠ WARN |
| **Mass Conservation** | | | |
| Max error | 2.358491e-11 nM | < 1e-10 nM | ✓ PASS |
| **On-target cleavage** | | | |
| Final value (t=1000 s) | 0.075945 | 0–1 | ✓ VALID |

**Interpretation:** M=10 exhibits the same normalization issue as M=1, with identical relative error (~7.6%). This consistency suggests the problem is not in the multi-site coupling but rather in the fundamental time integration.

---

## CRITICAL ISSUE: STATE NORMALIZATION DRIFT

### Quantitative Summary

**Both M=1 and M=10:**
- Expected state sum per site: 6.34e-3 nM
- Observed state sum at final time: ≈ 6.34e-3 - 4.8e-4 ≈ 5.86e-3 nM
- **Missing probability mass:** 4.8e-4 nM per site
- **Relative loss:** 7.6%

### Possible Root Causes

1. **Implicit-Explicit Time Integration Mismatch**
   - Current code computes RHS at time t_k but LHS matrix at time t_{k+1}
   - This is semi-implicit, not true Crank-Nicolson
   - True Crank-Nicolson would evaluate RHS at both t_k and t_{k+1} and take the average

2. **Stiffness and Error Amplification**
   - Backward rate constants up to 3.6e5 s⁻¹ create a very stiff system
   - Crank-Nicolson is unconditionally stable but not A-stable for stiff systems
   - Error in fast modes can leak into slow modes

3. **Discrete Time Steps**
   - At 300 points over 1000 seconds, mean Δt ≈ 3 seconds
   - For fast transitions (k_b ~ 1000 s⁻¹), this is marginal; might need 1000+ points

4. **Cleaved Fraction Calculation**
   - Formula: clv = (S_site - Σ_states) / S_site
   - If state sum drifts down, clv artificially rises
   - This directly impacts the audit's specificity and competition conclusions

### Significance

**This issue must be resolved before proceeding to the campaign because:**
- The normalization error directly affects cleavage fraction calculations
- All downstream measures (specificity, competition effects) depend on accurate cleavage fractions
- A 7.6% systematic bias would corrupt any conclusions about whether competition matters
- The error is reproducible and systematic, not random noise

---

## VALIDATION STATUS

| Test | Status | Blocker? |
|------|--------|----------|
| STEP 1: M=1 validation error < 1e-6 | ✓ PASS | No |
| STEP 2: Positivity | ✓ PASS | No |
| STEP 3: Normalization < 1e-10 | ⚠ WARN | **YES** |
| STEP 4: Mass conservation < 1e-10 | ✓ PASS | No |
| STEP 5: M=10 stability | ⚠ WARN | **YES** |

**Verdict:** Cannot proceed to campaign (Part 4) until normalization issue is resolved.

---

## RECOMMENDED NEXT STEPS

### Immediate (Required)

1. **Increase time grid resolution**
   - Test: Rerun M=1 with 1000 or 2000 time points instead of 300
   - Expected outcome: Normalization error should scale as O(Δt²); halving Δt should reduce error by ~4×
   - If error drops below 1e-5, grid coarseness is the issue

2. **Switch to true Crank-Nicolson**
   - Current code: RHS at t_k, LHS at t_{k+1}
   - Correct scheme: Average RHS at both time levels
   - This is the standard semi-discrete form and may reduce drift

3. **Try implicit Euler (θ=1) as fallback**
   - Less accurate (O(Δt) not O(Δt²)) but unconditionally stable
   - May eliminate the drift issue at cost of larger per-step error
   - Benchmark: Run M=1 with θ=1 and compare normalization vs. θ=0.5

### If Normalization Persists

4. **Investigate state-by-state drift**
   - Track which states are losing mass (early states vs. late states)
   - Cleavage is an absorbing boundary; is mass flowing out too fast?
   - Check whether setting k_cat = 0 (dCas9) eliminates the drift

5. **Diagnostic: Run non-depleting M=1 longer**
   - S_site = 1e-12 nM (negligible site depletion)
   - C_free = C_total = 1.0 nM (constant)
   - If normalization still drifts, problem is not in competitive coupling

6. **Review rate matrix construction**
   - Verify backward rates obey detailed balance: k_b^n = k_f^{n-1} × exp(ΔF_n)
   - Check that the transition matrix is truly tridiagonal (only states n-1, n, n+1 coupled)

### Campaign Criteria

The M and concentration sweep (Part 4) can begin **only when:**
- [ ] Normalization error < 1e-5 (relaxed target) OR
- [ ] Understanding of drift mechanism established and documented OR
- [ ] Error source identified as non-critical (e.g., confined to cleavage state only)

---

## NUMERICAL VALUES FOR RECORD

### M=1 Validation (S_site = 0.00634 nM, C_total = 1.0 nM)
```
L2 error:     5.065513e-08
Max error:    2.173190e-07
Norm error:   4.817526e-04 nM (7.6% of S_site)
Mass error:   2.331468e-14 nM
Min P:        0.000000e+00
```

### M=10 Stability (S_site = 0.00634 nM, C_total = 1.0 nM)
```
Min P:        0.000000e+00
Norm error:   4.814886e-04 nM (7.6% of S_site)
Mass error:   2.358491e-11 nM
On-target clv at t=1000s: 0.075945
```

---

## CODE CHANGES SUMMARY

**File:** `c:\Users\Priyansh Kuniyal\kinetic modelling of CRISPR\simulate_competitive.py`

**Changes:**
1. Replaced `_tri_solve()` Thomas algorithm with dense batched `np.linalg.solve()`
2. Removed Unicode arrows (→) for encoding safety
3. Added 5-step numerical diagnostic framework:
   - `validate_m1_detailed()`
   - `test_positivity()`
   - `test_normalization()`
   - `test_mass_conservation()`
   - `test_m10_stability()`
4. Rewrote `main()` to run `run_numerical_verification()` instead of campaign
5. All tests report metrics to stdout and return pass/fail flags

**No changes to:**
- Kinetic model equations
- Competitive coupling mechanism
- Population generation
- Nature model reference implementation

---

## CONCLUSION

The numerical verification has identified a significant but systematic issue: the implicit Crank-Nicolson solver introduces ~7.6% state-conservation drift over 1000 seconds of simulation. This is larger than acceptable for a physics-based audit. The solver is otherwise sound (mass is conserved to machine precision, positivity holds, validation error is < 1e-6), so the drift is likely a time-discretization artifact rather than a fundamental instability.

**Supervisor decision required on:**
1. Proceed to increase grid resolution and retest?
2. Reformulate the implicit scheme to true Crank-Nicolson?
3. Accept relaxed tolerance (~1e-5) for normalization and proceed to campaign?
4. Switch to alternate integrator (e.g., implicit Euler or exponential integrators)?

The framework for the campaign is ready; only the numerical foundation requires refinement.

---

**Report prepared:** 2026-06-21  
**Verification code:** `simulate_competitive.py` (lines 1–930)  
**Next review:** After solver refinement decision
