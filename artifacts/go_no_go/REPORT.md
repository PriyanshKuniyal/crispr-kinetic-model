# CRISPR Competitive-Pool Go/No-Go Report

**Final verdict: `PIVOT_EFFECT_ONLY_IN_EXTREME_REGIMES`**

## GPU execution

- Backend: CuPy 14.1.1, CUDA runtime 12090
- Device: NVIDIA GeForce RTX 5050 Laptop GPU (compute capability 12.0)
- Campaign ran on GPU (cupy arrays confirmed): True
- Peak GPU utilization during campaign: 49%  (peak mem 1021 MB)

## Controls

1. **M=1 reproduces Nature**: max |Δon| = 5.018e-02 pp — PASS
2. **S_site ×1e-4 → independent-site limit**: dev from Nature 3.202e-01 pp (normal) → 5.007e-02 pp (tiny); depletion 9.691e-01% → 9.785e-05% — PASS
3. **GPU vs CPU agreement**: M=1 max|Δ|=4.11e-16, M=10 max|Δ|=4.11e-16 — PASS
4. **Conservation/positivity/C_free≥0/finite**: cleaved∈[-3.65e-10,1.0000], min C_free=9.758e-01, NaN/Inf=False — PASS

## Headline differences (Competitive − Nature)

- **Largest on-target difference (any panel):** 17.289 pp at panel=DEFENSIBLE, M=10000, C=0.1 nM
- **Largest specificity change (any panel):** 10.67% at panel=DEFENSIBLE, M=10000, C=3.0 nM
- **Largest free-Cas9 depletion:** 37.95% at panel=DEFENSIBLE, M=10000, C=0.1 nM
- **Largest on-target diff under DEFENSIBLE:** 17.289 pp at M=10000, C=0.1 nM
- **Largest specificity change under DEFENSIBLE:** 10.67% at M=10000, C=3.0 nM

## Flagged conditions

| panel | M | C | dOnMax(pp) | spec%Δ | depl% | timing%Δ | rankInv | flags |
|---|---|---|---|---|---|---|---|---|
| DEFENSIBLE | 10000 | 0.1 | 17.289 | 8.44 | 37.95 | 59.09 | 0 | >5pp,depl,timing |
| DEFENSIBLE | 10000 | 0.3 | 17.205 | 8.20 | 37.86 | 60.21 | 0 | >5pp,depl,timing |
| DEFENSIBLE | 10000 | 1.0 | 16.910 | 10.17 | 37.58 | 58.96 | 0 | >5pp,spec,depl,timing |
| DEFENSIBLE | 10000 | 3.0 | 16.227 | 10.67 | 36.86 | 55.55 | 0 | >5pp,spec,depl,timing |
| DEFENSIBLE | 10000 | 10.0 | 14.390 | 8.61 | 34.85 | 46.77 | 0 | >5pp,depl,timing |
| CONTROLLED_STRESS | 10000 | 0.1 | 13.883 | 2.82 | 33.22 | 46.08 | 0 | >5pp,depl,timing |
| CONTROLLED_STRESS | 10000 | 0.3 | 13.795 | 1.75 | 32.99 | 46.23 | 0 | >5pp,depl,timing |
| CONTROLLED_STRESS | 10000 | 1.0 | 13.540 | 1.01 | 32.32 | 45.05 | 0 | >5pp,depl,timing |
| CONTROLLED_STRESS | 10000 | 3.0 | 12.867 | 1.11 | 30.92 | 42.10 | 0 | >5pp,depl,timing |
| DEFENSIBLE | 10000 | 30.0 | 11.420 | 5.69 | 30.97 | 33.49 | 0 | >5pp,depl,timing |
| CONTROLLED_STRESS | 10000 | 10.0 | 11.118 | 1.07 | 27.87 | 34.68 | 0 | >5pp,depl,timing |
| CONTROLLED_STRESS | 10000 | 30.0 | 8.562 | 0.62 | 23.46 | 24.02 | 0 | >5pp,depl,timing |
| DEFENSIBLE | 10000 | 100.0 | 6.856 | 2.62 | 23.33 | 15.58 | 0 | >5pp,depl,timing |
| CONTROLLED_STRESS | 10000 | 100.0 | 5.061 | 0.21 | 16.59 | 11.55 | 0 | >5pp,depl,timing |
| DEFENSIBLE | 1000 | 0.1 | 2.153 | 1.21 | 5.73 | 6.01 | 0 | >1pp,depl |
| DEFENSIBLE | 1000 | 0.3 | 2.137 | 1.39 | 5.70 | 5.91 | 0 | >1pp,depl |
| DEFENSIBLE | 1000 | 1.0 | 2.084 | 1.45 | 5.61 | 5.71 | 1 | >1pp,depl,rank |
| DEFENSIBLE | 1000 | 3.0 | 1.947 | 1.27 | 5.39 | 5.30 | 0 | >1pp,depl |
| DEFENSIBLE | 1000 | 10.0 | 1.641 | 0.99 | 4.85 | 4.39 | 0 | >1pp |
| CONTROLLED_STRESS | 1000 | 0.1 | 1.626 | 0.25 | 4.72 | 4.62 | 0 | >1pp |
| CONTROLLED_STRESS | 1000 | 0.3 | 1.613 | 0.16 | 4.65 | 4.54 | 0 | >1pp |
| CONTROLLED_STRESS | 1000 | 1.0 | 1.568 | 0.10 | 4.48 | 4.30 | 0 | >1pp |
| CONTROLLED_STRESS | 1000 | 3.0 | 1.461 | 0.12 | 4.14 | 4.00 | 0 | >1pp |
| DEFENSIBLE | 1000 | 30.0 | 1.222 | 0.63 | 3.98 | 2.97 | 0 | >1pp |
| CONTROLLED_STRESS | 1000 | 10.0 | 1.214 | 0.11 | 3.53 | 3.31 | 0 | >1pp |

## Full results

| panel | M | C | onA | onB | dOnMax(pp) | offΔmax(pp) | specΔ% | depl% | seq% | timing%Δ | rankInv |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CONTROLLED_STRESS | 1 | 0.1 | 0.9997 | 0.9997 | 0.000 | 0.000 | -0.00 | 0.00 | 0.00 | 0.00 | 0 |
| CONTROLLED_STRESS | 1 | 0.3 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| CONTROLLED_STRESS | 1 | 1.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| CONTROLLED_STRESS | 1 | 3.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| CONTROLLED_STRESS | 1 | 10.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| CONTROLLED_STRESS | 1 | 30.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| CONTROLLED_STRESS | 1 | 100.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| CONTROLLED_STRESS | 10 | 0.1 | 0.9997 | 0.9997 | 0.015 | 0.005 | 0.00 | 0.05 | 0.05 | 0.04 | 0 |
| CONTROLLED_STRESS | 10 | 0.3 | 1.0000 | 1.0000 | 0.015 | 0.005 | 0.00 | 0.05 | 0.05 | 0.04 | 0 |
| CONTROLLED_STRESS | 10 | 1.0 | 1.0000 | 1.0000 | 0.015 | 0.005 | 0.00 | 0.05 | 0.05 | 0.04 | 0 |
| CONTROLLED_STRESS | 10 | 3.0 | 1.0000 | 1.0000 | 0.015 | 0.005 | 0.00 | 0.04 | 0.04 | 0.04 | 0 |
| CONTROLLED_STRESS | 10 | 10.0 | 1.0000 | 1.0000 | 0.013 | 0.004 | 0.00 | 0.04 | 0.04 | 0.04 | 0 |
| CONTROLLED_STRESS | 10 | 30.0 | 1.0000 | 1.0000 | 0.010 | 0.003 | 0.00 | 0.03 | 0.03 | 0.03 | 0 |
| CONTROLLED_STRESS | 10 | 100.0 | 1.0000 | 1.0000 | 0.006 | 0.001 | 0.00 | 0.02 | 0.02 | 0.01 | 0 |
| CONTROLLED_STRESS | 100 | 0.1 | 0.9997 | 0.9997 | 0.124 | 0.062 | 0.02 | 0.37 | 0.37 | 0.35 | 0 |
| CONTROLLED_STRESS | 100 | 0.3 | 1.0000 | 1.0000 | 0.124 | 0.061 | 0.01 | 0.37 | 0.37 | 0.35 | 0 |
| CONTROLLED_STRESS | 100 | 1.0 | 1.0000 | 1.0000 | 0.122 | 0.060 | 0.01 | 0.36 | 0.36 | 0.33 | 0 |
| CONTROLLED_STRESS | 100 | 3.0 | 1.0000 | 1.0000 | 0.117 | 0.058 | 0.01 | 0.34 | 0.34 | 0.31 | 0 |
| CONTROLLED_STRESS | 100 | 10.0 | 1.0000 | 1.0000 | 0.103 | 0.051 | 0.01 | 0.30 | 0.30 | 0.29 | 0 |
| CONTROLLED_STRESS | 100 | 30.0 | 1.0000 | 1.0000 | 0.081 | 0.038 | 0.00 | 0.25 | 0.25 | 0.20 | 0 |
| CONTROLLED_STRESS | 100 | 100.0 | 1.0000 | 1.0000 | 0.047 | 0.020 | 0.00 | 0.17 | 0.17 | 0.10 | 0 |
| CONTROLLED_STRESS | 1000 | 0.1 | 0.9997 | 0.9996 | 1.626 | 0.797 | 0.25 | 4.72 | 4.72 | 4.62 | 0 |
| CONTROLLED_STRESS | 1000 | 0.3 | 1.0000 | 1.0000 | 1.613 | 0.790 | 0.16 | 4.65 | 4.65 | 4.54 | 0 |
| CONTROLLED_STRESS | 1000 | 1.0 | 1.0000 | 1.0000 | 1.568 | 0.768 | 0.10 | 4.48 | 4.48 | 4.30 | 0 |
| CONTROLLED_STRESS | 1000 | 3.0 | 1.0000 | 1.0000 | 1.461 | 0.714 | 0.12 | 4.14 | 4.14 | 4.00 | 0 |
| CONTROLLED_STRESS | 1000 | 10.0 | 1.0000 | 1.0000 | 1.214 | 0.585 | 0.11 | 3.53 | 3.53 | 3.31 | 0 |
| CONTROLLED_STRESS | 1000 | 30.0 | 1.0000 | 1.0000 | 0.890 | 0.410 | 0.06 | 2.77 | 2.77 | 2.21 | 0 |
| CONTROLLED_STRESS | 1000 | 100.0 | 1.0000 | 1.0000 | 0.503 | 0.206 | 0.02 | 1.79 | 1.79 | 1.09 | 0 |
| CONTROLLED_STRESS | 10000 | 0.1 | 0.9997 | 0.9973 | 13.883 | 6.694 | 2.82 | 33.22 | 33.22 | 46.08 | 0 |
| CONTROLLED_STRESS | 10000 | 0.3 | 1.0000 | 1.0000 | 13.795 | 6.658 | 1.75 | 32.99 | 32.99 | 46.23 | 0 |
| CONTROLLED_STRESS | 10000 | 1.0 | 1.0000 | 1.0000 | 13.540 | 6.525 | 1.01 | 32.32 | 32.32 | 45.05 | 0 |
| CONTROLLED_STRESS | 10000 | 3.0 | 1.0000 | 1.0000 | 12.867 | 6.172 | 1.11 | 30.92 | 30.92 | 42.10 | 0 |
| CONTROLLED_STRESS | 10000 | 10.0 | 1.0000 | 1.0000 | 11.118 | 5.279 | 1.07 | 27.87 | 27.87 | 34.68 | 0 |
| CONTROLLED_STRESS | 10000 | 30.0 | 1.0000 | 1.0000 | 8.562 | 3.889 | 0.62 | 23.46 | 23.46 | 24.02 | 0 |
| CONTROLLED_STRESS | 10000 | 100.0 | 1.0000 | 1.0000 | 5.061 | 2.051 | 0.21 | 16.59 | 16.59 | 11.55 | 0 |
| DEFENSIBLE | 1 | 0.1 | 0.9997 | 0.9997 | 0.000 | 0.000 | -0.00 | 0.00 | 0.00 | 0.00 | 0 |
| DEFENSIBLE | 1 | 0.3 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| DEFENSIBLE | 1 | 1.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| DEFENSIBLE | 1 | 3.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| DEFENSIBLE | 1 | 10.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| DEFENSIBLE | 1 | 30.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| DEFENSIBLE | 1 | 100.0 | 1.0000 | 1.0000 | 0.000 | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| DEFENSIBLE | 10 | 0.1 | 0.9997 | 0.9997 | 0.015 | 0.000 | 0.04 | 0.04 | 0.04 | 0.04 | 0 |
| DEFENSIBLE | 10 | 0.3 | 1.0000 | 1.0000 | 0.015 | 0.000 | 0.04 | 0.04 | 0.04 | 0.04 | 0 |
| DEFENSIBLE | 10 | 1.0 | 1.0000 | 1.0000 | 0.015 | 0.000 | 0.04 | 0.04 | 0.04 | 0.04 | 0 |
| DEFENSIBLE | 10 | 3.0 | 1.0000 | 1.0000 | 0.014 | 0.000 | 0.04 | 0.04 | 0.04 | 0.04 | 0 |
| DEFENSIBLE | 10 | 10.0 | 1.0000 | 1.0000 | 0.013 | 0.000 | 0.03 | 0.04 | 0.04 | 0.04 | 0 |
| DEFENSIBLE | 10 | 30.0 | 1.0000 | 1.0000 | 0.011 | 0.000 | 0.02 | 0.03 | 0.03 | 0.03 | 0 |
| DEFENSIBLE | 10 | 100.0 | 1.0000 | 1.0000 | 0.006 | 0.000 | 0.01 | 0.02 | 0.02 | 0.01 | 0 |
| DEFENSIBLE | 100 | 0.1 | 0.9997 | 0.9997 | 0.208 | 0.005 | 0.11 | 0.57 | 0.57 | 0.58 | 0 |
| DEFENSIBLE | 100 | 0.3 | 1.0000 | 1.0000 | 0.206 | 0.005 | 0.11 | 0.57 | 0.57 | 0.57 | 0 |
| DEFENSIBLE | 100 | 1.0 | 1.0000 | 1.0000 | 0.200 | 0.005 | 0.08 | 0.55 | 0.55 | 0.53 | 0 |
| DEFENSIBLE | 100 | 3.0 | 1.0000 | 1.0000 | 0.187 | 0.004 | 0.06 | 0.53 | 0.53 | 0.49 | 0 |
| DEFENSIBLE | 100 | 10.0 | 1.0000 | 1.0000 | 0.159 | 0.003 | 0.05 | 0.47 | 0.47 | 0.42 | 0 |
| DEFENSIBLE | 100 | 30.0 | 1.0000 | 1.0000 | 0.119 | 0.004 | 0.05 | 0.39 | 0.39 | 0.29 | 0 |
| DEFENSIBLE | 100 | 100.0 | 1.0000 | 1.0000 | 0.067 | 0.002 | 0.03 | 0.26 | 0.26 | 0.13 | 0 |
| DEFENSIBLE | 1000 | 0.1 | 0.9997 | 0.9995 | 2.153 | 0.049 | 1.21 | 5.73 | 5.73 | 6.01 | 0 |
| DEFENSIBLE | 1000 | 0.3 | 1.0000 | 1.0000 | 2.137 | 0.048 | 1.39 | 5.70 | 5.70 | 5.91 | 0 |
| DEFENSIBLE | 1000 | 1.0 | 1.0000 | 1.0000 | 2.084 | 0.067 | 1.45 | 5.61 | 5.61 | 5.71 | 1 |
| DEFENSIBLE | 1000 | 3.0 | 1.0000 | 1.0000 | 1.947 | 0.079 | 1.27 | 5.39 | 5.39 | 5.30 | 0 |
| DEFENSIBLE | 1000 | 10.0 | 1.0000 | 1.0000 | 1.641 | 0.082 | 0.99 | 4.85 | 4.85 | 4.39 | 0 |
| DEFENSIBLE | 1000 | 30.0 | 1.0000 | 1.0000 | 1.222 | 0.065 | 0.63 | 3.98 | 3.98 | 2.97 | 0 |
| DEFENSIBLE | 1000 | 100.0 | 1.0000 | 1.0000 | 0.682 | 0.033 | 0.27 | 2.67 | 2.67 | 1.42 | 0 |
| DEFENSIBLE | 10000 | 0.1 | 0.9997 | 0.9936 | 17.289 | 0.476 | 8.44 | 37.95 | 37.95 | 59.09 | 0 |
| DEFENSIBLE | 10000 | 0.3 | 1.0000 | 1.0000 | 17.205 | 0.473 | 8.20 | 37.86 | 37.86 | 60.21 | 0 |
| DEFENSIBLE | 10000 | 1.0 | 1.0000 | 1.0000 | 16.910 | 0.462 | 10.17 | 37.58 | 37.58 | 58.96 | 0 |
| DEFENSIBLE | 10000 | 3.0 | 1.0000 | 1.0000 | 16.227 | 0.605 | 10.67 | 36.86 | 36.86 | 55.55 | 0 |
| DEFENSIBLE | 10000 | 10.0 | 1.0000 | 1.0000 | 14.390 | 0.651 | 8.61 | 34.85 | 34.85 | 46.77 | 0 |
| DEFENSIBLE | 10000 | 30.0 | 1.0000 | 1.0000 | 11.420 | 0.540 | 5.69 | 30.97 | 30.97 | 33.49 | 0 |
| DEFENSIBLE | 10000 | 100.0 | 1.0000 | 1.0000 | 6.856 | 0.300 | 2.62 | 23.33 | 23.33 | 15.58 | 0 |

## Interpretation (adversarial)

**Mechanism.** The competitive model differs from the Nature independent-site
model only through the shared factor `C_free(t)`, which scales the PAM on-rate
of *every* site identically. Finite-pool competition therefore slows on- and
off-target binding *together*; it does not preferentially protect the
on-target. The result is a near-uniform kinetic slowdown, not a specificity
mechanism. This is visible directly in the data:

- **On-target differences are transient, not outcome changes.** The headline
  17.3 pp is the *maximum* gap during the approach to completion. Final
  on-target cleavage (`onB` vs `onA` columns) is ~1.0000 vs ~1.0000 for every
  condition except the single most extreme corner (M=10000, C=0.1: 0.9936 vs
  0.9997, a 0.6 pp final gap). Competition delays cleavage; it does not prevent
  it at realistic burdens.
- **Specificity is nearly invariant.** Max specificity change is 10.67% and it
  occurs only at M=10000. At every realistic burden (M<=1000) specificity moves
  <=1.5%.
- **No robust ranking inversions.** The one flagged inversion (DEFENSIBLE,
  M=1000, C=1.0) is a single near-degenerate tie flip with Spearman = 0.9999 —
  numerical noise, not a biological reordering. Off-target site rank order is
  preserved everywhere.
- **The real signature is enzyme starvation.** Free-Cas9 depletion (= sink
  sequestration) reaches 20-38% and t_half is delayed 30-60% — but only at
  M=10000 competing occupancy-holding loci with sub-nM to low-nM Cas9. That is
  textbook mass-action ("not enough enzyme to go round"), predictable without a
  new model.

**Regime dependence.** Every crossing of a *biological* threshold (>5 pp
on-target OR >10% specificity) occurs at M=10000. At M<=1000 the maximum effects
are ~2 pp transient on-target, ~1.5% specificity, ~5.7% depletion, no ranking
change, and identical final outcome. M=10000 distinct occupancy-holding
competitor loci for a single RNP species, combined with a starved (<=few nM)
Cas9 pool, is the extreme upper end of a genome-wide off-target burden — not a
typical genome-engineering regime.

## Verdict rationale

`PIVOT_EFFECT_ONLY_IN_EXTREME_REGIMES`. A genuine, reproducible, monotone,
numerically-certified difference exists (so it is not `ABORT`), but it crosses
biologically meaningful magnitudes only in the extreme-M / low-Cas9 corner, and
even there it is dominated by trivial enzyme depletion rather than any change in
discrimination or site ranking. The specificity hypothesis — that a shared
finite pool meaningfully reshapes on-vs-off-target selectivity — is **not**
supported at realistic competitor counts and Cas9 concentrations.

## Answers

1. **Largest difference from the Nature model:** 17.29 pp *transient* on-target
   cleavage (final gap only 0.6 pp); companion effects at that corner: 38%
   free-Cas9 depletion, ~59% delay in on-target t_half. Largest specificity
   change anywhere: 10.67%.
2. **Largest difference under the DEFENSIBLE panel:** identical to the global
   max — 17.29 pp transient on-target; 10.67% specificity (at C=3.0);
   37.95% depletion.
3. **Exact M and Cas9 concentration:** M = 10000 sites, C_total = 0.1 nM
   (on-target/depletion max); the specificity max is M = 10000, C_total = 3.0 nM.
4. **Biologically meaningful magnitude?** No, not under realistic conditions.
   At M<=1000 with Cas9 >= 1 nM the differences are transient (<=2 pp), leave
   final cleavage and site ranking unchanged, and move specificity <1.5%.
   Meaningful magnitudes appear only at M=10000 competing loci with sub-/low-nM
   Cas9, and are driven by enzyme starvation, not selectivity.
5. **Does the hypothesis survive?** As a general claim (competition meaningfully
   changes on/off specificity), **no** — kill it. As a narrow claim (finite-pool
   depletion measurably delays/attenuates cleavage under extreme off-target
   burden with scarce Cas9), it survives only in that corner and should be
   reframed as an enzyme-titration/timing effect, not a specificity mechanism.
