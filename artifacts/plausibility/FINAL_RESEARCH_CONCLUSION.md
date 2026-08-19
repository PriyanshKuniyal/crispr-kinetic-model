# Final Research Conclusion: Falsification of Finite-Pool Competition Models in CRISPR-Cas9
**Date:** July 11, 2026
**Investigator:** Priyansh Kuniyal / Claude Code

## Background
The research sought to extend the single-molecule kinetic model by Eslami-Mossallam et al. (Nature Communications, 2022) to determine if a shared, finite pool of SpCas9 creates a biologically meaningful competitive effect across thousands of genomic targets. Specifically, we tested if this competition mechanism provides an improvement in **on-vs-off target specificity** or reduces **editing yield** in plausible cellular regimes.

## Step 1: Solving the Numerical Drift (Solver Foundation)
Before running accurate diagnostic limits, the solver faced a systemic blocker: state normalization drifted by ~7.6% over a 1000s simulation due to artifacts in the `simulate_competitive.py` integration loop. 
- **Fix:** We replaced the unstable batched Thomas algorithm with CuPy's robust LAPACK batched dense solver (`cp.linalg.solve`). Furthermore, we corrected the semi-implicit Crank-Nicolson formulation to precisely propagate `c_free_prev` and `c_free_next`. 
- **Validation:** Mass conservation was successfully restored down to $10^{-14}$ precision limits, completely clearing the 7.6% error.

## Step 2: Testing Pivot 1 (Finite Cas9 Lifetime)
With the numerical foundation strictly verified, we ran diagnostics on our strongest alternative hypothesis (Pivot 1): If Cas9 is active only for a finite exposure window (e.g. transient delivery or RNP degradation), does the competitive burden induce a permanent deficit in the on-target editing yield?

We subjected this to the pre-registered mandatory kill-conditions (`pivot_audit.md`).

### Results:
1. **Defensible Burden Failure:** At a defensible burden ($M=300$, Cas9$=300$, $V=500$ fL) with physiological RNP exposure ($\tau = 24$ h), the yield deficit is only **0.008 percentage points (pp)**. This egregiously fails the required $\ge 2$ pp threshold.
2. **Implausible Conjunction Required:** To achieve an editing deficit greater than 2 pp, the model explicitly requires a triple-implausible conjunction: $M \ge 1000$ occupancy-competent loci **and** severe Cas9 scarcity ($N_c = 100$) **and** a short exposure window ($\tau \le 12$ h). For example, the `M1000/Cas100/250fL` corner showed a 4.4 pp deficit at 6 hours, but relaxing any single one of these variables collapses the deficit to $< 1$ pp.

## Conclusion: OFFICIALLY REJECTED
Both the original steady-state specificity hypothesis and the fallback finite-lifetime yield hypothesis are **rejected for practical biological relevance.**

While the math robustly confirms that finite exposure *does* convert transient competition delays into permanent yield deficits, the parameters required to produce an experimentally observable effect ($>2\%$) do not align with plausible biology. In virtually all realistic CRISPR editing scenarios, Cas9 competition does not tangibly influence specificity nor yield.
