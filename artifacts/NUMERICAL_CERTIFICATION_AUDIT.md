# Numerical Certification Audit (No New Features)

## Scope

This audit certifies the existing simulator mathematically. No solver, timestep, Crank-Nicolson scheme, competition model, or biology was changed.

## Probability Flow Traced from Implementation

From `simulate_competitive.py`:

1. **Free Cas9 -> PAM binding**  
   `cn_step()` scales only `kf_eff[:, 0]` by `c_free` (`kf_eff[:, 0] *= c_free`).
2. **PAM -> R-loop states**  
   Forward rates (`fwd[1:20]`) and backward rates (`bck[1:21]`) are applied in `_rhs_times_P()` and `_lhs_matrices()` as nearest-neighbor fluxes.
3. **R-loop state 20 -> cleavage sink**  
   `fwd[-1] = k_cat` is included as loss on the state-21 diagonal (`diag_K = -(kf_eff + kb)`), but there is no propagated product state receiving this flux.

So cleavage is an **intentional absorbing outflow** from the propagated 22-state vector.

## Correct Conservation Laws

For active Cas9 (`k_cat > 0`):

```text
d/dt sum(states) = -k_cat * P_20
d/dt cleaved     =  k_cat * P_20
=> sum(states) + cleaved_mass = S_site
```

For dCas9 (`forward_rates[-1] = 0`, so `k_cat = 0`):

```text
sum(states) = S_site
```

Therefore a decrease in `sum(states)` alone is expected for active Cas9 and is not, by itself, a numerical error.

## Numerical Verification

Conditions:

- `S_site = 6.34000000e-03 nM`
- `C_total = 1.00000000e+00 nM`
- final time = `1000 s`

Active Cas9:

- `sum(states_final) = 5.85824742e-03 nM`
- `cleaved_fraction_final = 7.59862110e-02`
- `cleaved_mass_final = 4.81752577e-04 nM`
- `max |S_site - sum(states) - cleaved_mass| = 2.71050543e-20 nM`

dCas9:

- `max |S_site - sum(states)| = 1.46584134e-16 nM`

Consistency of prior "7.6% drift":

- `1 - sum(states_final)/S_site = 7.59862110e-02`
- `| (1 - sum(states_final)/S_site) - cleaved_fraction_final | = 4.16333634e-17`

So the historical 7.6% "loss" equals cleavage flux, not numerical drift.

## Audit Finding

The previous 7.6% "normalization drift" is **expected physics** from irreversible cleavage leaving the transient propagated state vector.

## Numerical Certification Decision

**PASS**

The solver is mathematically validated against the correct conservation law.

Validation figure: `artifacts/conservation_audit.png`  
Structured metrics: `artifacts/numerical_certification_audit.json`
