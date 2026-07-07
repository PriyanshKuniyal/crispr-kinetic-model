# Strict Pivot Audit

**Verdict: `PIVOT_TO_FINITE_CAS9_LIFETIME`** (narrow, falsifiable; subsumes
transient/pulsed delivery, which is its experimentally-controllable twin).

Rule applied: a pivot is admitted only if it **changes the biological question**
(different outcome variable and/or a new physical parameter the killed hypothesis
ignored) and shows a **measurable effect in existing/diagnostic outputs** — not
if it merely re-parameterises the same steady-state competition.

The killed hypothesis compared models at t→∞ with **constant** Cas9, where the
on-target always completes (yield→1.0). Its fatal flaw: *final editing unchanged*.
Any credible pivot must break that specific flaw.

## Diagnostic evidence (step Cas9-exposure window, `finite_lifetime_diagnostic.json`)

Editing **yield deficit** (no-competition − with-burden) if active Cas9 exists
only for a window `T_cas`, read directly off the trusted competition dynamics:

| condition (plausibility) | ratio 2M/Nc | max transient | deficit @6 h | @12 h | @24 h | @48 h |
|---|---|---|---|---|---|---|
| M1000 / Cas100 / 250 fL (plausible corner) | 20 | 4.42 pp @5.5 h | **4.40 pp** | 2.99 | 0.68 | 0.02 |
| M1000 / Cas1000 / 500 fL (plausible, moderate Cas9) | 2 | 2.07 pp | 0.11 | 0.00 | ~0 | 0 |
| M300 / Cas300 / 500 fL (defensible burden) | 2 | 0.69 pp | 0.58 | 0.20 | 0.01 | ~0 |
| M3000 / Cas100 / 250 fL (IMPLAUSIBLE burden) | 60 | 11.6 pp | 11.63 | 8.77 | 2.52 | 0.11 |
| M10000 / Cas100 / 500 fL (IMPLAUSIBLE burden) | 200 | 17.7 pp @13.5 h | 13.85 | 17.61 | 14.32 | 4.94 |

Reading: finite exposure **does** convert the transient into a *permanent* yield
deficit — the one thing the killed hypothesis could not do. But in plausible
regimes it is ≤ ~4 pp and requires a **triple conjunction**: high burden
(M~1000) **and** scarce Cas9 (~100 molecules) **and** a short (~6–12 h) exposure
window. Relax any one (more Cas9, fewer competitors, ≥24 h window) and it decays
to <1 pp. Rigorous bound: yield deficit ≤ max transient dOn, already shown ≤4.39 pp
in the plausible envelope.

## Candidate-by-candidate audit

### 1. Finite Cas9 lifetime / degradation  →  **ADMITTED (rank 1)**
- **New question:** does off-target burden reduce on-target editing *yield* while
  Cas9 is active for a finite time τ — rather than steady-state selectivity?
- **Mathematically distinct?** Yes. Adds a competing timescale τ; `C_total`
  becomes `C_total·e^{−t/τ}` (or a window); the outcome variable changes from a
  steady ratio to a *finite-time yield* that can be < 1.0. Directly repairs the
  "final editing unchanged" failure.
- **Minimum code change:** ~3 lines — pass `C0·exp(−t_k/τ)` as the time-varying
  total to the existing root-find each step (bound-shedding when `C_total(t) <
  bound` handled by the C_free≥0 clamp already present). A step window is even
  simpler (read existing curves at `T_cas`, done here with zero solver change).
- **Effect already visible?** Yes (table above): 4.4 pp plausible-corner,
  10–18 pp at implausible burden.
- **Strongest reason it may also fail:** plausible-regime effect is small
  (≤~4 pp) and triple-gated; it is the *same* competition mechanism read at
  finite time, so it inherits the burden requirement. Only credible for
  promiscuous guides under short/low-dose Cas9 exposure; specificity still
  unchanged (a yield effect only).

### 2. Transient / pulsed / time-dependent delivery  →  merged into #1
- **Mathematically distinct from #1?** No — identical math (`C_total(t)`
  time-varying). It is the **experimentally controllable** instantiation (RNP
  bolus vs sustained expression; dose; degron-tagged Cas9), which is why #1 is
  worth pursuing at all: the exposure window is a real knob, unlike "10⁴
  competitors". Recommended as the *experimental handle* for testing #1.
- **Why not the standalone verdict:** the delivery knob alone (normal dose, low
  burden) gives ~0 effect; it only matters combined with the finite window of #1.

### 3. Sequential target exposure  →  REJECTED (weak, derivative)
- **New question:** if loci become accessible at staggered times (replication,
  chromatin), does order-of-exposure change on-target yield?
- **Distinct?** Partially (staggered initial conditions), but its worst case
  (on-target opens into an already-depleted pool) is bounded by the steady
  finite-pool depletion — shown small at plausible burden. It produces nothing
  new *alone*; it only bites when combined with #1's finite window.
- **Min change:** moderate (per-site activation times). **Effect in existing
  outputs:** none isolated. **Why it fails:** reduces to time-shifted competition
  under the same 2M/Nc bound.

### 4. Competitor turnover / release  →  REJECTED (anti-mechanism)
- **New question:** if competitors rapidly release Cas9 (fast turnover), is
  net sequestration reduced?
- **Distinct?** It is already partly in the model (reversible backward rates;
  cleaved sites release Cas9). **Increasing** turnover *reduces* sequestration,
  which **weakens** the competition effect — it works *against* the hypothesis.
- **Why it fails definitively:** it can only shrink an already-small effect.

## Ranking

1. **PIVOT_TO_FINITE_CAS9_LIFETIME** — only candidate that (a) is mathematically
   distinct, (b) changes the outcome variable to a permanent editing yield, and
   (c) already shows a measurable effect. Pursue via transient-delivery experiments.
2. Transient delivery — same math; the experimental knob for #1 (not standalone).
3. Sequential exposure — derivative; only relevant combined with #1.
4. Competitor turnover — anti-mechanism; discard.

## Why this is a pivot, not a rescue

It does **not** revive the specificity claim (specificity stays unchanged). It
asks a *different* question — editing **yield under finite Cas9 exposure** — using
a *real parameter* (RNP lifetime / delivery window) the killed model ignored, and
it changes the outcome from a washed-out transient to a permanent, observable
deficit. That is a change of question, not a cosmetic re-fit.

## Mandatory kill-conditions for the next phase (pre-registered)

Pursue #1 **only** with these falsification thresholds fixed in advance:
- **Kill** if, at defensible burden (M ≤ 300 occupancy-competent loci) and
  physiological RNP exposure (τ ≥ 24 h), the on-target yield deficit is < 2 pp.
- **Kill** if the effect requires simultaneously M > 1000 **and** Cas9 < 300
  molecules **and** τ < 12 h (triple-implausible conjunction).
- **Required first measurement (not simulation):** the number of
  occupancy-competent off-target loci per guide relative to active nuclear
  Cas9-RNP molecules (2M vs Nc). If unmeasurable, the pivot cannot be validated
  and defaults back to rejection.

Absent evidence that real guides + real delivery reach the M~1000 / short-τ
corner, this pivot is **narrow and provisional**, one honest experiment away from
`NO_CREDIBLE_PIVOT`.
