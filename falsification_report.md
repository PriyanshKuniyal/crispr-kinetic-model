# Falsification Report — Competitive Finite-Pool SpCas9 Hypothesis

**Status: REJECTED for practical biological relevance.** Frozen 2026-07-07.
Kinetics unchanged from the Eslami-Mossallam et al. (2022) Nature Communications
model throughout. No parameters were tuned to reach this conclusion.

## The original hypothesis

That modelling all genomic sites as drawing from a single *finite, shared* Cas9
pool — free concentration `C_free(t) = C_total − Σ bound` scaling every site's
PAM on-rate — produces a **biologically meaningful** difference from the original
independent-site, constant-Cas9 model: specifically an improvement in on-vs-off
target **specificity** driven by competition.

## Why it failed (three independent phases)

**Phase 1 — GPU go/no-go campaign** (`artifacts/go_no_go/`, verdict
`PIVOT_EFFECT_ONLY_IN_EXTREME_REGIMES`). Competition rescales *every* site's
on-rate by the same `C_free(t)`, so it slows on- and off-target binding together
and does not preferentially protect the on-target. Consequences measured across
M∈{1..10⁴}, C∈{0.1..100 nM}:
- Specificity change ≤ **10.67%** and only at M=10⁴; ≤1.5% at all realistic M.
- **No robust ranking inversion** anywhere (Spearman 0.9999; the lone flag was a
  tie flip within numerical tolerance).
- On-target "differences" are **transient delays**, not outcome changes: final
  cleavage is ~1.0000 in both models everywhere except a 0.6 pp gap at the single
  most extreme corner.

**Phase 2 — Panel-construction audit** (`artifacts/go_no_go/panel_audit_*`,
verdict `SIGNAL_ROBUST_TO_PANEL_COMPOSITION`). The surviving large-M signal is a
genuine mass-action *burden/titration* effect, not an artifact of injected strong
sinks: removing the near-cognate sinks retains 98.7–99.8% of the 17 pp signal;
the near-cognate class holds only 0.8% of sequestered Cas9 (the aggregate of
thousands of weak/relevant sites dominates). So the effect is real — but it is
merely enzyme titration by DNA, carrying no specificity information.

**Phase 3 — Molecule-count scaling** (`artifacts/plausibility/`, verdict
`PIVOT_NARROW_PLAUSIBLE_REGIME`). Replacing the universal `S_site = 0.00634 nM`
with per-volume copy-number accounting `S_site = 2/(N_A·V)` and expressing Cas9
as molecule counts shows the effect is governed by the ratio **2M/Nc**
(competitor target-copies per Cas9 molecule), essentially volume-independent.

## Strongest plausible-regime effect

Inside a generous plausibility envelope (occupancy-competent loci M ≤ 1000,
active Cas9 ≥ 100 molecules), the maximum effect anywhere (all nuclear volumes)
is **4.39 pp transient / 6.3% t50 delay / ~6% depletion**, at the very corner
(M = 1000, Cas9 = 100 molecules, 250 fL nucleus). At a defensible burden
(M ≤ 300) every plausible cell is **< 0.7 pp** and **< 2% t50 delay**.

## Threshold burden required to cross meaningful bars

| meaningful bar | minimum competitor loci M | at Cas9 / notes |
|---|---|---|
| **> 5 pp transient on-target** | **≥ 3000 occupancy-competent loci** | e.g. 10⁴ Cas9 molecules (66 nM), 250 fL; nothing below M=3000 reaches 5 pp at any Cas9/volume |
| **> 10% t50 delay** | **≥ 1000 loci** (first appears; ~M=1000 gives ~6%, robust ≥10% at M≥3000) | requires scarce Cas9 |
| **> 25% t50 delay** | **≥ 3000 loci** | e.g. 10³ Cas9 molecules (6.6 nM) |
| (1 pp transient, tripwire only) | ≥ 300 loci | — |

3000+ *occupancy-competent* off-target loci for a single guide has no empirical
support (GUIDE-seq/CIRCLE-seq detect tens–hundreds of cleavage off-targets).

## Why final editing remains essentially unchanged

The competitive coupling enters only through `C_free(t)`, a common scalar
multiplying the (slow, rate-limiting) PAM on-rate of every site. It changes *how
fast* the on-target reaches cleavage but not *whether* it does: as t→∞ the
on-target sink is absorbing and irreversible, so on-target cleaved fraction → 1.0
in both models for every condition tested (final gap ≤ 0.6 pp even at M=10⁴).
Competition redistributes *timing*, not *endpoint* — and specificity (the on/off
ratio) is nearly invariant because the same `C_free(t)` divides out.

## Why the extreme delays may be biologically moot

The large % delays live where Cas9 is scarce, where the *absolute* on-target t50
is already 10⁵–10⁷ s (hours to months). A "58% t50 delay" on a process that
takes weeks unfolds over a window far longer than: (i) Cas9-RNP active lifetime
(~hours–2 days), and (ii) the cell cycle (~1 day). A delay that only manifests
after the enzyme is gone and the cell has divided cannot change the observed
editing outcome. This is the crux that motivates the pivot audit: the transient
delay is only relevant *if* it acts within a finite Cas9-exposure window.

## Controlling unknown (now singular)

Relevance reduces to **one composite unmeasured quantity**: whether the number of
*occupancy-competent* off-target loci per guide can rival the number of *active
nuclear Cas9-RNP molecules* (2M ≳ Nc). Total NGG PAM sites (~4×10⁸) are NOT this
count. Until that ratio is measured (not assumed), the steady-state finite-pool
hypothesis stays rejected.

## Verdict

The steady-state finite-pool competition hypothesis is **falsified as a
specificity mechanism** and **rejected for practical relevance**: at realistic
molecule counts it produces only weak, outcome-neutral transient delays. See
`pivot_audit.md` for whether any mathematically-distinct phenomenon survives.
