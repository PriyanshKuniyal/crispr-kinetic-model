# Panel-Construction Audit of the M=1000-10000 Competitive Signal

**Verdict: `SIGNAL_ROBUST_TO_PANEL_COMPOSITION`**

## Baseline reproduction (fresh deterministic panel realization)

| condition | dOnMax(pp) | timing% | depletion% |
|---|---|---|---|
| DEFENSIBLE_M10000 | 17.15 | 58.8 | 37.76 |
| DEFENSIBLE_M1000 | 2.19 | 6.0 | 5.94 |

_Campaign headline was 17.29 pp / 59% / 38% at DEFENSIBLE M=10000, C=0.1._

## Sequestration fraction by class (who holds the scarce Cas9, C=0.1 nM)

### DEFENSIBLE_M10000
- class multiplicities: {'on_target': 1, 'relevant': 1800, 'sink': 200, 'weak': 7999}
- mismatch-count histogram: {'0': 1, '1': 110, '2': 90, '4': 3335, '3': 907, '5+': 5557}
- near-perfect (<=1 mm) sites: 110; distal-only near-cognates: 573
- **sequestration fraction by class:** {'on_target': 0.0, 'relevant': 0.403, 'sink': 0.008, 'weak': 0.589}

### DEFENSIBLE_M1000
- class multiplicities: {'on_target': 1, 'relevant': 180, 'sink': 20, 'weak': 799}
- mismatch-count histogram: {'0': 1, '2': 7, '1': 13, '4': 295, '3': 100, '5+': 584}
- near-perfect (<=1 mm) sites: 13; distal-only near-cognates: 61
- **sequestration fraction by class:** {'on_target': 0.0, 'relevant': 0.387, 'sink': 0.018, 'weak': 0.595}

### CONTROLLED_STRESS_M10000
- class multiplicities: {'distal': 2000, 'mid': 2000, 'near_perfect': 1999, 'on_target': 1, 'seed': 2000, 'weak': 2000}
- mismatch-count histogram: {'0': 1, '1': 5012, '2': 2987, '4': 975, '5+': 1025}
- near-perfect (<=1 mm) sites: 5012; distal-only near-cognates: 4022
- **sequestration fraction by class:** {'distal': 0.382, 'mid': 0.169, 'near_perfect': 0.003, 'on_target': 0.0, 'seed': 0.19, 'weak': 0.255}

### CONTROLLED_STRESS_M1000
- class multiplicities: {'distal': 200, 'mid': 200, 'near_perfect': 199, 'on_target': 1, 'seed': 200, 'weak': 200}
- mismatch-count histogram: {'0': 1, '1': 504, '2': 295, '5+': 107, '4': 93}
- near-perfect (<=1 mm) sites: 504; distal-only near-cognates: 402
- **sequestration fraction by class:** {'distal': 0.348, 'mid': 0.168, 'near_perfect': 0.003, 'on_target': 0.0, 'seed': 0.192, 'weak': 0.289}

## Ablation results (C=0.1 nM)

| panel | M | ablation | removed | M_eff | dOnMax(pp) | dOnFinal(pp) | timing% | depl% |
|---|---|---|---|---|---|---|---|---|
| DEFENSIBLE | 10000 | baseline | 0 | 10000 | 17.15 | -0.568 | 58.8 | 37.76 |
| DEFENSIBLE | 10000 | drop_near_perfect_le1mm | 110 | 9890 | 17.12 | -0.567 | 58.6 | 37.69 |
| DEFENSIBLE | 10000 | drop_1_2mm | 200 | 9800 | 16.92 | -0.558 | 57.7 | 37.27 |
| DEFENSIBLE | 10000 | drop_pam_distal_only | 573 | 9427 | 13.72 | -0.333 | 44.8 | 31.37 |
| DEFENSIBLE | 10000 | drop_top1pct_sink | 100 | 9900 | 16.19 | -0.487 | 54.6 | 36.07 |
| DEFENSIBLE | 10000 | drop_top5pct_sink | 500 | 9500 | 13.31 | -0.295 | 43.3 | 30.77 |
| DEFENSIBLE | 10000 | drop_top10pct_sink | 1000 | 9000 | 11.89 | -0.230 | 38.0 | 27.96 |
| DEFENSIBLE | 10000 | drop_random10pct_CTRL | 999 | 9001 | 15.72 | -0.451 | 52.6 | 35.22 |
| DEFENSIBLE | 1000 | baseline | 0 | 1000 | 2.19 | -0.002 | 6.0 | 5.94 |
| DEFENSIBLE | 1000 | drop_near_perfect_le1mm | 13 | 987 | 2.19 | -0.002 | 6.0 | 5.91 |
| DEFENSIBLE | 1000 | drop_1_2mm | 20 | 980 | 2.13 | -0.002 | 5.8 | 5.76 |
| DEFENSIBLE | 1000 | drop_pam_distal_only | 61 | 939 | 1.66 | 0.002 | 4.5 | 4.48 |
| DEFENSIBLE | 1000 | drop_top1pct_sink | 10 | 990 | 2.04 | -0.001 | 5.6 | 5.55 |
| DEFENSIBLE | 1000 | drop_top5pct_sink | 50 | 950 | 1.61 | 0.003 | 4.4 | 4.43 |
| DEFENSIBLE | 1000 | drop_top10pct_sink | 100 | 900 | 1.38 | 0.004 | 3.8 | 3.83 |
| DEFENSIBLE | 1000 | drop_random10pct_CTRL | 99 | 901 | 1.99 | -0.000 | 5.5 | 5.42 |
| CONTROLLED_STRESS | 10000 | baseline | 0 | 10000 | 13.96 | -0.220 | 46.5 | 33.51 |
| CONTROLLED_STRESS | 10000 | drop_near_perfect_le1mm | 5012 | 4988 | 11.60 | -0.169 | 37.3 | 28.15 |
| CONTROLLED_STRESS | 10000 | drop_1_2mm | 7999 | 2001 | 3.31 | -0.014 | 9.3 | 8.69 |
| CONTROLLED_STRESS | 10000 | drop_pam_distal_only | 4022 | 5978 | 7.92 | -0.079 | 24.0 | 20.03 |
| CONTROLLED_STRESS | 10000 | drop_top1pct_sink | 100 | 9900 | 12.99 | -0.174 | 42.8 | 31.71 |
| CONTROLLED_STRESS | 10000 | drop_top5pct_sink | 500 | 9500 | 10.07 | -0.088 | 31.9 | 26.03 |
| CONTROLLED_STRESS | 10000 | drop_top10pct_sink | 1000 | 9000 | 8.71 | -0.064 | 26.9 | 23.12 |
| CONTROLLED_STRESS | 10000 | drop_random10pct_CTRL | 999 | 9001 | 12.80 | -0.178 | 42.0 | 31.19 |

## Ranking-inversion audit (DEFENSIBLE, M=1000, C=1 nM)

- Spearman(off-target final cleavage, A vs B) = 0.999874
- inversions above margin 1e-9: 0
- max |A-B| final off-target cleavage = 1.975e-02
- numerical tolerance = 1e-06

_No inversions above a 1e-9 margin in this realization._

## Panel provenance (honesty check)

The `epsilon` energies and forward rates ARE literature-derived (fitted to the
Eslami-Mossallam et al. 2022 Nature Communications model). The competitor
**mismatch positions and class multiplicities are NOT literature-derived** — they
are algorithmically **sampled** from position ranges chosen to mimic seed / mid /
distal biology:

- CONTROLLED_STRESS multiplicities: 20/20/20/20/20% (near-perfect 1 mm @18-20 /
  seed 1-2 mm @1-8 / mid 1-2 mm @9-12 / distal 1-2 mm @13-20 / weak 4-5 mm). This
  panel DELIBERATELY over-represents strong binders and is an adversarial stress
  case, not a biological claim.
- DEFENSIBLE multiplicities: 2% sink (1-2 mm @15-20) / 18% relevant (3-4 mm
  @9-20, no seed mismatch) / 80% weak (seed + distal, or 5 spread mm). These
  fractions are the project's own `generate_population` heuristic, motivated by
  "most genome-wide off-targets carry many mismatches," but they are a **synthetic
  construction, not fitted to any specific guide/genome or GUIDE-seq/CIRCLE-seq
  dataset.**

Caveat on the label "DEFENSIBLE": the mismatch *positional biology* is defensible;
the *2%/18%/80% multiplicities and M itself* are assumptions, not evidence.

## Who holds the scarce Cas9 (C=0.1 nM)

Sequestration fraction by class, DEFENSIBLE M=10000: weak 58.9%, relevant 40.3%,
**near-cognate "sink" 0.8%**, on-target 0%. The Cas9 is soaked up by the
*aggregate* of ~9800 moderate/weak sites, NOT by the rare near-cognate sinks the
panel injects. Single-site residence times: on-target tau ~ 9.1e5 s, one distal
mismatch tau ~ 8.4e4 s, one seed mismatch tau ~ 1.5e3 s (a distal-mismatch site
holds Cas9 ~57x longer than a seed-mismatch site) — yet there are too few
near-cognates for them to matter in aggregate.

## Ablation verdict logic

DEFENSIBLE M=10000 baseline dOnMax = 17.15 pp (reproduces the campaign's 17.29 pp
on an independent random draw -> not seed-luck). Fraction of signal retained:

- remove near-perfect (<=1 mm): 0.998
- remove all 1-2 mm near-cognates: 0.987
- remove strongest 1% by occupancy: 0.944
- remove RANDOM 10%: 0.917
- remove strongest 10% by occupancy: 0.693
- remove PAM-distal-only near-cognates: 0.800

Removing the injected strong sinks does essentially nothing; random subsetting
reduces the signal ~proportionally to the count removed; only stripping the
top-occupancy 10% (a broad shoulder of ~1000 "relevant" sites) removes ~31%, and
69% still survives. This is a broad mass-action **burden/titration** effect,
robust to composition — the opposite of a rare-strong-sink artifact.

By contrast CONTROLLED_STRESS M=10000 collapses to 24% when all 1-2 mm sites are
removed (retain 0.237), because it is *built* from ~80% strong binders. That panel
IS composition-dependent by design and should not be used to argue biology.

## Ranking inversion audit (DEFENSIBLE, M=1000, C=1 nM)

Spearman(off-target final cleavage, independent vs competitive) = 0.9999; **zero**
inversions above a 1e-9 margin in the independent realization. The campaign's lone
flagged inversion (also Spearman 0.9999) is a near-degenerate tie flip within
numerical tolerance, NOT a robust reordering. Competition rescales all binding
uniformly and preserves site rank order.

## Verdict: SIGNAL_ROBUST_TO_PANEL_COMPOSITION

The surviving M=1000-10000 signal is a **genuine, composition-robust, mass-action
enzyme-titration effect**, not a construction artifact and not driven by rare
strong sinks. It reproduces on independent draws, survives removal of the injected
near-cognates, and scales with the *number* of competing sites.

## Project decision: PIVOT, do not kill

- **Do not kill on artifact grounds.** The effect is real and robust to how the
  panel is composed.
- **Pivot specifically to "transient kinetic delay caused by high off-target
  burden."** The signal is a delay (final cleavage essentially unchanged) driven
  by the aggregate DNA-site burden depleting a scarce Cas9 pool.
- **The remaining thing to justify is the burden MAGNITUDE, not the composition.**
  The effect needs M x S_site (63 nM of PAM-competent DNA at M=10000) to rival a
  sub-nM to low-nM Cas9 pool. Whether ~10^4 occupancy-competent off-target loci
  per guide and <=1 nM intracellular Cas9-RNP co-occur is the biological question
  that decides real-world relevance — and it is independent of this panel audit.
- **Drop the CONTROLLED_STRESS panel from any biological argument** (composition
  artifact by construction); keep it only as a labelled worst-case stress bound.
