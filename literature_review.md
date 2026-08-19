# Literature Review — Finite-Pool Cas9 Competition & the Specificity-Invariance Paper

**Compiled:** 2026-07-11 · Sources: PubMed, Consensus (bioRxiv/Semantic Scholar).
**Purpose:** Ground the planned paper (shared-scalar competition rescales *timing*, not
specificity; physiological bound via 2M/Nc). Each entry notes **how we use it**.

---

## Bucket A — Base kinetic model & its lineage (what we build on)

- **Eslami-Mossallam, Klein, van der Smagt, van der Sanden, Jones, Hawkins, Finkelstein, Depken (2022).**
  *A kinetic model predicts SpCas9 activity, improves off-target classification, and reveals the
  physical basis of targeting fidelity.* **Nat Commun 13:1367.**
  DOI: https://doi.org/10.1038/s41467-022-28994-2
  → **THE base model.** Our reconstruction (`biophysical_reconstruction.md`) reproduces this exact
  free-energy/master-equation framework and its parameters. Cite as the single-target model we extend.

- **Klein, Eslami-Mossallam, Gonzalez Arroyo, Depken (2018).** *Hybridization Kinetics Explains
  CRISPR-Cas Off-Targeting Rules.* **Cell Rep 22(6):1413.**
  DOI: https://doi.org/10.1016/j.celrep.2018.01.045
  → Predecessor 4-parameter kinetic model. Establishes the kinetic-origin-of-specificity framing our
  invariance argument extends. Same group as the base model.

- **Jones, Hawkins, Johnson, … Press, Finkelstein (2020).** *Massively parallel kinetic profiling of
  natural and engineered CRISPR nucleases (NucleaSeq).* **Nat Biotechnol 39(1):84.**
  DOI: https://doi.org/10.1038/s41587-020-0646-5
  → Source of the cleavage-kinetics training data behind the base model's parameters. Cite for
  parameter provenance; note "cleavage specificity ≠ binding specificity," which supports our
  timing-vs-endpoint distinction.

- **Fu, He, Dou, … Depken, Xu (2022).** *Systematic decomposition of sequence determinants governing
  CRISPR/Cas9 specificity (MOFF).* **Nat Commun 13:474.**
  DOI: https://doi.org/10.1038/s41467-022-28028-x
  → Independent multi-state kinetic model + "epistasis-like" multi-mismatch effects tied to the R-loop
  free-energy landscape. Cite to show the free-energy-landscape approach is field-standard, not ours alone.

- **Farasat & Salis (2016).** *A Biophysical Model of CRISPR/Cas9 Activity for Rational Design of
  Genome Editing and Gene Regulation.* **PLoS Comput Biol 12(1):e1004724.**
  DOI: https://doi.org/10.1371/journal.pcbi.1004724
  → ⚠️ **CLOSEST PRIOR ART / NOVELTY THREAT.** Already a "system-wide" biophysical model predicting
  genome-wide off-target *binding* and the effect of **Cas9/crRNA expression levels**, and explains why
  off-target activity can be high. We MUST distinguish: they compute (largely static, thermodynamic)
  binding occupancy per site; they do **not** analyze finite-pool depletion *dynamics* across a shared
  pool, nor prove the specificity-invariance under common-scalar coupling, nor derive the 2M/Nc bound.
  Position our contribution explicitly against this paper in the Intro and Discussion.

---

## Bucket B — "The hypothesis was live" (dose / titration / context-dependent fidelity)

- **Hsu, Scott, Weinstein, … Zhang (2013).** *DNA targeting specificity of RNA-guided Cas9 nucleases.*
  **Nat Biotechnol.** Consensus: https://consensus.app/papers/details/25aef8fd99e35ce5949e88734f250d3c/
  → **Key folklore citation.** States explicitly that "the dosage of SpCas9 and sgRNA can be titrated to
  minimize off-target modification." This is the intuitive belief our paper tests and bounds. Lead the
  Intro with it.

- **Slaymaker, Gao, … Zhang (2015).** *Rationally engineered Cas9 nucleases with improved specificity
  (eSpCas9).* **Science.** Consensus: https://consensus.app/papers/details/039c61e86b3b566f9ca14393d50275c0/
  → Shows specificity is treated as an engineerable enzyme property. Contrast with our finding that
  resource competition is *not* a lever on the on/off ratio.

- (Supporting) **Farasat & Salis 2016** (above) also explicitly frames Cas9 expression level as a
  fidelity knob — doubles as "hypothesis was live" evidence.

---

## Bucket C — Competition / sequestration / decoy precedent (direct)

- **Coelho, Bryan, … (2020).** *CRISPR GUARD protects off-target sites from Cas9 nuclease activity
  using short guide RNAs.* **Nat Commun.**
  Consensus: https://consensus.app/papers/details/ebbcd23eba57535ab1129de2c5d1ef41/
  → **Direct precedent that competition among sites is real and exploitable:** co-delivered short guides
  compete for off-target loci and reduce off-target editing while sparing on-target. Establishes that
  the *mechanism* we model (site competition for Cas9) genuinely operates — strengthens that our null is
  a meaningful, non-strawman result about *magnitude*, not existence.

- **Kuscu, Arslan, … Adli (2014).** *Genome-wide analysis reveals characteristics of off-target sites
  bound by the Cas9 endonuclease.* **Nat Biotechnol.**
  Consensus: https://consensus.app/papers/details/0b26004566815088b1b8282cbc7b4269/
  → ⚠️ **KEY EMPIRICAL ANCHOR AND THE PAPER'S SHARPEST TENSION.** dCas9 ChIP-seq shows the number of
  bound off-target sites ranges **~10 to >1000 depending on sgRNA.** Our whole falsification hinges on
  M (occupancy-competent loci). For promiscuous guides, *binding* can approach/exceed the M~1000 corner.
  Must address head-on: (i) most are weak/transient binding events, not sequestration-competent; (ii) our
  panel audit shows the aggregate of many weak sites is exactly what matters, so this is the crux, not a
  loophole; (iii) it makes the "measure 2M vs Nc" call concrete. Turn this threat into the central hook.

---

## Bucket D — Bounding M (genome-wide off-target counts by assay)

- **Kleinstiver, Prew, Tsai, … Joung (2015).** *Engineered CRISPR-Cas9 nucleases with altered PAM
  specificities.* **Nature 523:481.** DOI: https://doi.org/10.1038/nature14592
  → GUIDE-seq genome-wide specificity benchmark; cleavage off-targets number in the tens.

- **Malinin, Lee, Lazzarotto, … Tsai (2021).** *Defining genome-wide CRISPR-Cas genome-editing nuclease
  activity with GUIDE-seq.* **Nat Protoc 16:5592.** DOI: https://doi.org/10.1038/s41596-021-00626-x
  → Method reference; cellular off-target catalogs are typically tens–low hundreds of sites per guide.

- **Lazzarotto, Malinin, … Tsai (2020).** *CHANGE-seq reveals genetic and epigenetic effects on
  CRISPR-Cas9 genome-wide activity.* **Nat Biotechnol 38:1317.**
  DOI: https://doi.org/10.1038/s41587-020-0555-7
  → In vitro upper bound: 201,934 off-target sites across 110 guides (~1,800/guide) — but deproteinized,
  maximally sensitive; cellular activity concentrates near open chromatin/promoters and is far lower.
  Use to frame the gap between *in vitro detectable* vs *occupancy-competent in a nucleus*.

---

## Bucket E — Bounding Nc & RNP lifetime (residence time / single-turnover)

- **Richardson, Ray, DeWitt, Curie, Corn (2016).** *Enhancing homology-directed genome editing … using
  asymmetric donor DNA.* **Nat Biotechnol 34:339.** DOI: https://doi.org/10.1038/nbt.3481
  → **Cas9 dissociation from dsDNA is slow, lifetime ~6 h.** This is exactly the number the falsification
  uses to bound RNP active lifetime; it supports both (i) the sequestration premise and (ii) the
  finite-lifetime pivot's kill-condition (τ ~ hours). Central empirical peg for the pivot.

- **Kiernan et al. (2025).** *Visualization of a multi-turnover Cas9 after product release.* **Nat Commun.**
  Consensus: https://consensus.app/papers/details/9254b061f188561c8064544d1f39e202/
  → Confirms Cas9 is effectively single-turnover with long product residence (PAM-proximal product stays
  bound, blocks re-binding). Justifies the model's single-turnover / bound-post-cleavage assumption and
  the finite active-Cas9 premise.

---

## Bucket F — Single-molecule mechanism support (validates the model's internals)

- **Gong, Yu, Johnson, Taylor (2018).** *DNA Unwinding Is the Primary Determinant of CRISPR-Cas9
  Activity.* **Cell Rep 22:359.** DOI: https://doi.org/10.1016/j.celrep.2017.12.041
  → R-loop formation is rate-limiting and reversible — supports the slow PAM/on-rate step through which
  our shared C_free(t) scalar enters. Mechanistic backing for why competition acts on *timing*.

- **Ivanov, Wright, Cofsky, … Doudna, Bryant (2020).** *Cas9 interrogates DNA in discrete steps modulated
  by mismatches and supercoiling.* **PNAS 117:5853.** DOI: https://doi.org/10.1073/pnas.1913445117
  → Single-molecule evidence for the discrete intermediate R-loop state the base model predicts.
  Independent validation of the state structure we inherit.

---

## Strategic synthesis (for the paper)

1. **The base model and its provenance are rock-solid and precisely citable** (Eslami-Mossallam 2022 +
   NucleaSeq training data, same group). The reconstruction stands.

2. **"The hypothesis was live" is well-supported** — Hsu 2013 (titrate dose to cut off-targets), the
   high-fidelity-variant program (Slaymaker 2015), and Farasat & Salis 2016 (Cas9-level as a fidelity
   knob). The Intro can credibly claim resource competition was expected to matter.

3. **Two things must be handled explicitly or a reviewer will pounce:**
   - **Farasat & Salis 2016** — closest prior systems-level model. Differentiate on *dynamics +
     invariance + bound*, not just "we also modeled competition."
   - **Kuscu 2014 (10 to >1000 bound sites)** — the empirical number nearest our M~1000 threshold. This
     is simultaneously our biggest vulnerability and our sharpest hook: it makes the 2M/Nc measurement
     call concrete and non-hypothetical.

4. **The finite-lifetime pivot's kill-condition (τ ~ hours) is empirically grounded** by Richardson 2016
   (~6 h) and Kiernan 2025 (single-turnover), so rejecting it is defensible, not convenient.

## Gaps / still worth a targeted search before writing
- A clean cellular estimate of **active nuclear Cas9-RNP molecule counts (Nc)** per cell — needed to
  anchor the 2M/Nc physiological region quantitatively.
- Any prior **molecular-titration / decoy ("sponge") theory** from gene-regulation (ceRNA, TF-decoy)
  we can cite as the general principle our invariance specializes.
- Whether any group has **already reported a competition/finite-pool CRISPR model** post-2022 (protect
  the novelty claim).
