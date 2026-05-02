# crispr-kinetic-model
A time-resolved kinetic simulation framework for modeling CRISPR-Cas9 on-target and off-target dynamics.
# CRISPR Kinetic Model

A time-resolved simulation framework for modeling CRISPR-Cas9 binding, cleavage, and competition dynamics across genomic sites.

---

## Overview

Current CRISPR guide RNA design tools predict static efficiency and off-target scores. However, CRISPR editing is inherently a dynamic process involving repeated binding, unbinding, and cleavage events over time.

This project introduces a **kinetic modeling approach** to simulate:

* Time-dependent on-target editing
* Off-target accumulation
* Competition between genomic sites for Cas9 activity

Instead of predicting a single score, the model generates **time-resolved trajectories** of editing outcomes.

---

## Motivation

Existing tools (e.g., DeepCRISPR and similar models) treat CRISPR as a static prediction problem:

```
gRNA → efficiency score
```

In reality:

```
Cas9 + genome → dynamic interactions over time → editing outcomes
```

This project aims to bridge that gap by introducing **time as a core variable** in CRISPR modeling.

---

## Key Features

* ⏱️ Time-resolved simulation of CRISPR activity
* ⚔️ Competition between on-target and off-target sites
* 🔁 Repeated binding–unbinding–cleavage cycles
* 📊 Dynamic specificity metrics (on-target vs off-target over time)
* 🧠 Mechanistic modeling instead of purely statistical prediction

---

## Model Concept

Each genomic site is modeled using three key kinetic parameters:

* **k_on** → binding rate
* **k_off** → unbinding rate
* **k_cut** → cleavage rate

The system evolves over time as Cas9 molecules interact with multiple competing DNA sites.

---

## Expected Outputs

* On-target cleavage vs time
* Off-target cleavage vs time
* Time-dependent specificity curves
* Optimal exposure time estimation

---

## Why This Matters

CRISPR-based editing is often time-limited in real applications. Static models cannot capture how off-target effects accumulate over time.

This framework enables:

* Better guide RNA selection under time constraints
* Improved understanding of safety trade-offs
* Mechanistic insight into CRISPR dynamics

---

## Project Structure

```
/src        → core simulation code  
/models     → kinetic model implementations  
/experiments → simulation experiments  
/data       → input sequences / synthetic datasets  
/plots      → generated graphs and visualizations  
```

---

## Status

🚧 Initial development phase — building core kinetic simulation engine.

---

## Future Work

* Integrating sequence-based rate estimation
* Incorporating chromatin accessibility (epigenetics)
* Comparing with existing CRISPR prediction tools
* Extending to stochastic simulation models

---

## Author

Priyansh kuniyal

---

## License

MIT License
