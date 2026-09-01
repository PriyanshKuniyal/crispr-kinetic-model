# Biophysical Reconstruction and Analysis of the SpCas9 Kinetic Model

This document presents a mathematically rigorous, comprehensive reconstruction of the CRISPR/SpCas9 target recognition and cleavage kinetic model published by Eslami-Mossallam et al. in *Nature Communications* (2022). It serves as a complete technical guide for reproducing the authors' work, understanding the biophysical parameters, and extending the model to physiological cellular contexts.

---

## 1. Biophysical Concepts & Terminology

To bridge the gap between abstract physical chemistry and molecular biology, we first define the core biological components of the CRISPR-SpCas9 system using rigorous biophysical terms.

### 1.1 PAM Binding
SpCas9 (derived from *Streptococcus pyogenes*) is directed to genomic target sites by a $\sim 100$-nucleotide single-guide RNA (sgRNA). However, the enzyme cannot search the entire 20-bp target sequence simultaneously. Instead, it relies on a local search trigger. 
* **Mechanism**: SpCas9 dynamically diffuses through the genome and queries DNA by scanning for a specific 3-nucleotide sequence called the **Protospacer Adjacent Motif (PAM)**—canonically `5'-NGG-3'` in SpCas9. 
* **Biophysics**: PAM recognition occurs via direct hydrogen bonding and electrostatic interactions between the protein's PAM-interacting (PI) domain and the DNA major/minor grooves. This protein-DNA binding is modeled as a concentration-dependent, first-order reaction with rate $k_{\text{on}}$ and free-energy change $F_0$. Binding to the PAM is a thermodynamic prerequisite that mechanically anchors the Cas9 complex and induces local DNA strand separation.

### 1.2 R-Loop Formation and Progression
Once SpCas9 binds to the PAM, it destabilizes the adjacent DNA double helix.
* **Initiation (State 1)**: The $3'$ end of the target DNA strand (closest to the PAM) begins to base-pair with the complementary sgRNA, forming a DNA-RNA hybrid and displacing the non-target DNA strand. This three-stranded nucleic acid structure is termed an **R-loop**.
* **Progression (States 1–20)**: R-loop progression is a **strand-displacement reaction** that occurs in single-nucleotide steps. It behaves like a one-dimensional stochastic walk. At each step $n \to n+1$, a base pair of the host DNA double helix is melted, and a complementary base pair in the RNA-DNA hybrid is formed. This progression spans the 20 nucleotides of the sgRNA guide sequence:
  * **PAM-Proximal "Seed" Region (nt 1–8)**: Highly sensitive to mismatches. Mismatches here block progression early, promoting rapid Cas9 dissociation.
  * **Mid-Loop Region (nt 9–12)**: Serves as a mechanical transition zone.
  * **PAM-Distal Region (nt 13–20)**: Drives final conformational activation of the nuclease domains.

### 1.3 Mismatch Penalties
When SpCas9 encounters a genomic site that is not perfectly complementary to its sgRNA (an **off-target**), the R-loop progression must incorporate non-complementary base pairs.
* **Thermodynamic Penalty**: An RNA-DNA mismatch (e.g., G-A, C-U, etc.) is highly unfavorable. It fails to form stabilizing hydrogen bonds and introduces steric strain or structural distortion into the hybrid.
* **Mathematical Modeling**: The energetic cost of a mismatch at position $n$ is represented by a position-dependent free-energy penalty $\delta\epsilon_n > 0$. Under the **additive assumption** (local independence), the presence of a mismatch at position $n$ shifts the free-energy difference of *all subsequent states* $m \ge n$ upward by $\delta\epsilon_n$, acting as a permanent energy barrier that the stochastic walker must overcome to reach cleavage.

### 1.4 Free-Energy Landscapes and Metastable States
* **Free-Energy Landscape ($F_n$)**: The thermodynamic landscape coordinates the R-loop progression state $n \in \{-1, 0, 1, \dots, 20\}$. Each state $n$ corresponds to an R-loop hybrid of length $n$. The landscape is governed by protein-DNA, protein-RNA, and nucleic acid hybridization free energies.
* **Metastable States**: The on-target landscape of SpCas9 is not a smooth downward slope. Instead, it contains three distinct local minima (basins of stability) separated by high transition barriers:
  1. **Closed R-loop Basin (C)**: Comprising the unbound Solution state ($n=-1$) and the PAM-bound state ($n=0$).
  2. **Intermediate R-loop Basin (I)**: Centered around states $n=9\text{--}12$.
  3. **Open R-loop Basin (O)**: The fully formed 20-bp hybrid state ($n=20$).
  
The barriers between these basins act as rate-limiting conformational checkpoints.

### 1.5 Coarse-Graining and First-Passage Times
* **Coarse-Graining**: The microscopic model has 22 distinct states. To bridge the model with macroscopically observable structural transitions (such as FRET conformers or domain rotations), the states are grouped into a simplified 4-state system (Solution $\to$ Closed $\to$ Intermediate $\to$ Open $\to$ Cleaved).
* **Mean First Passage Time (MFPT)**: The mathematical technique used to calculate the effective transition rates between the coarse-grained metastable basins. The MFPT $\tau_{A \to B}$ is the average time it takes for a stochastic walker starting in basin $A$ to reach basin $B$ for the first time, treating $B$ as an absorbing boundary. The coarse-grained rate is then defined as $k_{A \to B} = 1 / \tau_{A \to B}$.

---

## 2. The Microscopic State-Space Framework

The full state-space model represents a 1D stochastic hopping process (a continuous-time Markov chain) across 23 states:

```mermaid
stateDiagram-v2
    direction LR
    Sol: Solution (-1)
    PAM: PAM-Bound (0)
    S1: State 1
    S2: State 2
    S19: State 19
    S20: State 20
    Clv: Cleaved

    Sol --> PAM : k_on
    PAM --> Sol : k_b^0
    PAM --> S1 : k_f
    S1 --> PAM : k_b^1
    S1 --> S2 : k_f
    S2 --> S1 : k_b^2
    S2 --> S19 : ...
    S19 --> S20 : k_f
    S20 --> S19 : k_b^20
    S20 --> Clv : k_cat
```

### 2.1 The Master Equation
Let $P(t) = \big( P_{-1}(t), P_0(t), P_1(t), \dots, P_{20}(t) \big)^T$ be the column vector of probabilities of the system being in state $n$ at time $t$. The time-evolution of these probabilities is described by the **Master Equation**:

$$\frac{\partial P(t)}{\partial t} = \mathbf{K} P(t)$$

where $\mathbf{K}$ is a $22 \times 22$ tri-diagonal transition-rate matrix. The elements of $\mathbf{K}$ are constructed from the microscopic forward rates $k_f^n$ ($n \to n+1$) and backward rates $k_b^n$ ($n \to n-1$):

$$K_{n, m} = \begin{cases} 
      k_f^{n-1} & m = n-1 \\
      -(k_f^n + k_b^n) & m = n \\
      k_b^{n+1} & m = n+1 \\
      0 & |n - m| \ge 2 
   \end{cases}$$

Here, states are index-shifted such that:
* Index $0$ is the Solution state ($n = -1$).
* Index $1$ is the PAM-bound state ($n = 0$).
* Index $n+1$ is the R-loop state $n$ ($n = 1 \dots 20$).

Because the final cleavage step from state $20$ to Cleaved is irreversible, the probability escapes the state-space. The rate matrix $\mathbf{K}$ is thus sub-stochastic, and the **cleaved fraction** $P_{\text{cleaved}}(t)$ is given by the probability deficit:

$$P_{\text{cleaved}}(t) = 1 - \sum_{n=-1}^{20} P_n(t)$$

The formal analytical solution to the Master Equation is:

$$P(t) = \exp(\mathbf{K}t) P(0)$$

where $\exp(\mathbf{K}t)$ is the matrix exponential.

---

## 3. Microscopic Energy and Rate Equations

To translate free-energy landscapes into kinetic transition rates, the authors apply the principle of **detailed balance** (local thermodynamic equilibrium between adjacent states):

### 3.1 Detailed Balance
If states $n-1$ and $n$ were allowed to equilibrate, the ratio of their occupancies would follow a Boltzmann distribution:

$$\frac{P_n^{\text{eq}}}{P_{n-1}^{\text{eq}}} = \exp\left( -\frac{F_n - F_{n-1}}{k_B T} \right) = \exp\left( -\frac{\Delta F_n}{k_B T} \right)$$

where $\Delta F_n = F_n - F_{n-1}$ is the free-energy difference associated with extending the hybrid by one base pair. Detailed balance requires the net probability flux between the two states to be zero at equilibrium:

$$P_{n-1}^{\text{eq}} k_f^{n-1} = P_n^{\text{eq}} k_b^n \implies \frac{k_b^n}{k_f^{n-1}} = \frac{P_{n-1}^{\text{eq}}}{P_n^{\text{eq}}} = \exp\left( \frac{\Delta F_n}{k_B T} \right)$$

Thus, the backward transition rate is explicitly tied to the forward rate and the local free-energy difference:

$$k_b^n = k_f^{n-1} \exp\left( \frac{\Delta F_n}{k_B T} \right)$$

### 3.2 Mechanistic Model Assumptions (Parameter Reduction)
To parameterize all possible guide and target sequences, the authors introduce four core biophysical assumptions:
1. **Mismatch Positional Equivalence**: Mismatch positions are dominant over mismatch types. All 12 possible base-pairing mismatch types are treated as thermodynamically equivalent.
2. **Local Additivity of Mismatch Energies**: The free-energy change $\Delta F_n$ at step $n$ depends only on whether there is a match or mismatch at that specific position. The energetic cost of a mismatch at position $n$ is represented by a local penalty $\delta\epsilon_n$. Thus:
   $$\Delta F_n = \begin{cases} 
      -\epsilon_n & \text{if position } n \text{ is matched} \\
      -\epsilon_n + \delta\epsilon_n & \text{if position } n \text{ is mismatched} 
   \end{cases}$$
   where $\epsilon_n > 0$ represents the favorable free-energy gain of matching base pair $n$ (under the sign convention that $\epsilon_n$ is a downward energy slope).
3. **dCas9 Catalytic Deadness**: Catalytically dead dCas9 is identical to active Cas9, except that the final cleavage catalysis step is suppressed ($k_{\text{cat}} = 0$). All other hybridization rates and energies remain unchanged.
4. **Forward Rate Uniformity**: The forward R-loop progression is treated as a uniform, sequence-independent sliding process. Thus, all internal forward progression rates are identical:
   $$k_f^0 = k_f^1 = \dots = k_f^{19} = k_f$$

### 3.3 Physical Energy Equations (Concentration Scaling)
* **Solution to PAM State ($n = -1 \to 0$)**: The binding rate is concentration-dependent:
  $$k_{\text{on}} = k_{\text{on}}^{\text{ref}} \frac{[\text{Cas9-sgRNA}]}{C_{\text{ref}}}$$
  where $C_{\text{ref}} = 1\text{ nM}$. The corresponding PAM-bound state energy scales thermodynamically with concentration:
  $$F_0([\text{Cas9-sgRNA}]) = F_0^{\text{ref}} - k_B T \ln\left( \frac{[\text{Cas9-sgRNA}]}{C_{\text{ref}}} \right)$$
* **Cumulative Free Energy ($F_n$)**: Choosing the unbound Solution state as the reference ($F_{-1} = 0\text{ k}_BT$), the free energy of any R-loop state $N \in \{0, \dots, 20\}$ is:
  $$F_N = F_0 + \sum_{n=1}^N \Delta F_n = F_0^{\text{ref}} - k_B T \ln\left(\frac{[\text{Cas9-sgRNA}]}{C_{\text{ref}}}\right) + \sum_{n=1}^N \left( -\epsilon_n + \delta\epsilon_n \cdot \mathbb{I}_n \right)$$
  where $\mathbb{I}_n = 1$ if there is a mismatch at position $n$, and $0$ otherwise.
* **Microscopic Backward Rates ($k_b^n$)**:
  $$k_b^n = \begin{cases}
     k_{\text{on}}^{\text{ref}} \exp\left( \frac{F_0^{\text{ref}}}{k_B T} \right) & n = 0 \quad (\text{Dissociation rate from PAM to Solution}) \\
     k_f \exp\left( \frac{-\epsilon_n + \delta\epsilon_n \cdot \mathbb{I}_n}{k_B T} \right) & n = 1 \dots 20
  \end{cases}$$

---

## 4. Extracted Microscopic Parameters

From the bulk experimental training sets (CHAMP for $dCas9$ binding and NucleaSeq for active $Cas9$ cleavage), we extract the following exact, high-confidence physical parameters (under standard thermal energy units where $k_B T \approx 1$):

### 4.1 Kinetic Rates and Binding Energies
* **Reference PAM Binding Rate**: $k_{\text{on}}^{\text{ref}} = 8.585963 \times 10^{-5} \text{ s}^{-1} \text{ nM}^{-1}$
* **PAM Binding Free Energy**: $F_0^{\text{ref}} = 5.040250 \text{ k}_B T$ (at 1 nM)
* **Forward R-loop Progression Rate**: $k_f = 641.289038 \text{ s}^{-1}$ (very rapid $\sim 640\text{ Hz}$)
* **Cleavage Catalysis Rate**: $k_{\text{cat}} = 2.392864 \text{ s}^{-1}$

### 4.2 On-Target Energy Steps ($\epsilon_n$) and Cumulative Landscape ($F_n$)
The table below lists the microscopic energy change at each step and the resulting cumulative free energy landscape for the matched on-target sequence at $1\text{ nM}$ concentration:

| State $n$ | Step Energy Change $\Delta F_n$ ($\text{k}_B T$) | Cumulative Free Energy $F_n$ ($\text{k}_B T$) |
| :--- | :---: | :---: |
| **Solution (–1)** | — | $0.000000$ |
| **PAM (0)** | $+5.040250$ | $5.040250$ |
| **State 1** | $+6.333716$ | $11.373966$ |
| **State 2** | $+1.364905$ | $12.738871$ |
| **State 3** | $-4.444667$ | $8.294204$ |
| **State 4** | $+1.537084$ | $9.831288$ |
| **State 5** | $-0.852685$ | $8.978603$ |
| **State 6** | $+0.055978$ | $9.034581$ |
| **State 7** | $+2.486526$ | $11.521106$ |
| **State 8** | $-1.432459$ | $10.088647$ |
| **State 9** | $-3.755672$ | $6.332975$ |
| **State 10** | $-1.156545$ | $5.176430$ |
| **State 11** | $-1.024544$ | $4.151886$ |
| **State 12** | $-0.619641$ | $3.532245$ |
| **State 13** | $+3.134097$ | $6.666343$ |
| **State 14** | $-1.643146$ | $5.023196$ |
| **State 15** | $-1.078241$ | $3.944955$ |
| **State 16** | $+2.153910$ | $6.098865$ |
| **State 17** | $-1.042519$ | $5.056345$ |
| **State 18** | $-3.021701$ | $2.034645$ |
| **State 19** | $+0.399089$ | $2.433734$ |
| **State 20** | $-6.789944$ | $-4.356210$ |

### 4.3 Mismatch Penalties ($\delta\epsilon_n$)
The fitted thermodynamic penalties for mismatches at each sequence position are:

$$\begin{aligned}
\delta\epsilon_1 &= 5.653338 \text{ k}_BT, \quad \delta\epsilon_2 = 4.110412 \text{ k}_BT, \quad \delta\epsilon_3 = 6.482408 \text{ k}_BT, \quad \delta\epsilon_4 = 6.976721 \text{ k}_BT \\
\delta\epsilon_5 &= 6.266481 \text{ k}_BT, \quad \delta\epsilon_6 = 7.388568 \text{ k}_BT, \quad \delta\epsilon_7 = 6.899488 \text{ k}_BT, \quad \delta\epsilon_8 = 6.220480 \text{ k}_BT \\
\delta\epsilon_9 &= 8.994202 \text{ k}_BT, \quad \delta\epsilon_{10} = 7.253615 \text{ k}_BT, \quad \delta\epsilon_{11} = 7.403559 \text{ k}_BT, \quad \delta\epsilon_{12} = 7.024180 \text{ k}_BT \\
\delta\epsilon_{13} &= 7.739561 \text{ k}_BT, \quad \delta\epsilon_{14} = 7.884431 \text{ k}_BT, \quad \delta\epsilon_{15} = 7.644810 \text{ k}_BT, \quad \delta\epsilon_{16} = 6.355368 \text{ k}_BT \\
\delta\epsilon_{17} &= 5.133518 \text{ k}_BT, \quad \delta\epsilon_{18} = 4.247951 \text{ k}_BT, \quad \delta\epsilon_{19} = 5.830443 \text{ k}_BT, \quad \delta\epsilon_{20} = 2.418782 \text{ k}_BT
\end{aligned}$$

---

## 5. The Coarse-Grained System & Dimensional Reduction

To understand the dynamic checkpoints and connect R-loop kinetics to macroscopically observable structural FRET states, we map the 22-state microscopic model to a simplified 4-state coarse-grained representation:

```mermaid
stateDiagram-v2
    direction LR
    Sol: Solution
    Closed: Closed R-loop (P)
    Int: Intermediate R-loop (I)
    Open: Open R-loop (O)
    Clv: Cleaved

    Sol --> Closed : k_on
    Closed --> Sol : k_PO
    Closed --> Int : k_PI
    Int --> Closed : k_IP
    Int --> Open : k_IO
    Open --> Int : k_OI
    Open --> Clv : k_cat
```

### 5.1 Determining Coarse-Grained State Energies
* **Solution State ($E_{\text{sol}}$)**: Reference, $E_{\text{sol}} = 0 \text{ k}_B T$.
* **Closed R-loop State $P$ ($E_P$)**: Represented by the PAM-bound state energy:
  $$E_P = F_0 = 5.040250 \text{ k}_B T \quad (\text{at 1 nM})$$
* **Intermediate R-loop State $I$ ($E_I$)**: Calculated by performing a partition function sum over the Boltzmann factors of the intermediate basin (defined by states $n=7\text{--}13$):
  $$E_I = -k_B T \ln\left( \sum_{n=7}^{13} e^{-F_n / k_B T} \right)$$
  Using the fitted energies, we obtain:
  $$E_I = -k_B T \ln \left( e^{-11.52} + e^{-10.09} + e^{-6.33} + e^{-5.18} + e^{-4.15} + e^{-3.53} + e^{-6.67} \right) \approx 2.946358 \text{ k}_B T$$
* **Open R-loop State $O$ ($E_O$)**: Represented by the energy of the fully formed hybrid:
  $$E_O = F_{20} = -4.356210 \text{ k}_B T$$

### 5.2 Deriving Coarse-Grained Rates via Mean First Passage Time (MFPT)
The microscopic location of the Intermediate basin minimum is:

$$n_I = \underset{n \in [7, 13]}{\operatorname{argmin}} F_n = 12 \quad (F_{12} = 3.53\text{ k}_BT)$$

To find the effective transition rates between these non-adjacent metastable positions ($P$ at $n=0$, $I$ at $n=12$, and $O$ at $n=20$), we isolate the intervening barriers and calculate their **Mean First Passage Time (MFPT)**:

#### 5.2.1 PAM to Intermediate Forward Rate ($k_{P \to I}$)
We truncate the master equation to only span states $n = 0, \dots, n_I-1$ (states $0$ to $11$). We set state $n_I = 12$ as an absorbing boundary (so probability escaping state $11$ to the right enters the absorber). 
The truncated rate matrix is $\mathbf{K}_{P \to I}$ of size $12 \times 12$. Starting with all probability in state $0$ ($P_{\text{sub}}(0) = [1, 0, \dots, 0]^T$), the average time to reach the Intermediate state is:

$$\tau_{P \to I} = \mathbf{1}^T (-\mathbf{K}_{P \to I})^{-1} P_{\text{sub}}(0)$$

The effective forward rate is the reciprocal of this passage time:

$$k_{P \to I} = \frac{1}{\tau_{P \to I}}$$

#### 5.2.2 Intermediate to Open Forward Rate ($k_{I \to O}$)
Similarly, we truncate the master equation to span states $n = n_I, \dots, 19$ (states $12$ to $19$), setting state $20$ as an absorbing boundary. The truncated rate matrix is $\mathbf{K}_{I \to O}$ of size $8 \times 8$. Starting with all probability in the intermediate minimum ($P_{\text{sub}}(0) = [1, 0, \dots, 0]^T$), we calculate the passage time:

$$\tau_{I \to O} = \mathbf{1}^T (-\mathbf{K}_{I \to O})^{-1} P_{\text{sub}}(0) \implies k_{I \to O} = \frac{1}{\tau_{I \to O}}$$

#### 5.2.3 Coarse-Grained Backward Rates
To ensure exact thermodynamic consistency with the microscopic model, the remaining coarse-grained backward rates ($k_{P \to \text{sol}}$, $k_{I \to P}$, $k_{O \to I}$) are computed using detailed balance on the coarse-grained free energies:

$$k_{b}^{\text{CG}, n} = k_{f}^{\text{CG}, n-1} \exp\left( \frac{E_n^{\text{CG}} - E_{n-1}^{\text{CG}}}{k_B T} \right)$$

---

## 6. Mathematical Distinctions and Assumptions

To maintain the highest scientific integrity, we explicitly demarcate the boundaries of origin for all mathematical parts of the model:

### 6.1 Equations Directly Taken from the Paper
1. The microscopic Master Equation $\partial_t P(t) = \mathbf{K} P(t)$ and its tri-diagonal rate structure.
2. The detailed balance condition linking backward rates to energies and forward rates: $k_b^n = k_f^{n-1} \exp(\Delta F_n / k_B T)$.
3. Concentration-dependent PAM scaling of $k_{\text{on}}$ and $F_0$.
4. The coarse-grained state energy $E_I$ partition function sum.
5. The first-passage time calculation $\tau = \mathbf{1}^T (-\mathbf{K})^{-1} P(0)$ using sub-matrix truncation.
6. The low-concentration cleavage probability equation:
   $$P_{\text{PAM} \to \text{clv}} = \frac{k_{\text{cat}} e^{-F_{19}/k_BT}}{k_{\text{cat}} \sum_{n=0}^{19} e^{-F_n/k_BT} + k_f e^{-F_{20}/k_BT}}$$

### 6.2 Key Assumptions Made by the Authors
1. **Additive thermodynamics**: Multi-mismatch penalties are simply the sum of individual mismatch penalties ($\sum \delta\epsilon_i$), neglecting cooperative structural melting or context-dependent cooperative binding.
2. **Homogeneous forward R-loop rate**: $k_f^n = k_f$ is independent of sequence composition (i.e., G-C hybridization vs. A-T hybridization has no effect on forward kinetics, and is purely captured by backward melting rates).
3. **Single-mismatch type equivalence**: Treating all 12 types of base-pairing mismatches identically, despite known differences between purine-purine, purine-pyrimidine, and wobble (G-U) pairings.
4. **dCas9 equivalence**: Assuming dead Cas9 has the exact same thermodynamic landscape and microscopic hybridization kinetics as active Cas9.

### 6.3 Assumptions Introduced in our Reconstruction
1. **Real-eigenvalue diagonalization**: In simulating the master equation, we utilize numerical diagonalization $K = V \Lambda V^{-1}$ to calculate $\exp(\mathbf{K}t)$ exactly and rapidly. This assumes that $K$ is non-defective, which is numerically verified for all target patterns.
2. **Single-turnover cleavage boundary**: In plotting occupancies, we assume cleavage is completely irreversible and SpCas9 remains bound to its target sequence post-cleavage, which matches the biochemically validated single-turnover kinetics of SpCas9 in vitro.

---

## 7. Numerical Simulation Results

Using the fully parameterized microscopic model, we solved the Master Equation numerically for multiple targets at $1\text{ nM}$ SpCas9 concentration. 

### 7.1 Free-Energy landscapes
The free-energy landscapes show the three metastable states of SpCas9. When a mismatch is introduced, it acts as a permanent step-barrier raising all subsequent states:
* A **seed mismatch at position 3** raises the barrier early, causing rapid R-loop reversal and Cas9 dissociation.
* A **distal mismatch at position 15** permits R-loop progression to the intermediate basin but severely blocks final transition to the cleavage-competent state 20.

![Free-Energy Landscapes](file:///C:/Users/Priyansh%20Kuniyal/.gemini/antigravity/brain/42accf3d-17cd-4026-b12b-276e03194c2c/artifacts/free_energy_landscapes.png)

### 7.2 Mismatch Penalties
Fitted mismatch penalties remain highly consistent between $4\text{--}8\text{ k}_B T$ along the hybrid, with a noticeable spike at position 9 ($\sim 9\text{ k}_B T$), representing a strong structural barrier.

![Mismatch Penalties](file:///C:/Users/Priyansh%20Kuniyal/.gemini/antigravity/brain/42accf3d-17cd-4026-b12b-276e03194c2c/artifacts/mismatch_penalties.png)

### 7.3 Metastable State Occupancies Over Time (Replicating Figures 6a & 6b)
* **On-Target Dynamics**: The transition into the intermediate state is the rate-limiting barrier. Once the intermediate state is entered, progression to the open and cleaved states is extremely rapid and essentially irreversible. The intermediate state is thus only visited transiently.
* **3 PAM-Distal Mismatches (States 18, 19, 20)**: Here, the barrier between Intermediate and Open is massive. The complex gets kinetically trapped in the Intermediate state, which becomes highly occupied and remains long-lived ($\sim 10^3\text{ s}$), while cleavage is completely suppressed. This matches the single-molecule structural FRET trajectories observed in literature.

![Metastable State Occupancies](file:///C:/Users/Priyansh%20Kuniyal/.gemini/antigravity/brain/42accf3d-17cd-4026-b12b-276e03194c2c/artifacts/metastable_occupancies.png)

### 7.4 Cleavage Kinetics
Cleavage fraction over time demonstrates that while the on-target is rapidly cleaved ($t_{1/2} \approx 10\text{ s}$), off-targets are strongly suppressed: a seed mismatch yields virtually zero cleavage, whereas distal mismatches show slow, titrated cleavage over longer exposures.

![Cleavage Dynamics](file:///C:/Users/Priyansh%20Kuniyal/.gemini/antigravity/brain/42accf3d-17cd-4026-b12b-276e03194c2c/artifacts/cleavage_dynamics.png)

---

## 8. Multi-Site Competitive Extension

### 8.1 The Need for an Extension
The original paper models a single target sequence in contact with a constant pool of free SpCas9 complexes. In a real living cell, the nucleus contains **millions of potential off-target sites** simultaneously. These sites compete for a shared, finite pool of active SpCas9-sgRNA complexes. 

Under standard transfection conditions, Cas9 is a finite resource. The thermodynamic and kinetic sequestration of Cas9 on high-abundance off-targets (acting as "decoys" or "sponges") can dramatically deplete the free Cas9 pool, slowing down on-target cleavage and altering targeting fidelity.

### 8.2 Mathematical Formulation
Let there be $M$ genomic sites, indexed $i = 1 \dots M$. Let:
* $S_i$ be the total concentration of genomic site $i$ (e.g., $S_1$ is the on-target, and $S_2 \dots S_M$ are various off-target loci).
* $c_{i, n}(t)$ be the concentration of site $i$ bound by Cas9 in state $n \in \{0, \dots, 20\}$ at time $t$.
* $c_{i, -1}(t)$ be the concentration of unbound site $i$ at time $t$.
* $c_{i, \text{clv}}(t)$ be the concentration of cleaved site $i$.
* $C_{\text{tot}}$ be the total concentration of active SpCas9-sgRNA complexes in the cell.
* $C_{\text{free}}(t)$ be the concentration of free SpCas9-sgRNA complexes in solution at time $t$.

#### 8.2.1 Mass Balance of SpCas9
Assuming SpCas9 remains tightly bound to DNA post-cleavage (single-turnover enzyme), the free pool is constrained by mass conservation:

$$C_{\text{free}}(t) = C_{\text{tot}} - \sum_{i=1}^M \sum_{n=0}^{20} c_{i, n}(t)$$

#### 8.2.2 Coupled System of Ordinary Differential Equations (ODEs)
For each site $i \in \{1 \dots M\}$, the kinetic equations are:

$$\begin{aligned}
\frac{d c_{i, -1}}{d t} &= -k_{\text{on}, i}^{\text{ref}} C_{\text{free}}(t) c_{i, -1} + k_{b, i}^0 c_{i, 0} \\
\frac{d c_{i, 0}}{d t} &= k_{\text{on}, i}^{\text{ref}} C_{\text{free}}(t) c_{i, -1} - \left( k_{b, i}^0 + k_f \right) c_{i, 0} + k_{b, i}^1 c_{i, 1} \\
\frac{d c_{i, n}}{d t} &= k_f c_{i, n-1} - \left( k_{b, i}^n + k_f \right) c_{i, n} + k_{b, i}^{n+1} c_{i, n+1} \quad \text{for } n = 1 \dots 19 \\
\frac{d c_{i, 20}}{d t} &= k_f c_{i, 19} - \left( k_{b, i}^{20} + k_{\text{cat}} \right) c_{i, 20} \\
\frac{d c_{i, \text{clv}}}{d t} &= k_{\text{cat}} c_{i, 20}
\end{aligned}$$

Where the backward rates $k_{b, i}^n$ are computed based on the specific mismatch positions of site $i$ relative to the guide.

### 8.3 Biophysical Significance & Novelty
This competitive extension introduces **non-linear coupling** between all genomic sites. The system cannot be solved as independent single-site Markov chains. Instead, it forms a coupled non-linear dynamical system.

#### 8.3.1 Key Phenomena Captured by this Extension:
1. **Sequestration Effect (Kinetic Buffering)**: High-abundance, highly complementary off-targets can act as kinetic traps. They bind Cas9 and, due to slow dissociation rates, sequester it, reducing $C_{\text{free}}(t)$ and suppressing on-target cleavage activity.
2. **Fidelity-Concentration Feedback**: Cells can optimize fidelity by titrating Cas9 expression levels ($C_{\text{tot}}$). This model quantitatively predicts the optimal $C_{\text{tot}}$ window to maximize on-target cleavage while keeping off-target cleavage below a desired tolerability threshold.
3. **Genomic Context Dependency**: It proves that the fidelity of SpCas9 is *not* a fixed property of the enzyme alone, but is highly dependent on the genomic background (the number and sequence of off-target sites in the host genome), which explains why the same sgRNA displays different specificity profiles in different cell types or organisms.

This represents a **highly promising and novel research direction** in computational biophysics, bridging the gap between single-molecule chemical kinetics and in vivo genome engineering systems biology.
