# 🌌 Space Radiobiology Pipeline
### *In Silico Structural Vulnerability Analysis of Human Stress-Response Proteins Under Deep Space Conditions*

[![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=flat&logo=python&logoColor=white)](https://www.python.org/)
[![AlphaFold](https://img.shields.io/badge/AlphaFold-EBI%20API-005B94?style=flat)](https://alphafold.ebi.ac.uk/)
[![UniProt](https://img.shields.io/badge/UniProt-REST%20API-003366?style=flat)](https://www.uniprot.org/)
[![WebGL](https://img.shields.io/badge/3Dmol.js-Interactive%20WebGL-F34F29?style=flat&logo=webgl&logoColor=white)](https://3dmol.csb.pitt.edu/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 📌 Abstract & Overview

During deep-space exploration, biological systems face a harsh combination of cosmic ionizing radiation and microgravity. While cellular signaling networks rapidly **upregulate emergency defense and repair genes**, the translated proteins often possess structural motifs highly susceptible to **radiolytic reactive oxygen species (ROS) and free-radical attack**.

This repository provides an automated, end-to-end computational pipeline that:
1. **Performs differential gene expression analysis** on microarray spaceflight datasets (1G Earth Control vs. Deep Space environment).
2. **Quantifies biophysical vulnerability** by profiling radiolytically sensitive residues (*Cysteine, Methionine, and Aromatic side chains*).
3. **Computes an empirical Paradox Score** to isolate highly expressed yet structurally fragile targets.
4. **Fetches 3D atomic structures dynamically** via UniProt and AlphaFold REST APIs.
5. **Renders an interactive 3D WebGL multi-panel dashboard (`model.html`)** highlighting vulnerable atomic sites in real time.

---

## 🖥️ 3D Dashboard Architecture (3x2 Grid)

* **Top Row:** `CCL2` | `TMEM176A` | `HDC`
* **Bottom Row:** `PPA1` | `ADAMDEC1` | `UTS2`

Each panel provides independent 3D camera controls (orbit, zoom, pan), real-time atom rendering, and contextual stress-response annotations.

---

## 🧮 Mathematical Model & Scoring Pipeline

### 1. Differential Expression ($\Delta_{\text{Space}}$)
Given the baseline control group $C$ (Earth Normal) and the deep-space group $D$ (Microgravity + Radiation):
$$\Delta_{\text{Space}} = \bar{X}_{\text{Deep Space}} - \bar{X}_{\text{Earth Normal}}$$
Candidates are filtered using a two-sample Welch's $t$-test ($p < 0.05$) and a fold-change threshold ($\Delta_{\text{Space}} > 1.0$).

### 2. Radical Vulnerability Score ($\text{Weakness}$)
Water radiolysis generates reactive hydroxyl radicals ($^\bullet\text{OH}$), primarily targeting sulfur-containing and aromatic residues:
* **Cysteine (C) [Weight: 2.0]:** Disruption of structural disulfide bonds ($-\text{S}-\text{S}-$).
* **Methionine (M) [Weight: 1.0]:** High susceptibility to thioether oxidation.
* **Tryptophan, Tyrosine, Histidine (W, Y, H) [Weight: 0.5]:** Oxidation-prone aromatic rings.

$$\text{Weakness (\%)} = \frac{2 \cdot N_{\text{C}} + 1 \cdot N_{\text{M}} + 0.5 \cdot (N_{\text{W}} + N_{\text{Y}} + N_{\text{H}})}{\text{Total Sequence Length } (L)} \times 100$$

### 3. Paradox Score
$$\text{Paradox} = \text{Weakness} \times 2^{\Delta_{\text{Space}}}$$

---

## 🧬 Top Identified Targets (Deep Space)

| Gene Symbol | Full Gene Name | $\Delta_{\text{Space}}$ | Weakness (%) | Paradox Score | Primary Radiolytic Vulnerability |
|:---|:---|:---:|:---:|:---:|:---|
| **CCL2** | C-C motif chemokine ligand 2 | 2.31 | 14.14 | **0.70** | Disulfide bridge cleavage prevents immune cell recruitment |
| **TMEM176A** | Transmembrane protein 176A | 1.33 | 12.34 | **0.31** | Membrane core degradation disrupts ion homeostasis |
| **HDC** | Histidine decarboxylase | 1.21 | 12.08 | **0.28** | Catalytic active-site breakdown halts histamine defense |
| **PPA1** | Inorganic pyrophosphatase 1 | 1.27 | 11.59 | **0.28** | Oxidative inactivation starves DNA repair of energy |
| **ADAMDEC1** | ADAM like decysin 1 | 1.04 | 13.40 | **0.28** | High cysteine density leads to scaffold misfolding |
| **UTS2** | Urotensin 2 | 1.06 | 10.89 | **0.23** | Cyclic ring cleavage disables vasoconstrictive signaling |

---

## 🎨 3D Visualization Palette

The pipeline maps chemical vulnerabilities onto AlphaFold models using color-coded styles:

* 🔲 **Grey Cartoon:** Native secondary structure backbone.
* 🔴 **Red Sticks & Spheres:** **Cysteine (CYS)** — fragile disulfide bridges.
* 🟠 **Orange Sticks & Spheres:** **Methionine (MET)** — oxidation-sensitive thioethers.
* 🟡 **Yellow Sticks:** **Aromatic residues (TRP, TYR, HIS)** — oxidation-prone ring networks.

---

## 🚀 Quick Start

### 1. Installation
Clone the repository and install dependencies:
```bash
git clone [https://github.com/](https://github.com/)<your-username>/<your-repo-name>.git
cd <your-repo-name>
pip install -r requirements.txt
```

### 2. Data Preparation
Ensure your microarray expression dataset is placed as `cumulative.csv` in the root directory.

### 3. Execution
Run the automated pipeline:
```bash
python main.py
```
The script will process the expression matrix, query UniProt & AlphaFold APIs, generate `model.html`, and automatically launch the interactive 3D viewer in your default browser.

---

## 📦 Dependencies

* `numpy>=1.24.0`
* `pandas>=2.0.0`
* `scipy>=1.10.0`
* `3Dmol.js` *(embedded via CDN in generated HTML)*

---

## 🙏 Acknowledgments & Data Sources

* **Microarray Expression Data:** Publicly available spaceflight transcriptomic datasets from NCBI GEO / NASA GeneLab.
* **Structural Coordinates & Sequences:** [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) (EMBL-EBI) & [UniProt Consortium](https://www.uniprot.org/).
* **3D Molecular Rendering:** Powered by [3Dmol.js](https://3dmol.csb.pitt.edu/).
* **Development Assistance:** AI-assisted mathematical structuring and documentation refinement.

---

## 📄 License

Distributed under the **MIT License**. See `LICENSE` for more information.
