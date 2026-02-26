# Galactic Kinematics Pipeline & Statistical Evaluation (UHST v15)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18749640.svg)](https://doi.org/10.5281/zenodo.18749640)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)

## 📌 Executive Summary
This repository contains the data engineering pipeline, statistical modeling, and reproducibility framework for the **Unified Holographic Superfluid Theory (UHST v15)**. 
It demonstrates the application of rigorous Data Science techniques to astrophysical open data (SPARC database) to solve the galactic missing mass problem without relying on dark matter halo tuning.

## 🚀 Quick Start (Local Reproducibility)
To reproduce the entire validation pipeline locally using the provided manifest:

`python -m venv .venv`
`source .venv/bin/activate`          *(Windows: .venv\Scripts\activate)*
`pip install -r requirements.txt`

`python src/validate_data.py --data_dir data/sparc_dat`
`python src/run_reproduce.py`


## 📂 Repository Structure & Contents
This repository is a clean, GitHub-oriented mirror of the UHST v15.1 Full Package, reorganized into a standard Data Engineering structure:

* `data/sparc_dat/` — Raw SPARC kinematic datasets (`*_rotmod.dat`).
* `data/params/` — Global parameter configurations (`v15_params.json`), manifests, and `SHA256SUMS.txt`.
* `src/` — Python execution scripts, including the verbatim upstream `reproduce_v15.py` and data validators.
* `results/out_v15/` — Precomputed output artifacts (CSV summaries, generated PNG plots).
* `paper/` — The comprehensive manuscript (`UHST_v15_Manuscript.pdf`) and LaTeX placeholder.

## 🔬 Scientific Context & Validation
Standard astrophysical models require 2-3 free parameters *per galaxy*, leading to massive overfitting. This project implements the UHST framework using only **6 global parameters** shared across the entire sample of 171 galaxies (3,373 data points).

**Key Results:**
* **Elimination of Tilt:** Residual systematic radial bias is bounded near zero.
* **Out-of-Sample Stability:** Verified via rigorous cross-validation methodologies.
* **Information Criterion:** $\Delta$BIC $\approx$ 2800 improvement over standard tuned halo models.

## 🔗 Upstream Package Provenance
For full data provenance and exact checksums verifying the integrity of this scientific package, please see:
- `docs/UHST_PACKAGE_README.md` (if applicable)
- `data/params/SHA256SUMS.txt`
- `data/params/github_relayout_manifest.json`

## 🤖 AI Transparency & Methodology
I believe in full transparency regarding the tools used in modern research. 
* **Conceptual Framework:** The physics logic, mathematical derivations, validation strategies (e.g., K-Fold / Monte Carlo CV), and data interpretation were developed entirely independently.
* **Code Implementation:** I am not a traditional software developer. The Python boilerplate, `pandas` data pipelines, and `scipy` optimization scripts were written with the assistance of advanced Large Language Models acting as my coding assistants. All AI-generated code was strictly directed and mathematically verified against the raw SPARC database to ensure absolute deterministic accuracy.

---
**Author / Maintainer:** Krzysztof Tylec (Independent Researcher & Data Analyst)