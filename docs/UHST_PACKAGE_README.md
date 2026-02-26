# UHST v15 — Full Reproducibility Package (SPARC Rotmod_LTG)

**Author:** Krzysztof Tylec  
**Date:** 2026-02-22  

This repository contains a **self-contained, full** reproducibility package for the UHST v15 model evaluation on SPARC Rotmod late-type galaxies.

## What is included (FULL)
- `UHST_v15_Manuscript.pdf` — manuscript (PDF)
- `reproduce_v15.py` — standalone reproducer script (CLI)
- `v15_params.json` — global parameters (and fixed M/L values)
- `v15_manifest.txt` — list of **171** galaxy IDs included in the final run
- `Rotmod_LTG.zip` — SPARC rotation-curve tables (rotmod) used as input (public dataset)
- `out_v15/` — **precomputed outputs** generated from the above inputs:
  - `global_summary.json`
  - `per_gal_results.csv`
  - `point_residuals_all.csv`
  - `failed_list.csv`
  - `run_metadata.json`
  - `figures/*.png` (4 showcase galaxies)

## Key pass/fail checks (for this exact package)
- `n_galaxies_processed` = **171**
- `n_points_total` = **3373**
- `chi2_total` = **422606.658570**

> Note: if your manuscript currently reports different totals (e.g. different χ² or point count), that means **a different SPARC release / different filtering / different parameters** were used. In that case, update the manuscript or update the package so they match **exactly**.

## How to reproduce from scratch (recommended)
1. Unzip SPARC rotmod tables:
   ```bash
   unzip -q Rotmod_LTG.zip -d Rotmod_LTG
   ```
2. Install dependencies:
   ```bash
   python -m pip install --upgrade numpy matplotlib
   ```
3. Run (deterministic, via manifest):
   ```bash
   python reproduce_v15.py \
     --sparc_rotmod_dir Rotmod_LTG \
     --params_json v15_params.json \
     --manifest v15_manifest.txt \
     --out_dir out_v15_reproduced
   ```

Compare `out_v15_reproduced/global_summary.json` with `out_v15/global_summary.json`.

## Statistical convention
The script computes:
- `chi2 = Σ ((Vobs - Vmodel)/σ)^2` with **σᵢ := Verr (raw SPARC errors)**  
No per-galaxy tuning and no intrinsic-scatter inflation are applied.

## Data provenance
The input archive `Rotmod_LTG.zip` contains rotmod tables for the SPARC sample of late-type galaxies (175 total in the full dataset). Please cite the original SPARC paper/database when using the data.

