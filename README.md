# MANHATTAN — UHST v15.1 (FULL PACKAGE) — GitHub-ready mirror

This repository is a **clean GitHub-oriented mirror** of the **UHST v15.1 FULL PACKAGE**, reorganized into a standard structure:

- `data/` → raw SPARC inputs + parameter files and provenance
- `src/` → Python scripts (including the upstream `reproduce_v15.py`)
- `results/` → output artifacts (`out_v15/`, CSV/PNG, etc.)
- `paper/` → manuscript PDF + a LaTeX placeholder directory (optional)

## Quickstart

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python src/validate_data.py --data_dir data/sparc_dat
python src/run_reproduce.py
```

## Contents

- `data/sparc_dat/` — `*_rotmod.dat` files (SPARC Rotmod_LTG)
- `data/params/` — `v15_params.json`, `v15_manifest.txt`, `SHA256SUMS.txt` + `github_relayout_manifest.json`
- `src/reproduce_v15.py` — upstream reproducibility script (kept **verbatim**)
- `results/out_v15/` — upstream output bundle included in the full package
- `paper/UHST_v15_Manuscript.pdf` — manuscript PDF from the full package

## Upstream package provenance

See:
- `docs/UHST_PACKAGE_README.md`
- `data/params/SHA256SUMS.txt`
- `data/params/github_relayout_manifest.json`

Maintainer: Krzysztof Tylec
