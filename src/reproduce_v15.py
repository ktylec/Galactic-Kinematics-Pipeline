#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""reproduce_v15.py

Reproducer "Open Science" dla modelu UHST v15 (Two-Phase Superfluid Vacuum)
na bazie plików rotmod SPARC (*.dat).
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

PARSEC_M = 3.0856775814913673e16
KPC_TO_M = 1.0e3 * PARSEC_M

@dataclass(frozen=True)
class UHSTParams:
    y0_mid: float = -0.6
    w_mid: float = 1.6
    s_mid: float = 0.3
    y0_out: float = -1.2
    w_out: float = 1.6
    s_out: float = 0.4
    a_crit: float = 1.2e-10
    upsilon_disk: float = 0.5
    upsilon_bulge: float = 0.7
    outer_frac_r_over_rmax: float = 0.6

def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')

def stable_logistic_gate(y: np.ndarray, y0: float, w: float, exponent_sign: int) -> np.ndarray:
    if w == 0:
        raise ValueError("w (szerokość bramki) nie może być 0")
    x = (y - y0) / w
    expo = exponent_sign * x
    expo = np.clip(expo, -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(expo))

def load_params_json(path: Optional[Path]) -> UHSTParams:
    base = UHSTParams()
    if path is None or not path.exists():
        return base
    with path.open('r', encoding='utf-8') as f:
        raw = json.load(f)
    p = raw.get('v15_params', raw)
    ups = raw.get('upsilon', {})
    outer = raw.get('outer_definition', {})

    def pick(key: str, default: float) -> float:
        val = p.get(key, default)
        try:
            return float(val)
        except Exception:
            return default

    return UHSTParams(
        y0_mid=pick('y0_mid', base.y0_mid), w_mid=pick('w_mid', base.w_mid), s_mid=pick('s_mid', base.s_mid),
        y0_out=pick('y0_out', base.y0_out), w_out=pick('w_out', base.w_out), s_out=pick('s_out', base.s_out),
        a_crit=pick('a_crit', base.a_crit),
        upsilon_disk=float(ups.get('disk', base.upsilon_disk)) if ups else base.upsilon_disk,
        upsilon_bulge=float(ups.get('bulge', base.upsilon_bulge)) if ups else base.upsilon_bulge,
        outer_frac_r_over_rmax=float(outer.get('outer_frac_r_over_rmax', base.outer_frac_r_over_rmax)) if outer else base.outer_frac_r_over_rmax,
    )

def parse_manifest(path: Optional[Path]) -> Optional[set]:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(f"manifest nie istnieje: {path}")
    names: List[str] = []
    if path.suffix.lower() == '.txt':
        for ln in path.read_text(encoding='utf-8').splitlines():
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            names.append(ln)
    else:
        with path.open('r', encoding='utf-8', newline='') as f:
            reader = csv.reader(f)
            rows = list(reader)
        if not rows:
            return set()
        header = [c.strip().lower() for c in rows[0]]
        if 'name' in header:
            idx = header.index('name')
            for r in rows[1:]:
                if len(r) > idx and r[idx].strip():
                    names.append(r[idx].strip())
        else:
            for r in rows[1:] if any(any(ch.isalpha() for ch in c) for c in rows[0]) else rows:
                if r and r[0].strip():
                    names.append(r[0].strip())
    return set(names)

def infer_galaxy_name_from_filename(path: Path) -> str:
    stem = path.stem
    for suf in ['_rotmod', '-rotmod', '.rotmod']:
        if stem.endswith(suf):
            stem = stem[:-len(suf)]
    return stem

def load_rotmod_table(path: Path) -> np.ndarray:
    with path.open('r', encoding='utf-8', errors='replace') as f:
        raw_lines = [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith('#')]
    if not raw_lines:
        raise ValueError("pusty plik lub same komentarze")
    first = raw_lines[0]
    if any(ch.isalpha() for ch in first):
        raw_lines = raw_lines[1:]
        if not raw_lines:
            raise ValueError("nagłówek bez danych")
    rows = []
    for ln in raw_lines:
        try:
            rows.append([float(p) for p in ln.split()])
        except ValueError as e:
            raise ValueError(f"nie-numeryczny token w linii: {ln[:120]}") from e
    arr = np.asarray(rows, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 5:
        raise ValueError(f"zbyt mało kolumn ({arr.shape[1]}). Oczekuję >= 5")
    return arr

def compute_baryonic_velocity_sq(vgas: np.ndarray, vdisk: np.ndarray, vbul: np.ndarray, upsilon_disk: float, upsilon_bulge: float) -> np.ndarray:
    return (np.abs(vgas) * vgas + upsilon_disk * np.abs(vdisk) * vdisk + upsilon_bulge * np.abs(vbul) * vbul)

def uhst_v15_predict(r_kpc: np.ndarray, vgas_kms: np.ndarray, vdisk_kms: np.ndarray, vbul_kms: np.ndarray, params: UHSTParams) -> Dict[str, np.ndarray]:
    r_m = r_kpc * KPC_TO_M
    vbar2 = compute_baryonic_velocity_sq(vgas_kms, vdisk_kms, vbul_kms, params.upsilon_disk, params.upsilon_bulge)
    a_bar = (vbar2 * 1.0e6) / r_m
    ratio = a_bar / params.a_crit
    y = np.full_like(ratio, np.nan, dtype=float)
    mask_pos = ratio > 0
    y[mask_pos] = np.log10(ratio[mask_pos])

    G_mid = stable_logistic_gate(y, params.y0_mid, params.w_mid, exponent_sign=-1)
    G_out = stable_logistic_gate(y, params.y0_out, params.w_out, exponent_sign=+1)

    a_vac_mid = params.s_mid * params.a_crit * np.power(ratio, 0.5) * G_mid
    a_vac_out = params.s_out * params.a_crit * np.power(ratio, 0.25) * G_out
    a_tot = a_bar + a_vac_mid + a_vac_out

    Vbar_kms = np.sqrt(np.maximum(a_bar * r_m, 0.0)) / 1000.0
    Vmid_kms = np.sqrt(np.maximum(a_vac_mid * r_m, 0.0)) / 1000.0
    Vout_kms = np.sqrt(np.maximum(a_vac_out * r_m, 0.0)) / 1000.0
    Vmodel_kms = np.sqrt(np.maximum(a_tot * r_m, 0.0)) / 1000.0

    return {
        'Vbar_kms': Vbar_kms, 'Vmid_kms': Vmid_kms, 'Vout_kms': Vout_kms, 'Vmodel_kms': Vmodel_kms,
        'a_bar': a_bar, 'a_vac_mid': a_vac_mid, 'a_vac_out': a_vac_out, 'a_tot': a_tot,
        'y': y, 'G_mid': G_mid, 'G_out': G_out,
    }

@dataclass
class GalaxyResult:
    name: str
    file: str
    n_points: int
    chi2: float
    chi2_over_n: float
    Vmax_kms: float
    Rmax_kpc: float

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def write_csv(path: Path, header: List[str], rows: Iterable[Iterable]) -> None:
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(list(r))

def plot_showcase_rotation_curve(out_png: Path, galaxy: str, r_kpc: np.ndarray, vobs: np.ndarray, verr: np.ndarray, vmodel: np.ndarray, vbar: np.ndarray, vmid: np.ndarray, vout: np.ndarray, chi2: float, params: UHSTParams) -> None:
    if plt is None:
        return
    fig = plt.figure(figsize=(7.5, 4.8))
    ax = fig.gca()
    ax.errorbar(r_kpc, vobs, yerr=verr, fmt='o', markersize=3, capsize=2, label='Vobs (SPARC)')
    ax.plot(r_kpc, vmodel, label='Vmodel (UHST v15)')
    ax.plot(r_kpc, vbar, label='Vbar (baryony)')
    ax.plot(r_kpc, vmid, label='Vvac,mid')
    ax.plot(r_kpc, vout, label='Vvac,out')
    ax.set_xlabel('R [kpc]')
    ax.set_ylabel('V [km/s]')
    ax.set_title(f"{galaxy} | n={len(r_kpc)} chi2={chi2:.1f} | y0_mid={params.y0_mid}, w_mid={params.w_mid}, s_mid={params.s_mid} y0_out={params.y0_out}, w_out={params.w_out}, s_out={params.s_out}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog='reproduce_v15.py')
    ap.add_argument('--sparc_rotmod_dir', required=True, type=Path)
    ap.add_argument('--out_dir', default=Path('out_v15'), type=Path)
    ap.add_argument('--params_json', default=None, type=Path)
    ap.add_argument('--manifest', default=None, type=Path)
    ap.add_argument('--no_plots', action='store_true')
    ap.add_argument('--showcase', nargs='*', default=['F579-V1', 'UGC00128', 'UGC02487', 'UGC07577'])
    ap.add_argument('--k_params', type=int, default=6)
    ap.add_argument('--expected_chi2', type=float, default=416493.0)
    ap.add_argument('--tolerance', type=float, default=50.0)
    args = ap.parse_args(argv)

    if not args.sparc_rotmod_dir.exists():
        print(f"[ERROR] katalog nie istnieje: {args.sparc_rotmod_dir}", file=sys.stderr)
        return 2

    params = load_params_json(args.params_json)
    manifest = parse_manifest(args.manifest)
    ensure_dir(args.out_dir)
    ensure_dir(args.out_dir / 'figures')

    meta = {'timestamp_utc': _utcnow_iso(), 'python': sys.version.split()[0], 'sparc_rotmod_dir': str(args.sparc_rotmod_dir), 'params': params.__dict__}
    (args.out_dir / 'run_metadata.json').write_text(json.dumps(meta, indent=2, sort_keys=True), encoding='utf-8')

    files = sorted([p for p in args.sparc_rotmod_dir.glob('*.dat')])
    if not files:
        print(f"[ERROR] brak plików *.dat w {args.sparc_rotmod_dir}", file=sys.stderr)
        return 2

    per_gal, failed, point_rows = [], [], []
    showcase_cache, showcase_obs = {}, {}
    total_chi2, total_n = 0.0, 0

    for f in files:
        name = infer_galaxy_name_from_filename(f)
        if manifest is not None and name not in manifest:
            continue
        try:
            arr = load_rotmod_table(f)
            r, vobs, verr, vgas, vdisk = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]
            vbul = arr[:, 5] if arr.shape[1] >= 6 else np.zeros_like(r)

            m = np.isfinite(r) & np.isfinite(vobs) & np.isfinite(verr) & np.isfinite(vgas) & np.isfinite(vdisk) & np.isfinite(vbul) & (r > 0) & (verr > 0)
            r, vobs, verr, vgas, vdisk, vbul = r[m], vobs[m], verr[m], vgas[m], vdisk[m], vbul[m]

            idx = np.argsort(r)
            r, vobs, verr, vgas, vdisk, vbul = r[idx], vobs[idx], verr[idx], vgas[idx], vdisk[idx], vbul[idx]

            pred = uhst_v15_predict(r, vgas, vdisk, vbul, params)
            good = np.isfinite(pred['y']) & np.isfinite(pred['Vmodel_kms']) & np.isfinite(vobs) & np.isfinite(verr) & (verr > 0)
            
            r_g, vobs_g, verr_g = r[good], vobs[good], verr[good]
            if r_g.size == 0:
                raise ValueError('wszystkie punkty odrzucone')

            vmodel_g = pred['Vmodel_kms'][good]
            resid = vobs_g - vmodel_g
            z = resid / verr_g
            chi2, n = float(np.sum(z ** 2)), int(r_g.size)
            total_chi2 += chi2
            total_n += n

            per_gal.append(GalaxyResult(name, str(f.name), n, chi2, chi2 / n, float(np.nanmax(vobs_g)), float(np.nanmax(r_g))))

            for i in range(n):
                point_rows.append([name, float(r_g[i]), float(vobs_g[i]), float(verr_g[i]), float(vmodel_g[i]), float(resid[i]), float(z[i]), float(pred['y'][good][i]), float(pred['a_bar'][good][i]), float(pred['a_vac_mid'][good][i]), float(pred['a_vac_out'][good][i])])

            if (not args.no_plots) and (name in args.showcase):
                showcase_cache[name] = {'r_kpc': r_g, 'Vmodel_kms': vmodel_g, 'Vbar_kms': pred['Vbar_kms'][good], 'Vmid_kms': pred['Vmid_kms'][good], 'Vout_kms': pred['Vout_kms'][good]}
                showcase_obs[name] = {'Vobs_kms': vobs_g, 'Verr_kms': verr_g, 'chi2': np.array([chi2], dtype=float)}
        except Exception as e:
            failed.append((str(f.name), str(e)))

    per_gal.sort(key=lambda g: g.chi2, reverse=True)
    write_csv(args.out_dir / 'failed_list.csv', ['file', 'reason'], failed)
    write_csv(args.out_dir / 'per_gal_results.csv', ['name', 'file', 'n_points', 'chi2', 'chi2_over_n', 'Vmax_kms', 'Rmax_kpc'], [[g.name, g.file, g.n_points, f"{g.chi2:.6f}", f"{g.chi2_over_n:.6f}", f"{g.Vmax_kms:.6f}", f"{g.Rmax_kpc:.6f}"] for g in per_gal])
    write_csv(args.out_dir / 'point_residuals_all.csv', ['galaxy', 'R_kpc', 'Vobs_kms', 'Verr_kms', 'Vmodel_kms', 'resid_kms', 'z', 'y', 'a_bar_mps2', 'a_vac_mid_mps2', 'a_vac_out_mps2'], point_rows)

    k, n, chi2 = int(args.k_params), int(total_n), float(total_chi2)
    bic_known = chi2 + k * math.log(max(n, 1))
    summary = {'n_galaxies_processed': len(per_gal), 'n_points_total': n, 'chi2_total': chi2, 'k_params': k, 'BIC_known_sigma': bic_known, 'expected_chi2': float(args.expected_chi2), 'abs_diff_chi2': abs(chi2 - float(args.expected_chi2)), 'tolerance': float(args.tolerance)}
    (args.out_dir / 'global_summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True), encoding='utf-8')

    print("=" * 78)
    print("UHST v15 — reprodukcja")
    print(f"Galaktyk przetworz.:   {len(per_gal)}")
    print(f"Punktów łącznie:       {n}")
    print(f"chi^2_total:           {chi2:.6f}")
    print(f"|diff|:                {abs(chi2 - args.expected_chi2):.6f}")
    print("=" * 78)

    if not args.no_plots and plt is not None:
        for gname in args.showcase:
            if gname in showcase_cache:
                plot_showcase_rotation_curve(args.out_dir / 'figures' / f"{gname}_Vcurves.png", gname, showcase_cache[gname]['r_kpc'], showcase_obs[gname]['Vobs_kms'], showcase_obs[gname]['Verr_kms'], showcase_cache[gname]['Vmodel_kms'], showcase_cache[gname]['Vbar_kms'], showcase_cache[gname]['Vmid_kms'], showcase_cache[gname]['Vout_kms'], float(showcase_obs[gname]['chi2'][0]), params)

    return 0

if __name__ == '__main__':
    raise SystemExit(main())
