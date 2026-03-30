#!/usr/bin/env python3
"""
A-9: Energy decomposition analysis for PME QM/MM trajectory.

Parses GROMACS energy output (.xvg) and md.log embedding diagnostics to produce
publication-quality data files for pgfplots energy decomposition figure.

Usage:
    python tools/analyze_energy_decomposition.py \
        --xvg docs/session/20260213/a9_decomposition/energy_components.xvg \
        --log tests/data/b5_pme_drift/md.log \
        --outdir docs/session/20260213/a9_decomposition \
        --tmax 10.0

Outputs:
    total_energy.dat       - Total energy time series (ΔE from mean)
    components.dat         - Component fluctuations (ΔE from mean)
    embedding_diag.dat     - E_embedding_check from md.log (kJ/mol)
    decomposition_stats.json - Summary statistics and verification
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np


def parse_xvg(xvg_path):
    """Parse GROMACS .xvg file, return time array and column data dict."""
    times = []
    data_rows = []
    with open(xvg_path) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 10:
                times.append(float(parts[0]))
                data_rows.append([float(x) for x in parts[1:]])
    times = np.array(times)
    data = np.array(data_rows)
    # Column mapping (0-indexed within data array):
    # 0: Connect Bonds, 1: LJ (SR), 2: Disper. corr., 3: Coulomb (SR),
    # 4: Coul. recip., 5: Quantum En., 6: Potential, 7: Kinetic En., 8: Total Energy
    columns = {
        'bonds': data[:, 0],
        'lj_sr': data[:, 1],
        'disp_corr': data[:, 2],
        'coul_sr': data[:, 3],
        'coul_recip': data[:, 4],
        'quantum_en': data[:, 5],
        'potential': data[:, 6],
        'kinetic_en': data[:, 7],
        'total_energy': data[:, 8],
    }
    return times, columns


def parse_md_log_embedding(log_path, tmax_step=None, dt=0.0005):
    """Parse E_embedding_check and E_correction from md.log."""
    pattern = re.compile(
        r'Step\s+(\d+):\s+E_embedding_check=([-\d.eE+]+)\s+Ha,\s+'
        r'E_correction=([-\d.eE+]+)\s+Ha\s+\(([-\d.eE+]+)\s+kJ/mol\)'
    )
    steps, times, e_emb, e_corr, e_corr_kjmol = [], [], [], [], []
    # Hartree to kJ/mol
    ha2kjmol = 2625.499639  # CODATA 2022
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                step = int(m.group(1))
                t = step * dt
                if tmax_step is not None and step > tmax_step:
                    break
                steps.append(step)
                times.append(t)
                e_emb.append(float(m.group(2)) * ha2kjmol)  # Convert to kJ/mol
                e_corr.append(float(m.group(3)) * ha2kjmol)
                e_corr_kjmol.append(float(m.group(4)))
    return {
        'steps': np.array(steps),
        'times': np.array(times),
        'e_embedding_check': np.array(e_emb),
        'e_correction': np.array(e_corr),
    }


def main():
    parser = argparse.ArgumentParser(description='A-9 Energy decomposition analysis')
    parser.add_argument('--xvg', required=True, help='GROMACS energy .xvg file')
    parser.add_argument('--log', required=True, help='GROMACS md.log file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--tmax', type=float, default=10.0, help='Max time (ps)')
    parser.add_argument('--dt', type=float, default=0.0005, help='Timestep (ps)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse energy data
    times, cols = parse_xvg(args.xvg)
    mask = times <= args.tmax + 1e-6
    times = times[mask]
    for k in cols:
        cols[k] = cols[k][mask]
    nframes = len(times)

    # Compute derived quantities
    coul_total = cols['coul_sr'] + cols['coul_recip']
    classical_nonelec = cols['bonds'] + cols['lj_sr'] + cols['disp_corr']

    # Verify energy sum
    e_sum = (cols['bonds'] + cols['lj_sr'] + cols['disp_corr'] +
             cols['coul_sr'] + cols['coul_recip'] + cols['quantum_en'] +
             cols['kinetic_en'])
    max_sum_err = np.max(np.abs(e_sum - cols['total_energy']))
    print(f"Max energy sum error: {max_sum_err:.6e} kJ/mol (should be ~0)")

    # Component means and fluctuations
    means = {k: np.mean(v) for k, v in cols.items()}
    means['coul_total'] = np.mean(coul_total)
    means['classical_nonelec'] = np.mean(classical_nonelec)

    # Statistics
    stats = {
        'tmax_ps': args.tmax,
        'nframes': nframes,
        'dt_ps': times[1] - times[0] if nframes > 1 else args.dt,
        'max_energy_sum_error_kjmol': float(max_sum_err),
        'means_kjmol': {k: float(v) for k, v in means.items()},
        'std_kjmol': {k: float(np.std(cols[k])) for k in cols},
    }
    stats['std_kjmol']['coul_total'] = float(np.std(coul_total))

    # --- Write data files for pgfplots ---

    # 1. Total energy fluctuation (subsample for manageable file size)
    stride = max(1, nframes // 2000)
    idx = np.arange(0, nframes, stride)
    with open(outdir / 'total_energy.dat', 'w') as f:
        f.write('# A-9 Energy Decomposition: Total energy fluctuation\n')
        f.write(f'# Source: {args.xvg}\n')
        f.write(f'# Mean Total Energy: {means["total_energy"]:.4f} kJ/mol\n')
        f.write('# time_ps  delta_E_total_kjmol\n')
        for i in idx:
            de = cols['total_energy'][i] - means['total_energy']
            f.write(f'{times[i]:.6f}  {de:.6f}\n')

    # 2. Component fluctuations
    components = [
        ('quantum_en', 'E_QM'),
        ('kinetic_en', 'E_kin'),
        ('coul_sr', 'E_Coul_SR'),
        ('coul_recip', 'E_Coul_recip'),
        ('lj_sr', 'E_LJ'),
    ]
    with open(outdir / 'components.dat', 'w') as f:
        f.write('# A-9 Energy Decomposition: Component fluctuations around mean\n')
        f.write(f'# Source: {args.xvg}\n')
        hdr = 'time_ps  ' + '  '.join(f'delta_{name}' for _, name in components)
        f.write(f'# {hdr}\n')
        for i in idx:
            vals = [cols[key][i] - means[key] for key, _ in components]
            line = f'{times[i]:.6f}  ' + '  '.join(f'{v:.4f}' for v in vals)
            f.write(line + '\n')

    # 3. Embedding diagnostic from md.log
    tmax_step = int(args.tmax / args.dt)
    emb_data = parse_md_log_embedding(args.log, tmax_step=tmax_step, dt=args.dt)
    with open(outdir / 'embedding_diag.dat', 'w') as f:
        f.write('# A-9 Energy Decomposition: QM-MM embedding energy diagnostic\n')
        f.write(f'# Source: {args.log}\n')
        f.write('# time_ps  E_embedding_check_kjmol  E_correction_kjmol\n')
        for i in range(len(emb_data['times'])):
            f.write(f'{emb_data["times"][i]:.6f}  '
                    f'{emb_data["e_embedding_check"][i]:.4f}  '
                    f'{emb_data["e_correction"][i]:.4f}\n')

    # 4. Correlation matrix (for discussion)
    comp_arrays = np.column_stack([
        cols['quantum_en'][idx], cols['kinetic_en'][idx],
        cols['coul_sr'][idx], cols['coul_recip'][idx], cols['lj_sr'][idx],
    ])
    corr = np.corrcoef(comp_arrays.T)
    comp_names = ['E_QM', 'E_kin', 'E_Coul_SR', 'E_Coul_recip', 'E_LJ']
    stats['correlation_matrix'] = {
        'labels': comp_names,
        'values': corr.tolist()
    }

    # Add embedding stats
    if len(emb_data['e_embedding_check']) > 0:
        stats['embedding_check'] = {
            'npoints': len(emb_data['e_embedding_check']),
            'mean_kjmol': float(np.mean(emb_data['e_embedding_check'])),
            'std_kjmol': float(np.std(emb_data['e_embedding_check'])),
            'min_kjmol': float(np.min(emb_data['e_embedding_check'])),
            'max_kjmol': float(np.max(emb_data['e_embedding_check'])),
            'e_correction_all_zero': bool(np.all(emb_data['e_correction'] == 0)),
        }

    # Print summary
    print(f"\n=== A-9 Energy Decomposition Summary ({args.tmax} ps) ===")
    print(f"Frames: {nframes} (stride {stride} → {len(idx)} points in output)")
    print(f"\nComponent means (kJ/mol):")
    for key in ['total_energy', 'quantum_en', 'coul_sr', 'coul_recip',
                'lj_sr', 'disp_corr', 'kinetic_en']:
        print(f"  {key:20s}: {means[key]:12.2f} ± {np.std(cols[key]):8.2f}")
    print(f"\nVerification: sum of components = total?")
    print(f"  Max |sum - total| = {max_sum_err:.2e} kJ/mol")

    if len(emb_data['e_embedding_check']) > 0:
        print(f"\nEmbedding diagnostic ({len(emb_data['e_embedding_check'])} points):")
        print(f"  E_embedding_check: {stats['embedding_check']['mean_kjmol']:.2f} "
              f"± {stats['embedding_check']['std_kjmol']:.2f} kJ/mol")
        print(f"  E_correction all zero: {stats['embedding_check']['e_correction_all_zero']}")

    # Key correlations
    print(f"\nKey correlations:")
    for i, ni in enumerate(comp_names):
        for j, nj in enumerate(comp_names):
            if j > i and abs(corr[i, j]) > 0.3:
                print(f"  {ni} vs {nj}: r = {corr[i, j]:.3f}")

    # Save stats
    with open(outdir / 'decomposition_stats.json', 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"\nOutput files in {outdir}/")
    print(f"  total_energy.dat, components.dat, embedding_diag.dat, decomposition_stats.json")


if __name__ == '__main__':
    main()
