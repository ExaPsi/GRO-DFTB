#!/usr/bin/env python3
"""Analyze B6 multi-seed NVE drift data.

Computes drift rate for each seed and reports mean, std, 95% CI.
Also computes B6/B5 drift ratio with error propagation.

Usage:
    python3 tools/analyze_b6_multiseed.py docs/session/20260213/t2_b6_multiseed/
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import numpy as np
from datetime import date

GMX = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "external", "gromacs", "_install_dftb", "bin", "gmx"
)

SEEDS = [42, 137, 2718, 31415, 98765]
NATOM = 3163  # B6 total atoms
N_QM = 12     # 10 real + 2 link H
THRESHOLD = 0.005  # kJ/mol/ps/atom

# B5 PME 10-seed reference
B5_DRIFT_MEAN = 6.1e-4
B5_DRIFT_STD = 1.0e-4


def extract_total_energy(edr_path):
    """Extract total energy time series from .edr file.
    Uses 'Total-Energy' term (15 for B6 with bonded terms)."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xvg', delete=False) as f:
        xvg_path = f.name
    try:
        # Use term name selection (Total-Energy) — works regardless of numbering
        proc = subprocess.run(
            [GMX, "energy", "-f", edr_path, "-o", xvg_path],
            input="Total-Energy\n\n", capture_output=True, text=True, timeout=60
        )
        data = []
        with open(xvg_path) as f:
            for line in f:
                if not line.startswith(('#', '@')):
                    parts = line.split()
                    if len(parts) >= 2:
                        data.append((float(parts[0]), float(parts[1])))
        return np.array(data)
    finally:
        os.unlink(xvg_path)


def compute_drift(t, E, natom):
    """Compute drift rate from linear regression."""
    coeffs = np.polyfit(t, E, 1)
    drift = coeffs[0] / natom
    E_fit = np.polyval(coeffs, t)
    r2 = 1 - np.sum((E - E_fit)**2) / np.sum((E - np.mean(E))**2)

    n2 = len(E) // 2
    d1 = np.polyfit(t[:n2], E[:n2], 1)[0] / natom
    d2 = np.polyfit(t[n2:], E[n2:], 1)[0] / natom
    block_ratio = max(abs(d1), abs(d2)) / min(abs(d1), abs(d2)) \
        if min(abs(d1), abs(d2)) > 0 else float('inf')

    return {
        'drift': float(drift),
        'r_squared': float(r2),
        'block_ratio': float(block_ratio),
        'e_mean': float(np.mean(E)),
        'e_std': float(np.std(E)),
        'duration_ps': float(t[-1]),
        'n_frames': len(t),
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze B6 multi-seed NVE')
    parser.add_argument('datadir', help='Directory with seed_* subdirectories')
    parser.add_argument('-o', '--output', default=None, help='Output JSON')
    args = parser.parse_args()

    seed_results = {}
    drifts = []

    for seed in SEEDS:
        sdir = os.path.join(args.datadir, f'seed_{seed}')
        edr = os.path.join(sdir, 'nve.edr')

        if not os.path.exists(edr):
            print(f"  Seed {seed}: MISSING (no edr file)")
            continue

        try:
            data = extract_total_energy(edr)
            t, E = data[:, 0], data[:, 1]
            result = compute_drift(t, E, NATOM)
            result['seed'] = seed
            result['pass'] = abs(result['drift']) < THRESHOLD
            result['margin'] = THRESHOLD / abs(result['drift'])
            seed_results[str(seed)] = result
            drifts.append(result['drift'])
            print(f"  Seed {seed}: drift = {result['drift']:.4e}, "
                  f"R² = {result['r_squared']:.4f}, "
                  f"block = {result['block_ratio']:.2f}, "
                  f"margin = {result['margin']:.1f}x, "
                  f"{'PASS' if result['pass'] else 'FAIL'}")
        except Exception as e:
            print(f"  Seed {seed}: ERROR - {e}")

    if len(drifts) < 2:
        print(f"\nOnly {len(drifts)} seeds completed. Need at least 2 for statistics.")
        return

    drifts = np.array(drifts)
    n = len(drifts)
    mean = float(np.mean(drifts))
    std = float(np.std(drifts, ddof=1))
    stderr = std / np.sqrt(n)

    # 95% CI using t-distribution
    from scipy.stats import t as t_dist
    t_crit = t_dist.ppf(0.975, df=n-1)
    ci_low = mean - t_crit * stderr
    ci_high = mean + t_crit * stderr

    # B6/B5 ratio with error propagation
    ratio = abs(mean) / B5_DRIFT_MEAN
    # σ_ratio/ratio = sqrt((σ_B6/μ_B6)² + (σ_B5/μ_B5)²)
    ratio_err = ratio * np.sqrt((std/abs(mean))**2 + (B5_DRIFT_STD/B5_DRIFT_MEAN)**2)

    print(f"\n{'='*60}")
    print(f"B6 Multi-Seed NVE Drift Summary ({n} seeds)")
    print(f"{'='*60}")
    print(f"  Mean drift: {mean:.4e} kJ/mol/ps/atom")
    print(f"  Std:        {std:.4e}")
    print(f"  Stderr:     {stderr:.4e}")
    print(f"  95% CI:     [{ci_low:.4e}, {ci_high:.4e}]")
    print(f"  Mean margin: {THRESHOLD/abs(mean):.1f}x below {THRESHOLD}")
    print(f"  All pass:    {all(abs(d) < THRESHOLD for d in drifts)}")
    print(f"  B6/B5 ratio: {ratio:.2f} ± {ratio_err:.2f}")

    summary = {
        'description': 'B6 solvated alanine dipeptide, 5-seed NVE drift analysis',
        'system': 'B6 (12 QM atoms, 3163 total, PME, sigma=0.1 nm)',
        'n_seeds': n,
        'n_total_atoms': NATOM,
        'n_qm_atoms': N_QM,
        'threshold': THRESHOLD,
        'mean_drift': mean,
        'std_drift': std,
        'stderr_drift': stderr,
        'ci_95_low': ci_low,
        'ci_95_high': ci_high,
        'mean_margin': float(THRESHOLD / abs(mean)),
        'all_pass': bool(all(abs(d) < THRESHOLD for d in drifts)),
        'b6_b5_ratio': float(ratio),
        'b6_b5_ratio_err': float(ratio_err),
        'b5_reference': {
            'mean_drift': B5_DRIFT_MEAN,
            'std_drift': B5_DRIFT_STD,
            'source': '10-seed PME NVE (docs/session/20260212/a2_statistical/)'
        },
        'seeds': seed_results,
        'provenance': {
            'date': str(date.today()),
            'script': 'tools/analyze_b6_multiseed.py',
            'gmx': GMX,
        }
    }

    outfile = args.output or os.path.join(args.datadir, 'multiseed_results.json')
    with open(outfile, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to: {outfile}")


if __name__ == '__main__':
    main()
