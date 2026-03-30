#!/usr/bin/env python3
"""Analyze A-2 statistical NVE drift data from multiple independent seeds.

Computes drift rate (via linear regression of total energy vs time) for each
independent trajectory, then reports mean ± std across seeds for cutoff and
PME embedding modes.

Usage:
    python3 tools/analyze_a2_statistical.py docs/session/20260212/a2_statistical/
    python3 tools/analyze_a2_statistical.py docs/session/20260212/a2_statistical/ --gmx /path/to/gmx
    python3 tools/analyze_a2_statistical.py docs/session/20260212/a2_statistical/ -o results.json
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from datetime import date

VERSION = "0.1.0"
GMX_DEFAULT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "external", "gromacs", "_install_dftb", "bin", "gmx"
)
DRIFT_THRESHOLD = 0.01  # kJ/mol/ps/atom (specs.md section 21.1)


def extract_energy_xvg(gmx, edr_path, terms="Total-Energy"):
    """Extract energy terms from .edr using gmx energy, return (time[], val[])."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xvg', delete=False) as f:
        xvg_path = f.name
    try:
        # Map term names to selection numbers
        proc = subprocess.run(
            [gmx, "energy", "-f", edr_path],
            input="0\n", capture_output=True, text=True, timeout=30
        )
        lines = proc.stdout.splitlines() + proc.stderr.splitlines()
        term_num = None
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == terms:
                term_num = parts[0]
                break
        if term_num is None:
            raise RuntimeError(f"Could not find term '{terms}' in {edr_path}")

        proc = subprocess.run(
            [gmx, "energy", "-f", edr_path, "-o", xvg_path],
            input=f"{term_num}\n\n", capture_output=True, text=True, timeout=60
        )
        if proc.returncode != 0:
            raise RuntimeError(f"gmx energy failed: {proc.stderr}")

        times, vals = [], []
        with open(xvg_path) as fh:
            for line in fh:
                if line.startswith(('#', '@')):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    times.append(float(parts[0]))
                    vals.append(float(parts[1]))
        return times, vals
    finally:
        if os.path.exists(xvg_path):
            os.unlink(xvg_path)


def linear_regression(x, y):
    """Simple linear regression. Returns slope, intercept, R², stderr, p-value."""
    n = len(x)
    if n < 3:
        return None
    sx = sum(x)
    sy = sum(y)
    sxx = sum(xi * xi for xi in x)
    sxy = sum(xi * yi for xi, yi in zip(x, y))
    syy = sum(yi * yi for yi in y)

    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return None
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n

    ss_res = sum((yi - slope * xi - intercept) ** 2 for xi, yi in zip(x, y))
    ss_tot = syy - sy * sy / n
    r2 = 1.0 - ss_res / ss_tot if abs(ss_tot) > 1e-30 else 0.0

    se_slope = math.sqrt(ss_res / ((n - 2) * (sxx - sx * sx / n))) if n > 2 else float('inf')

    # t-statistic and p-value (two-tailed)
    if se_slope > 0:
        t_stat = slope / se_slope
        # Approximate p-value using student's t with n-2 df
        # For the purpose of this analysis, we use scipy if available
        try:
            from scipy import stats as sp_stats
            p_value = sp_stats.t.sf(abs(t_stat), n - 2) * 2
        except ImportError:
            p_value = None
    else:
        t_stat = float('inf')
        p_value = 0.0

    return {
        'slope': slope,
        'intercept': intercept,
        'r2': r2,
        'stderr': se_slope,
        't_stat': t_stat,
        'p_value': p_value,
        'n_points': n,
    }


def analyze_single_run(gmx, run_dir):
    """Analyze a single NVE run directory. Returns drift info dict or None."""
    edr_path = os.path.join(run_dir, "nve.edr")
    if not os.path.exists(edr_path):
        return None

    try:
        times, energies = extract_energy_xvg(gmx, edr_path, "Total-Energy")
    except Exception as e:
        return {"error": str(e)}

    if len(times) < 100:
        return {"error": f"Too few data points: {len(times)}"}

    # Get natoms from log
    log_path = os.path.join(run_dir, "nve.log")
    natoms = None
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                if "number of atoms" in line.lower() or "Atoms" in line:
                    parts = line.split()
                    for i, p in enumerate(parts):
                        if p == "Atoms" and i + 1 < len(parts):
                            try:
                                natoms = int(parts[i + 1])
                            except ValueError:
                                pass
        # Fallback: try from topology
        if natoms is None:
            for line in open(log_path):
                if "System" in line and "atoms" in line:
                    parts = line.split()
                    for p in parts:
                        try:
                            n = int(p)
                            if n > 100:  # reasonable atom count
                                natoms = n
                                break
                        except ValueError:
                            pass
    if natoms is None:
        natoms = 2652  # B5 default

    reg = linear_regression(times, energies)
    if reg is None:
        return {"error": "Linear regression failed"}

    drift_per_atom = reg['slope'] / natoms  # kJ/mol/ps/atom
    duration_ps = times[-1] - times[0]

    # Temperature
    try:
        t_times, temps = extract_energy_xvg(gmx, edr_path, "Temperature")
        temp_mean = sum(temps) / len(temps)
        temp_std = math.sqrt(sum((t - temp_mean) ** 2 for t in temps) / len(temps))
    except Exception:
        temp_mean = None
        temp_std = None

    # QM energy
    try:
        qm_times, qm_energies = extract_energy_xvg(gmx, edr_path, "Quantum-En.")
        qm_mean = sum(qm_energies) / len(qm_energies)
        qm_std = math.sqrt(sum((e - qm_mean) ** 2 for e in qm_energies) / len(qm_energies))
    except Exception:
        qm_mean = None
        qm_std = None

    return {
        'duration_ps': duration_ps,
        'n_frames': len(times),
        'natoms': natoms,
        'drift_slope': reg['slope'],
        'drift_per_atom': drift_per_atom,
        'drift_stderr': reg['stderr'] / natoms,
        'drift_r2': reg['r2'],
        'drift_p_value': reg.get('p_value'),
        'energy_mean': sum(energies) / len(energies),
        'energy_std': math.sqrt(sum((e - sum(energies)/len(energies))**2 for e in energies) / len(energies)),
        'temp_mean': temp_mean,
        'temp_std': temp_std,
        'qm_mean': qm_mean,
        'qm_std': qm_std,
        'pass': abs(drift_per_atom) < DRIFT_THRESHOLD,
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze A-2 statistical NVE drift data")
    parser.add_argument("basedir", help="Base directory with seed_N/{cutoff,pme}/ subdirs")
    parser.add_argument("--gmx", default=GMX_DEFAULT, help="Path to gmx binary")
    parser.add_argument("-o", "--output", help="Output JSON file")
    parser.add_argument("--seeds", default="1,2,3,4,5", help="Comma-separated seed numbers")
    args = parser.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    modes = ["cutoff", "pme"]

    results = {"version": VERSION, "date": str(date.today()), "seeds": seeds, "modes": {}}

    for mode in modes:
        mode_results = []
        for seed in seeds:
            run_dir = os.path.join(args.basedir, f"seed_{seed}", mode)
            print(f"  Analyzing {mode} seed {seed}...", end=" ", flush=True)
            r = analyze_single_run(args.gmx, run_dir)
            if r is None:
                print("MISSING")
                continue
            if "error" in r:
                print(f"ERROR: {r['error']}")
                mode_results.append({"seed": seed, "error": r["error"]})
                continue
            status = "PASS" if r['pass'] else "FAIL"
            print(f"drift={r['drift_per_atom']:.4e} kJ/mol/ps/atom [{status}]")
            r['seed'] = seed
            mode_results.append(r)

        # Compute statistics across seeds
        valid = [r for r in mode_results if 'drift_per_atom' in r]
        if valid:
            drifts = [r['drift_per_atom'] for r in valid]
            n = len(drifts)
            mean_drift = sum(drifts) / n
            if n > 1:
                std_drift = math.sqrt(sum((d - mean_drift) ** 2 for d in drifts) / (n - 1))
                sem_drift = std_drift / math.sqrt(n)
            else:
                std_drift = 0.0
                sem_drift = 0.0

            stats = {
                "n_seeds": n,
                "mean_drift": mean_drift,
                "std_drift": std_drift,
                "sem_drift": sem_drift,
                "min_drift": min(drifts),
                "max_drift": max(drifts),
                "all_pass": all(r.get('pass', False) for r in valid),
                "margin_over_threshold": DRIFT_THRESHOLD / abs(mean_drift) if abs(mean_drift) > 0 else float('inf'),
            }
        else:
            stats = {"n_seeds": 0, "error": "No valid results"}

        results["modes"][mode] = {"runs": mode_results, "statistics": stats}

    # Print summary
    print("\n" + "=" * 70)
    print("A-2 Statistical NVE Drift Summary")
    print("=" * 70)
    for mode in modes:
        s = results["modes"][mode]["statistics"]
        if "error" in s:
            print(f"\n{mode.upper()}: {s['error']}")
            continue
        print(f"\n{mode.upper()} embedding ({s['n_seeds']} seeds):")
        print(f"  Mean drift:  {s['mean_drift']:+.4e} kJ/mol/ps/atom")
        print(f"  Std drift:   {s['std_drift']:.4e} kJ/mol/ps/atom")
        print(f"  SEM:         {s['sem_drift']:.4e} kJ/mol/ps/atom")
        print(f"  Range:       [{s['min_drift']:.4e}, {s['max_drift']:.4e}]")
        print(f"  Margin:      {s['margin_over_threshold']:.1f}× below threshold ({DRIFT_THRESHOLD})")
        print(f"  All pass:    {'YES' if s['all_pass'] else 'NO'}")

    # PME/cutoff comparison
    if all(m in results["modes"] for m in modes):
        cutoff_s = results["modes"]["cutoff"]["statistics"]
        pme_s = results["modes"]["pme"]["statistics"]
        if "error" not in cutoff_s and "error" not in pme_s:
            ratio = abs(pme_s['mean_drift']) / abs(cutoff_s['mean_drift']) if abs(cutoff_s['mean_drift']) > 0 else float('inf')
            print(f"\nPME/cutoff drift ratio: {ratio:.2f}")
            results["pme_cutoff_ratio"] = ratio

    print("=" * 70)

    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults written to {args.output}")

    # Return exit code based on all_pass
    all_modes_pass = all(
        results["modes"][m]["statistics"].get("all_pass", False)
        for m in modes if "error" not in results["modes"][m]["statistics"]
    )
    return 0 if all_modes_pass else 1


if __name__ == "__main__":
    sys.exit(main())
