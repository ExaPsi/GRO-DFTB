#!/usr/bin/env python3
"""Analyze sigma (damping parameter) sensitivity sweep for A-7.

Runs drift analysis on each sigma directory and compiles a comparison.

Usage:
    python3 tools/analyze_sigma_sweep.py docs/session/20260213/a7_sigma_sweep/
    python3 tools/analyze_sigma_sweep.py docs/session/20260213/a7_sigma_sweep/ -o results.json
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from datetime import date

DRIFT_THRESHOLD = 0.005  # kJ/mol/ps/atom (reconciled GROMACS standard)
QM_ENERGY_REF = -10743.5  # kJ/mol (PME reference)
QM_ENERGY_TOL = 5.0  # kJ/mol


def extract_energy_data(gmx_binary, edr_file, tmpdir):
    """Extract total and QM energy from .edr file using gmx energy."""
    total_xvg = os.path.join(tmpdir, "total.xvg")
    qm_xvg = os.path.join(tmpdir, "qm.xvg")

    # Total energy
    proc = subprocess.run(
        [gmx_binary, "energy", "-f", edr_file, "-o", total_xvg],
        input="Total-Energy\n0\n",
        capture_output=True, text=True, timeout=60
    )
    if proc.returncode != 0:
        raise RuntimeError(f"gmx energy (total) failed: {proc.stderr}")

    # QM energy
    proc = subprocess.run(
        [gmx_binary, "energy", "-f", edr_file, "-o", qm_xvg],
        input="Quantum-En.\n0\n",
        capture_output=True, text=True, timeout=60
    )
    if proc.returncode != 0:
        raise RuntimeError(f"gmx energy (QM) failed: {proc.stderr}")

    total_t, total_e = parse_xvg(total_xvg)
    qm_t, qm_e = parse_xvg(qm_xvg)
    return total_t, total_e, qm_t, qm_e


def parse_xvg(xvg_file):
    """Parse XVG file into time and value arrays."""
    times, values = [], []
    with open(xvg_file) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                values.append(float(parts[1]))
    return times, values


def linear_regression(x, y):
    """Simple linear regression: y = a + b*x. Returns slope, intercept, R², stderr."""
    n = len(x)
    sx = sum(x)
    sy = sum(y)
    sxx = sum(xi * xi for xi in x)
    sxy = sum(xi * yi for xi, yi in zip(x, y))

    denom = n * sxx - sx * sx
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n

    y_mean = sy / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    ss_res = sum((yi - a - b * xi) ** 2 for xi, yi in zip(x, y))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    if n > 2:
        mse = ss_res / (n - 2)
        se_b = math.sqrt(mse / (sxx - sx * sx / n))
    else:
        se_b = 0.0

    return b, a, r_squared, se_b


def count_scc_failures(log_file):
    """Count SCC convergence failures from GROMACS log."""
    count = 0
    if not os.path.exists(log_file):
        return -1
    with open(log_file) as f:
        for line in f:
            if "SCC did not converge" in line or "SCCNotConverged" in line:
                count += 1
    return count


def count_consistency_warnings(log_file):
    """Count embedding consistency warnings from GROMACS log."""
    count = 0
    if not os.path.exists(log_file):
        return -1
    with open(log_file) as f:
        for line in f:
            if "EMBEDDING CONSISTENCY WARNING" in line:
                count += 1
    return count


def analyze_one_sigma(gmx_binary, sigma_dir, n_atoms, tmpdir):
    """Analyze a single sigma run."""
    edr_file = os.path.join(sigma_dir, "nve_pme.edr")
    log_file = os.path.join(sigma_dir, "nve_pme.log")

    if not os.path.exists(edr_file):
        return {"completed": False, "error": f"Missing {edr_file}"}

    total_t, total_e, qm_t, qm_e = extract_energy_data(gmx_binary, edr_file, tmpdir)

    duration_ps = total_t[-1] - total_t[0] if total_t else 0
    n_steps = len(total_t)

    # Drift analysis
    slope, intercept, r2, se_slope = linear_regression(total_t, total_e)
    drift_per_atom = slope / n_atoms

    # Drift standard error per atom
    drift_stderr = se_slope / n_atoms

    # QM energy statistics
    qm_mean = sum(qm_e) / len(qm_e)
    qm_std = math.sqrt(sum((q - qm_mean) ** 2 for q in qm_e) / len(qm_e))
    qm_min = min(qm_e)
    qm_max = max(qm_e)

    # Total energy statistics
    e_mean = sum(total_e) / len(total_e)
    e_std = math.sqrt(sum((e - e_mean) ** 2 for e in total_e) / len(total_e))

    # SCC failures
    scc_failures = count_scc_failures(log_file)
    consistency_warnings = count_consistency_warnings(log_file)

    # Block analysis (two halves)
    mid = len(total_t) // 2
    slope1, _, _, _ = linear_regression(total_t[:mid], total_e[:mid])
    slope2, _, _, _ = linear_regression(total_t[mid:], total_e[mid:])
    drift1 = slope1 / n_atoms
    drift2 = slope2 / n_atoms
    block_ratio = max(abs(drift1), abs(drift2)) / min(abs(drift1), abs(drift2)) if min(abs(drift1), abs(drift2)) > 0 else float('inf')

    # Energy spikes (3-sigma)
    spike_count = sum(1 for e in total_e if abs(e - e_mean) > 3 * e_std)

    return {
        "completed": True,
        "duration_ps": duration_ps,
        "n_frames": n_steps,
        "drift_kj_mol_ps_atom": drift_per_atom,
        "drift_slope_kj_mol_ps": slope,
        "drift_r_squared": r2,
        "drift_stderr": drift_stderr,
        "drift_pass": abs(drift_per_atom) < DRIFT_THRESHOLD,
        "drift_margin": DRIFT_THRESHOLD / abs(drift_per_atom) if drift_per_atom != 0 else float('inf'),
        "block_ratio": block_ratio,
        "block_drift_1st_half": drift1,
        "block_drift_2nd_half": drift2,
        "qm_energy_mean_kj_mol": qm_mean,
        "qm_energy_std_kj_mol": qm_std,
        "qm_energy_min_kj_mol": qm_min,
        "qm_energy_max_kj_mol": qm_max,
        "qm_energy_pass": abs(qm_mean - QM_ENERGY_REF) < QM_ENERGY_TOL,
        "total_energy_mean_kj_mol": e_mean,
        "total_energy_std_kj_mol": e_std,
        "scc_failures": scc_failures,
        "consistency_warnings": consistency_warnings,
        "energy_spikes_3sigma": spike_count,
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze sigma sensitivity sweep")
    parser.add_argument("base_dir", help="Base directory containing sigma_* subdirs")
    parser.add_argument("-o", "--output", help="Output JSON file")
    parser.add_argument("--gmx", default=None, help="Path to gmx binary")
    parser.add_argument("--n-atoms", type=int, default=2652, help="Number of atoms (default: 2652 for B5 3.0nm)")
    args = parser.parse_args()

    gmx = args.gmx
    if gmx is None:
        gmx = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "external", "gromacs", "_install_dftb", "bin", "gmx")
    if not os.path.exists(gmx):
        print(f"Error: gmx binary not found at {gmx}")
        sys.exit(1)

    sigma_values = [0.05, 0.10, 0.15, 0.20]
    sigma_dirs = {
        f"{s:.2f}": os.path.join(args.base_dir, f"sigma_{f'{s:.2f}'.replace('.','p')}")
        for s in sigma_values
    }

    results = {
        "story": "A-7",
        "date": str(date.today()),
        "gmx_binary": gmx,
        "base_dir": os.path.abspath(args.base_dir),
        "n_atoms": args.n_atoms,
        "thresholds": {
            "drift_kj_mol_ps_atom": DRIFT_THRESHOLD,
            "qm_energy_ref_kj_mol": QM_ENERGY_REF,
            "qm_energy_tol_kj_mol": QM_ENERGY_TOL,
        },
        "sigma_results": {},
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        for sigma_str, sigma_dir in sorted(sigma_dirs.items()):
            print(f"\n=== Analyzing sigma = {sigma_str} nm ===")
            if not os.path.isdir(sigma_dir):
                print(f"  Directory not found: {sigma_dir}")
                results["sigma_results"][sigma_str] = {"completed": False, "error": "Directory not found"}
                continue

            try:
                result = analyze_one_sigma(gmx, sigma_dir, args.n_atoms, tmpdir)
                results["sigma_results"][sigma_str] = result
                if result["completed"]:
                    print(f"  Drift: {result['drift_kj_mol_ps_atom']:.4e} kJ/mol/ps/atom "
                          f"(R²={result['drift_r_squared']:.4f}, "
                          f"margin={result['drift_margin']:.1f}×)")
                    print(f"  QM energy: {result['qm_energy_mean_kj_mol']:.1f} ± "
                          f"{result['qm_energy_std_kj_mol']:.2f} kJ/mol")
                    print(f"  SCC failures: {result['scc_failures']}")
                    print(f"  Pass: {'YES' if result['drift_pass'] else 'NO'}")
                else:
                    print(f"  Error: {result.get('error', 'Unknown')}")
            except Exception as e:
                print(f"  Analysis failed: {e}")
                results["sigma_results"][sigma_str] = {"completed": False, "error": str(e)}

    # Summary statistics
    completed = {k: v for k, v in results["sigma_results"].items() if v.get("completed")}
    if len(completed) >= 2:
        drifts = {k: v["drift_kj_mol_ps_atom"] for k, v in completed.items()}
        qm_energies = {k: v["qm_energy_mean_kj_mol"] for k, v in completed.items()}

        results["summary"] = {
            "n_completed": len(completed),
            "all_pass_drift": all(v["drift_pass"] for v in completed.values()),
            "all_pass_qm_energy": all(v["qm_energy_pass"] for v in completed.values()),
            "drift_range": [min(abs(d) for d in drifts.values()), max(abs(d) for d in drifts.values())],
            "drift_ratio": max(abs(d) for d in drifts.values()) / min(abs(d) for d in drifts.values()) if min(abs(d) for d in drifts.values()) > 0 else float('inf'),
            "qm_energy_spread_kj_mol": max(qm_energies.values()) - min(qm_energies.values()),
            "total_scc_failures": sum(v["scc_failures"] for v in completed.values()),
        }

        print("\n" + "=" * 60)
        print("SIGMA SWEEP SUMMARY")
        print("=" * 60)
        print(f"{'sigma (nm)':<12} {'Drift (kJ/mol/ps/atom)':<25} {'R²':<8} {'Margin':<8} {'QM E (kJ/mol)':<18} {'SCC fail':<10} {'Pass'}")
        print("-" * 100)
        for sigma_str in sorted(completed.keys()):
            v = completed[sigma_str]
            print(f"{sigma_str:<12} {v['drift_kj_mol_ps_atom']:>+.4e}{'':>12} "
                  f"{v['drift_r_squared']:.4f}  {v['drift_margin']:>6.1f}×  "
                  f"{v['qm_energy_mean_kj_mol']:>10.1f} ± {v['qm_energy_std_kj_mol']:>5.2f}  "
                  f"{v['scc_failures']:>6d}    {'YES' if v['drift_pass'] else 'NO'}")
        print("-" * 100)
        print(f"Drift ratio (max/min): {results['summary']['drift_ratio']:.2f}")
        print(f"QM energy spread: {results['summary']['qm_energy_spread_kj_mol']:.2f} kJ/mol")
        print(f"Overall: {'ALL PASS' if results['summary']['all_pass_drift'] else 'SOME FAIL'}")

    output_file = args.output or os.path.join(args.base_dir, "analysis_results.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults written to: {output_file}")


if __name__ == "__main__":
    main()
