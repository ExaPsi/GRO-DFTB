#!/usr/bin/env python3
"""Analyze long-duration PME NVE energy conservation for US-050.

Checks all 7 acceptance criteria:
  AC-1: PME NVE trajectory completes 100+ ps without crash
  AC-2: NVE drift < 0.01 kJ/mol/ps/atom over full trajectory
  AC-3: Zero SCC failures over entire trajectory
  AC-4: US-048 consistency check passes (0 warnings) throughout
  AC-5: Mean QM energy within +/-5 kJ/mol of PME reference (-10743.5 kJ/mol)
  AC-6: No late-onset instabilities: block drift ratio <= 3
  AC-7: Temperature physically reasonable (T = 300 +/- 10 K, sigma_T < 8 K)

Additional diagnostics:
  - Standard error of drift rate (T-US-050-2.1d)
  - Signal-to-noise ratio (T-US-050-2.1f)
  - Energy jump detection (T-US-050-H.1)
  - PME/cutoff drift comparison (T-US-050-3.1)

SDD:specs.md:S4 -- tools/ helper scripts
ThDD:T-US-050-2.1 -- linear regression for drift measurement
ThDD:T-US-050-2.2 -- two-block analysis for late-onset detection

Usage:
    python3 tools/analyze_pme_drift.py tests/data/b5_pme_drift/
    python3 tools/analyze_pme_drift.py tests/data/b5_pme_drift/ --gmx /path/to/gmx
    python3 tools/analyze_pme_drift.py tests/data/b5_pme_drift/ -o results.json
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

VERSION = "0.1.0"

# Acceptance criteria thresholds (from specs.md and US-050 story)
DRIFT_THRESHOLD = 0.01          # kJ/mol/ps/atom (specs.md section 21.1)
QM_ENERGY_REF = -10743.5        # kJ/mol (PME reference from US-047 verification)
QM_ENERGY_TOL = 5.0             # kJ/mol (US-050 AC-5)
BLOCK_DRIFT_RATIO_THRESHOLD = 3.0  # max ratio first/second half (US-050 AC-6)
TEMP_MEAN_LOW = 290.0           # K (US-050 AC-7)
TEMP_MEAN_HIGH = 310.0          # K (US-050 AC-7)
TEMP_STD_MAX = 8.0              # K (US-050 AC-7)
ENERGY_JUMP_SIGMA = 5.0         # sigma threshold for energy jumps (T-US-050-H.1)
MIN_DURATION_PS = 100.0         # minimum trajectory length (US-050 AC-1)

# Reference values from archived prior runs (Golden Rule: real calculations only)
CUTOFF_DRIFT = -4.53e-05        # kJ/mol/ps/atom (US-043b, tests/data/b5/drift_analysis_damped_500ps.json)
CUTOFF_QM_ENERGY = -10706.14    # kJ/mol (US-043b, tests/data/b5/drift_analysis_damped_500ps.json)
PME_CUTOFF_RATIO_MIN = 2.0      # Expected range (T-US-050-3.1)
PME_CUTOFF_RATIO_MAX = 20.0     # Expected range (T-US-050-3.1)


def find_file(directory, extensions, preferred_names=None):
    """Find a file in a directory by extension, preferring specific names."""
    if preferred_names is None:
        preferred_names = []

    for name in preferred_names:
        path = os.path.join(directory, name)
        if os.path.isfile(path):
            return path

    for f in sorted(os.listdir(directory)):
        if any(f.endswith(ext) for ext in extensions):
            return os.path.join(directory, f)
    return None


def find_edr_file(data_dir):
    """Find the .edr file."""
    return find_file(data_dir, [".edr"], ["md.edr", "nve_pme.edr"])


def find_log_file(data_dir):
    """Find the .log file."""
    return find_file(data_dir, [".log"], ["md.log", "nve_pme.log"])


def find_gro_file(data_dir):
    """Find a .gro file."""
    return find_file(data_dir, [".gro"], ["md.gro", "confout.gro", "initial.gro"])


def get_n_atoms(data_dir, log_path):
    """Get atom count from GRO file or log file."""
    gro_path = find_gro_file(data_dir)
    if gro_path:
        try:
            with open(gro_path) as f:
                f.readline()  # title
                return int(f.readline().strip())
        except (ValueError, IOError):
            pass

    if log_path:
        with open(log_path) as f:
            for line in f:
                m = re.search(r"number of atoms\s*=\s*(\d+)", line, re.IGNORECASE)
                if m:
                    return int(m.group(1))
                m = re.search(r"^\s*Atoms:\s+(\d+)", line)
                if m:
                    return int(m.group(1))
    return None


def extract_energy_xvg(gmx_bin, edr_path, energy_term):
    """Extract energy time series from .edr using gmx energy.

    Returns list of (time_ps, energy_kj_mol) tuples, or None on failure.

    ThDD:T-US-050-U.2 -- Time in ps (GROMACS native)
    ThDD:T-US-050-U.1 -- Energy in kJ/mol (GROMACS native)
    """
    with tempfile.NamedTemporaryFile(suffix=".xvg", delete=False) as tmp:
        xvg_path = tmp.name

    try:
        stdin_text = f"{energy_term}\n\n"
        result = subprocess.run(
            [gmx_bin, "energy", "-f", edr_path, "-o", xvg_path],
            input=stdin_text,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            print(f"  WARNING: gmx energy failed for {energy_term}:", file=sys.stderr)
            print(f"  stderr: {result.stderr[:500]}", file=sys.stderr)
            return None

        return parse_xvg(xvg_path)
    except subprocess.TimeoutExpired:
        print(f"  WARNING: gmx energy timed out for {energy_term}", file=sys.stderr)
        return None
    except FileNotFoundError:
        print(f"  ERROR: gmx binary not found: {gmx_bin}", file=sys.stderr)
        return None
    finally:
        if os.path.isfile(xvg_path):
            os.unlink(xvg_path)


def parse_xvg(xvg_path):
    """Parse XVG file, returning list of (time, value) tuples."""
    data = []
    with open(xvg_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    t = float(parts[0])
                    v = float(parts[1])
                    data.append((t, v))
                except ValueError:
                    continue
    return data if data else None


def compute_statistics(values):
    """Compute mean, std, min, max of a list of floats."""
    n = len(values)
    if n == 0:
        return {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0, "n": 0}

    mean = sum(values) / n
    if n > 1:
        variance = sum((v - mean) ** 2 for v in values) / (n - 1)
        std = variance ** 0.5
    else:
        std = 0.0

    return {
        "mean": mean,
        "std": std,
        "min": min(values),
        "max": max(values),
        "n": n,
    }


def compute_drift_with_stderr(data):
    """Compute drift rate with standard error via least-squares linear regression.

    ThDD:T-US-050-2.1  -- E(t) = E_0 + alpha * t + epsilon(t)
    ThDD:T-US-050-2.1c -- alpha = S_tE / S_tt (least-squares estimator)
    ThDD:T-US-050-2.1d -- sigma_alpha = sqrt(SSR / ((N-2) * S_tt))

    data: list of (time_ps, energy_kj_mol) tuples.
    Returns dict with slope, intercept, r_squared, stderr, n_points.
    """
    n = len(data)
    if n < 3:
        return {
            "slope": 0.0, "intercept": 0.0, "r_squared": 0.0,
            "stderr": float("inf"), "n_points": n,
        }

    # Compute sums for least-squares
    sum_t = sum(d[0] for d in data)
    sum_e = sum(d[1] for d in data)
    sum_tt = sum(d[0] * d[0] for d in data)
    sum_te = sum(d[0] * d[1] for d in data)

    mean_t = sum_t / n
    mean_e = sum_e / n

    # S_tt = sum (t_i - t_bar)^2
    s_tt = sum_tt - n * mean_t * mean_t
    # S_tE = sum (t_i - t_bar)(E_i - E_bar)
    s_te = sum_te - n * mean_t * mean_e

    if abs(s_tt) < 1e-30:
        return {
            "slope": 0.0, "intercept": mean_e, "r_squared": 0.0,
            "stderr": float("inf"), "n_points": n,
        }

    # Slope (drift rate) and intercept  (T-US-050-2.1c)
    slope = s_te / s_tt
    intercept = mean_e - slope * mean_t

    # R-squared
    ss_tot = sum((d[1] - mean_e) ** 2 for d in data)
    ss_res = sum((d[1] - (slope * d[0] + intercept)) ** 2 for d in data)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # Standard error of slope  (T-US-050-2.1d)
    # sigma_alpha = sqrt(SSR / ((N-2) * S_tt))
    if n > 2 and s_tt > 0:
        stderr = math.sqrt(ss_res / ((n - 2) * s_tt))
    else:
        stderr = float("inf")

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "stderr": stderr,
        "n_points": n,
    }


def block_analysis(data):
    """Two-block analysis: split trajectory at midpoint and compare drift rates.

    ThDD:T-US-050-2.2  -- alpha_1 = drift(0 -> T/2), alpha_2 = drift(T/2 -> T)
    ThDD:T-US-050-2.2b -- R = |alpha_2| / |alpha_1| <= 3

    Returns dict with first_half, second_half regression results and ratio.
    """
    n = len(data)
    mid = n // 2

    if mid < 10:
        return {
            "first_half": None,
            "second_half": None,
            "ratio": None,
            "pass": False,
            "error": "Insufficient data for block analysis",
        }

    first_half = compute_drift_with_stderr(data[:mid])
    second_half = compute_drift_with_stderr(data[mid:])

    # Compute ratio |alpha_2| / |alpha_1|
    abs_alpha1 = abs(first_half["slope"])
    abs_alpha2 = abs(second_half["slope"])

    if abs_alpha1 > 1e-30:
        ratio = abs_alpha2 / abs_alpha1
    else:
        # First half drift is essentially zero — any second half drift is finite/zero
        ratio = 0.0 if abs_alpha2 < 1e-30 else float("inf")

    # Special case: both halves have very small drift (below noise floor)
    # In this case, the ratio is meaningless — pass automatically
    both_tiny = abs_alpha1 < 1e-8 and abs_alpha2 < 1e-8

    return {
        "first_half": first_half,
        "second_half": second_half,
        "ratio": ratio,
        "pass": ratio <= BLOCK_DRIFT_RATIO_THRESHOLD or both_tiny,
    }


def detect_energy_jumps(data, sigma_threshold=ENERGY_JUMP_SIGMA):
    """Detect energy discontinuities larger than sigma_threshold * sigma.

    ThDD:T-US-050-H.1 -- Shadow Hamiltonian: jumps indicate non-symplectic behavior.

    Returns dict with count, locations, and expected count.
    """
    if len(data) < 10:
        return {"count": 0, "locations_ps": [], "expected": 0, "n_points": len(data)}

    energies = [d[1] for d in data]
    stats = compute_statistics(energies)
    n = stats["n"]

    if stats["std"] < 1e-30:
        return {"count": 0, "locations_ps": [], "expected": 0, "n_points": n}

    # Count jumps beyond sigma_threshold standard deviations
    locations = []
    for t, e in data:
        if abs(e - stats["mean"]) > sigma_threshold * stats["std"]:
            locations.append(t)

    # Expected count from Gaussian statistics
    # P(|X| > k*sigma) ~ erfc(k/sqrt(2))
    from math import erfc
    p_tail = erfc(sigma_threshold / math.sqrt(2.0))
    expected = n * p_tail

    return {
        "count": len(locations),
        "locations_ps": locations[:20],  # first 20 for brevity
        "expected": expected,
        "n_points": n,
    }


def qm_energy_statistics(data):
    """Compute QM energy statistics and compare with PME reference.

    ThDD:T-US-050-4.1 -- |<E_QM> - E_ref| < 5 kJ/mol

    Returns dict with mean, std, min, max, deviation from ref, pass.
    """
    if not data:
        return None

    values = [d[1] for d in data]
    stats = compute_statistics(values)
    deviation = abs(stats["mean"] - QM_ENERGY_REF)

    return {
        "mean": stats["mean"],
        "std": stats["std"],
        "min": stats["min"],
        "max": stats["max"],
        "n": stats["n"],
        "reference": QM_ENERGY_REF,
        "deviation": deviation,
        "pass": deviation < QM_ENERGY_TOL,
    }


def check_scc_convergence(log_path):
    """Parse md.log for SCC non-convergence events.

    ThDD:T-US-050-1.23 -- Gaussian damping bounds potential to Phi_max = 5.97 Ha/e
    ThDD:T-US-050-1.24 -- No SCC failure expected with damping

    Returns dict with failure count.
    """
    count = 0
    if log_path and os.path.isfile(log_path):
        with open(log_path) as f:
            for line in f:
                if "SCC is NOT converged" in line:
                    count += 1
    return {"failures": count, "pass": count == 0}


def check_consistency(log_path):
    """Parse md.log for US-048 consistency check failures.

    ThDD:T-US-048-3.2/3.3/4.3 -- Three algebraic identity checks

    Returns dict with warning count.
    """
    count = 0
    if log_path and os.path.isfile(log_path):
        with open(log_path) as f:
            for line in f:
                if "FAILED" in line or "EXCEEDED" in line:
                    if any(kw in line.lower() for kw in
                           ("consistency", "relation", "embedding", "check", "tolerance")):
                        count += 1
    return {"warnings": count, "pass": count == 0}


def temperature_statistics(data):
    """Compute temperature statistics and verify microcanonical ensemble.

    ThDD:T-US-050-5.1 -- <T> = 2<E_kin> / ((3N-6) k_B)
    ThDD:T-US-050-5.2 -- sigma_T = T * sqrt(2/(3N-6)) ~ 4.8 K for N=2652

    Returns dict with mean, std, pass.
    """
    if not data:
        return None

    values = [d[1] for d in data]
    stats = compute_statistics(values)

    mean_pass = TEMP_MEAN_LOW <= stats["mean"] <= TEMP_MEAN_HIGH
    std_pass = stats["std"] < TEMP_STD_MAX
    overall = mean_pass and std_pass

    return {
        "mean": stats["mean"],
        "std": stats["std"],
        "min": stats["min"],
        "max": stats["max"],
        "n": stats["n"],
        "mean_pass": mean_pass,
        "std_pass": std_pass,
        "pass": overall,
        "theory_std": 4.8,  # T * sqrt(2/(3*2652-6)) for T=300K, N=2652
    }


def compare_with_cutoff(pme_drift):
    """Compare PME drift with cutoff drift from US-043b.

    ThDD:T-US-050-3.1 -- drift_PME / drift_cutoff in [2, 20]

    Returns dict with ratio and assessment.
    """
    if pme_drift is None or abs(CUTOFF_DRIFT) < 1e-30:
        return {"ratio": None, "in_range": False, "error": "Missing data"}

    ratio = abs(pme_drift) / abs(CUTOFF_DRIFT)

    return {
        "pme_drift": pme_drift,
        "cutoff_drift": CUTOFF_DRIFT,
        "ratio": ratio,
        "expected_range": [PME_CUTOFF_RATIO_MIN, PME_CUTOFF_RATIO_MAX],
        "in_range": PME_CUTOFF_RATIO_MIN <= ratio <= PME_CUTOFF_RATIO_MAX,
    }


def compute_snr(drift_slope, duration_ps, energy_std):
    """Compute signal-to-noise ratio for drift detection.

    ThDD:T-US-050-2.1f -- SNR = |alpha| * T / sigma_E

    Returns dict with SNR and reliability assessment.
    """
    if energy_std < 1e-30 or duration_ps < 1e-30:
        return {"snr": 0.0, "reliable": False}

    snr = abs(drift_slope) * duration_ps / energy_std

    if snr > 10:
        reliability = "Excellent"
    elif snr > 3:
        reliability = "Good"
    elif snr > 1:
        reliability = "Marginal"
    else:
        reliability = "Unreliable (noise-dominated)"

    return {
        "snr": snr,
        "reliability": reliability,
        "reliable": snr > 1.0,
    }


def check_log_completion(log_path):
    """Check if GROMACS log file indicates successful completion."""
    if not log_path or not os.path.isfile(log_path):
        return False
    with open(log_path) as f:
        content = f.read()
    if "Finished mdrun" in content:
        return True
    if "Performance:" in content and "Wall t" in content:
        return True
    return False


def get_duration_and_steps(log_path):
    """Extract simulation duration and step count from log."""
    nsteps = None
    dt = None
    if log_path and os.path.isfile(log_path):
        with open(log_path) as f:
            for line in f:
                m = re.search(r"nsteps\s*=\s*(\d+)", line)
                if m:
                    nsteps = int(m.group(1))
                m = re.search(r"\bdt\s*=\s*([\d.eE+-]+)", line)
                if m:
                    dt = float(m.group(1))
    duration = None
    if nsteps and dt:
        duration = nsteps * dt
    return nsteps, duration


def analyze(data_dir, gmx_bin):
    """Run full analysis on the PME drift data directory.

    Returns dict with all results.
    """
    results = {
        "story": "US-050",
        "date": str(date.today()),
        "data_dir": os.path.abspath(data_dir),
        "gmx_binary": gmx_bin,
    }

    # Find files
    edr_path = find_edr_file(data_dir)
    log_path = find_log_file(data_dir)

    if edr_path is None:
        results["error"] = "No .edr file found"
        return results
    if log_path is None:
        results["error"] = "No .log file found"
        return results

    print(f"Data directory: {data_dir}")
    print(f"EDR file: {os.path.basename(edr_path)}")
    print(f"Log file: {os.path.basename(log_path)}")

    # Get atom count
    n_atoms = get_n_atoms(data_dir, log_path)
    if n_atoms is None:
        results["error"] = "Cannot determine atom count"
        return results
    results["n_atoms"] = n_atoms
    print(f"Atoms: {n_atoms}")

    # Check completion (AC-1)
    completed = check_log_completion(log_path)
    nsteps, expected_duration = get_duration_and_steps(log_path)
    results["completed"] = completed
    results["nsteps"] = nsteps
    results["expected_duration_ps"] = expected_duration
    print(f"Completed: {completed}")

    # SCC convergence (AC-3)
    scc = check_scc_convergence(log_path)
    results["scc"] = scc
    print(f"SCC failures: {scc['failures']} ({'PASS' if scc['pass'] else 'FAIL'})")

    # Consistency check (AC-4)
    consistency = check_consistency(log_path)
    results["consistency"] = consistency
    print(f"Consistency warnings: {consistency['warnings']} "
          f"({'PASS' if consistency['pass'] else 'FAIL'})")

    # --- Extract energy time series ---

    # Total energy (AC-2, AC-6)
    print("\nExtracting Total-Energy...")
    total_energy_data = extract_energy_xvg(gmx_bin, edr_path, "Total-Energy")

    if total_energy_data and len(total_energy_data) > 1:
        duration_ps = total_energy_data[-1][0] - total_energy_data[0][0]
        results["duration_ps"] = duration_ps
        results["n_energy_frames"] = len(total_energy_data)
        print(f"  Duration: {duration_ps:.3f} ps ({len(total_energy_data)} frames)")

        # Total energy statistics
        total_energies = [d[1] for d in total_energy_data]
        results["total_energy"] = compute_statistics(total_energies)

        # AC-2: Drift measurement with standard error
        drift_result = compute_drift_with_stderr(total_energy_data)
        drift_per_atom = drift_result["slope"] / n_atoms
        drift_stderr_per_atom = drift_result["stderr"] / n_atoms
        drift_pass = abs(drift_per_atom) < DRIFT_THRESHOLD

        results["drift"] = {
            "slope_kj_mol_ps": drift_result["slope"],
            "intercept_kj_mol": drift_result["intercept"],
            "r_squared": drift_result["r_squared"],
            "stderr_kj_mol_ps": drift_result["stderr"],
            "drift_kj_mol_ps_atom": drift_per_atom,
            "drift_stderr_kj_mol_ps_atom": drift_stderr_per_atom,
            "margin": DRIFT_THRESHOLD / abs(drift_per_atom) if abs(drift_per_atom) > 1e-30 else float("inf"),
            "pass": drift_pass,
        }

        print(f"  Drift: {drift_per_atom:.6e} +/- {drift_stderr_per_atom:.2e} kJ/mol/ps/atom")
        print(f"  Drift pass: {'PASS' if drift_pass else 'FAIL'} "
              f"(margin: {results['drift']['margin']:.1f}x)")
        print(f"  R²: {drift_result['r_squared']:.6f}")

        # AC-6: Block analysis
        block = block_analysis(total_energy_data)
        results["block_analysis"] = {
            "first_half_slope": block["first_half"]["slope"] if block["first_half"] else None,
            "second_half_slope": block["second_half"]["slope"] if block["second_half"] else None,
            "first_half_r_squared": block["first_half"]["r_squared"] if block["first_half"] else None,
            "second_half_r_squared": block["second_half"]["r_squared"] if block["second_half"] else None,
            "ratio": block["ratio"],
            "pass": block["pass"],
        }

        if block["first_half"] and block["second_half"]:
            fh_per_atom = block["first_half"]["slope"] / n_atoms
            sh_per_atom = block["second_half"]["slope"] / n_atoms
            print(f"\n  Block analysis:")
            print(f"    First half drift:  {fh_per_atom:.6e} kJ/mol/ps/atom "
                  f"(R² = {block['first_half']['r_squared']:.4f})")
            print(f"    Second half drift: {sh_per_atom:.6e} kJ/mol/ps/atom "
                  f"(R² = {block['second_half']['r_squared']:.4f})")
            print(f"    Ratio: {block['ratio']:.2f} "
                  f"({'PASS' if block['pass'] else 'FAIL'}, threshold: {BLOCK_DRIFT_RATIO_THRESHOLD})")

        # Energy jump detection
        jumps = detect_energy_jumps(total_energy_data)
        results["energy_jumps"] = jumps
        print(f"\n  Energy jumps (>{ENERGY_JUMP_SIGMA}σ): {jumps['count']} "
              f"(expected ~{jumps['expected']:.1f} from statistics)")

        # SNR
        snr = compute_snr(drift_result["slope"], duration_ps, results["total_energy"]["std"])
        results["snr"] = snr
        print(f"  SNR: {snr['snr']:.2f} ({snr['reliability']})")

        # PME/cutoff comparison
        comparison = compare_with_cutoff(drift_per_atom)
        results["pme_cutoff_comparison"] = comparison
        if comparison["ratio"] is not None:
            print(f"\n  PME/cutoff drift ratio: {comparison['ratio']:.2f} "
                  f"(expected {PME_CUTOFF_RATIO_MIN}–{PME_CUTOFF_RATIO_MAX})")
    else:
        print("  WARNING: Could not extract Total-Energy")
        results["duration_ps"] = 0.0

    # QM energy (AC-5)
    print("\nExtracting Quantum-En....")
    qm_energy_data = extract_energy_xvg(gmx_bin, edr_path, "Quantum-En.")
    if qm_energy_data:
        qm_stats = qm_energy_statistics(qm_energy_data)
        results["qm_energy"] = qm_stats
        print(f"  QM energy: {qm_stats['mean']:.2f} +/- {qm_stats['std']:.2f} kJ/mol")
        print(f"  Deviation from ref: {qm_stats['deviation']:.2f} kJ/mol "
              f"({'PASS' if qm_stats['pass'] else 'FAIL'})")
    else:
        print("  WARNING: Could not extract Quantum-En.")

    # Temperature (AC-7)
    print("\nExtracting Temperature...")
    temp_data = extract_energy_xvg(gmx_bin, edr_path, "Temperature")
    if temp_data:
        temp_stats = temperature_statistics(temp_data)
        results["temperature"] = temp_stats
        print(f"  Temperature: {temp_stats['mean']:.2f} +/- {temp_stats['std']:.2f} K")
        print(f"  Mean in [290, 310] K: {'PASS' if temp_stats['mean_pass'] else 'FAIL'}")
        print(f"  Std < 8 K: {'PASS' if temp_stats['std_pass'] else 'FAIL'}")
    else:
        print("  WARNING: Could not extract Temperature")

    # --- Evaluate acceptance criteria ---
    ac = evaluate_acceptance_criteria(results)
    results["acceptance_criteria"] = ac

    overall = all(v["pass"] for v in ac.values())
    results["overall_pass"] = overall

    return results


def evaluate_acceptance_criteria(results):
    """Evaluate all 7 acceptance criteria from the analysis results."""
    ac = {}

    # AC-1: Completes 100+ ps without crash
    duration = results.get("duration_ps", 0.0)
    ac["AC-1"] = {
        "description": "PME NVE trajectory completes 100+ ps without crash",
        "pass": results.get("completed", False) and duration >= MIN_DURATION_PS,
        "duration_ps": duration,
        "completed": results.get("completed", False),
    }

    # AC-2: NVE drift < 0.01 kJ/mol/ps/atom
    drift = results.get("drift", {})
    ac["AC-2"] = {
        "description": "NVE drift < 0.01 kJ/mol/ps/atom",
        "pass": drift.get("pass", False),
        "drift_kj_mol_ps_atom": drift.get("drift_kj_mol_ps_atom"),
        "drift_stderr": drift.get("drift_stderr_kj_mol_ps_atom"),
        "margin": drift.get("margin"),
    }

    # AC-3: Zero SCC failures
    scc = results.get("scc", {})
    ac["AC-3"] = {
        "description": "Zero SCC failures over entire trajectory",
        "pass": scc.get("pass", False),
        "failures": scc.get("failures"),
    }

    # AC-4: US-048 consistency check passes
    consistency = results.get("consistency", {})
    ac["AC-4"] = {
        "description": "US-048 consistency check passes (0 warnings)",
        "pass": consistency.get("pass", False),
        "warnings": consistency.get("warnings"),
    }

    # AC-5: QM energy within +/-5 kJ/mol of reference
    qm = results.get("qm_energy", {})
    ac["AC-5"] = {
        "description": "Mean QM energy within +/-5 kJ/mol of -10743.5 kJ/mol",
        "pass": qm.get("pass", False) if qm else False,
        "mean": qm.get("mean") if qm else None,
        "deviation": qm.get("deviation") if qm else None,
    }

    # AC-6: Block drift ratio <= 3
    block = results.get("block_analysis", {})
    ac["AC-6"] = {
        "description": "No late-onset instabilities: block drift ratio <= 3",
        "pass": block.get("pass", False),
        "ratio": block.get("ratio"),
    }

    # AC-7: Temperature physically reasonable
    temp = results.get("temperature", {})
    ac["AC-7"] = {
        "description": "Temperature physically reasonable (300 +/- 10 K, sigma < 8 K)",
        "pass": temp.get("pass", False) if temp else False,
        "mean": temp.get("mean") if temp else None,
        "std": temp.get("std") if temp else None,
    }

    return ac


def print_summary(results):
    """Print human-readable summary."""
    print("\n" + "=" * 78)
    print("US-050: Long-Duration PME NVE Energy Conservation Test")
    print("=" * 78)

    ac = results.get("acceptance_criteria", {})

    print(f"\n{'AC':>5} | {'Result':>6} | Description")
    print("-" * 70)
    for ac_id in sorted(ac.keys()):
        v = ac[ac_id]
        status = "PASS" if v["pass"] else "FAIL"
        desc = v["description"]

        # Add key metric
        detail = ""
        if ac_id == "AC-1":
            detail = f" [{v.get('duration_ps', 0):.1f} ps]"
        elif ac_id == "AC-2" and v.get("drift_kj_mol_ps_atom") is not None:
            detail = f" [{v['drift_kj_mol_ps_atom']:.2e}, {v.get('margin', 0):.0f}x margin]"
        elif ac_id == "AC-3":
            detail = f" [{v.get('failures', '?')} failures]"
        elif ac_id == "AC-4":
            detail = f" [{v.get('warnings', '?')} warnings]"
        elif ac_id == "AC-5" and v.get("deviation") is not None:
            detail = f" [delta={v['deviation']:.2f} kJ/mol]"
        elif ac_id == "AC-6" and v.get("ratio") is not None:
            detail = f" [ratio={v['ratio']:.2f}]"
        elif ac_id == "AC-7" and v.get("mean") is not None:
            detail = f" [{v['mean']:.1f} +/- {v.get('std', 0):.1f} K]"

        print(f"{ac_id:>5} | {status:>6} | {desc}{detail}")

    overall = results.get("overall_pass", False)
    print(f"\n{'':>5} | {'PASS' if overall else 'FAIL':>6} | OVERALL")
    print("=" * 78)

    # M2b gate criterion
    print(f"\nM2b-D4 (Energy conservation under PBC): {'SATISFIED' if overall else 'NOT SATISFIED'}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze long-duration PME NVE energy conservation for US-050.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s tests/data/b5_pme_drift/\n"
            "  %(prog)s tests/data/b5_pme_drift/ --gmx /path/to/gmx\n"
            "  %(prog)s tests/data/b5_pme_drift/ -o custom_results.json\n"
        ),
    )
    parser.add_argument(
        "data_dir",
        help="Directory containing md.edr, md.log (or nve_pme.*)",
    )
    parser.add_argument(
        "--gmx",
        default="external/gromacs/_install_dftb/bin/gmx",
        help="Path to gmx binary (default: external/gromacs/_install_dftb/bin/gmx)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output JSON file (default: <data_dir>/analysis_results.json)",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}",
    )

    args = parser.parse_args()

    data_dir = os.path.abspath(args.data_dir)
    if not os.path.isdir(data_dir):
        print(f"ERROR: Data directory not found: {data_dir}", file=sys.stderr)
        sys.exit(1)

    gmx_bin = args.gmx
    if not os.path.isabs(gmx_bin):
        gmx_bin = os.path.abspath(gmx_bin)
    if not os.path.isfile(gmx_bin):
        print(f"WARNING: gmx binary not found at {gmx_bin}", file=sys.stderr)

    output_path = args.output or os.path.join(data_dir, "analysis_results.json")

    # Run full analysis
    results = analyze(data_dir, gmx_bin)

    # Print summary
    print_summary(results)

    # Serialize results to JSON
    # Handle inf/nan for JSON serialization
    def sanitize_for_json(obj):
        if isinstance(obj, float):
            if math.isinf(obj) or math.isnan(obj):
                return str(obj)
            return obj
        elif isinstance(obj, dict):
            return {k: sanitize_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [sanitize_for_json(v) for v in obj]
        return obj

    output_data = sanitize_for_json(results)

    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"\nResults written to: {output_path}")

    sys.exit(0 if results.get("overall_pass", False) else 1)


if __name__ == "__main__":
    main()
