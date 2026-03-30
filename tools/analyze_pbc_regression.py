#!/usr/bin/env python3
"""Analyze PBC regression test results across 3 box sizes for PME QM/MM embedding.

Checks all 7 acceptance criteria for US-049:
  AC-1: All 3 boxes complete without crash
  AC-2: NVE drift < 0.01 kJ/mol/ps/atom per box
  AC-3: US-048 consistency check: 0 warnings
  AC-4: Drift rates comparable: max/min ratio < 10
  AC-5: Mean QM energy within +/-5 kJ/mol of -10743.5 kJ/mol (PME reference)
  AC-6: Zero SCC failures
  AC-7: No edge artifacts at 2.5 nm (no energy spikes > 3 sigma)

SDD:specs.md:S4 -- tools/ helper scripts

Usage:
    python3 tools/analyze_pbc_regression.py tests/data/b5_pbc_regression/
    python3 tools/analyze_pbc_regression.py tests/data/b5_pbc_regression/ --gmx /path/to/gmx
    python3 tools/analyze_pbc_regression.py tests/data/b5_pbc_regression/ -o results.json
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from datetime import date

VERSION = "0.1.0"

# Subdirectory names and their box sizes in nm
BOX_DIRS = {
    "b5_2p5nm": 2.5,
    "b5_3p0nm": 3.0,
    "b5_4p0nm": 4.0,
}

# Acceptance criteria thresholds (from specs.md and US-049 story)
DRIFT_THRESHOLD = 0.01          # kJ/mol/ps/atom (specs.md section 21.1)
DRIFT_RATIO_THRESHOLD = 10.0    # max/min ratio (US-049 AC-4)
QM_ENERGY_REF = -10743.5        # kJ/mol (PME reference from US-047 verification, B5 3.0nm box)
QM_ENERGY_TOL = 5.0             # kJ/mol (US-049 AC-5)
SPIKE_SIGMA = 3.0               # sigma threshold for edge artifacts (US-049 AC-7)


def find_edr_file(box_dir):
    """Find the .edr file in a box directory. Tries common names."""
    for name in ("nve_pme.edr", "md.edr", "nve.edr"):
        path = os.path.join(box_dir, name)
        if os.path.isfile(path):
            return path
    # Fallback: any .edr file
    for f in os.listdir(box_dir):
        if f.endswith(".edr"):
            return os.path.join(box_dir, f)
    return None


def find_log_file(box_dir):
    """Find the .log file in a box directory. Tries common names."""
    for name in ("nve_pme.log", "md.log", "nve.log"):
        path = os.path.join(box_dir, name)
        if os.path.isfile(path):
            return path
    for f in os.listdir(box_dir):
        if f.endswith(".log"):
            return os.path.join(box_dir, f)
    return None


def find_gro_file(box_dir):
    """Find a .gro file to read atom count from."""
    for name in ("initial.gro", "nve_pme.gro", "confout.gro", "md.gro"):
        path = os.path.join(box_dir, name)
        if os.path.isfile(path):
            return path
    for f in os.listdir(box_dir):
        if f.endswith(".gro"):
            return os.path.join(box_dir, f)
    return None


def get_n_atoms_from_gro(gro_path):
    """Read atom count from second line of GRO file."""
    with open(gro_path) as f:
        f.readline()  # title
        n_atoms = int(f.readline().strip())
    return n_atoms


def get_n_atoms_from_log(log_path):
    """Read atom count from GROMACS log file."""
    with open(log_path) as f:
        for line in f:
            # "System total charge:" line comes after atom listing
            m = re.search(r"number of atoms\s*=\s*(\d+)", line, re.IGNORECASE)
            if m:
                return int(m.group(1))
            # Alternative: "Atoms:" in topology summary
            m = re.search(r"^\s*Atoms:\s+(\d+)", line)
            if m:
                return int(m.group(1))
    return None


def extract_energy_xvg(gmx_bin, edr_path, energy_term):
    """Extract energy time series from .edr using gmx energy.

    Returns list of (time_ps, energy_kj_mol) tuples, or None on failure.
    """
    with tempfile.NamedTemporaryFile(suffix=".xvg", delete=False) as tmp:
        xvg_path = tmp.name

    try:
        # Pipe the energy term name to gmx energy stdin
        # Two newlines: first selects the term, second confirms
        stdin_text = f"{energy_term}\n\n"
        result = subprocess.run(
            [gmx_bin, "energy", "-f", edr_path, "-o", xvg_path],
            input=stdin_text,
            capture_output=True,
            text=True,
            timeout=60,
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


def linear_regression(data):
    """Compute slope and intercept via least-squares.

    data: list of (x, y) tuples.
    Returns (slope, intercept, r_squared).
    """
    n = len(data)
    if n < 2:
        return 0.0, 0.0, 0.0

    sum_x = sum(d[0] for d in data)
    sum_y = sum(d[1] for d in data)
    sum_xx = sum(d[0] * d[0] for d in data)
    sum_xy = sum(d[0] * d[1] for d in data)

    denom = n * sum_xx - sum_x * sum_x
    if abs(denom) < 1e-30:
        return 0.0, sum_y / n, 0.0

    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n

    # R-squared
    mean_y = sum_y / n
    ss_tot = sum((d[1] - mean_y) ** 2 for d in data)
    ss_res = sum((d[1] - (slope * d[0] + intercept)) ** 2 for d in data)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return slope, intercept, r_squared


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


def count_energy_spikes(values, sigma_threshold):
    """Count values more than sigma_threshold standard deviations from mean."""
    stats = compute_statistics(values)
    if stats["n"] < 3 or stats["std"] < 1e-30:
        return 0
    return sum(
        1 for v in values if abs(v - stats["mean"]) > sigma_threshold * stats["std"]
    )


def check_log_completion(log_path):
    """Check if GROMACS log file indicates successful completion."""
    with open(log_path) as f:
        content = f.read()

    # GROMACS prints this line on successful completion
    if "Finished mdrun" in content:
        return True
    # Alternative completion markers
    if "Performance:" in content and "Wall t" in content:
        return True
    return False


def count_scc_failures(log_path):
    """Count SCC non-convergence events in GROMACS log."""
    count = 0
    with open(log_path) as f:
        for line in f:
            if "SCC is NOT converged" in line:
                count += 1
    return count


def count_consistency_warnings(log_path):
    """Count US-048 consistency check failures/warnings in GROMACS log."""
    count = 0
    with open(log_path) as f:
        for line in f:
            if "FAILED" in line or "EXCEEDED" in line:
                # Filter to consistency-check lines (avoid false matches)
                if "consistency" in line.lower() or "relation" in line.lower() \
                        or "embedding" in line.lower() or "check" in line.lower() \
                        or "tolerance" in line.lower():
                    count += 1
    return count


def parse_ewald_alpha(log_path):
    """Extract Ewald splitting parameter alpha from GROMACS log."""
    with open(log_path) as f:
        for line in f:
            # Line like: "Using a Gaussian width (1/beta) of ... nm for Ewald"
            m = re.search(r"Ewald\s+(?:parameter|rtol|alpha)\s*[=:]\s*([\d.eE+-]+)", line, re.IGNORECASE)
            if m:
                return float(m.group(1))
            # Alternative: "beta = ..."
            m = re.search(r"Using.*1/beta.*of\s+([\d.eE+-]+)\s*nm", line)
            if m:
                # beta = 1/width, but we want alpha which is usually reported differently
                pass
            # "Ewald: rtol = ..., alpha = ..."
            m = re.search(r"alpha\s*=\s*([\d.eE+-]+)", line)
            if m:
                return float(m.group(1))
    return None


def parse_pme_grid(log_path):
    """Extract PME grid dimensions from GROMACS log."""
    with open(log_path) as f:
        for line in f:
            m = re.search(r"PME\s+(?:mesh|grid)\s*[=:]?\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)", line, re.IGNORECASE)
            if m:
                return (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    return None


def get_n_steps_from_log(log_path):
    """Extract number of steps from GROMACS log."""
    with open(log_path) as f:
        for line in f:
            m = re.search(r"nsteps\s*=\s*(\d+)", line)
            if m:
                return int(m.group(1))
    return None


def analyze_box(box_dir, box_label, box_size_nm, gmx_bin):
    """Analyze a single box size. Returns dict of results."""
    result = {
        "box_size_nm": box_size_nm,
        "n_atoms": None,
        "duration_ps": None,
        "n_steps": None,
        "drift_kj_mol_ps_atom": None,
        "drift_slope_kj_mol_ps": None,
        "drift_r_squared": None,
        "drift_pass": False,
        "qm_energy_mean_kj_mol": None,
        "qm_energy_std_kj_mol": None,
        "qm_energy_min_kj_mol": None,
        "qm_energy_max_kj_mol": None,
        "qm_energy_pass": False,
        "scc_failures": None,
        "consistency_warnings": None,
        "energy_spikes_3sigma": None,
        "ewald_alpha": None,
        "pme_grid": None,
        "completed": False,
        "error": None,
    }

    if not os.path.isdir(box_dir):
        result["error"] = f"Directory not found: {box_dir}"
        return result

    # Find files
    edr_path = find_edr_file(box_dir)
    log_path = find_log_file(box_dir)
    gro_path = find_gro_file(box_dir)

    if edr_path is None:
        result["error"] = "No .edr file found"
        return result
    if log_path is None:
        result["error"] = "No .log file found"
        return result

    print(f"\n--- Analyzing {box_label} ({box_size_nm} nm) ---")
    print(f"  EDR: {os.path.basename(edr_path)}")
    print(f"  LOG: {os.path.basename(log_path)}")

    # Get atom count
    n_atoms = None
    if gro_path:
        try:
            n_atoms = get_n_atoms_from_gro(gro_path)
        except (ValueError, IOError):
            pass
    if n_atoms is None:
        n_atoms = get_n_atoms_from_log(log_path)
    if n_atoms is None:
        result["error"] = "Cannot determine atom count"
        return result
    result["n_atoms"] = n_atoms
    print(f"  Atoms: {n_atoms}")

    # Check completion
    result["completed"] = check_log_completion(log_path)
    print(f"  Completed: {result['completed']}")

    # Parse Ewald parameters
    result["ewald_alpha"] = parse_ewald_alpha(log_path)
    result["pme_grid"] = parse_pme_grid(log_path)
    if result["ewald_alpha"]:
        print(f"  Ewald alpha: {result['ewald_alpha']:.4f}")
    if result["pme_grid"]:
        print(f"  PME grid: {result['pme_grid'][0]}x{result['pme_grid'][1]}x{result['pme_grid'][2]}")

    # Parse nsteps
    result["n_steps"] = get_n_steps_from_log(log_path)

    # SCC failures (AC-6)
    result["scc_failures"] = count_scc_failures(log_path)
    print(f"  SCC failures: {result['scc_failures']}")

    # Consistency warnings (AC-3)
    result["consistency_warnings"] = count_consistency_warnings(log_path)
    print(f"  Consistency warnings: {result['consistency_warnings']}")

    # Extract Total-Energy time series (AC-2)
    print(f"  Extracting Total-Energy from EDR...")
    total_energy_data = extract_energy_xvg(gmx_bin, edr_path, "Total-Energy")
    if total_energy_data and len(total_energy_data) > 1:
        result["duration_ps"] = total_energy_data[-1][0] - total_energy_data[0][0]

        # Linear regression for drift
        slope, intercept, r_sq = linear_regression(total_energy_data)
        drift_per_atom = slope / n_atoms
        result["drift_slope_kj_mol_ps"] = slope
        result["drift_kj_mol_ps_atom"] = drift_per_atom
        result["drift_r_squared"] = r_sq
        result["drift_pass"] = abs(drift_per_atom) < DRIFT_THRESHOLD
        print(f"  Total-Energy frames: {len(total_energy_data)}")
        print(f"  Duration: {result['duration_ps']:.3f} ps")
        print(f"  Drift: {drift_per_atom:.6e} kJ/mol/ps/atom "
              f"({'PASS' if result['drift_pass'] else 'FAIL'})")

        # Energy spike detection (AC-7)
        total_energies = [d[1] for d in total_energy_data]
        result["energy_spikes_3sigma"] = count_energy_spikes(total_energies, SPIKE_SIGMA)
        print(f"  Energy spikes (>3sigma): {result['energy_spikes_3sigma']}")
    else:
        print(f"  WARNING: Could not extract Total-Energy")

    # Extract Quantum-En. time series (AC-5)
    print(f"  Extracting Quantum-En. from EDR...")
    qm_energy_data = extract_energy_xvg(gmx_bin, edr_path, "Quantum-En.")
    if qm_energy_data and len(qm_energy_data) > 0:
        qm_values = [d[1] for d in qm_energy_data]
        stats = compute_statistics(qm_values)
        result["qm_energy_mean_kj_mol"] = stats["mean"]
        result["qm_energy_std_kj_mol"] = stats["std"]
        result["qm_energy_min_kj_mol"] = stats["min"]
        result["qm_energy_max_kj_mol"] = stats["max"]
        result["qm_energy_pass"] = abs(stats["mean"] - QM_ENERGY_REF) < QM_ENERGY_TOL
        print(f"  QM energy: {stats['mean']:.2f} +/- {stats['std']:.2f} kJ/mol "
              f"({'PASS' if result['qm_energy_pass'] else 'FAIL'})")
    else:
        print(f"  WARNING: Could not extract Quantum-En.")

    return result


def evaluate_acceptance_criteria(box_results):
    """Evaluate all 7 acceptance criteria from the per-box results.

    Returns dict of AC results and overall pass/fail.
    """
    ac = {}

    # Collect valid boxes (those with actual results)
    valid = {k: v for k, v in box_results.items() if v.get("error") is None}

    # AC-1: All 3 boxes complete without crash
    all_complete = (len(valid) == 3 and
                    all(v["completed"] for v in valid.values()))
    ac["AC-1"] = {
        "description": "All boxes complete without crash",
        "pass": all_complete,
        "detail": f"{sum(1 for v in valid.values() if v['completed'])}/3 completed",
    }

    # AC-2: NVE drift < 0.01 kJ/mol/ps/atom per box
    drift_results = {}
    all_drift_pass = True
    for label, v in valid.items():
        d = v.get("drift_kj_mol_ps_atom")
        if d is not None:
            drift_results[label] = d
            if not v["drift_pass"]:
                all_drift_pass = False
        else:
            all_drift_pass = False
    ac["AC-2"] = {
        "description": "NVE drift < 0.01 kJ/mol/ps/atom per box",
        "pass": all_drift_pass and len(drift_results) == 3,
        "drift_values": drift_results,
    }

    # AC-3: US-048 consistency check: 0 warnings
    all_consistent = all(
        v.get("consistency_warnings", -1) == 0 for v in valid.values()
    ) and len(valid) == 3
    ac["AC-3"] = {
        "description": "Consistency check: 0 warnings",
        "pass": all_consistent,
        "warnings": {k: v.get("consistency_warnings") for k, v in valid.items()},
    }

    # AC-4: Drift rates comparable: max/min ratio < 10
    drift_ratio = None
    ac4_pass = False
    abs_drifts = [abs(d) for d in drift_results.values() if d is not None]
    if len(abs_drifts) == 3:
        min_drift = min(abs_drifts)
        max_drift = max(abs_drifts)
        if min_drift > 0:
            drift_ratio = max_drift / min_drift
            ac4_pass = drift_ratio < DRIFT_RATIO_THRESHOLD
        else:
            # If minimum drift is essentially zero, check if max is also tiny
            # Use 1e-10 floor to avoid division by zero
            drift_ratio = max_drift / max(min_drift, 1e-10)
            ac4_pass = max_drift < DRIFT_THRESHOLD
    ac["AC-4"] = {
        "description": "Drift ratio < 10",
        "pass": ac4_pass,
        "drift_ratio": drift_ratio,
    }

    # AC-5: Mean QM energy within +/-5 kJ/mol of -10706 kJ/mol
    all_qm_pass = all(
        v.get("qm_energy_pass", False) for v in valid.values()
    ) and len(valid) == 3
    ac["AC-5"] = {
        "description": "QM energy within +/-5 kJ/mol of reference",
        "pass": all_qm_pass,
        "reference_kj_mol": QM_ENERGY_REF,
        "tolerance_kj_mol": QM_ENERGY_TOL,
        "means": {k: v.get("qm_energy_mean_kj_mol") for k, v in valid.items()},
    }

    # AC-6: Zero SCC failures
    all_scc_clean = all(
        v.get("scc_failures", -1) == 0 for v in valid.values()
    ) and len(valid) == 3
    ac["AC-6"] = {
        "description": "Zero SCC failures",
        "pass": all_scc_clean,
        "counts": {k: v.get("scc_failures") for k, v in valid.items()},
    }

    # AC-7: No edge artifacts at 2.5 nm
    box_2p5 = valid.get("2.5nm")
    ac7_pass = False
    spikes_2p5 = None
    if box_2p5:
        spikes_2p5 = box_2p5.get("energy_spikes_3sigma", -1)
        # AC-7 requires: no spikes AND drift pass AND consistency pass for 2.5 nm
        ac7_pass = (
            spikes_2p5 == 0
            and box_2p5.get("drift_pass", False)
            and box_2p5.get("consistency_warnings", -1) == 0
        )
    ac["AC-7"] = {
        "description": "No edge artifacts at 2.5 nm",
        "pass": ac7_pass,
        "energy_spikes_3sigma": spikes_2p5,
    }

    overall = all(v["pass"] for v in ac.values())

    return ac, drift_ratio, overall


def print_summary(box_results, ac, drift_ratio, overall):
    """Print human-readable summary table."""
    print("\n" + "=" * 78)
    print("US-049: PBC Regression Test Results")
    print("=" * 78)

    # Per-box table
    print(f"\n{'Box':>8} | {'N_atoms':>7} | {'Drift (kJ/mol/ps/at)':>20} | "
          f"{'QM E (kJ/mol)':>15} | {'SCC':>3} | {'Cons':>4} | {'Spikes':>6}")
    print("-" * 78)
    for label in sorted(box_results.keys()):
        v = box_results[label]
        if v.get("error"):
            print(f"{label:>8} | {'ERROR':>7} | {v['error']}")
            continue

        n_str = str(v["n_atoms"]) if v["n_atoms"] else "?"
        drift_str = (f"{v['drift_kj_mol_ps_atom']:.2e}"
                     if v["drift_kj_mol_ps_atom"] is not None else "N/A")
        qm_str = (f"{v['qm_energy_mean_kj_mol']:.1f}"
                  if v["qm_energy_mean_kj_mol"] is not None else "N/A")
        scc_str = str(v["scc_failures"]) if v["scc_failures"] is not None else "?"
        cons_str = (str(v["consistency_warnings"])
                    if v["consistency_warnings"] is not None else "?")
        spike_str = (str(v["energy_spikes_3sigma"])
                     if v["energy_spikes_3sigma"] is not None else "?")

        print(f"{label:>8} | {n_str:>7} | {drift_str:>20} | "
              f"{qm_str:>15} | {scc_str:>3} | {cons_str:>4} | {spike_str:>6}")

    if drift_ratio is not None:
        print(f"\nDrift ratio (max/min): {drift_ratio:.2f} (threshold: {DRIFT_RATIO_THRESHOLD})")

    # AC summary
    print(f"\n{'AC':>5} | {'Result':>6} | Description")
    print("-" * 60)
    for ac_id in sorted(ac.keys()):
        v = ac[ac_id]
        status = "PASS" if v["pass"] else "FAIL"
        print(f"{ac_id:>5} | {status:>6} | {v['description']}")

    print(f"\n{'OVERALL':>5} | {'PASS' if overall else 'FAIL':>6} |")
    print("=" * 78)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze PBC regression test results for US-049.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s tests/data/b5_pbc_regression/\n"
            "  %(prog)s tests/data/b5_pbc_regression/ --gmx /usr/local/bin/gmx\n"
            "  %(prog)s tests/data/b5_pbc_regression/ -o custom_results.json\n"
        ),
    )
    parser.add_argument(
        "base_dir",
        help="Base directory containing b5_2p5nm/, b5_3p0nm/, b5_4p0nm/ subdirectories",
    )
    parser.add_argument(
        "--gmx",
        default="external/gromacs/_install_dftb/bin/gmx",
        help="Path to gmx binary (default: external/gromacs/_install_dftb/bin/gmx)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output JSON file (default: <base_dir>/analysis_results.json)",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}",
    )

    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    if not os.path.isdir(base_dir):
        print(f"ERROR: Base directory not found: {base_dir}", file=sys.stderr)
        sys.exit(1)

    gmx_bin = args.gmx
    # Resolve relative gmx path from CWD (not from base_dir)
    if not os.path.isabs(gmx_bin):
        gmx_bin = os.path.abspath(gmx_bin)
    if not os.path.isfile(gmx_bin):
        print(f"WARNING: gmx binary not found at {gmx_bin}", file=sys.stderr)
        print("Energy extraction will fail. Use --gmx to specify correct path.",
              file=sys.stderr)

    output_path = args.output or os.path.join(base_dir, "analysis_results.json")

    print(f"Base directory: {base_dir}")
    print(f"GMX binary: {gmx_bin}")

    # Analyze each box
    box_results = {}
    for subdir, box_nm in sorted(BOX_DIRS.items(), key=lambda x: x[1]):
        box_dir = os.path.join(base_dir, subdir)
        label = f"{box_nm}nm"
        box_results[label] = analyze_box(box_dir, label, box_nm, gmx_bin)

    # Evaluate acceptance criteria
    ac, drift_ratio, overall = evaluate_acceptance_criteria(box_results)

    # Print summary
    print_summary(box_results, ac, drift_ratio, overall)

    # Build output JSON
    output = {
        "story": "US-049",
        "date": str(date.today()),
        "gmx_binary": gmx_bin,
        "base_dir": base_dir,
        "thresholds": {
            "drift_kj_mol_ps_atom": DRIFT_THRESHOLD,
            "drift_ratio": DRIFT_RATIO_THRESHOLD,
            "qm_energy_ref_kj_mol": QM_ENERGY_REF,
            "qm_energy_tol_kj_mol": QM_ENERGY_TOL,
            "spike_sigma": SPIKE_SIGMA,
        },
        "boxes": {},
        "drift_ratio": drift_ratio,
        "acceptance_criteria": {},
        "overall_pass": overall,
    }

    for label, v in box_results.items():
        # Serialize box results (convert tuples to lists for JSON)
        box_out = dict(v)
        if box_out.get("pme_grid"):
            box_out["pme_grid"] = list(box_out["pme_grid"])
        output["boxes"][label] = box_out

    for ac_id, v in sorted(ac.items()):
        output["acceptance_criteria"][ac_id] = v

    # Write JSON
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults written to: {output_path}")

    # Exit with non-zero code if any AC failed
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()
