#!/usr/bin/env python3
"""Analyze nstlist sweep experiment for PME drift attribution.

Determines whether charge-lag (nstlist-synchronized QM charge insertion)
is the dominant source of NVE energy drift in PME-embedded QM/MM.

Theory: QM Mulliken charges are inserted into GROMACS chargeA[] only at
pair-search steps (every nstlist steps). Between insertions, PME uses stale
charges. The charge-lag error scales linearly with nstlist:

    drift_charge_lag(nstlist) = a * nstlist

where a is the per-step charge-lag contribution. A linear fit of measured
drift vs nstlist identifies:
  - slope a: charge-lag contribution per nstlist unit
  - intercept b: non-charge-lag drift (integration error, rounding, etc.)

If R^2 of the linear fit is high and a > 0, charge-lag is the dominant source.

SDD:specs.md:S4 -- tools/ helper scripts
ThDD:T-US-050-3.1 -- PME/cutoff drift comparison
ThDD:T-US-047-1.2 -- Charge-lag from nstlist-synchronized insertion

Usage:
    python3 tools/analyze_nstlist_sweep.py docs/session/20260212/nstlist_sweep/ \\
        --gmx external/gromacs/_install_dftb/bin/gmx \\
        -o docs/session/20260212/nstlist_sweep/analysis_results.json \\
        --plot docs/Manuscript/M1/Manuscript/figures/nstlist_drift.tex
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

# System parameters (B5: QM water in TIP3P box, 3.0 nm)
N_ATOMS_B5 = 2652

# Cutoff drift reference (Golden Rule: from real calculation)
# Provenance: tests/data/b5/drift_analysis_damped_500ps.json, US-043b
# total_energy.drift_per_atom_kj_mol_ps = -4.532983071206687e-05
CUTOFF_DRIFT_PER_ATOM = -4.53e-05  # kJ/mol/ps/atom

# NVE drift threshold (reconciled to 0.005 per GROMACS standard, see C-4)
DRIFT_THRESHOLD = 0.005  # kJ/mol/ps/atom

# Expected nstlist subdirectory pattern
NSTLIST_DIR_PATTERN = re.compile(r"nstlist_(\d+)")


def find_nstlist_dirs(sweep_dir):
    """Find all nstlist_XX subdirectories and return sorted (nstlist, path) pairs."""
    results = []
    for entry in sorted(os.listdir(sweep_dir)):
        m = NSTLIST_DIR_PATTERN.match(entry)
        if m:
            nstlist = int(m.group(1))
            path = os.path.join(sweep_dir, entry)
            if os.path.isdir(path):
                results.append((nstlist, path))
    return sorted(results, key=lambda x: x[0])


def extract_energy_xvg(gmx_bin, edr_path, energy_term):
    """Extract energy time series from .edr using gmx energy.

    Returns list of (time_ps, energy_kj_mol) tuples, or None on failure.
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

    sum_t = sum(d[0] for d in data)
    sum_e = sum(d[1] for d in data)
    sum_tt = sum(d[0] * d[0] for d in data)
    sum_te = sum(d[0] * d[1] for d in data)

    mean_t = sum_t / n
    mean_e = sum_e / n

    s_tt = sum_tt - n * mean_t * mean_t
    s_te = sum_te - n * mean_t * mean_e

    if abs(s_tt) < 1e-30:
        return {
            "slope": 0.0, "intercept": mean_e, "r_squared": 0.0,
            "stderr": float("inf"), "n_points": n,
        }

    slope = s_te / s_tt
    intercept = mean_e - slope * mean_t

    ss_tot = sum((d[1] - mean_e) ** 2 for d in data)
    ss_res = sum((d[1] - (slope * d[0] + intercept)) ** 2 for d in data)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

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

    abs_alpha1 = abs(first_half["slope"])
    abs_alpha2 = abs(second_half["slope"])

    if abs_alpha1 > 1e-30:
        ratio = abs_alpha2 / abs_alpha1
    else:
        ratio = 0.0 if abs_alpha2 < 1e-30 else float("inf")

    both_tiny = abs_alpha1 < 1e-8 and abs_alpha2 < 1e-8

    return {
        "first_half": first_half,
        "second_half": second_half,
        "ratio": ratio,
        "pass": ratio <= 3.0 or both_tiny,
    }


def check_scc_convergence(log_path):
    """Parse md.log for SCC non-convergence events.

    Returns dict with failure count.
    """
    count = 0
    if log_path and os.path.isfile(log_path):
        with open(log_path) as f:
            for line in f:
                if "SCC is NOT converged" in line:
                    count += 1
    return {"failures": count, "pass": count == 0}


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


def get_n_atoms(data_dir, log_path):
    """Get atom count from GRO file or log file."""
    # Try .gro file first
    for name in ["initial.gro", "confout.gro", "md.gro"]:
        gro_path = os.path.join(data_dir, name)
        if os.path.isfile(gro_path):
            # Resolve symlinks
            real_path = os.path.realpath(gro_path)
            if os.path.isfile(real_path):
                try:
                    with open(real_path) as f:
                        f.readline()  # title
                        return int(f.readline().strip())
                except (ValueError, IOError):
                    pass

    # Fall back to log file
    if log_path and os.path.isfile(log_path):
        with open(log_path) as f:
            for line in f:
                m = re.search(r"number of atoms\s*=\s*(\d+)", line, re.IGNORECASE)
                if m:
                    return int(m.group(1))
    return None


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


def linear_fit(x_vals, y_vals, y_errs=None):
    """Weighted or unweighted linear regression: y = a*x + b.

    If y_errs provided and non-zero, uses weighted least squares.
    Returns dict with slope (a), intercept (b), r_squared, slope_stderr, intercept_stderr.
    """
    n = len(x_vals)
    if n < 2:
        return {
            "slope": 0.0, "intercept": 0.0, "r_squared": 0.0,
            "slope_stderr": float("inf"), "intercept_stderr": float("inf"),
            "n_points": n,
        }

    # Use unweighted fit (more robust for small N)
    sum_x = sum(x_vals)
    sum_y = sum(y_vals)
    sum_xx = sum(x * x for x in x_vals)
    sum_xy = sum(x * y for x, y in zip(x_vals, y_vals))

    mean_x = sum_x / n
    mean_y = sum_y / n

    s_xx = sum_xx - n * mean_x * mean_x
    s_xy = sum_xy - n * mean_x * mean_y

    if abs(s_xx) < 1e-30:
        return {
            "slope": 0.0, "intercept": mean_y, "r_squared": 0.0,
            "slope_stderr": float("inf"), "intercept_stderr": float("inf"),
            "n_points": n,
        }

    slope = s_xy / s_xx
    intercept = mean_y - slope * mean_x

    ss_tot = sum((y - mean_y) ** 2 for y in y_vals)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(x_vals, y_vals))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    if n > 2 and s_xx > 0:
        mse = ss_res / (n - 2)
        slope_stderr = math.sqrt(mse / s_xx)
        intercept_stderr = math.sqrt(mse * sum_xx / (n * s_xx))
    else:
        slope_stderr = float("inf")
        intercept_stderr = float("inf")

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "slope_stderr": slope_stderr,
        "intercept_stderr": intercept_stderr,
        "n_points": n,
    }


def analyze_single_run(nstlist, run_dir, gmx_bin):
    """Analyze a single nstlist run directory.

    Returns dict with drift results for this nstlist value.
    """
    result = {
        "nstlist": nstlist,
        "run_dir": run_dir,
    }

    # Find files
    edr_path = os.path.join(run_dir, "md.edr")
    log_path = os.path.join(run_dir, "md.log")

    if not os.path.isfile(edr_path):
        result["error"] = f"No md.edr found in {run_dir}"
        return result
    if not os.path.isfile(log_path):
        result["error"] = f"No md.log found in {run_dir}"
        return result

    # Check completion
    result["completed"] = check_log_completion(log_path)

    # Get atom count
    n_atoms = get_n_atoms(run_dir, log_path)
    if n_atoms is None:
        n_atoms = N_ATOMS_B5  # Fallback for B5 system
    result["n_atoms"] = n_atoms

    # SCC failures
    scc = check_scc_convergence(log_path)
    result["scc_failures"] = scc["failures"]

    # Extract total energy
    total_energy_data = extract_energy_xvg(gmx_bin, edr_path, "Total-Energy")
    if not total_energy_data or len(total_energy_data) < 3:
        result["error"] = "Insufficient total energy data"
        return result

    duration_ps = total_energy_data[-1][0] - total_energy_data[0][0]
    result["duration_ps"] = duration_ps
    result["n_frames"] = len(total_energy_data)

    # Total energy statistics
    energies = [d[1] for d in total_energy_data]
    stats = compute_statistics(energies)
    result["total_energy_mean"] = stats["mean"]
    result["total_energy_std"] = stats["std"]

    # Drift via linear regression
    drift = compute_drift_with_stderr(total_energy_data)
    drift_per_atom = drift["slope"] / n_atoms
    drift_stderr_per_atom = drift["stderr"] / n_atoms

    result["drift_kj_mol_ps"] = drift["slope"]
    result["drift_per_atom"] = drift_per_atom
    result["drift_stderr_per_atom"] = drift_stderr_per_atom
    result["r_squared"] = drift["r_squared"]
    result["drift_pass"] = abs(drift_per_atom) < DRIFT_THRESHOLD

    # Block analysis
    block = block_analysis(total_energy_data)
    result["block_ratio"] = block["ratio"]
    result["block_pass"] = block["pass"]

    # Extract QM energy
    qm_data = extract_energy_xvg(gmx_bin, edr_path, "Quantum-En.")
    if qm_data:
        qm_vals = [d[1] for d in qm_data]
        qm_stats = compute_statistics(qm_vals)
        result["qm_energy_mean"] = qm_stats["mean"]
        result["qm_energy_std"] = qm_stats["std"]

    # Extract temperature
    temp_data = extract_energy_xvg(gmx_bin, edr_path, "Temperature")
    if temp_data:
        temp_vals = [d[1] for d in temp_data]
        temp_stats = compute_statistics(temp_vals)
        result["temperature_mean"] = temp_stats["mean"]
        result["temperature_std"] = temp_stats["std"]

    return result


def generate_pgfplots_tex(plot_path, nstlist_values, drift_values, drift_errors,
                          fit_result, cutoff_drift):
    """Generate publication-quality pgfplots .tex file.

    ThDD:T-US-047-1.2 -- Charge-lag from nstlist-synchronized insertion
    """
    # Build data points string
    data_points = "\n".join(
        f"        {n} {d:.8e} {e:.8e}"
        for n, d, e in zip(nstlist_values, drift_values, drift_errors)
    )

    # Linear fit line: two endpoints
    x_min = 0
    x_max = max(nstlist_values) + 2
    y_fit_min = fit_result["slope"] * x_min + fit_result["intercept"]
    y_fit_max = fit_result["slope"] * x_max + fit_result["intercept"]

    # Format fit equation
    slope_sci = f"{fit_result['slope']:.2e}"
    intercept_sci = f"{fit_result['intercept']:.2e}"
    r2_str = f"{fit_result['r_squared']:.3f}"

    # Cutoff drift value (absolute)
    cutoff_abs = abs(cutoff_drift)

    tex_content = r"""\begin{tikzpicture}
\begin{axis}[
    width=0.85\columnwidth,
    height=0.65\columnwidth,
    xlabel={$n_\mathrm{stlist}$},
    ylabel={Drift rate (kJ\,mol$^{-1}$\,ps$^{-1}$\,atom$^{-1}$)},
    xmin=0, xmax=""" + str(x_max) + r""",
    ymin=0,
    grid=major,
    grid style={gray!30},
    legend pos=north west,
    legend style={font=\footnotesize, draw=none, fill=white, fill opacity=0.8},
    tick label style={font=\footnotesize},
    label style={font=\small},
    every axis plot/.append style={thick},
]

% Data points with error bars
\addplot[
    only marks,
    mark=*,
    mark size=2.5pt,
    color=blue!80!black,
    error bars/.cd,
    y dir=both,
    y explicit,
] table [x=nstlist, y=drift, y error=err] {
    nstlist drift err
""" + data_points + r"""
};
\addlegendentry{PME NVE drift}

% Linear fit
\addplot[
    no marks,
    color=blue!80!black,
    dashed,
    domain=0:""" + str(x_max) + r""",
] {""" + f"{fit_result['slope']:.10e}" + r"""*x + """ + f"{fit_result['intercept']:.10e}" + r"""};
\addlegendentry{Linear fit ($R^2=""" + r2_str + r"""$)}

% Cutoff drift reference (horizontal line)
\addplot[
    no marks,
    color=red!70!black,
    densely dotted,
    thick,
] coordinates {(0, """ + f"{cutoff_abs:.8e}" + r""") (""" + str(x_max) + r""", """ + f"{cutoff_abs:.8e}" + r""")};
\addlegendentry{Cutoff drift (US-043b)}

% Fit annotation
\node[anchor=south west, font=\scriptsize, text=blue!80!black] at (rel axis cs:0.05,0.70) {%
    $\alpha = """ + slope_sci + r""" \times n_\mathrm{stlist} + """ + intercept_sci + r"""$};

\end{axis}
\end{tikzpicture}
"""
    # Ensure parent directory exists
    plot_dir = os.path.dirname(plot_path)
    if plot_dir and not os.path.isdir(plot_dir):
        os.makedirs(plot_dir, exist_ok=True)

    with open(plot_path, "w") as f:
        f.write(tex_content)

    print(f"\nPlot written to: {plot_path}")


def sanitize_for_json(obj):
    """Handle inf/nan for JSON serialization."""
    if isinstance(obj, float):
        if math.isinf(obj) or math.isnan(obj):
            return str(obj)
        return obj
    elif isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [sanitize_for_json(v) for v in obj]
    return obj


def main():
    parser = argparse.ArgumentParser(
        description="Analyze nstlist sweep experiment for PME drift attribution.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s docs/session/20260212/nstlist_sweep/\n"
            "  %(prog)s docs/session/20260212/nstlist_sweep/ --gmx /path/to/gmx\n"
            "  %(prog)s docs/session/20260212/nstlist_sweep/ -o results.json "
            "--plot figures/nstlist_drift.tex\n"
        ),
    )
    parser.add_argument(
        "sweep_dir",
        help="Directory containing nstlist_XX subdirectories",
    )
    parser.add_argument(
        "--gmx",
        default="external/gromacs/_install_dftb/bin/gmx",
        help="Path to gmx binary (default: external/gromacs/_install_dftb/bin/gmx)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output JSON file (default: <sweep_dir>/analysis_results.json)",
    )
    parser.add_argument(
        "--plot",
        default=None,
        help="Output pgfplots .tex file for drift vs nstlist figure",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}",
    )

    args = parser.parse_args()

    sweep_dir = os.path.abspath(args.sweep_dir)
    if not os.path.isdir(sweep_dir):
        print(f"ERROR: Sweep directory not found: {sweep_dir}", file=sys.stderr)
        sys.exit(1)

    gmx_bin = args.gmx
    if not os.path.isabs(gmx_bin):
        gmx_bin = os.path.abspath(gmx_bin)
    if not os.path.isfile(gmx_bin):
        print(f"WARNING: gmx binary not found at {gmx_bin}", file=sys.stderr)

    output_path = args.output or os.path.join(sweep_dir, "analysis_results.json")

    # Find nstlist directories
    nstlist_dirs = find_nstlist_dirs(sweep_dir)
    if not nstlist_dirs:
        print(f"ERROR: No nstlist_XX directories found in {sweep_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"nstlist sweep directory: {sweep_dir}")
    print(f"Found {len(nstlist_dirs)} nstlist values: "
          f"{', '.join(str(n) for n, _ in nstlist_dirs)}")
    print()

    # Analyze each run
    run_results = []
    for nstlist, run_dir in nstlist_dirs:
        print(f"{'='*60}")
        print(f"Analyzing nstlist = {nstlist}  ({run_dir})")
        print(f"{'='*60}")

        result = analyze_single_run(nstlist, run_dir, gmx_bin)
        run_results.append(result)

        if "error" in result:
            print(f"  ERROR: {result['error']}")
        else:
            print(f"  Duration: {result.get('duration_ps', 0):.3f} ps "
                  f"({result.get('n_frames', 0)} frames)")
            print(f"  Completed: {result.get('completed', False)}")
            print(f"  SCC failures: {result.get('scc_failures', '?')}")
            drift = result.get("drift_per_atom", 0)
            stderr = result.get("drift_stderr_per_atom", 0)
            print(f"  Drift: {drift:.6e} +/- {stderr:.2e} kJ/mol/ps/atom")
            print(f"  R^2: {result.get('r_squared', 0):.6f}")
            print(f"  Block ratio: {result.get('block_ratio', '?')}")
            print(f"  Drift pass (<{DRIFT_THRESHOLD}): {result.get('drift_pass', False)}")

            if "temperature_mean" in result:
                print(f"  Temperature: {result['temperature_mean']:.1f} +/- "
                      f"{result['temperature_std']:.1f} K")
            if "qm_energy_mean" in result:
                print(f"  QM energy: {result['qm_energy_mean']:.2f} +/- "
                      f"{result['qm_energy_std']:.2f} kJ/mol")
        print()

    # Collect valid results for linear fit
    valid = [r for r in run_results if "error" not in r and r.get("drift_per_atom") is not None]

    if len(valid) < 2:
        print("ERROR: Fewer than 2 valid runs — cannot fit linear model.", file=sys.stderr)
        results = {
            "experiment": "nstlist_sweep",
            "date": str(date.today()),
            "sweep_dir": sweep_dir,
            "runs": run_results,
            "error": "Insufficient valid runs for linear fit",
        }
        with open(output_path, "w") as f:
            json.dump(sanitize_for_json(results), f, indent=2)
        print(f"\nResults written to: {output_path}")
        sys.exit(1)

    nstlist_values = [r["nstlist"] for r in valid]
    drift_values = [abs(r["drift_per_atom"]) for r in valid]
    drift_errors = [r["drift_stderr_per_atom"] for r in valid]

    # Linear fit: |drift| = a * nstlist + b
    fit = linear_fit(nstlist_values, drift_values, drift_errors)

    print(f"\n{'='*60}")
    print("LINEAR FIT: |drift| = a * nstlist + b")
    print(f"{'='*60}")
    print(f"  a (slope) = {fit['slope']:.6e} +/- {fit['slope_stderr']:.2e} "
          f"kJ/mol/ps/atom per nstlist unit")
    print(f"  b (intercept) = {fit['intercept']:.6e} +/- {fit['intercept_stderr']:.2e} "
          f"kJ/mol/ps/atom")
    print(f"  R^2 = {fit['r_squared']:.6f}")
    print(f"  N points = {fit['n_points']}")

    # Interpretation
    if fit["r_squared"] > 0.9 and fit["slope"] > 0:
        interpretation = "Strong linear scaling confirms charge-lag as dominant drift source"
    elif fit["r_squared"] > 0.7 and fit["slope"] > 0:
        interpretation = "Moderate linear scaling suggests charge-lag is a significant drift source"
    elif fit["slope"] <= 0:
        interpretation = "Non-positive slope — charge-lag does NOT explain drift"
    else:
        interpretation = "Weak correlation — charge-lag may not be the dominant drift source"
    print(f"  Interpretation: {interpretation}")

    # Charge-lag fraction at nstlist=10 (our production value)
    if fit["intercept"] > 0 and 10 in nstlist_values:
        idx_10 = nstlist_values.index(10)
        total_drift_10 = drift_values[idx_10]
        charge_lag_10 = fit["slope"] * 10
        non_lag_10 = fit["intercept"]
        fraction_lag = charge_lag_10 / (charge_lag_10 + non_lag_10) if (charge_lag_10 + non_lag_10) > 0 else 0
        print(f"\n  At nstlist=10 (production):")
        print(f"    Measured |drift|:   {total_drift_10:.6e} kJ/mol/ps/atom")
        print(f"    Charge-lag part:    {charge_lag_10:.6e} ({fraction_lag*100:.1f}%)")
        print(f"    Non-lag part:       {non_lag_10:.6e} ({(1-fraction_lag)*100:.1f}%)")

    # Compare with cutoff reference
    cutoff_abs = abs(CUTOFF_DRIFT_PER_ATOM)
    print(f"\n  Cutoff drift reference: {CUTOFF_DRIFT_PER_ATOM:.6e} kJ/mol/ps/atom")
    print(f"    (from tests/data/b5/drift_analysis_damped_500ps.json, US-043b)")

    for r in valid:
        ratio = abs(r["drift_per_atom"]) / cutoff_abs if cutoff_abs > 0 else float("inf")
        print(f"    nstlist={r['nstlist']:2d}: |drift|/|cutoff_drift| = {ratio:.2f}")

    # nstlist=1 is the irreducible baseline (no charge lag)
    nstlist_1_result = next((r for r in valid if r["nstlist"] == 1), None)
    if nstlist_1_result:
        baseline = abs(nstlist_1_result["drift_per_atom"])
        print(f"\n  nstlist=1 baseline (no charge lag): {baseline:.6e} kJ/mol/ps/atom")
        for r in valid:
            if r["nstlist"] > 1:
                excess = abs(r["drift_per_atom"]) - baseline
                print(f"    nstlist={r['nstlist']:2d}: excess over baseline = {excess:.6e}")

    # Build full results dict
    results = {
        "experiment": "nstlist_sweep",
        "purpose": "Determine if charge-lag is the dominant PME drift source",
        "date": str(date.today()),
        "sweep_dir": sweep_dir,
        "gmx_binary": gmx_bin,
        "system": "B5: 1 QM water in TIP3P box (3.0 nm)",
        "n_atoms": valid[0]["n_atoms"] if valid else N_ATOMS_B5,
        "cutoff_drift_reference": {
            "value_kj_mol_ps_atom": CUTOFF_DRIFT_PER_ATOM,
            "provenance": "tests/data/b5/drift_analysis_damped_500ps.json (US-043b)",
        },
        "runs": run_results,
        "linear_fit": {
            "model": "|drift| = a * nstlist + b",
            "slope_a": fit["slope"],
            "slope_a_stderr": fit["slope_stderr"],
            "intercept_b": fit["intercept"],
            "intercept_b_stderr": fit["intercept_stderr"],
            "r_squared": fit["r_squared"],
            "n_points": fit["n_points"],
            "interpretation": interpretation,
        },
        "nstlist_values": nstlist_values,
        "drift_values_abs": drift_values,
        "drift_errors": drift_errors,
    }

    # Generate pgfplots .tex if requested
    if args.plot:
        plot_path = args.plot
        if not os.path.isabs(plot_path):
            plot_path = os.path.abspath(plot_path)
        generate_pgfplots_tex(
            plot_path, nstlist_values, drift_values, drift_errors,
            fit, CUTOFF_DRIFT_PER_ATOM,
        )
        results["plot_file"] = plot_path

    # Write JSON
    with open(output_path, "w") as f:
        json.dump(sanitize_for_json(results), f, indent=2)
    print(f"\nResults written to: {output_path}")

    # Summary table
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'nstlist':>8} | {'Drift (kJ/mol/ps/atom)':>22} | {'R^2':>8} | {'SCC':>4} | {'Pass':>4}")
    print(f"{'-'*60}")
    for r in run_results:
        if "error" in r:
            print(f"{r['nstlist']:>8} | {'ERROR':>22} | {'--':>8} | {'--':>4} | {'--':>4}")
        else:
            d = r.get("drift_per_atom", 0)
            e = r.get("drift_stderr_per_atom", 0)
            r2 = r.get("r_squared", 0)
            scc = r.get("scc_failures", 0)
            p = "PASS" if r.get("drift_pass", False) else "FAIL"
            print(f"{r['nstlist']:>8} | {d:>12.4e} +/- {e:.1e} | {r2:>8.4f} | {scc:>4} | {p:>4}")
    print(f"{'-'*60}")
    print(f"{'Fit':>8} | a={fit['slope']:.4e}, b={fit['intercept']:.4e}, R^2={fit['r_squared']:.4f}")
    print(f"{'Cutoff':>8} | {CUTOFF_DRIFT_PER_ATOM:>12.4e} (US-043b reference)")
    print(f"{'='*60}")

    # Exit code: 0 if all runs pass drift threshold
    all_pass = all(r.get("drift_pass", False) for r in run_results if "error" not in r)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
