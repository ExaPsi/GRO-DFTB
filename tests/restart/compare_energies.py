#!/usr/bin/env python3
"""
US-052 Restart Validation: Energy Comparison Script

Compares QM energies between a continuous and a restarted GROMACS QM/MM
trajectory to verify that checkpoint/restart produces bit-identical (or
near-identical within SCC convergence bounds) results.

Theory references:
  ThDD:Eq. T-US-052-N-1.1 -- trajectory comparison protocol
  ThDD:Eq. T-US-052-N-2.1 -- max absolute energy deviation
  ThDD:Eq. T-US-052-U-1.2 -- kJ/mol to Ha conversion
  ThDD:Eq. T-US-052-V-1.1 -- 1e-6 Ha tolerance (cutoff)
  ThDD:Eq. T-US-052-V-1.2 -- 1e-6 Ha tolerance (PME)
  ThDD:Eq. T-US-052-V-2.2 -- energy discontinuity at restart

LGPL-3.0-or-later
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile


# ThDD:Eq. T-US-052-U-1.2 -- conversion factor kJ/mol -> Hartree
# 1 Ha = 2625.5 kJ/mol (CODATA 2022)
KJ_MOL_PER_HARTREE = 2625.5

# ThDD:Eq. T-US-052-V-1.1 -- tolerance in Hartree
TOLERANCE_HA = 1.0e-6


def run_gmx_energy(gmx_binary, edr_file, term, output_xvg):
    """Extract an energy term from .edr to .xvg using gmx energy.

    Args:
        gmx_binary: Path to gmx executable.
        edr_file: Path to .edr file.
        term: Energy term name (e.g., "Quantum-En." or "Total-Energy").
        output_xvg: Path to write .xvg output.
    """
    cmd = [gmx_binary, "energy", "-f", edr_file, "-o", output_xvg, "-dp"]
    proc = subprocess.run(
        cmd,
        input=f"{term}\n",
        capture_output=True,
        text=True,
        timeout=60,
    )
    if proc.returncode != 0:
        print(f"ERROR: gmx energy failed for {edr_file}", file=sys.stderr)
        print(f"stdout: {proc.stdout}", file=sys.stderr)
        print(f"stderr: {proc.stderr}", file=sys.stderr)
        sys.exit(1)


def parse_xvg(xvg_file):
    """Parse an XVG file, returning list of (time, value) tuples.

    Skips lines starting with # or @.
    """
    data = []
    with open(xvg_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                t = float(parts[0])
                v = float(parts[1])
                data.append((t, v))
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Compare QM energies between continuous and restarted trajectories"
    )
    parser.add_argument(
        "--gmx", required=True, help="Path to gmx binary"
    )
    parser.add_argument(
        "--cont-edr", required=True, help="Continuous run .edr file"
    )
    parser.add_argument(
        "--restart-edr", required=True, help="Restarted run .edr file"
    )
    parser.add_argument(
        "--restart-step", type=int, default=21,
        help="First post-restart step (default: 21)"
    )
    parser.add_argument(
        "--output", required=True, help="Output JSON file path"
    )
    parser.add_argument(
        "--mode", choices=["cutoff", "pme"], required=True,
        help="Embedding mode (cutoff or pme)"
    )
    args = parser.parse_args()

    # Validate inputs
    for path in [args.gmx, args.cont_edr, args.restart_edr]:
        if not os.path.exists(path):
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract Quantum-En. from both .edr files
        cont_qm_xvg = os.path.join(tmpdir, "cont_qm.xvg")
        restart_qm_xvg = os.path.join(tmpdir, "restart_qm.xvg")

        print("Extracting Quantum-En. from continuous run...")
        run_gmx_energy(args.gmx, args.cont_edr, "Quantum-En.", cont_qm_xvg)

        print("Extracting Quantum-En. from restarted run...")
        run_gmx_energy(args.gmx, args.restart_edr, "Quantum-En.", restart_qm_xvg)

        cont_qm = parse_xvg(cont_qm_xvg)
        restart_qm = parse_xvg(restart_qm_xvg)

        # Extract Total-Energy for discontinuity check (AC-5)
        cont_etot_xvg = os.path.join(tmpdir, "cont_etot.xvg")
        restart_etot_xvg = os.path.join(tmpdir, "restart_etot.xvg")

        print("Extracting Total-Energy from continuous run...")
        run_gmx_energy(args.gmx, args.cont_edr, "Total-Energy", cont_etot_xvg)

        print("Extracting Total-Energy from restarted run...")
        run_gmx_energy(args.gmx, args.restart_edr, "Total-Energy", restart_etot_xvg)

        cont_etot = parse_xvg(cont_etot_xvg)
        restart_etot = parse_xvg(restart_etot_xvg)

    # Build time-indexed dictionaries for continuous data
    cont_qm_dict = {round(t, 6): v for t, v in cont_qm}
    cont_etot_dict = {round(t, 6): v for t, v in cont_etot}

    # Compute dt from continuous data
    if len(cont_qm) < 2:
        print("ERROR: Continuous run has fewer than 2 energy data points", file=sys.stderr)
        sys.exit(1)
    dt = cont_qm[1][0] - cont_qm[0][0]  # ps

    # Restart time = restart_step * dt (step is 0-indexed in GROMACS)
    # The first post-restart step is args.restart_step
    restart_time = round((args.restart_step - 1) * dt, 6)  # time of restart point (step 20)

    print(f"\ndt = {dt} ps")
    print(f"Restart point time = {restart_time} ps (step {args.restart_step - 1})")
    print(f"First post-restart step = {args.restart_step}")
    print(f"Continuous run: {len(cont_qm)} QM energy data points")
    print(f"Restarted run: {len(restart_qm)} QM energy data points")

    # ThDD:Eq. T-US-052-N-1.1 -- trajectory comparison
    # Compare QM energies at each post-restart step
    deltas_ha = []
    deltas_detail = []

    for t_restart, e_restart in restart_qm:
        t_key = round(t_restart, 6)
        # Only compare post-restart steps
        if t_key < round(args.restart_step * dt, 6) - 1e-9:
            continue

        if t_key not in cont_qm_dict:
            print(f"WARNING: Time {t_key} ps not found in continuous run", file=sys.stderr)
            continue

        e_cont = cont_qm_dict[t_key]

        # ThDD:Eq. T-US-052-U-1.2 -- kJ/mol to Hartree
        delta_kjmol = e_restart - e_cont
        delta_ha = delta_kjmol / KJ_MOL_PER_HARTREE

        deltas_ha.append(delta_ha)
        deltas_detail.append({
            "time_ps": t_key,
            "step": round(t_key / dt),
            "E_cont_kJmol": e_cont,
            "E_restart_kJmol": e_restart,
            "deltaE_kJmol": delta_kjmol,
            "deltaE_Ha": delta_ha,
        })

    if not deltas_ha:
        print("ERROR: No overlapping post-restart data points found", file=sys.stderr)
        sys.exit(1)

    # ThDD:Eq. T-US-052-N-2.1 -- max absolute energy deviation
    abs_deltas = [abs(d) for d in deltas_ha]
    max_delta_ha = max(abs_deltas)
    max_delta_step = deltas_detail[abs_deltas.index(max_delta_ha)]["step"]

    # ThDD:Eq. T-US-052-N-2.2 -- RMS energy deviation
    rms_delta_ha = math.sqrt(sum(d**2 for d in deltas_ha) / len(deltas_ha))

    # First post-restart step comparison
    first_delta_ha = abs(deltas_ha[0]) if deltas_ha else float("nan")
    first_step = deltas_detail[0]["step"] if deltas_detail else -1

    # ThDD:Eq. T-US-052-V-2.2 -- energy discontinuity at restart point
    # Check total energy jump at restart vs normal step-to-step fluctuations
    etot_jump = None
    etot_sigma = None
    etot_jump_pass = None

    if cont_etot:
        # Compute step-to-step total energy differences for pre-restart region
        pre_restart_diffs = []
        for i in range(1, len(cont_etot)):
            t = round(cont_etot[i][0], 6)
            if t <= restart_time:
                diff = cont_etot[i][1] - cont_etot[i - 1][1]
                pre_restart_diffs.append(diff)

        if pre_restart_diffs:
            etot_sigma = math.sqrt(
                sum(d**2 for d in pre_restart_diffs) / len(pre_restart_diffs)
            )

            # Find total energy at restart step and step after in restarted run
            restart_etot_dict = {round(t, 6): v for t, v in restart_etot}
            t_restart_point = round(restart_time, 6)
            t_next = round(restart_time + dt, 6)

            if t_restart_point in restart_etot_dict and t_next in restart_etot_dict:
                etot_jump = restart_etot_dict[t_next] - restart_etot_dict[t_restart_point]
                # ThDD:Eq. T-US-052-V-2.2 -- within 3-sigma
                if etot_sigma > 0:
                    etot_jump_pass = abs(etot_jump) <= 3.0 * etot_sigma
                else:
                    etot_jump_pass = True  # sigma = 0 means no fluctuation

    # ThDD:Eq. T-US-052-V-1.1 / V-1.2 -- tolerance check
    tolerance_pass = max_delta_ha < TOLERANCE_HA

    # Print results
    print("\n" + "=" * 70)
    print(f"US-052 Restart Validation: {args.mode.upper()} Mode")
    print("=" * 70)
    print(f"  Compared steps: {len(deltas_ha)} post-restart steps")
    print(f"  First post-restart step: {first_step}")
    print(f"  |deltaE| at first step: {first_delta_ha:.2e} Ha ({first_delta_ha * KJ_MOL_PER_HARTREE:.2e} kJ/mol)")
    print(f"  max |deltaE|: {max_delta_ha:.2e} Ha (at step {max_delta_step})")
    print(f"  RMS deltaE: {rms_delta_ha:.2e} Ha")
    print(f"  Tolerance: {TOLERANCE_HA:.0e} Ha")
    print(f"  PASS: {'YES' if tolerance_pass else 'NO'}")

    if etot_jump is not None:
        print(f"\n  Energy discontinuity check:")
        print(f"    Etot jump at restart: {etot_jump:.4f} kJ/mol")
        print(f"    Normal sigma(dEtot): {etot_sigma:.4f} kJ/mol")
        if etot_sigma > 0:
            print(f"    Jump / sigma: {abs(etot_jump) / etot_sigma:.2f}")
        print(f"    Within 3-sigma: {'YES' if etot_jump_pass else 'NO'}")
    print("=" * 70)

    # Build results dictionary
    results = {
        "mode": args.mode,
        "cont_edr": os.path.abspath(args.cont_edr),
        "restart_edr": os.path.abspath(args.restart_edr),
        "restart_step": args.restart_step,
        "n_compared_steps": len(deltas_ha),
        "first_postrestart_step": first_step,
        "first_deltaE_Ha": first_delta_ha,
        "max_deltaE_Ha": max_delta_ha,
        "max_deltaE_step": int(max_delta_step),
        "rms_deltaE_Ha": rms_delta_ha,
        "tolerance_Ha": TOLERANCE_HA,
        "tolerance_pass": tolerance_pass,
        "etot_jump_kJmol": etot_jump,
        "etot_sigma_kJmol": etot_sigma,
        "etot_jump_within_3sigma": etot_jump_pass,
        "per_step_details": deltas_detail,
        "theory_refs": {
            "comparison_protocol": "ThDD:Eq. T-US-052-N-1.1",
            "max_deviation": "ThDD:Eq. T-US-052-N-2.1",
            "rms_deviation": "ThDD:Eq. T-US-052-N-2.2",
            "unit_conversion": "ThDD:Eq. T-US-052-U-1.2",
            "tolerance_cutoff": "ThDD:Eq. T-US-052-V-1.1",
            "tolerance_pme": "ThDD:Eq. T-US-052-V-1.2",
            "discontinuity": "ThDD:Eq. T-US-052-V-2.2",
        },
    }

    # Write JSON output
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to: {args.output}")

    # Exit code
    if tolerance_pass:
        print("\nPASS: Restart validation successful")
        sys.exit(0)
    else:
        print("\nFAIL: Energy deviation exceeds tolerance")
        sys.exit(1)


if __name__ == "__main__":
    main()
