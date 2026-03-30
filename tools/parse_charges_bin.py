#!/usr/bin/env python3
"""Parse DFTB+ charges.bin binary restart file and display Mulliken charges.

Reads the binary (stream) format written by DFTB+ (sccinit.F90, writeQToFile).
Supports format versions 4-8. Requires species orbital information from either
an HSD input file (auto-detected) or explicit --species argument.

SDD:specs.md:S4 -- tools/ helper scripts

Usage examples:
    # Auto-detect HSD in same directory
    python3 tools/parse_charges_bin.py tests/data/b5/charges.bin

    # Explicit species specification
    python3 tools/parse_charges_bin.py charges.bin --species "O:p H:s"

    # JSON output for scripting
    python3 tools/parse_charges_bin.py charges.bin --format json

    # Per-orbital breakdown
    python3 tools/parse_charges_bin.py charges.bin -v
"""

import argparse
import json
import os
import re
import struct
import sys

VERSION = "0.1.0"

# Angular momentum letter -> l quantum number
ANG_MOM_MAP = {"s": 0, "p": 1, "d": 2, "f": 3}

# nOrb = (l_max + 1)^2 for each max angular momentum
NORB_FROM_LMAX = {0: 1, 1: 4, 2: 9, 3: 16}

# Default valence electron counts (mio-1-1 / 3ob conventions)
VALENCE_ELECTRONS = {
    "H": 1, "C": 4, "N": 5, "O": 6, "F": 7,
    "P": 5, "S": 6, "Cl": 7, "Br": 7, "I": 7,
    "Na": 1, "Mg": 2, "K": 1, "Ca": 2, "Zn": 2,
    "Ti": 4, "Fe": 8,
}

# Orbital shell labels for verbose output
SHELL_LABELS = {
    0: ["s"],
    1: ["s", "py", "pz", "px"],
    2: ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2-y2"],
    3: ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2-y2",
        "fy(3x2-y2)", "fx2+y2+z2", "fyz2", "fz3", "fxz2", "fz(x2-y2)", "fx(x2-3y2)"],
}


# ---------------------------------------------------------------------------
# Binary parsing
# ---------------------------------------------------------------------------

def parse_header(data):
    """Parse charges.bin header. Returns dict with metadata and header_size."""
    if len(data) < 4:
        raise ValueError("File too small to contain a valid header")

    off = 0
    fmt_ver = struct.unpack_from("<i", data, off)[0]
    off += 4

    if fmt_ver not in (4, 5, 6, 7, 8):
        raise ValueError(f"Unsupported format version: {fmt_ver} (expected 4-8)")

    flags = {}
    if fmt_ver in (4, 5):
        vals = struct.unpack_from("<3i", data, off)
        off += 12
        flags["block"] = bool(vals[0])
        flags["iblock"] = bool(vals[1])
        flags["rho"] = bool(vals[2])
        flags["kpoint"] = False
        flags["multipolar"] = False
    elif fmt_ver == 6:
        vals = struct.unpack_from("<4i", data, off)
        off += 16
        flags["block"] = bool(vals[0])
        flags["iblock"] = bool(vals[1])
        flags["rho"] = bool(vals[2])
        flags["kpoint"] = False
        flags["multipolar"] = bool(vals[3])
    else:  # 7, 8
        vals = struct.unpack_from("<5i", data, off)
        off += 20
        flags["block"] = bool(vals[0])
        flags["iblock"] = bool(vals[1])
        flags["rho"] = bool(vals[2])
        flags["kpoint"] = bool(vals[3])
        flags["multipolar"] = bool(vals[4])

    if fmt_ver == 4:
        n_atom = None  # not stored in v4
        n_spin = struct.unpack_from("<i", data, off)[0]
        off += 4
    else:
        n_atom, n_spin = struct.unpack_from("<2i", data, off)
        off += 8

    checksums = list(struct.unpack_from(f"<{n_spin}d", data, off))
    off += n_spin * 8

    return {
        "format_version": fmt_ver,
        "flags": flags,
        "n_atom": n_atom,
        "n_spin": n_spin,
        "checksums": checksums,
        "header_size": off,
    }


def parse_orbital_charges(data, offset, n_atom, n_spin, n_orb_per_atom):
    """Read orbital charge arrays. Returns (charges_by_spin, new_offset).

    charges_by_spin[spin][atom] = list of orbital populations.
    """
    charges = []
    off = offset
    for _ in range(n_spin):
        spin_charges = []
        for iatom in range(n_atom):
            norb = n_orb_per_atom[iatom]
            vals = list(struct.unpack_from(f"<{norb}d", data, off))
            off += norb * 8
            spin_charges.append(vals)
        charges.append(spin_charges)
    return charges, off


def skip_optional_sections(data, offset, header, n_orb_per_atom):
    """Skip over optional data blocks after orbital charges. Returns new offset."""
    off = offset
    n_atom = header["n_atom"]
    n_spin = header["n_spin"]
    flags = header["flags"]

    # Multipoles
    if flags["multipolar"]:
        n_dip, n_quad = struct.unpack_from("<2i", data, off)
        off += 8
        # dipoles: n_spin * n_atom * n_dip doubles
        off += n_spin * n_atom * n_dip * 8
        # quadrupoles: n_spin * n_atom * n_quad doubles
        off += n_spin * n_atom * n_quad * 8

    # Block Mulliken
    if flags["block"]:
        for _ in range(n_spin):
            for iatom in range(n_atom):
                norb = n_orb_per_atom[iatom]
                off += norb * norb * 8

    # Imaginary block Mulliken
    if flags["iblock"]:
        for _ in range(n_spin):
            for iatom in range(n_atom):
                norb = n_orb_per_atom[iatom]
                off += norb * norb * 8

    return off


def validate_checksum(charges_by_spin, checksums, tolerance=1e-8):
    """Validate stored checksum against computed sum of orbital populations."""
    results = []
    for ispin, (spin_charges, stored) in enumerate(zip(charges_by_spin, checksums)):
        computed = sum(sum(orbs) for orbs in spin_charges)
        ok = abs(computed - stored) < tolerance
        results.append({"spin": ispin, "stored": stored, "computed": computed, "ok": ok})
    return results


# ---------------------------------------------------------------------------
# HSD / geometry parsing
# ---------------------------------------------------------------------------

def parse_hsd_max_angular_momentum(hsd_text):
    """Extract MaxAngularMomentum from HSD text. Returns dict like {"O": "p", "H": "s"}."""
    pattern = r"MaxAngularMomentum\s*=?\s*\{([^}]+)\}"
    match = re.search(pattern, hsd_text)
    if not match:
        return None
    block = match.group(1)
    pairs = re.findall(r"(\w+)\s*=\s*\"?([spdf])\"?", block)
    if not pairs:
        return None
    return {species: ang for species, ang in pairs}


def parse_gen_format(text):
    """Parse DFTB+ GenFormat geometry.

    Returns (species_names, atom_species_indices) where indices are 0-based.
    """
    lines = [l.strip() for l in text.strip().splitlines() if l.strip()]
    if len(lines) < 3:
        raise ValueError("GenFormat text too short")

    parts = lines[0].split()
    n_atom = int(parts[0])
    species_names = lines[1].split()
    atom_species = []
    for i in range(2, 2 + n_atom):
        if i >= len(lines):
            raise ValueError(f"GenFormat: expected {n_atom} atoms, got {i - 2}")
        atom_parts = lines[i].split()
        sp_idx = int(atom_parts[1]) - 1  # convert 1-based to 0-based
        atom_species.append(sp_idx)

    return species_names, atom_species


def parse_xyz_format(text):
    """Parse XYZ format geometry.

    Returns (species_names, atom_species_indices) where indices are 0-based.
    """
    lines = [l.strip() for l in text.strip().splitlines() if l.strip()]
    if len(lines) < 3:
        raise ValueError("XYZ text too short")

    n_atom = int(lines[0])
    # line 1 is comment
    species_names = []
    atom_species = []
    for i in range(2, 2 + n_atom):
        if i >= len(lines):
            raise ValueError(f"XYZ: expected {n_atom} atoms, got {i - 2}")
        elem = lines[i].split()[0]
        if elem not in species_names:
            species_names.append(elem)
        atom_species.append(species_names.index(elem))

    return species_names, atom_species


def parse_hsd_geometry(hsd_text, hsd_dir):
    """Extract species list and atom-species mapping from HSD Geometry block.

    Returns (species_names, atom_species_indices) or (None, None) on failure.
    """
    # Try GenFormat
    gen_match = re.search(
        r"Geometry\s*=?\s*GenFormat\s*\{(.*?)\}", hsd_text, re.DOTALL
    )
    if gen_match:
        block = gen_match.group(1).strip()
        # Check for external file reference
        file_ref = re.search(r'<<<\s*"([^"]+)"', block)
        if file_ref:
            gen_path = os.path.join(hsd_dir, file_ref.group(1))
            if os.path.isfile(gen_path):
                with open(gen_path) as f:
                    return parse_gen_format(f.read())
            return None, None
        else:
            return parse_gen_format(block)

    # Try xyzFormat
    xyz_match = re.search(
        r"Geometry\s*=?\s*xyzFormat\s*\{(.*?)\}", hsd_text, re.DOTALL
    )
    if xyz_match:
        block = xyz_match.group(1).strip()
        file_ref = re.search(r'<<<\s*"([^"]+)"', block)
        if file_ref:
            xyz_path = os.path.join(hsd_dir, file_ref.group(1))
            if os.path.isfile(xyz_path):
                with open(xyz_path) as f:
                    return parse_xyz_format(f.read())
            return None, None
        else:
            return parse_xyz_format(block)

    return None, None


def resolve_hsd_file(charges_path):
    """Auto-detect dftb_pin.hsd or dftb_in.hsd in same directory as charges.bin."""
    charges_dir = os.path.dirname(os.path.abspath(charges_path))
    for name in ("dftb_pin.hsd", "dftb_in.hsd"):
        candidate = os.path.join(charges_dir, name)
        if os.path.isfile(candidate):
            return candidate
    return None


def build_orbital_map(species_names, atom_species_indices, max_ang_mom):
    """Build per-atom orbital count list.

    Returns (n_orb_per_atom, atom_labels) where atom_labels[i] is species name.
    """
    n_orb_per_atom = []
    atom_labels = []
    for sp_idx in atom_species_indices:
        sp_name = species_names[sp_idx]
        atom_labels.append(sp_name)
        ang = max_ang_mom.get(sp_name)
        if ang is None:
            raise ValueError(
                f"Species '{sp_name}' not found in MaxAngularMomentum. "
                f"Available: {list(max_ang_mom.keys())}"
            )
        lmax = ANG_MOM_MAP[ang]
        n_orb_per_atom.append(NORB_FROM_LMAX[lmax])
    return n_orb_per_atom, atom_labels


# ---------------------------------------------------------------------------
# Species specification parsing
# ---------------------------------------------------------------------------

def parse_species_arg(spec_str):
    """Parse --species argument like 'O:p H:s C:p'.

    Returns (species_names, max_ang_mom) where species_names is ordered list
    and max_ang_mom maps name -> letter.
    """
    max_ang_mom = {}
    species_names = []
    for token in spec_str.split():
        if ":" not in token:
            raise ValueError(f"Invalid species spec '{token}', expected 'Element:angmom' (e.g. 'O:p')")
        name, ang = token.split(":", 1)
        if ang not in ANG_MOM_MAP:
            raise ValueError(f"Unknown angular momentum '{ang}' for {name}, expected s/p/d/f")
        species_names.append(name)
        max_ang_mom[name] = ang
    return species_names, max_ang_mom


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def get_atom_population(orb_charges):
    """Sum orbital populations for one atom (handles OrbitalEquiv_expand packing)."""
    return sum(orb_charges)


def format_human(header, charges_by_spin, n_orb_per_atom, atom_labels,
                 verbose=False):
    """Format charges as human-readable table."""
    lines = []
    n_atom = header["n_atom"]
    n_spin = header["n_spin"]
    fmt_ver = header["format_version"]
    flags = header["flags"]

    flag_str = " ".join(
        f"{k}={'Yes' if v else 'No'}" for k, v in flags.items()
    )
    lines.append(f"DFTB+ charges.bin -- v{fmt_ver}, {n_atom} atoms, {n_spin} spin")
    lines.append(f"Flags: {flag_str}")
    lines.append("")

    for ispin in range(n_spin):
        if n_spin > 1:
            lines.append(f"--- Spin channel {ispin + 1} ---")

        if atom_labels:
            lines.append(f" {'Atom':>4}  {'Species':>7}  {'nOrb':>4}  {'Population':>10}  {'Mulliken dq':>11}")
            lines.append(f" {'----':>4}  {'-------':>7}  {'----':>4}  {'----------':>10}  {'-----------':>11}")
        else:
            lines.append(f" {'Atom':>4}  {'nOrb':>4}  {'Population':>10}")
            lines.append(f" {'----':>4}  {'----':>4}  {'----------':>10}")

        total_pop = 0.0
        total_ref = 0.0
        has_valence = True

        for iatom in range(n_atom):
            orbs = charges_by_spin[ispin][iatom]
            pop = get_atom_population(orbs)
            total_pop += pop
            norb = n_orb_per_atom[iatom]
            label = atom_labels[iatom] if atom_labels else "?"

            q0 = VALENCE_ELECTRONS.get(label)
            if q0 is not None:
                total_ref += q0
                dq = q0 - pop
                dq_str = f"{dq:>+11.6f}"
            else:
                has_valence = False
                dq_str = f"{'N/A':>11}"

            if atom_labels:
                lines.append(
                    f" {iatom + 1:>4}  {label:>7}  {norb:>4}  {pop:>10.6f}  {dq_str}"
                )
            else:
                lines.append(f" {iatom + 1:>4}  {norb:>4}  {pop:>10.6f}")

            if verbose and norb > 1:
                # Show per-orbital breakdown
                lmax = {1: 0, 4: 1, 9: 2, 16: 3}.get(norb, -1)
                labels = SHELL_LABELS.get(lmax, [f"orb{i}" for i in range(norb)])
                for iorb, val in enumerate(orbs):
                    if iorb < len(labels):
                        orb_label = labels[iorb]
                    else:
                        orb_label = f"orb{iorb}"
                    lines.append(f"        {orb_label:>8}: {val:>12.8f}")

        lines.append("")

        # Checksum
        chk = header["checksums"][ispin]
        chk_ok = abs(chk - total_pop) < 1e-8
        chk_status = "OK" if chk_ok else f"MISMATCH (stored={chk:.8f})"
        lines.append(f"Total electrons: {total_pop:.6f} (checksum {chk_status})")

        if has_valence and atom_labels:
            net_charge = total_ref - total_pop
            lines.append(f"Net charge: {net_charge:+.6f}")

    return "\n".join(lines)


def format_json(header, charges_by_spin, n_orb_per_atom, atom_labels):
    """Format charges as JSON."""
    n_atom = header["n_atom"]
    n_spin = header["n_spin"]

    chk_results = validate_checksum(charges_by_spin, header["checksums"])

    result = {
        "format_version": header["format_version"],
        "n_atom": n_atom,
        "n_spin": n_spin,
        "flags": header["flags"],
        "checksum_stored": header["checksums"],
        "checksum_computed": [r["computed"] for r in chk_results],
        "checksum_valid": all(r["ok"] for r in chk_results),
        "spins": [],
    }

    for ispin in range(n_spin):
        atoms = []
        for iatom in range(n_atom):
            orbs = charges_by_spin[ispin][iatom]
            pop = get_atom_population(orbs)
            label = atom_labels[iatom] if atom_labels else None
            q0 = VALENCE_ELECTRONS.get(label) if label else None

            entry = {
                "index": iatom + 1,
                "n_orbitals": n_orb_per_atom[iatom],
                "population": pop,
                "orbital_populations": orbs,
            }
            if label:
                entry["species"] = label
            if q0 is not None:
                entry["valence_electrons"] = q0
                entry["mulliken_charge"] = q0 - pop
            atoms.append(entry)
        result["spins"].append({"channel": ispin + 1, "atoms": atoms})

    return json.dumps(result, indent=2)


def format_raw(header, charges_by_spin, n_orb_per_atom):
    """Format raw orbital populations."""
    lines = []
    for ispin, spin_charges in enumerate(charges_by_spin):
        if header["n_spin"] > 1:
            lines.append(f"# spin {ispin + 1}")
        for iatom, orbs in enumerate(spin_charges):
            vals = " ".join(f"{v:.10f}" for v in orbs)
            lines.append(f"{iatom + 1} {vals}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Parse DFTB+ charges.bin binary restart file and display Mulliken charges.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s tests/data/b5/charges.bin\n"
            "  %(prog)s charges.bin --species 'O:p H:s'\n"
            "  %(prog)s charges.bin --format json\n"
            "  %(prog)s charges.bin -v\n"
        ),
    )
    parser.add_argument(
        "charges_file", nargs="?", default="./charges.bin",
        help="Path to charges.bin (default: ./charges.bin)",
    )
    parser.add_argument(
        "--hsd", metavar="FILE",
        help="Path to dftb_in.hsd or dftb_pin.hsd (default: auto-detect)",
    )
    parser.add_argument(
        "--species", metavar="SPEC",
        help="Manual species spec, e.g. 'O:p H:s C:p' (overrides HSD)",
    )
    parser.add_argument(
        "--format", choices=["human", "json", "raw"], default="human",
        dest="output_format",
        help="Output format (default: human)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show per-orbital breakdown",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}",
    )

    args = parser.parse_args()

    # Read binary file
    if not os.path.isfile(args.charges_file):
        print(f"Error: file not found: {args.charges_file}", file=sys.stderr)
        sys.exit(1)

    with open(args.charges_file, "rb") as f:
        data = f.read()

    # Parse header
    try:
        header = parse_header(data)
    except ValueError as e:
        print(f"Error parsing header: {e}", file=sys.stderr)
        sys.exit(1)

    n_atom = header["n_atom"]
    if n_atom is None:
        print("Error: format v4 does not store atom count; use --species and provide atom count", file=sys.stderr)
        sys.exit(1)

    # Resolve species information
    species_names = None
    atom_species_indices = None
    max_ang_mom = None
    atom_labels = None
    n_orb_per_atom = None
    hsd_source = None

    if args.species:
        # Manual species specification
        try:
            species_names, max_ang_mom = parse_species_arg(args.species)
        except ValueError as e:
            print(f"Error in --species: {e}", file=sys.stderr)
            sys.exit(1)

        # Try to get atom-species mapping from HSD
        hsd_path = args.hsd or resolve_hsd_file(args.charges_file)
        if hsd_path and os.path.isfile(hsd_path):
            with open(hsd_path) as f:
                hsd_text = f.read()
            hsd_dir = os.path.dirname(os.path.abspath(hsd_path))
            _, atom_species_indices = parse_hsd_geometry(hsd_text, hsd_dir)

        if atom_species_indices is None:
            # Infer: if only one species, all atoms are that species
            if len(species_names) == 1:
                atom_species_indices = [0] * n_atom
            else:
                print(
                    f"Warning: cannot determine atom-species mapping for {n_atom} atoms "
                    f"with {len(species_names)} species. Need HSD file or single-species system.",
                    file=sys.stderr,
                )
                # Fallback: use uniform nOrb based on max across species
                max_norb = max(
                    NORB_FROM_LMAX[ANG_MOM_MAP[max_ang_mom[sp]]]
                    for sp in species_names
                )
                n_orb_per_atom = [max_norb] * n_atom
                atom_labels = [None] * n_atom

        if n_orb_per_atom is None and atom_species_indices is not None:
            try:
                n_orb_per_atom, atom_labels = build_orbital_map(
                    species_names, atom_species_indices, max_ang_mom
                )
            except ValueError as e:
                print(f"Error building orbital map: {e}", file=sys.stderr)
                sys.exit(1)
            hsd_source = "--species argument"
    else:
        # Auto-detect from HSD
        hsd_path = args.hsd or resolve_hsd_file(args.charges_file)
        if hsd_path and os.path.isfile(hsd_path):
            with open(hsd_path) as f:
                hsd_text = f.read()
            hsd_dir = os.path.dirname(os.path.abspath(hsd_path))

            max_ang_mom = parse_hsd_max_angular_momentum(hsd_text)
            if max_ang_mom is None:
                print(f"Warning: could not parse MaxAngularMomentum from {hsd_path}", file=sys.stderr)
            else:
                species_names_geo, atom_species_indices = parse_hsd_geometry(hsd_text, hsd_dir)
                if atom_species_indices is None:
                    print(f"Warning: could not parse geometry from {hsd_path}", file=sys.stderr)
                else:
                    species_names = species_names_geo
                    try:
                        n_orb_per_atom, atom_labels = build_orbital_map(
                            species_names, atom_species_indices, max_ang_mom
                        )
                    except ValueError as e:
                        print(f"Error building orbital map: {e}", file=sys.stderr)
                        sys.exit(1)
                    hsd_source = os.path.basename(hsd_path)

    if n_orb_per_atom is None:
        print(
            "Error: no species information available.\n"
            "Provide --species 'Element:angmom ...' or ensure dftb_pin.hsd / dftb_in.hsd "
            "is in the same directory as charges.bin.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Validate atom count
    if len(n_orb_per_atom) != n_atom:
        print(
            f"Error: atom count mismatch -- charges.bin has {n_atom} atoms, "
            f"but species info gives {len(n_orb_per_atom)} atoms.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Parse orbital charges
    try:
        charges_by_spin, off = parse_orbital_charges(
            data, header["header_size"], n_atom, header["n_spin"], n_orb_per_atom
        )
    except struct.error as e:
        print(f"Error: file truncated while reading orbital charges: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate checksum
    chk_results = validate_checksum(charges_by_spin, header["checksums"])
    for r in chk_results:
        if not r["ok"]:
            print(
                f"Warning: checksum mismatch for spin {r['spin'] + 1}: "
                f"stored={r['stored']:.10f}, computed={r['computed']:.10f}",
                file=sys.stderr,
            )

    # Format output
    if atom_labels is None:
        atom_labels = [None] * n_atom

    if args.output_format == "human":
        if hsd_source:
            print(f"Source: {hsd_source}")
            print()
        print(format_human(header, charges_by_spin, n_orb_per_atom, atom_labels,
                           verbose=args.verbose))
    elif args.output_format == "json":
        print(format_json(header, charges_by_spin, n_orb_per_atom, atom_labels))
    elif args.output_format == "raw":
        print(format_raw(header, charges_by_spin, n_orb_per_atom))


if __name__ == "__main__":
    main()
