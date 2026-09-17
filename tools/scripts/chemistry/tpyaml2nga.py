#!/usr/bin/env python3
"""
Generate Fortran module from a two-phase (gas + liquid) thermo YAML.

Input YAML has no reactions. Liquid-phase species are identified by names
starting with ``L`` (e.g. LH2O); all others are treated as gas-phase.

Output (default ``chem_data_tp.f90``) includes:
  - Species order: all gas species first, then all liquid species (YAML order within each group)
  - nS, nG, nL, nP=2, iGmin, iGmax, iLmin, iLmax
  - thermo_data(nS, 15): col1 = T_switch (K), cols 2--8 = low-T NASA7 a1..a7, cols 9--15 = high-T
    (a copy of the low-T row for species with a single temperature range)
  - Element matrix elem_mat(nS, nA) — species index first, atom index second
  - phse_mat(nS, nP): column 1 = gas (1 if species is gas), column 2 = liquid

NASA7 polynomials (same as Cantera): Cp/R, H/(RT), S/R from the seven coefficients
below/above T_switch (standard NASA piecewise formulas).

Usage:
  python tpyaml2nga.py <mechanism.yaml> [output.f90]
  python tpyaml2nga.py kinetics/tpmech_sample.yaml src/chem_data_tp.f90
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# Shared Fortran emitters / YAML loader (YAML 1.2 scalars, chunked constructors, 132-column check)
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))
from yaml2nga import (  # noqa: E402
    MechanismError,
    canonical_element,
    check_fortran_limits,
    f90_2d_array_parameter as _f90_2d_array_parameter,
    f90_array_parameter as _f90_array_parameter,
    f90_int,
    f90_real,
    load_yaml,
    species_composition,
    species_thermo,
)


def f90_array_parameter(name, values, max_per_line=4):
    return _f90_array_parameter(name, values, per_line=max_per_line)


def _fort_dim_expr(fort_dim, nrows, ncols):
    """Dimension/shape expression for Fortran: use parameter names if fort_dim given, else numeric."""
    if fort_dim:
        return ", ".join(p.strip() for p in fort_dim.split(","))
    return f"{nrows}, {ncols}"


def f90_real_2d_parameter(name, arr, nrows, ncols, max_per_line=4, fort_dim=None):
    """Fortran real parameter matrix, column-major; chunked so no statement exceeds the continuation limit."""
    return _f90_2d_array_parameter(name, arr, nrows, ncols, per_line=max_per_line,
                                   dim_expr=_fort_dim_expr(fort_dim, nrows, ncols))


def f90_int_2d_parameter(name, arr, nrows, ncols, max_per_line=8, fort_dim=None):
    """Integer 2D parameter matrix, column-major; chunked like the real version."""
    return _f90_2d_array_parameter(name, arr, nrows, ncols, kind="integer", fmt=f90_int, per_line=max_per_line,
                                   dim_expr=_fort_dim_expr(fort_dim, nrows, ncols))


def canon_species_name(n):
    return str(n)


def is_liquid_name(name: str) -> bool:
    """Liquid species: name starts with 'L' (e.g. LH2O)."""
    return len(name) > 0 and name[0] == "L"


def extract_tp_thermo(spec: dict) -> tuple[float, list[float], list[float]]:
    """
    Return (T_switch, coeffs_low[7], coeffs_high[7]) using yaml2nga's validated NASA7 reader (exactly
    7 coefficients per row, increasing temperature bounds, cp > 0 over the range, T_mid discontinuity
    warning). A single-range species repeats its row above T_switch = T_high.
    """
    low, high, _, T_mid, _ = species_thermo(spec, str(spec.get("name")))
    return T_mid, list(low), list(high)


def fortran_safe_symbol(sym: str) -> str:
    s = re.sub(r"[^0-9A-Za-z_]", "_", str(sym))
    if s and s[0].isdigit():
        s = "e_" + s
    return s or "X"


def load_tp_species(path: Path):
    data = load_yaml(path)

    raw = data.get("species", [])
    species_dicts = [s for s in raw if isinstance(s, dict)]

    gas_sp = []
    liq_sp = []
    for s in species_dicts:
        nm = canon_species_name(s.get("name", "X"))
        if is_liquid_name(nm):
            liq_sp.append(s)
        else:
            gas_sp.append(s)

    ordered = gas_sp + liq_sp
    names = [canon_species_name(s.get("name", "X")) for s in ordered]
    nG = len(gas_sp)
    nL = len(liq_sp)
    nS = len(ordered)

    # Atoms: sorted union of the (validated, canonicalized) composition keys; every species needs one
    compositions = [species_composition(s, nm) for s, nm in zip(ordered, names)]
    atom_names = sorted({e for comp in compositions for e in comp}, key=lambda e: canonical_element(e))
    nA = len(atom_names)
    atom_index = {a: j for j, a in enumerate(atom_names)}

    elem_mat = np.zeros((nS, nA), dtype=int)
    for i, comp in enumerate(compositions):
        for el, cnt in comp.items():
            elem_mat[i, atom_index[el]] = cnt

    nP = 2
    phse_mat = np.zeros((nS, nP), dtype=int)
    for i in range(nS):
        if i < nG:
            phse_mat[i, 0] = 1  # gas
        else:
            phse_mat[i, 1] = 1  # liquid

    iGmin = 1
    iGmax = nG
    iLmin = nG + 1
    iLmax = nS

    T_switch = np.zeros(nS)
    thermo_low = np.zeros((nS, 7))
    thermo_high = np.zeros((nS, 7))
    for i, s in enumerate(ordered):
        ts, lo, hi = extract_tp_thermo(s)
        T_switch[i] = ts
        thermo_low[i, :] = lo
        thermo_high[i, :] = hi

    return {
        "names": names,
        "nS": nS,
        "nG": nG,
        "nL": nL,
        "nA": nA,
        "nP": nP,
        "iGmin": iGmin,
        "iGmax": iGmax,
        "iLmin": iLmin,
        "iLmax": iLmax,
        "atom_names": atom_names,
        "elem_mat": elem_mat,
        "phse_mat": phse_mat,
        "T_switch": T_switch,
        "thermo_low": thermo_low,
        "thermo_high": thermo_high,
        "yaml_source": path.name,
    }


def emit_fortran(mech: dict, out_path: Path) -> None:
    names = mech["names"]
    nS, nG, nL, nA, nP = mech["nS"], mech["nG"], mech["nL"], mech["nA"], mech["nP"]

    lines = []
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append(f"!  FILE {out_path.name}")
    lines.append("!  Two-phase thermodynamic data (gas + liquid species)")
    lines.append("!  Generated by tpyaml2nga.py. Do not edit manually.")
    lines.append(f"!  Source: {mech['yaml_source']}")
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append("!  Thermodynamics: NASA7 piecewise polynomials (same convention as Cantera / GRI).")
    lines.append("!    For species s: T_switch = thermo_data(s,1); low coeffs thermo_data(s,2:8); high thermo_data(s,9:15).")
    lines.append("!    For T <= T_switch(s): use columns 2..8; for T > T_switch(s): use columns 9..15 (Cantera convention).")
    lines.append("!    Coefficients a1..a7 (per row) define (R = universal gas constant):")
    lines.append("!      Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4")
    lines.append("!      H/(R*T) = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T")
    lines.append("!      S/R = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7")
    lines.append("!    Species with a single temperature range repeat the low-T row in cols 9..15 (T_switch = T_high).")
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append("!  Layout summary:")
    lines.append("!    nS, nG, nL, nA, nP     : counts (total species, gas, liquid, atom types, phases=2).")
    lines.append("!    iGmin..iGmax           : Fortran species index range for gas (always iGmin=1).")
    lines.append("!    iLmin..iLmax           : Fortran species index range for liquid (iLmin = iGmax+1).")
    lines.append("!    s<Name>                : named alias for species index 1..nS (gas then liquid).")
    lines.append("!    thermo_data(s,1:15)    : col1 = T_switch (K); cols 2--8 = low NASA7 a1..a7; cols 9--15 = high a1..a7.")
    lines.append("!    elem_mat(s,a)          : stoichiometric count of atom a in species s (column order below).")
    lines.append("!    phse_mat(s,p)          : p=1 gas, p=2 liquid; 1 if species s exists in that phase, else 0.")
    lines.append("!--------------------------------------------------------------------------------------------------")
    atom_order = ", ".join(mech["atom_names"])
    lines.append(f"!  Atom column order for elem_mat(:,a): a=1..nA => {atom_order}")
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append("")
    lines.append("module chem_two_phase")
    lines.append("  use precision")
    lines.append("  implicit none")
    lines.append("")
    lines.append("  ! --- Species and phase counts (gas species listed first in all species-indexed arrays) ---")
    lines.append("  integer, parameter :: nS = " + str(nS) + "  ! total number of species (nG + nL)")
    lines.append("  integer, parameter :: nG = " + str(nG) + "  ! number of gas-phase species (indices 1..nG)")
    lines.append("  integer, parameter :: nL = " + str(nL) + "  ! number of liquid-phase species (indices nG+1..nS)")
    lines.append("  integer, parameter :: nA = " + str(nA) + "  ! number of chemical elements (atom types) in elem_mat")
    lines.append("  integer, parameter :: nP = " + str(nP) + "  ! number of phases represented in phse_mat (gas, liquid)")
    lines.append("  integer, parameter :: nThermoDataCols = 1 + 7 + 7  ! thermo_data: T_switch + low NASA7 + high NASA7")
    lines.append("")
    lines.append("  ! --- Index ranges for querying species by phase (1-based Fortran species index) ---")
    lines.append("  integer, parameter :: iGmin = " + str(mech["iGmin"]) + "  ! first gas species index (always 1)")
    lines.append("  integer, parameter :: iGmax = " + str(mech["iGmax"]) + "  ! last gas species index (= nG)")
    lines.append("  integer, parameter :: iLmin = " + str(mech["iLmin"]) + "  ! first liquid species index (= nG+1)")
    lines.append(
        "  integer, parameter :: iLmax = "
        + str(mech["iLmax"])
        + "  ! last liquid species index (= nS)"
    )
    lines.append("")

    # Species index parameters
    lines.append("  ! --- Named species indices: use s<Name> as shorthand for species index 1..nS ---")
    for i, nm in enumerate(names):
        sym = fortran_safe_symbol(nm)
        lines.append(f"  integer, parameter :: s{sym} = {i + 1}  ! species '{nm}'")
    lines.append("")

    # Max species name length for tp_species_name result (no species_name_* parameters)
    max_len = max(len(nm) for nm in names) if names else 1
    max_len = max(max_len, 1)

    thermo_data = np.zeros((nS, 15), dtype=float)
    thermo_data[:, 0] = mech["T_switch"]
    thermo_data[:, 1:8] = mech["thermo_low"]
    thermo_data[:, 8:15] = mech["thermo_high"]
    lines.append(
        "  ! --- thermo_data(s,j): j=1 T_switch (K); j=2..8 low-T NASA7 a1..a7; j=9..15 high-T a1..a7 (see file header) ---"
    )
    lines.append(
        "  !     dimension(nS, nThermoDataCols) with nThermoDataCols = 1+7+7; high-T cols = low-T row for single-range species ---"
    )
    lines.extend(
        f90_real_2d_parameter("thermo_data", thermo_data, nS, 15, fort_dim="nS, nThermoDataCols")
    )
    lines.append("")
    lines.append(
        "  ! --- Elemental composition: elem_mat(s,a) = count of atom type a in species s ---"
    )
    lines.append(
        "  !     Column a matches atom order: "
        + atom_order
        + " (see file header line 'Atom column order')."
    )
    lines.extend(
        f90_int_2d_parameter("elem_mat", mech["elem_mat"], nS, nA, fort_dim="nS, nA")
    )
    lines.append("")
    lines.append(
        "  ! --- Phase membership: phse_mat(s,1)=1 if gas species s, phse_mat(s,2)=1 if liquid species s ---"
    )
    lines.append("  !     Each species is in exactly one phase here (gas XOR liquid).")
    lines.extend(
        f90_int_2d_parameter("phse_mat", mech["phse_mat"], nS, nP, fort_dim="nS, nP")
    )
    lines.append("")
    lines.append("contains")
    lines.append("")
    lines.append("  ! Return species name for index i (1..nS); '?' if out of range.")
    lines.append("  pure function tp_species_name(i) result(nm)")
    lines.append("    integer, intent(in) :: i  ! species index")
    lines.append(f"    character(len={max_len}) :: nm")
    lines.append("    select case (i)")
    for i, nm in enumerate(names):
        lines.append(f"    case ({i + 1}); nm = '{nm}'")
    lines.append("    case default; nm = '?'")
    lines.append("    end select")
    lines.append("  end function tp_species_name")
    lines.append("")
    atom_len = max(len(an) for an in mech["atom_names"]) if mech["atom_names"] else 3
    atom_len = max(atom_len, 3)
    lines.append(f"  ! Return chemical symbol for atom column index a (1..nA); '?' if out of range.")
    lines.append("  pure function tp_atom_name(a) result(nm)")
    lines.append("    integer, intent(in) :: a  ! atom column index, same as elem_mat(:,a)")
    lines.append(f"    character(len={atom_len}) :: nm")
    lines.append("    select case (a)")
    for j, an in enumerate(mech["atom_names"]):
        # pad to len 3 for Fortran char
        q = an.replace("'", "''")
        lines.append(f"    case ({j + 1}); nm = '{q}'")
    lines.append("    case default; nm = '?'")
    lines.append("    end select")
    lines.append("  end function tp_atom_name")
    lines.append("")
    lines.append("end module chem_two_phase")
    lines.append("")

    check_fortran_limits(lines)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate chem_data_tp.f90 from two-phase thermo YAML")
    parser.add_argument("input", type=Path, help="Input YAML (species only, no reactions)")
    parser.add_argument("output", nargs="?", type=Path, default=Path("chem_data_tp.f90"), help="Output Fortran file")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: input not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    try:
        mech = load_tp_species(args.input)
        print(f"Species: nS={mech['nS']} (nG={mech['nG']}, nL={mech['nL']}), nA={mech['nA']}")
        print("Order:", ", ".join(mech["names"]))
        emit_fortran(mech, args.output)
    except MechanismError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
