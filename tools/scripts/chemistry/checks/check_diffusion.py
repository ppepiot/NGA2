#!/usr/bin/env python3
"""
Validate the mixture-averaged transport data emitted by yaml2nga against Cantera and tabulated
literature values.

The binary diffusion coefficients are evaluated with the SAME precomputed arrays and collision
integral fit that yaml2nga writes into the generated Fortran (build_transport_arrays / omega_D are
imported from yaml2nga, not re-derived here), so a PASS means the shipped data are right.
Pure-species viscosities are checked the same way (mucoeff / omega_mu).

Cantera applies Stockmayer (dipole) corrections for polar species that the Lennard-Jones/Neufeld
model used here does not have, so polar species are compared with a looser tolerance.

Usage:
  python check_diffusion.py kinetics/gri30.yaml
  python check_diffusion.py kinetics/gri30.yaml --temperatures 298.15 300 1000 --pairs H2:N2 O2:N2
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from yaml2nga import MechanismError, build_transport_arrays, load_mechanism, omega_D, omega_mu  # noqa: E402

# Binary D at 298.15 K, 1 atm (m^2/s): dilute-gas tabulations (Fuller-Schettler-Giddings, Perry's, Incropera).
NIST_BINARY_D_298K = {
    ("H2", "N2"): 7.58e-5,
    ("N2", "O2"): 2.21e-5,
    ("H2O", "N2"): 2.42e-5,
    ("CH4", "N2"): 2.24e-5,
}

DEFAULT_PAIRS = [("H2", "N2"), ("O2", "N2"), ("H2O", "N2"), ("CH4", "N2")]
DEFAULT_TEMPERATURES = (298.15, 300.0, 1000.0)
DEFAULT_PRESSURE = 101325.0
DEFAULT_CANTERA_RTOL = 0.035  # non-polar species / pairs
DEFAULT_POLAR_RTOL = 0.30  # pairs or species involving a polar molecule (no Stockmayer correction here)
DEFAULT_NIST_RTOL = 0.08


def binary_diffusion_coeff(tr: dict, i: int, j: int, T: float, P: float) -> float:
    """D_ij (m^2/s) exactly as fcmech_get_invDij evaluates it."""
    tp_term = P / (T * np.sqrt(T))
    inv_d = tp_term * omega_D(T * tr["Ocoeffs"][i, j]) / tr["Dcoeffs"][i, j]
    return 1.0 / inv_d


def species_viscosity(tr: dict, i: int, T: float) -> float:
    """mu_i (Pa s) exactly as fcmech_get_viscosity evaluates it."""
    return tr["mucoeff"][i] * np.sqrt(T) / omega_mu(T * tr["koveps"][i])


def check_diffusion(yaml_path, *, pairs=None, temperatures=DEFAULT_TEMPERATURES, pressure=DEFAULT_PRESSURE,
                    cantera_rtol=DEFAULT_CANTERA_RTOL, polar_rtol=DEFAULT_POLAR_RTOL,
                    nist_rtol=DEFAULT_NIST_RTOL, verbose=True) -> bool:
    """Compare the generated transport data with Cantera (and NIST at 298.15 K). True if all pass."""
    try:
        import cantera as ct
    except ImportError:
        print("Error: Cantera is required for transport checks (pip install cantera).", file=sys.stderr)
        return False

    yaml_path = Path(yaml_path)
    mech = load_mechanism(yaml_path)
    tr = build_transport_arrays(mech)
    if tr is None:
        print("Error: mechanism has species without transport data; nothing to check.", file=sys.stderr)
        return False
    species_names = mech["species_names"]
    name_to_idx = {name: idx for idx, name in enumerate(species_names)}
    pairs = pairs or DEFAULT_PAIRS

    gas = ct.Solution(str(yaml_path), transport_model="mixture-averaged")
    if gas.species_names != species_names:
        print("Error: species order differs between yaml2nga and Cantera.", file=sys.stderr)
        return False
    polar = {nm for nm in species_names if gas.species(nm).transport.dipole != 0.0}

    all_ok = True
    n_checked = 0
    lines: list[str] = []
    missing = sorted({s for p in pairs for s in p if s not in name_to_idx})
    if missing:
        lines.append(f"Note: species not in this mechanism, pairs involving them are skipped: {', '.join(missing)}")

    for T in temperatures:
        gas.TP = T, pressure
        if verbose:
            lines += ["", f"T = {T:.2f} K, P = {pressure:.0f} Pa",
                      f"{'Pair':<12} {'NGA2 D (m2/s)':>14} {'Cantera':>14} {'rel err':>10} {'tol':>6} {'status':>8}"]
        for a, b in pairs:
            if a not in name_to_idx or b not in name_to_idx:
                continue
            n_checked += 1
            i, j = name_to_idx[a], name_to_idx[b]
            d_nga = binary_diffusion_coeff(tr, i, j, T, pressure)
            d_ct = float(gas.binary_diff_coeffs[i, j])
            rel = (d_nga - d_ct) / d_ct
            tol = polar_rtol if (a in polar or b in polar) else cantera_rtol
            ok = abs(rel) <= tol
            all_ok &= ok
            if verbose:
                lines.append(f"  {a + '-' + b:<10} {d_nga:14.6e} {d_ct:14.6e} {rel:+10.2%} {tol:6.0%} {'PASS' if ok else 'FAIL':>8}")

        if verbose:
            lines += ["", f"{'Species':<12} {'NGA2 mu (Pa s)':>14} {'Cantera':>14} {'rel err':>10} {'tol':>6} {'status':>8}"]
        mu_ct = gas.species_viscosities
        for nm in sorted({s for p in pairs for s in p if s in name_to_idx}):
            i = name_to_idx[nm]
            mu_nga = species_viscosity(tr, i, T)
            rel = (mu_nga - mu_ct[i]) / mu_ct[i]
            tol = polar_rtol if nm in polar else cantera_rtol
            ok = abs(rel) <= tol
            all_ok &= ok
            if verbose:
                lines.append(f"  {nm:<10} {mu_nga:14.6e} {mu_ct[i]:14.6e} {rel:+10.2%} {tol:6.0%} {'PASS' if ok else 'FAIL':>8}")

    nist_T = 298.15
    if any(abs(nist_T - T) < 1.0 for T in temperatures):
        if verbose:
            lines += ["", f"Literature tabulation at {nist_T:.2f} K, 1 atm (secondary check)",
                      f"{'Pair':<12} {'NGA2 D (m2/s)':>14} {'NIST':>14} {'rel err':>10} {'status':>8}"]
        for a, b in pairs:
            key = tuple(sorted((a, b)))
            if key not in NIST_BINARY_D_298K or a not in name_to_idx or b not in name_to_idx:
                continue
            d_nga = binary_diffusion_coeff(tr, name_to_idx[a], name_to_idx[b], nist_T, DEFAULT_PRESSURE)
            d_ref = NIST_BINARY_D_298K[key]
            rel = (d_nga - d_ref) / d_ref
            ok = abs(rel) <= nist_rtol
            all_ok &= ok
            if verbose:
                lines.append(f"  {a + '-' + b:<10} {d_nga:14.6e} {d_ref:14.6e} {rel:+10.2%} {'PASS' if ok else 'FAIL':>8}")

    if n_checked == 0:
        lines.append("Error: none of the requested species pairs exists in this mechanism (use --pairs)")
        all_ok = False
    if verbose:
        print("\n".join(lines))
        print("\nTransport check:", "PASS" if all_ok else "FAIL")
    return all_ok


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate yaml2nga transport data against Cantera and NIST references.")
    parser.add_argument("yaml_file", help="Mechanism YAML file (e.g. kinetics/gri30.yaml)")
    parser.add_argument("--temperatures", nargs="+", type=float, default=list(DEFAULT_TEMPERATURES),
                        help="Temperatures in K (default: 298.15 300 1000)")
    parser.add_argument("--pressure", type=float, default=DEFAULT_PRESSURE, help="Pressure in Pa (default: 101325)")
    parser.add_argument("--pairs", nargs="+", default=None, help="Species pairs as A:B (default: H2:N2 O2:N2 H2O:N2 CH4:N2)")
    parser.add_argument("--cantera-rtol", type=float, default=DEFAULT_CANTERA_RTOL,
                        help=f"Relative tolerance vs Cantera for non-polar species (default: {DEFAULT_CANTERA_RTOL})")
    parser.add_argument("--polar-rtol", type=float, default=DEFAULT_POLAR_RTOL,
                        help=f"Relative tolerance vs Cantera when a polar species is involved (default: {DEFAULT_POLAR_RTOL})")
    parser.add_argument("--nist-rtol", type=float, default=DEFAULT_NIST_RTOL,
                        help=f"Relative tolerance vs tabulated values at 298.15 K (default: {DEFAULT_NIST_RTOL})")
    args = parser.parse_args()

    yaml_path = Path(args.yaml_file)
    if not yaml_path.is_file():
        print(f"Error: file not found: {yaml_path}", file=sys.stderr)
        sys.exit(1)
    pairs = None
    if args.pairs:
        try:
            pairs = [tuple(p.split(":", 1)) for p in args.pairs]
            assert all(len(p) == 2 for p in pairs)
        except (AssertionError, ValueError):
            print("Error: --pairs entries must look like A:B", file=sys.stderr)
            sys.exit(1)

    try:
        ok = check_diffusion(yaml_path, pairs=pairs, temperatures=tuple(args.temperatures), pressure=args.pressure,
                             cantera_rtol=args.cantera_rtol, polar_rtol=args.polar_rtol, nist_rtol=args.nist_rtol)
    except MechanismError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
