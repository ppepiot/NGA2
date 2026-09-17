#!/usr/bin/env python3
"""
Generate chem_data_fc.f90 (Fortran module `fcmech`) from a Cantera-format YAML mechanism.

Usage: python yaml2nga.py <mechanism.yaml> [output.f90]
       Default output: chem_data_fc.f90 (a reaction index file <stem>_reactions.txt is written next to it)

Supported reaction types
  - elementary (Arrhenius), including `negative-A: true`
  - three-body (`+ M`, or `(+X)` with a specific collider), with `efficiencies` / `default-efficiency`
  - falloff (Lindemann / Troe)
  - pressure-dependent-Arrhenius (PLOG; duplicate pressures are summed, exact table hits are honored)
Reverse rates of reversible reactions are, by default, evaluated at run time from detailed balance,
k_r = k_f / K_c(T), with K_c from the NASA7 table stored in the module (exactly consistent with the
thermodynamics). With --fit-reverse, k_f/K_c is instead fitted to an Arrhenius form over a temperature
range and the reverse rates are evaluated exactly like forward ones (cheaper: no Gibbs energies at run
time; the fit error per reaction is reported and written to the reaction index file).

Everything that is unsupported, ambiguous, or inconsistent (unknown units, unknown reaction types or
keys, unresolved species, unbalanced elements, non-physical thermo, ...) raises a MechanismError instead
of being silently defaulted. Mixture-averaged transport data are emitted only if every species has
transport data; otherwise the transport routines are omitted with a warning.

Generated Fortran is kept within the F2008 free-form limits (132 columns, <=255 continuation lines).
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print("Error: PyYAML required. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Error: NumPy required. Install with: pip install numpy", file=sys.stderr)
    sys.exit(1)


class MechanismError(ValueError):
    """Raised for any mechanism content that cannot be converted faithfully."""


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------------------------------
R_UNIV = 8.314462618  # J/(mol K)
K_BOLTZMANN = 1.380649e-23  # J/K
N_AVOGADRO = 6.02214076e23  # 1/mol
ATM_TO_PA = 101325.0
P_REF = ATM_TO_PA  # standard-state pressure of NASA7 thermo data (1 atm, as in Cantera)

# Standard atomic weights (kg/mol): IUPAC 2021 conventional values, as used by Cantera.
ATOMIC_MASSES_KG = {
    "H": 1.008e-3, "He": 4.002602e-3, "Li": 6.94e-3, "Be": 9.0121831e-3, "B": 10.81e-3,
    "C": 12.011e-3, "N": 14.007e-3, "O": 15.999e-3, "F": 18.998403163e-3, "Ne": 20.1797e-3,
    "Na": 22.98976928e-3, "Mg": 24.305e-3, "Al": 26.9815384e-3, "Si": 28.085e-3,
    "P": 30.973761998e-3, "S": 32.06e-3, "Cl": 35.45e-3, "Ar": 39.95e-3, "K": 39.0983e-3,
    "Ca": 40.078e-3, "Ti": 47.867e-3, "Cr": 51.9961e-3, "Mn": 54.938043e-3, "Fe": 55.845e-3,
    "Ni": 58.6934e-3, "Cu": 63.546e-3, "Zn": 65.38e-3, "Br": 79.904e-3, "Kr": 83.798e-3,
    "I": 126.90447e-3, "Xe": 131.293e-3,
    "D": 2.0141017781e-3, "Tr": 3.0160492820e-3,  # deuterium, tritium (Cantera symbols)
    "E": 5.48579909065e-7,  # electron
}
_ELEMENT_BY_UPPER = {k.upper(): k for k in ATOMIC_MASSES_KG}


def canonical_element(symbol) -> str:
    """Canonical element symbol (e.g. 'AR' -> 'Ar'); raises for unknown elements."""
    key = str(symbol).strip().upper()
    if key not in _ELEMENT_BY_UPPER:
        raise MechanismError(
            f"Unknown element '{symbol}'. Add it to ATOMIC_MASSES_KG in yaml2nga.py "
            f"(known: {', '.join(ATOMIC_MASSES_KG)})"
        )
    return _ELEMENT_BY_UPPER[key]


# Neufeld et al. (1972) rational fits of the Lennard-Jones collision integrals (same coefficients in
# the emitted Fortran and in the Python helpers below, so check scripts test the shipped formula).
NEUFELD_OMEGA_MU = [3.3530622607, 2.53272006, 2.9024238575, 0.11186138893, 0.8662326188,
                    1.3913958626, 3.158490576, 0.18973411754, 0.00018682962894]
NEUFELD_OMEGA_D = [6.8728271691, 9.4122316321, 7.7442359037, 0.23424661229, 1.45337701568,
                   5.2269794238, 9.7108519575, 0.46539437353, 0.00041908394781]


def _neufeld_rational(m, T_star):
    dr = m[8]
    for k in range(1, 5):
        dr = m[8 - k] + T_star * dr
    nr = m[3]
    for k in range(1, 4):
        nr = m[3 - k] + T_star * nr
    return nr / dr


def omega_mu(T_star: float) -> float:
    """Reduced collision integral Omega^(2,2)*(T*) for viscosity (Neufeld fit)."""
    return _neufeld_rational(NEUFELD_OMEGA_MU, T_star)


def omega_D(T_star: float) -> float:
    """Reduced collision integral Omega^(1,1)*(T*) for diffusion (Neufeld fit)."""
    return _neufeld_rational(NEUFELD_OMEGA_D, T_star)


# ---------------------------------------------------------------------------------------------------
# Units (Cantera conventions). Values convert TO SI: m, s, mol, J/mol, Pa.
# ---------------------------------------------------------------------------------------------------
CANTERA_DEFAULT_UNITS = {
    "length": "m", "time": "s", "quantity": "kmol", "energy": "J", "pressure": "Pa",
}
_LENGTH_UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "dm": 1e-1, "km": 1e3}
_ENERGY_UNITS = {"J": 1.0, "kJ": 1e3, "cal": 4.184, "kcal": 4184.0, "eV": 1.602176634e-19, "erg": 1e-7}
_TIME_UNITS = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "min": 60.0, "hr": 3600.0}
_QUANTITY_UNITS = {"mol": 1.0, "kmol": 1e3, "molec": 1.0 / N_AVOGADRO, "molecule": 1.0 / N_AVOGADRO}
_EA_UNITS = {
    "J/mol": 1.0, "kJ/mol": 1e3, "J/kmol": 1e-3, "kJ/kmol": 1.0,
    "cal/mol": 4.184, "kcal/mol": 4184.0, "cal/kmol": 4.184e-3, "kcal/kmol": 4.184,
    "K": R_UNIV, "eV": 96485.33212,
}
_PRESSURE_UNITS = {
    "Pa": 1.0, "kPa": 1e3, "MPa": 1e6, "bar": 1e5, "mbar": 1e2, "atm": ATM_TO_PA,
    "atmosphere": ATM_TO_PA, "torr": ATM_TO_PA / 760.0, "mmHg": ATM_TO_PA / 760.0,
    "psi": 6894.757293168, "dyn/cm^2": 0.1, "dyn/cm2": 0.1,
}


def _lookup_unit(table: dict, unit, what: str) -> float:
    key = str(unit).strip()
    if key in table:
        return table[key]
    lower = {k.lower(): v for k, v in table.items()}
    if key.lower() in lower:
        return lower[key.lower()]
    raise MechanismError(f"Unsupported {what} unit '{unit}' (supported: {', '.join(table)})")


def parse_units(data: dict) -> dict:
    """Resolve the mechanism-level `units:` block to SI conversion factors."""
    block = data.get("units")
    if block is None:
        print("  No 'units' block: using Cantera defaults (m, s, kmol, J/kmol, Pa)")
        block = {}
    if not isinstance(block, dict):
        raise MechanismError("'units' must be a mapping")
    known = set(CANTERA_DEFAULT_UNITS) | {"mass", "activation-energy"}  # 'mass' affects nothing used here
    unknown = set(block) - known
    if unknown:
        raise MechanismError(f"Unknown key(s) in 'units' block: {sorted(unknown)}")
    u = dict(CANTERA_DEFAULT_UNITS)
    u.update({k: v for k, v in block.items() if k in CANTERA_DEFAULT_UNITS})
    quantity_mol = _lookup_unit(_QUANTITY_UNITS, u["quantity"], "quantity")
    if "activation-energy" in block:
        u["activation-energy"] = block["activation-energy"]
        Ea_J_per_mol = _lookup_unit(_EA_UNITS, u["activation-energy"], "activation-energy")
    else:
        # Cantera: activation energies default to energy/quantity (J/kmol with the default units)
        u["activation-energy"] = f"{u['energy']}/{u['quantity']}"
        Ea_J_per_mol = _lookup_unit(_ENERGY_UNITS, u["energy"], "energy") / quantity_mol
    return {
        "length_m": _lookup_unit(_LENGTH_UNITS, u["length"], "length"),
        "time_s": _lookup_unit(_TIME_UNITS, u["time"], "time"),
        "quantity_mol": quantity_mol,
        "Ea_J_per_mol": Ea_J_per_mol,
        "pressure_Pa": _lookup_unit(_PRESSURE_UNITS, u["pressure"], "pressure"),
        "names": u,
    }


_NUMBER_UNIT_RE = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*([^\s].*?)?\s*$")


def value_with_unit(val, table: dict, default_factor: float, what: str) -> float:
    """Number (mechanism units) or 'number unit' string -> SI."""
    if isinstance(val, bool) or val is None:
        raise MechanismError(f"{what}: missing or invalid value '{val}'")
    if isinstance(val, (int, float)):
        return float(val) * default_factor
    m = _NUMBER_UNIT_RE.match(str(val))
    if not m:
        raise MechanismError(f"{what}: cannot parse value '{val}'")
    num = float(m.group(1))
    unit = m.group(2)
    if not unit:
        return num * default_factor
    return num * _lookup_unit(table, unit, what)


def parse_pressure(val, units: dict) -> float:
    """PLOG pressure entry -> Pa (bare numbers use the mechanism pressure unit)."""
    return value_with_unit(val, _PRESSURE_UNITS, units["pressure_Pa"], "pressure")


def A_to_si(A: float, order: float, units: dict) -> float:
    """Pre-exponential factor from mechanism units to SI ((m^3/mol)^(order-1)/s)."""
    return A * (units["length_m"] ** 3 / units["quantity_mol"]) ** (order - 1.0) / units["time_s"]


# ---------------------------------------------------------------------------------------------------
# YAML loading. PyYAML implements YAML 1.1 (NO/yes/on/off are booleans, 012 is octal, 12:30 is
# sexagesimal, `1e5` and `3.5e13` are strings); Cantera files are YAML 1.2. Replace the bool/int/float
# implicit resolvers with the YAML 1.2 core schema so species named NO keep their name and every
# numeric spelling Cantera accepts parses as a number.
# ---------------------------------------------------------------------------------------------------
class _Yaml12Loader(yaml.SafeLoader):
    pass


_YAML12_REPLACED = {"tag:yaml.org,2002:bool", "tag:yaml.org,2002:int", "tag:yaml.org,2002:float"}
_Yaml12Loader.yaml_implicit_resolvers = {
    ch: [(tag, rx) for tag, rx in resolvers if tag not in _YAML12_REPLACED]
    for ch, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
_Yaml12Loader.add_implicit_resolver(
    "tag:yaml.org,2002:bool", re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$"), list("tTfF")
)
_Yaml12Loader.add_implicit_resolver(
    "tag:yaml.org,2002:int", re.compile(r"^(?:[-+]?[0-9]+|0o[0-7]+|0x[0-9a-fA-F]+)$"), list("-+0123456789")
)
_Yaml12Loader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(r"^(?:[-+]?(?:\.[0-9]+|[0-9]+(?:\.[0-9]*)?)(?:[eE][-+]?[0-9]+)?|[-+]?\.(?:inf|Inf|INF)|\.(?:nan|NaN|NAN))$"),
    list("-+.0123456789"),
)


def _construct_int12(loader, node):
    """YAML 1.2 integers: decimal (a leading 0 is NOT octal), 0o octal, 0x hexadecimal."""
    v = loader.construct_scalar(node)
    if v[:2] in ("0o", "0O"):
        return int(v[2:], 8)
    if v[:2] in ("0x", "0X"):
        return int(v[2:], 16)
    return int(v, 10)


_Yaml12Loader.add_constructor("tag:yaml.org,2002:int", _construct_int12)


def load_yaml(path) -> dict:
    with open(path) as f:
        data = yaml.load(f, Loader=_Yaml12Loader)
    if not isinstance(data, dict):
        raise MechanismError(f"{path}: top level must be a mapping")
    return data


# ---------------------------------------------------------------------------------------------------
# Fortran formatting helpers
# ---------------------------------------------------------------------------------------------------
F90_MAX_LINE = 132  # free-form line length limit (F2008 3.3.2.1)
F90_MAX_CONTINUATION = 255  # continuation lines per statement (F2008 3.3.2.6)
_MAX_CONSTRUCTOR_VALS = 152  # values per array constructor (4/line -> 38 continuation lines, F90-safe)
_WRAP_WIDTH = 100  # target width when wrapping expressions
_MAX_LINES_PER_STMT = 40  # statement length cap for the production-rate sums


def f90_real(x) -> str:
    """Shortest round-trip literal for a double, e.g. 3.298677_WP, 1.2e+17_WP."""
    x = float(x)
    if not math.isfinite(x):
        raise MechanismError(f"Non-finite value {x!r} in generated data")
    return repr(x) + "_WP"


def f90_int(x) -> str:
    return str(int(x))


def _constructor_body(values, fmt, per_line, comments=None, indent="       "):
    """Inner lines of an array constructor: `v1, v2, ..., &` with `]` on the last line."""
    lines = []
    n = len(values)
    for i in range(0, n, per_line):
        chunk = values[i : i + per_line]
        text = indent + ", ".join(fmt(v) for v in chunk)
        text += " ]" if i + per_line >= n else ", &"
        if comments is not None:
            comment = "  ! " + ", ".join(str(c) for c in comments[i : i + per_line])
            if len(text) + len(comment) > F90_MAX_LINE:
                comment = comment[: max(0, F90_MAX_LINE - len(text) - 3)] + "..."
            text += comment
        lines.append(text)
    return lines


def f90_array_parameter(name, values, kind="real(WP)", fmt=f90_real, per_line=4, comments=None,
                        dim_expr=None, shape=None):
    """
    `kind, parameter, dimension(dim) :: name = [ ... ]` (or reshape([...], shape=[...]) if shape given).
    Long arrays are split into chunk parameters name_c1, name_c2, ... and concatenated, so no single
    statement exceeds the continuation-line limit.
    """
    values = list(values)
    n = len(values)
    if n == 0:
        raise ValueError(f"Refusing to emit empty parameter array '{name}'")
    dim = dim_expr if dim_expr is not None else str(n)
    if shape is not None:
        head_open, tail = "reshape([ &", f", shape=[{shape}])"
    else:
        head_open, tail = "[ &", ""
    if n <= _MAX_CONSTRUCTOR_VALS:
        lines = [f"  {kind}, parameter, dimension({dim}) :: {name} = {head_open}"]
        lines.extend(_constructor_body(values, fmt, per_line, comments))
        lines[-1] = lines[-1].replace(" ]", " ]" + tail, 1) if tail else lines[-1]
        return lines
    lines = []
    chunk_names = []
    for c, start in enumerate(range(0, n, _MAX_CONSTRUCTOR_VALS)):
        chunk = values[start : start + _MAX_CONSTRUCTOR_VALS]
        cname = f"{name}_c{c + 1}"
        chunk_names.append(cname)
        lines.append(f"  {kind}, parameter, dimension({len(chunk)}) :: {cname} = [ &")
        lines.extend(_constructor_body(chunk, fmt, per_line,
                                       comments[start : start + len(chunk)] if comments else None))
    lines.append(f"  {kind}, parameter, dimension({dim}) :: {name} = {head_open}")
    body = _constructor_body(chunk_names, str, 4)
    if tail:
        body[-1] = body[-1].replace(" ]", " ]" + tail, 1)
    lines.extend(body)
    return lines


def f90_2d_array_parameter(name, arr, nrows, ncols, kind="real(WP)", fmt=f90_real, per_line=4, dim_expr=None):
    """Column-major (Fortran order) 2D parameter array via reshape."""
    flat = np.asarray(arr).flatten(order="F")
    dim = dim_expr if dim_expr is not None else f"{nrows},{ncols}"
    return f90_array_parameter(name, flat, kind=kind, fmt=fmt, per_line=per_line, dim_expr=dim, shape=dim)


def f90_2d_data_variable(name, arr, nrows, ncols, kind="real(WP)", fmt=f90_real, per_line=4, dim_expr=None):
    """
    Column-major 2D module variable initialized with DATA statements (one per column chunk). Unlike a
    parameter array constructor this has no 65535-element limit, so it is used for the nS x nS tables.
    The variable is PROTECTED (read-only outside the module) and statically initialized (no run-time cost).
    """
    a = np.asarray(arr)
    dim = dim_expr if dim_expr is not None else f"{nrows},{ncols}"
    lines = [f"  {kind}, protected, dimension({dim}) :: {name}"]
    for j in range(ncols):
        col = list(a[:, j])
        for start in range(0, nrows, _MAX_CONSTRUCTOR_VALS):
            chunk = col[start : start + _MAX_CONSTRUCTOR_VALS]
            lines.append(f"  data {name}({start + 1}:{start + len(chunk)},{j + 1}) / &")
            body = _constructor_body(chunk, fmt, per_line)
            body[-1] = body[-1].replace(" ]", " /", 1)
            lines.extend(body)
    return lines


def f90_1d_data_variable(name, values, kind="real(WP)", fmt=f90_real, per_line=4, dim_expr=None, attrs="save"):
    """1D module variable (modifiable, e.g. `save`) statically initialized with chunked DATA statements."""
    values = list(values)
    n = len(values)
    dim = dim_expr if dim_expr is not None else str(n)
    lines = [f"  {kind}, {attrs}, dimension({dim}) :: {name}"]
    for start in range(0, n, _MAX_CONSTRUCTOR_VALS):
        chunk = values[start : start + _MAX_CONSTRUCTOR_VALS]
        lines.append(f"  data {name}({start + 1}:{start + len(chunk)}) / &")
        body = _constructor_body(chunk, fmt, per_line)
        body[-1] = body[-1].replace(" ]", " /", 1)
        lines.extend(body)
    return lines


def f90_wrap(head: str, terms, indent: str, tail: str = "", width: int = _WRAP_WIDTH):
    """
    Join `terms` (each carrying its own leading operator) after `head` into one Fortran statement,
    breaking into continuation lines so each stays below `width` columns.
    """
    lines = []
    cur = head
    for k, t in enumerate(terms):
        piece = t if (k == 0 and head.endswith(("(", "= "))) else " " + t
        if k > 0 and len(cur) + len(piece) > width:
            lines.append(cur + " &")
            cur = indent + t
        else:
            cur += piece
    lines.append(cur + tail)
    return lines


def f90_sum_statements(lhs: str, terms, indent: str):
    """
    `lhs = t1 t2 ...` where terms carry their sign ('+ x' / '- x'); the first term drops a leading '+'.
    Splits into several `lhs = lhs ...` statements when the wrapped statement would exceed the
    continuation-line cap.
    """
    lines = []
    stmt_lines = []
    cur = None
    first_stmt = True

    def flush():
        nonlocal stmt_lines, cur, first_stmt
        if cur is not None:
            stmt_lines.append(cur)
        lines.extend(stmt_lines)
        stmt_lines = []
        cur = None
        first_stmt = False

    for t in terms:
        if cur is None:
            head = f"{lhs} = " if first_stmt else f"{lhs} = {lhs.strip()} "
            cur = head + (t[2:] if (first_stmt and t.startswith("+ ")) else t)
            continue
        if len(cur) + 1 + len(t) > _WRAP_WIDTH:
            stmt_lines.append(cur + " &")
            cur = indent + t
            if len(stmt_lines) >= _MAX_LINES_PER_STMT:
                flush()
        else:
            cur += " " + t
    flush()
    return lines


def check_fortran_limits(lines) -> None:
    """Defensive check that the generated source respects the free-form length limits."""
    cont = 0
    for ln, text in enumerate(lines, 1):
        if len(text) > F90_MAX_LINE:
            raise RuntimeError(f"Generated line {ln} has {len(text)} > {F90_MAX_LINE} characters:\n{text}")
        code = text.split("!", 1)[0] if "'" not in text else text
        if code.rstrip().endswith("&"):
            cont += 1
            if cont > F90_MAX_CONTINUATION:
                raise RuntimeError(f"Generated statement ending near line {ln} exceeds "
                                   f"{F90_MAX_CONTINUATION} continuation lines")
        else:
            cont = 0


# ---------------------------------------------------------------------------------------------------
# Species: names, composition, thermo, transport
# ---------------------------------------------------------------------------------------------------
_RESERVED_F90 = {
    "ns", "nr", "na", "rcst", "p_ref", "w_sp", "koveps", "mucoeff", "dcoeffs", "ocoeffs", "cpsp", "hsp",
    "t_mid", "thermo_coeffs", "atom_names", "atom_masses", "comp_matrix", "n_arr", "n_arr_fwd",
    "n_arr_rev", "n_tb", "n_tb_rev", "n_fo", "n_fo_rev", "n_plog", "n_plog_rev", "n_plog_pts", "n_m",
    "ln_a_tb_b", "b_tb_b", "e_r_tb_b", "a_sign_tb_b", "rev_lna_fo", "rev_b_fo", "rev_e_r_fo", "fo_rev_src",
    "rev_lna_plog", "rev_b_plog", "rev_e_r_plog", "plog_rev_src", "ln_p_plog", "n_rev", "n_kcinv", "rev_fwd_w",
    "rev_bwd_w", "fo_m_idx", "t_nasa_cache", "cp_r_cache", "h_rt_cache", "s_r_cache", "t_rate_cache", "p_rate_cache",
    "k_t_cache", "k0_cache", "kinf_cache", "fc_cache", "kcinv_cache", "wn", "fcmech_vexp", "update_rate_cache",
    "update_plog_cache", "cdes", "wdes", "ddot", "get_destruction_rates", "fcmech_get_ydot_ddot",
    "n_eff", "eff_species_idx", "eff_matrix_reduced", "default_eff", "ln_a_arr", "b_arr", "e_r_arr",
    "a_sign_arr", "ln_a_tb", "b_tb", "e_r_tb", "a_sign_tb", "ln_a0_fo", "b0_fo", "e_r0_fo",
    "ln_ainf_fo", "binf_fo", "e_rinf_fo", "troe_a", "troe_t3", "troe_t1", "troe_t2", "p_plog",
    "ln_a_plog", "b_plog", "e_r_plog", "plog_reac_start", "plog_reac_len", "i_arr_end", "i_tb_end",
    "i_tb_bwd_end", "i_fo_end", "i_fo_bwd_end", "i_plog_start", "i_plog_end", "i_plog_bwd_start",
    "i_plog_bwd_end", "precision", "string", "wp",
}


def build_f90_species_names(species_names) -> dict:
    """Unique, valid Fortran identifiers sXXX for every species (case-insensitive uniqueness)."""
    names = {}
    used = set(_RESERVED_F90)
    for nm in species_names:
        base = "s" + re.sub(r"[^A-Za-z0-9_]", "_", str(nm))
        base = base[:40]  # keeps every emitted statement head well inside 132 columns
        cand = base
        k = 2
        while cand.lower() in used:
            cand = f"{base}_{k}"
            k += 1
        used.add(cand.lower())
        names[nm] = cand
    return names


def species_composition(spec: dict, name: str) -> dict:
    comp = spec.get("composition")
    if not isinstance(comp, dict) or not comp:
        raise MechanismError(f"Species '{name}': missing or empty 'composition'")
    out = {}
    for elem, count in comp.items():
        sym = canonical_element(elem)
        try:
            cnt = float(count)
        except (TypeError, ValueError):
            raise MechanismError(f"Species '{name}': composition count '{count}' for {sym} is not a number")
        if not cnt.is_integer() or cnt < 0:
            raise MechanismError(f"Species '{name}': composition count {count} for {sym} must be a non-negative integer")
        out[sym] = out.get(sym, 0) + int(cnt)
    return out


def species_molar_mass(comp: dict) -> float:
    """Molar mass (kg/mol) from an element->count composition."""
    return sum(ATOMIC_MASSES_KG[e] * n for e, n in comp.items())


def _nasa7_eval(a, T):
    """cp/R, h/(RT), s/R from one 7-coefficient NASA7 row at temperature T."""
    T2, T3, T4 = T * T, T ** 3, T ** 4
    cp_R = a[0] + a[1] * T + a[2] * T2 + a[3] * T3 + a[4] * T4
    h_RT = a[0] + a[1] * T / 2 + a[2] * T2 / 3 + a[3] * T3 / 4 + a[4] * T4 / 5 + a[5] / T
    s_R = a[0] * math.log(T) + a[1] * T + a[2] * T2 / 2 + a[3] * T3 / 3 + a[4] * T4 / 4 + a[6]
    return cp_R, h_RT, s_R


def species_thermo(spec: dict, name: str):
    """
    Validated NASA7 data -> (low[7], high[7], T_low, T_mid, T_high).
    Single-range species use the same row for both and T_mid = T_high.
    """
    thermo = spec.get("thermo")
    if isinstance(thermo, list):
        thermo = thermo[0] if thermo else None
    if not isinstance(thermo, dict):
        raise MechanismError(f"Species '{name}': missing 'thermo' block")
    model = str(thermo.get("model", "")).strip().upper()
    if model != "NASA7":
        raise MechanismError(f"Species '{name}': thermo model '{thermo.get('model')}' is not supported (NASA7 only)")
    data = thermo.get("data")
    tr = thermo.get("temperature-ranges")
    if not isinstance(data, list) or not data or not isinstance(tr, list):
        raise MechanismError(f"Species '{name}': NASA7 thermo needs 'temperature-ranges' and 'data'")
    if isinstance(data[0], (int, float)):
        data = [data]
    rows = []
    for row in data:
        if not isinstance(row, (list, tuple)) or len(row) != 7:
            raise MechanismError(f"Species '{name}': each NASA7 data row must have exactly 7 coefficients")
        rows.append([float(x) for x in row])
    tr = [float(x) for x in tr]
    if len(rows) == 1 and len(tr) == 2:
        low = high = rows[0]
        T_low, T_high = tr
        T_mid = T_high
    elif len(rows) == 2 and len(tr) == 3:
        low, high = rows
        T_low, T_mid, T_high = tr
    else:
        raise MechanismError(f"Species '{name}': {len(rows)} NASA7 row(s) with {len(tr)} temperature bounds; "
                             f"expected 1 row + 2 bounds or 2 rows + 3 bounds")
    if not (T_low < T_mid <= T_high):
        raise MechanismError(f"Species '{name}': temperature-ranges {tr} are not increasing")

    # Physical sanity: cp must stay positive over the declared range. A negative cp is the
    # signature of swapped low/high rows (CHEMKIN card order) or corrupted coefficients.
    for T in np.linspace(T_low, T_high, 25):
        a = low if T <= T_mid else high
        cp_R = _nasa7_eval(a, T)[0]
        if not cp_R > 0.0:
            raise MechanismError(f"Species '{name}': NASA7 gives cp/R = {cp_R:.4g} <= 0 at T = {T:.0f} K; "
                                 f"check the row order (low-T range first) and the coefficients")
    if len(rows) == 2:
        cp_lo, h_lo, s_lo = _nasa7_eval(low, T_mid)
        cp_hi, h_hi, s_hi = _nasa7_eval(high, T_mid)
        jump = max(abs(cp_hi - cp_lo) / abs(cp_lo), abs(h_hi - h_lo) / max(abs(h_lo), 1.0),
                   abs(s_hi - s_lo) / max(abs(s_lo), 1.0))
        if jump > 1e-2:
            warn(f"Species '{name}': NASA7 low/high polynomials differ by {100 * jump:.1f}% at T_mid = {T_mid:.0f} K")
    return low, high, T_low, T_mid, T_high


_GEOMETRIES = {"atom", "linear", "nonlinear"}


def species_transport(spec: dict, name: str):
    """Validated gas transport data -> dict(diameter [A], well_depth [K], geometry) or None if absent."""
    trans = spec.get("transport")
    if isinstance(trans, list):
        trans = trans[0] if trans else None
    if trans is None:
        return None
    if not isinstance(trans, dict):
        raise MechanismError(f"Species '{name}': 'transport' must be a mapping")
    model = str(trans.get("model", "gas")).strip().lower()
    if model != "gas":
        raise MechanismError(f"Species '{name}': transport model '{trans.get('model')}' is not supported (gas only)")
    geometry = str(trans.get("geometry", "")).strip().lower()
    if geometry not in _GEOMETRIES:
        raise MechanismError(f"Species '{name}': transport geometry '{trans.get('geometry')}' must be one of {sorted(_GEOMETRIES)}")
    try:
        d = float(trans["diameter"])
        eps = float(trans["well-depth"])
    except (KeyError, TypeError, ValueError):
        raise MechanismError(f"Species '{name}': transport needs numeric 'diameter' and 'well-depth'")
    if d <= 0.0 or eps <= 0.0:
        raise MechanismError(f"Species '{name}': transport diameter and well-depth must be positive")
    return {"diameter": d, "well_depth": eps, "geometry": geometry}


# ---------------------------------------------------------------------------------------------------
# Mechanism loading
# ---------------------------------------------------------------------------------------------------
def load_mechanism(path):
    """Load and validate a YAML mechanism; returns the normalized mechanism dict."""
    data = load_yaml(path)
    units = parse_units(data)

    all_species = data.get("species")
    if not isinstance(all_species, list) or not all_species:
        raise MechanismError("No 'species' list found in mechanism")
    by_name = {}
    for s in all_species:
        if not isinstance(s, dict) or "name" not in s:
            raise MechanismError("Every entry of 'species' must be a mapping with a 'name'")
        nm = str(s["name"])
        if nm in by_name:
            raise MechanismError(f"Duplicate species definition '{nm}'")
        by_name[nm] = s

    # Phase (first one): species order, elements, reaction selection, third-body policy
    phase = None
    phases = data.get("phases")
    if isinstance(phases, list) and phases:
        phase = phases[0]
        if not isinstance(phase, dict):
            raise MechanismError("Entries of 'phases' must be mappings")
        if len(phases) > 1:
            print(f"  Using first phase '{phase.get('name', '?')}' of {len(phases)}")
    species_names = []
    if phase is not None and isinstance(phase.get("species"), list):
        for n in phase["species"]:
            nm = str(n)
            if nm not in by_name:
                raise MechanismError(f"Phase species '{nm}' has no definition in the 'species' list")
            species_names.append(nm)
    else:
        species_names = list(by_name)
    if len(set(species_names)) != len(species_names):
        raise MechanismError("Phase species list contains duplicates")
    species_list = [by_name[nm] for nm in species_names]
    name_to_idx = {nm: i for i, nm in enumerate(species_names)}  # 0-based

    reactions_mode = "all"
    skip_undeclared_tb = False
    if phase is not None:
        rmode = phase.get("reactions", "all")
        if isinstance(rmode, str) and rmode in ("all", "none", "declared-species"):
            reactions_mode = rmode
        elif rmode is None:
            reactions_mode = "none"
        else:
            raise MechanismError(f"Phase 'reactions: {rmode}' is not supported (use all, none, or declared-species)")
        skip_undeclared_tb = bool(phase.get("skip-undeclared-third-bodies", False))

    # Composition, molar mass, thermo, transport
    compositions = []
    W_sp = []
    thermo = []
    transport = []
    for nm, sp in zip(species_names, species_list):
        comp = species_composition(sp, nm)
        compositions.append(comp)
        W_sp.append(species_molar_mass(comp))
        thermo.append(species_thermo(sp, nm))
        transport.append(species_transport(sp, nm))

    # Elements: phase list (canonicalized) or union of compositions in a conventional order
    used_elements = set()
    for comp in compositions:
        used_elements.update(comp)
    if phase is not None and isinstance(phase.get("elements"), list):
        atom_names = [canonical_element(e) for e in phase["elements"]]
        missing = used_elements - set(atom_names)
        if missing:
            raise MechanismError(f"Elements {sorted(missing)} used in species compositions are not in the phase 'elements' list")
    else:
        atom_names = [e for e in ("O", "H", "C", "N", "Ar", "He") if e in used_elements]
        atom_names += sorted(used_elements - set(atom_names))
    nA = len(atom_names)
    nS = len(species_names)
    comp_matrix = np.zeros((nA, nS), dtype=int)
    for j, comp in enumerate(compositions):
        for e, n in comp.items():
            comp_matrix[atom_names.index(e), j] = n
    atom_masses = np.array([ATOMIC_MASSES_KG[a] for a in atom_names])

    reactions = data.get("reactions", []) or []
    if reactions_mode == "none":
        reactions = []
    if not isinstance(reactions, list):
        raise MechanismError("'reactions' must be a list")

    return {
        "species": species_list,
        "species_names": species_names,
        "name_to_idx": name_to_idx,
        "compositions": compositions,
        "W_sp": np.array(W_sp),
        "thermo": thermo,
        "transport": transport,
        "reactions": reactions,
        "reactions_mode": reactions_mode,
        "skip_undeclared_third_bodies": skip_undeclared_tb,
        "units": units,
        "atom_names": atom_names,
        "atom_masses": atom_masses,
        "comp_matrix": comp_matrix,
        "nA": nA,
        "_resolve_warned": set(),
    }


def extract_yaml_labels(yaml_path) -> dict:
    """
    Trailing '# ...' comments on the `equation:` lines of the top-level `reactions:` block (YAML strips
    them), keyed by reaction number = position in that block. Items are counted by their list dash, so a
    reaction whose mapping starts with another key (e.g. `- type: falloff`) is numbered correctly.
    """
    labels = {}
    with open(yaml_path) as f:
        lines = f.read().splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if re.match(r"^reactions\s*:\s*(#.*)?$", ln))
    except StopIteration:
        return labels
    rid = 0
    item_indent = None
    for ln in lines[start + 1:]:
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        indent = len(ln) - len(ln.lstrip())
        is_item = re.match(r"^\s*-\s", ln) is not None
        if indent == 0 and not is_item:
            break  # next top-level key (sequence items themselves may sit at column 0)
        if is_item and (item_indent is None or indent == item_indent):
            item_indent = indent
            rid += 1
        if re.match(r"^\s*(-\s*)?equation\s*:", ln) and "#" in ln:
            labels[rid] = "# " + ln.split("#", 1)[1].strip()
    return labels


# ---------------------------------------------------------------------------------------------------
# Reaction equations
# ---------------------------------------------------------------------------------------------------
_M_TOKEN_RE = re.compile(r"^M'*$")  # M (Cantera) or M', M'' (FlameMaster collider sets)
_PAREN_COLLIDER_RE = re.compile(r"\(\+\s*([^()\s]+)\s*\)")
_TERM_COLLIDER_RE = re.compile(r"^(.*?)\s*\(\+([^()\s]+)\)$")
_COEFF_SPECIES_RE = re.compile(r"^([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s+(\S+)$")


def parse_equation(eq: str):
    """
    'a A + b B (+M) <=> c C (+M)' -> (left, right, reversible, collider)
    left/right: [(species_token, coeff)]; collider: None, or {'name': 'M' | species, 'paren': bool}.
    Species tokens are not resolved here (see resolve_side). '=' alone (CHEMKIN) means reversible.
    """
    s = str(eq).strip()
    if "<=>" in s:
        rev, (lhs, rhs) = True, s.split("<=>", 1)
    elif "=>" in s:
        rev, (lhs, rhs) = False, s.split("=>", 1)
    elif "<=" in s:
        raise MechanismError(f"Equation '{eq}': arrow '<=' is not supported")
    elif "=" in s:
        rev, (lhs, rhs) = True, s.split("=", 1)
    else:
        raise MechanismError(f"Equation '{eq}': no reaction arrow found")
    if any(ch in lhs + rhs for ch in "=<>"):
        raise MechanismError(f"Equation '{eq}': more than one reaction arrow")

    colliders = []

    def parse_side(text):
        text = _PAREN_COLLIDER_RE.sub(lambda m: f" (+{m.group(1)}) ", text)
        terms = re.split(r"\s\+\s", " " + text.strip() + " ")
        out = []
        for t in terms:
            t = t.strip()
            if not t:
                continue
            m = _TERM_COLLIDER_RE.match(t)
            if m:
                colliders.append({"name": m.group(2), "paren": True})
                t = m.group(1).strip()
                if not t:
                    continue
            m = _COEFF_SPECIES_RE.match(t)
            if m:
                coeff, name = float(m.group(1)), m.group(2)
            else:
                coeff, name = 1.0, t
            if _M_TOKEN_RE.match(name):
                if coeff != 1.0:
                    raise MechanismError(f"Equation '{eq}': third body M cannot carry a coefficient")
                colliders.append({"name": "M", "paren": False})
                continue
            if coeff <= 0.0:
                raise MechanismError(f"Equation '{eq}': non-positive stoichiometric coefficient for '{name}'")
            out.append((name, coeff))
        return out

    left = parse_side(lhs)
    right = parse_side(rhs)
    if not left or not right:
        raise MechanismError(f"Equation '{eq}': both sides must contain at least one species")
    collider = None
    if colliders:
        if len(colliders) != 2 or colliders[0] != colliders[1]:
            raise MechanismError(f"Equation '{eq}': the third body must appear exactly once on each side, identically")
        collider = colliders[0]
        if collider["paren"] and _M_TOKEN_RE.match(collider["name"]):
            collider["name"] = "M"
    return left, right, rev, collider


def _resolve_species_name(token: str, mech: dict):
    """Exact match, else case-insensitive match, else 'NAME-suffix' -> NAME (each with a warning)."""
    names = mech["name_to_idx"]
    if token in names:
        return token
    warned = mech["_resolve_warned"]
    upper = token.upper()
    for s in names:
        if s.upper() == upper:
            if token not in warned:
                warn(f"Species '{token}' in an equation resolved case-insensitively to '{s}'")
                warned.add(token)
            return s
    for s in sorted(names, key=len, reverse=True):
        if token.startswith(s + "-") or upper.startswith(s.upper() + "-"):
            if token not in warned:
                warn(f"Species '{token}' in an equation resolved to '{s}' by dropping its suffix")
                warned.add(token)
            return s
    return None


def resolve_side(side, mech: dict, ctx: str):
    """Resolve tokens to mechanism species names, merge repeated species (O + O -> 2 O)."""
    out = {}
    order = []
    for token, coeff in side:
        name = _resolve_species_name(token, mech)
        if name is None and coeff == 1.0:
            # CHEMKIN-style coefficient glued to the name ('2OH'): only if the remainder is a species
            m = re.fullmatch(r"([0-9]+)(\D\S*)", token)
            if m and _resolve_species_name(m.group(2), mech) is not None:
                name = _resolve_species_name(m.group(2), mech)
                coeff = float(m.group(1))
                if token not in mech["_resolve_warned"]:
                    warn(f"{ctx}: '{token}' read as {coeff:g} {name}")
                    mech["_resolve_warned"].add(token)
        if name is None:
            raise MechanismError(f"{ctx}: unknown species '{token}'")
        if name not in out:
            out[name] = 0.0
            order.append(name)
        out[name] += coeff
    return [(n, out[n]) for n in order]


def format_side(side, collider=None) -> str:
    """'2 O + M' style rendering of a reaction side for comments and the index file."""
    parts = []
    for sp, c in side:
        parts.append(sp if c == 1.0 else f"{c:g} {sp}")
    text = " + ".join(parts)
    if collider is not None:
        text += f" (+{collider['name']})" if collider["paren"] else f" + {collider['name']}"
    return text


# ---------------------------------------------------------------------------------------------------
# Reactions
# ---------------------------------------------------------------------------------------------------
class Reaction:
    """One YAML reaction, fully validated, rate parameters in SI (m, mol, s, J/mol)."""

    def __init__(self, rid: int, eq: str):
        self.rid = rid  # 1-based position in the YAML reactions list
        self.eq = eq
        self.left = []  # [(species, coeff)]
        self.right = []
        self.rev = False
        self.collider = None  # None | {'name', 'paren'}
        self.kind = None  # 'arrhenius' | 'three-body' | 'falloff' | 'plog'
        self.A = self.b = self.Ea = None  # arrhenius / three-body, A in SI
        self.sign = 1.0  # sign of A (negative-A: true)
        self.eff = {}  # species -> third-body efficiency
        self.default_eff = 1.0
        self.A0 = self.b0 = self.Ea0 = None  # falloff low-pressure limit
        self.Ainf = self.binf = self.Eainf = None  # falloff high-pressure limit
        self.troe = None  # (a, T3, T1, T2) with T2 = 1e30 when absent (Lindemann: a=1, T3=T1=1e30)
        self.plog = []  # [(P_Pa, A_SI, b, Ea)] sorted by P
        self.nu = {}  # species index (0-based) -> net stoichiometric coefficient
        self.dnu = 0.0
        # filled by the generator: 0-based positions in the k/w arrays and third-body index
        self.w_fwd = None
        self.w_rev = None
        self.m_idx = None

    def label(self) -> str:
        arrow = "<=>" if self.rev else "=>"
        return f"{format_side(self.left, self.collider)} {arrow} {format_side(self.right, self.collider)}"

    @property
    def order_fwd(self) -> float:
        return sum(c for _, c in self.left)

    @property
    def order_rev(self) -> float:
        return sum(c for _, c in self.right)


_COMMON_KEYS = {"equation", "type", "duplicate", "note", "negative-A"}
_KIND_KEYS = {
    "arrhenius": {"rate-constant"},
    "three-body": {"rate-constant", "efficiencies", "default-efficiency"},
    "falloff": {"low-P-rate-constant", "high-P-rate-constant", "Troe", "efficiencies", "default-efficiency"},
    "plog": {"rate-constants"},
}
_UNSUPPORTED_KEYS = {
    "orders": "explicit reaction orders",
    "negative-orders": "explicit reaction orders",
    "nonreactant-orders": "explicit reaction orders",
    "SRI": "SRI falloff",
    "units": "per-reaction units",
    "reverse-rate-constant": "explicit reverse rate constants (write two irreversible reactions)",
    "Chebyshev": "Chebyshev rates",
}
_SUPPORTED_TYPES = {
    None: "arrhenius", "elementary": "arrhenius", "Arrhenius": "arrhenius",
    "three-body": "three-body", "three-body-Arrhenius": "three-body",
    "falloff": "falloff", "pressure-dependent-Arrhenius": "plog", "plog": "plog",
}


def _parse_arrhenius(rc, units: dict, ctx: str, allow_negative_A: bool = False):
    """{A, b, Ea} in mechanism units -> (A, b, Ea_J_per_mol); A is NOT yet converted to SI."""
    if not isinstance(rc, dict):
        raise MechanismError(f"{ctx}: rate constant must be a mapping {{A, b, Ea}}")
    unknown = set(rc) - {"A", "b", "Ea"}
    if unknown:
        raise MechanismError(f"{ctx}: unknown rate-constant key(s) {sorted(unknown)}")
    if "A" not in rc:
        raise MechanismError(f"{ctx}: rate constant is missing 'A'")
    A = rc["A"]
    if isinstance(A, str):
        try:
            A = float(A)
        except ValueError:
            raise MechanismError(f"{ctx}: pre-exponential factor with an explicit unit string ('{A}') is not supported")
    A = float(A)
    b = float(rc.get("b", 0.0))
    Ea = value_with_unit(rc.get("Ea", 0.0), _EA_UNITS, units["Ea_J_per_mol"], f"{ctx}: Ea")
    if A == 0.0:
        raise MechanismError(f"{ctx}: A = 0 (identically zero rate) is not supported; remove the reaction")
    if A < 0.0 and not allow_negative_A:
        raise MechanismError(f"{ctx}: negative pre-exponential factor requires 'negative-A: true'")
    return A, b, Ea


def _parse_troe(t, ctx: str):
    """Troe: {A, T3, T1[, T2]} or [a, T3, T1[, T2]]. T2 absent or 0 -> term omitted (Cantera behaviour)."""
    if isinstance(t, dict):
        unknown = set(t) - {"A", "a", "T3", "T1", "T2"}
        if unknown:
            raise MechanismError(f"{ctx}: unknown Troe key(s) {sorted(unknown)}")
        if ("A" not in t and "a" not in t) or "T3" not in t or "T1" not in t:
            raise MechanismError(f"{ctx}: Troe needs A, T3 and T1")
        a = float(t["A"] if "A" in t else t["a"])
        T3, T1 = float(t["T3"]), float(t["T1"])
        T2 = float(t["T2"]) if t.get("T2") is not None else 0.0
    elif isinstance(t, (list, tuple)) and len(t) in (3, 4):
        a, T3, T1 = (float(x) for x in t[:3])
        T2 = float(t[3]) if len(t) == 4 and t[3] is not None else 0.0
    else:
        raise MechanismError(f"{ctx}: Troe must be a mapping {{A, T3, T1[, T2]}} or a list of 3 or 4 values")
    if T3 == 0.0 or T1 == 0.0:
        raise MechanismError(f"{ctx}: Troe T3 and T1 must be non-zero")
    return a, T3, T1, (T2 if T2 != 0.0 else 1.0e30)


def parse_reactions(mech: dict) -> list:
    """Validate and normalize all reactions; returns a list of Reaction (SI units, resolved species)."""
    units = mech["units"]
    name_to_idx = mech["name_to_idx"]
    species_set = set(mech["species_names"])
    out = []
    n_skipped = 0
    for i, r in enumerate(mech["reactions"]):
        rid = i + 1
        ctx = f"Reaction {rid}"
        if not isinstance(r, dict):
            raise MechanismError(f"{ctx}: entry is not a mapping")
        eq = r.get("equation")
        if not eq:
            raise MechanismError(f"{ctx}: missing 'equation'")
        ctx = f"Reaction {rid} '{eq}'"
        for key, what in _UNSUPPORTED_KEYS.items():
            if key in r:
                raise MechanismError(f"{ctx}: {what} ('{key}') are not supported")

        left_tok, right_tok, rev, collider = parse_equation(eq)
        if mech["reactions_mode"] == "declared-species":
            tokens = [t for t, _ in left_tok + right_tok]
            if any(_resolve_species_name(t, mech) is None for t in tokens):
                n_skipped += 1
                continue
        rx = Reaction(rid, str(eq))
        rx.left = resolve_side(left_tok, mech, ctx)
        rx.right = resolve_side(right_tok, mech, ctx)
        rx.rev = rev
        rx.collider = collider

        # Reaction type: from `type:`, cross-checked against the third-body syntax in the equation
        rtype = r.get("type")
        if rtype not in _SUPPORTED_TYPES:
            raise MechanismError(f"{ctx}: reaction type '{rtype}' is not supported "
                                 f"(supported: elementary, three-body, falloff, pressure-dependent-Arrhenius)")
        kind = _SUPPORTED_TYPES[rtype]
        if kind == "arrhenius" and collider is not None:
            if collider["paren"]:
                raise MechanismError(f"{ctx}: '(+{collider['name']})' requires 'type: falloff'")
            if rtype is not None:
                raise MechanismError(f"{ctx}: 'type: {rtype}' contradicts the '+ M' third body in the equation")
            kind = "three-body"  # Cantera >= 3.0 infers three-body from '+ M'
        if kind == "three-body" and (collider is None or collider["paren"]):
            raise MechanismError(f"{ctx}: 'type: three-body' requires '+ M' (or '+ X' for a specific collider) "
                                 f"on both sides; write explicit colliders as ordinary reactants instead")
        if kind == "falloff" and (collider is None or not collider["paren"]):
            raise MechanismError(f"{ctx}: 'type: falloff' requires '(+M)' or '(+X)' on both sides")
        if kind == "plog" and collider is not None:
            raise MechanismError(f"{ctx}: pressure-dependent-Arrhenius reactions cannot have a third body")
        rx.kind = kind

        allowed = _COMMON_KEYS | _KIND_KEYS[kind]
        unknown = set(r) - allowed
        if unknown:
            raise MechanismError(f"{ctx}: unknown key(s) {sorted(unknown)} for a {kind} reaction")
        negative_A = bool(r.get("negative-A", False))
        if negative_A and kind not in ("arrhenius", "three-body"):
            raise MechanismError(f"{ctx}: 'negative-A' is only supported for elementary and three-body reactions")

        if kind in ("arrhenius", "three-body"):
            if "rate-constant" not in r:
                raise MechanismError(f"{ctx}: missing 'rate-constant'")
            A, b, Ea = _parse_arrhenius(r["rate-constant"], units, ctx, allow_negative_A=negative_A)
            order = rx.order_fwd + (1.0 if kind == "three-body" else 0.0)
            rx.sign = -1.0 if A < 0.0 else 1.0
            rx.A, rx.b, rx.Ea = A_to_si(abs(A), order, units), b, Ea
        elif kind == "falloff":
            for key in ("low-P-rate-constant", "high-P-rate-constant"):
                if key not in r:
                    raise MechanismError(f"{ctx}: missing '{key}'")
            A0, b0, Ea0 = _parse_arrhenius(r["low-P-rate-constant"], units, f"{ctx} low-P")
            Ainf, binf, Eainf = _parse_arrhenius(r["high-P-rate-constant"], units, f"{ctx} high-P")
            rx.A0, rx.b0, rx.Ea0 = A_to_si(A0, rx.order_fwd + 1.0, units), b0, Ea0
            rx.Ainf, rx.binf, rx.Eainf = A_to_si(Ainf, rx.order_fwd, units), binf, Eainf
            rx.troe = _parse_troe(r["Troe"], ctx) if "Troe" in r else (1.0, 1.0e30, 1.0e30, 1.0e30)
        elif kind == "plog":
            rcs = r.get("rate-constants")
            if not isinstance(rcs, list) or not rcs:
                raise MechanismError(f"{ctx}: 'rate-constants' must be a non-empty list of {{P, A, b, Ea}}")
            pts = []
            for rd in rcs:
                if not isinstance(rd, dict) or "P" not in rd:
                    raise MechanismError(f"{ctx}: each PLOG entry needs 'P' and 'A'")
                P = parse_pressure(rd["P"], units)
                A, b, Ea = _parse_arrhenius({k: v for k, v in rd.items() if k != "P"}, units, f"{ctx} P={rd['P']}")
                pts.append((P, A_to_si(A, rx.order_fwd, units), b, Ea))
            pts.sort(key=lambda p: p[0])
            rx.plog = pts

        # Third-body efficiencies
        if kind in ("three-body", "falloff"):
            eff = r.get("efficiencies", {}) or {}
            if not isinstance(eff, dict):
                raise MechanismError(f"{ctx}: 'efficiencies' must be a mapping species: value")
            if collider["name"] != "M":
                if eff or "default-efficiency" in r:
                    raise MechanismError(f"{ctx}: explicit collider '(+{collider['name']})' cannot be combined with efficiencies")
                if collider["name"] not in species_set:
                    raise MechanismError(f"{ctx}: collider '{collider['name']}' is not a species of the mechanism")
                rx.eff = {collider["name"]: 1.0}
                rx.default_eff = 0.0
            else:
                rx.default_eff = float(r.get("default-efficiency", 1.0))
                for sp, val in eff.items():
                    sp = str(sp)
                    if sp not in species_set:
                        if mech["skip_undeclared_third_bodies"]:
                            continue
                        raise MechanismError(f"{ctx}: third-body efficiency for undeclared species '{sp}' "
                                             f"(set 'skip-undeclared-third-bodies: true' in the phase to drop these)")
                    rx.eff[sp] = float(val)

        # Net stoichiometry and element balance
        nu = {}
        for sp, c in rx.left:
            nu[name_to_idx[sp]] = nu.get(name_to_idx[sp], 0.0) - c
        for sp, c in rx.right:
            nu[name_to_idx[sp]] = nu.get(name_to_idx[sp], 0.0) + c
        rx.nu = {k: v for k, v in nu.items() if v != 0.0}
        rx.dnu = sum(rx.nu.values())
        imbalance = mech["comp_matrix"].astype(float) @ np.array([nu.get(j, 0.0) for j in range(len(species_set))])
        bad = [f"{a}: {d:+g}" for a, d in zip(mech["atom_names"], imbalance) if abs(d) > 1e-8]
        if bad:
            raise MechanismError(f"{ctx}: element balance violated ({', '.join(bad)})")
        out.append(rx)
    if n_skipped:
        print(f"  Skipped {n_skipped} reaction(s) involving undeclared species (reactions: declared-species)")
    return out


# ---------------------------------------------------------------------------------------------------
# Transport precomputation (shared with check_diffusion.py)
# ---------------------------------------------------------------------------------------------------
def build_transport_arrays(mech: dict):
    """
    Chapman-Enskog / Lennard-Jones precomputation, or None if any species lacks transport data.
      koveps(i)     = 1/(eps_i/k)                        [1/K]; T* = T*koveps
      mucoeff(i)    : mu_i = mucoeff(i)*sqrt(T)/Omega_mu(T*koveps(i))                   [Pa s]
      Dcoeffs(i,j)  : D_ij = Dcoeffs(i,j)*T^1.5/(P*Omega_D(T*Ocoeffs(i,j)))             [m^2/s], P in Pa
      Ocoeffs(i,j)  = sqrt(koveps(i)*koveps(j))  (geometric-mean well depth; no polar correction)
    """
    transport = mech["transport"]
    if any(t is None for t in transport):
        return None
    W = np.asarray(mech["W_sp"], dtype=float)  # kg/mol
    nS = len(W)
    sigma = np.array([t["diameter"] * 1e-10 for t in transport])  # m
    eps_k = np.array([t["well_depth"] for t in transport])  # K
    m = W / N_AVOGADRO  # kg per molecule
    koveps = 1.0 / eps_k
    mucoeff = 5.0 / 16.0 * np.sqrt(math.pi * m * K_BOLTZMANN) / (math.pi * sigma ** 2)
    Dcoeffs = np.ones((nS, nS))
    Ocoeffs = np.ones((nS, nS))
    for i in range(nS):
        for j in range(nS):
            if i == j:
                continue
            sigma_ij = 0.5 * (sigma[i] + sigma[j])
            m_ij = m[i] * m[j] / (m[i] + m[j])
            Dcoeffs[i, j] = 3.0 / 16.0 * math.sqrt(2.0 * math.pi * K_BOLTZMANN ** 3 / m_ij) / (math.pi * sigma_ij ** 2)
            Ocoeffs[i, j] = math.sqrt(koveps[i] * koveps[j])
    return {"koveps": koveps, "mucoeff": mucoeff, "Dcoeffs": Dcoeffs, "Ocoeffs": Ocoeffs}


# ---------------------------------------------------------------------------------------------------
# Code generation
# ---------------------------------------------------------------------------------------------------
def _coef(c: float) -> str:
    return f90_real(abs(c))


def _signed_term(c: float, expr: str) -> str:
    """'+ expr', '- expr', '+ 2.0_WP*expr', ..."""
    sgn = "+" if c > 0 else "-"
    if abs(c) == 1.0:
        return f"{sgn} {expr}"
    return f"{sgn} {_coef(c)}*{expr}"


def _reactant_terms(side, f90n: dict):
    """['c(sA)', '* c(sB)**2', ...] for the concentration product of one side."""
    terms = []
    for k, (sp, c) in enumerate(side):
        base = f"c({f90n[sp]})"
        if c != 1.0:
            base += f"**{int(c)}" if float(c).is_integer() else f"**{f90_real(c)}"
        terms.append(base if k == 0 else "* " + base)
    return terms


def _kc_terms(rx: Reaction, f90n: dict, species_names) -> list:
    """Terms of ln(k_r/k_f) = sum_i nu_i g_i/(RT) - dnu ln(P_ref/(RT)); products first, then reactants."""
    items = sorted(rx.nu.items(), key=lambda kv: (-np.sign(kv[1]), kv[0]))
    terms = [_signed_term(nu, f"g_RT({f90n[species_names[idx]]})") for idx, nu in items]
    if rx.dnu != 0.0:
        terms.append(_signed_term(-rx.dnu, "ln_P0RT"))
    return terms


def ln_kc_table(rx: Reaction, T_arr, mech: dict):
    """ln K_c(T) from the mechanism NASA7 data, with the same low/high switching as the generated Fortran."""
    thermo = mech["thermo"]
    out = np.zeros(len(T_arr))
    for k, T in enumerate(T_arr):
        s = 0.0
        for idx, nu in rx.nu.items():
            low, high, _, T_mid, _ = thermo[idx]
            _, h_RT, s_R = _nasa7_eval(low if T <= T_mid else high, T)
            s += nu * (h_RT - s_R)
        out[k] = -s + rx.dnu * math.log(P_REF / (R_UNIV * T))
    return out


def fit_arrhenius(T_arr, ln_k, n_iter: int = 60):
    """
    Minimax fit of ln k = ln A + b ln T - (E/R)/T to tabulated ln k(T): the maximum |ln k_fit - ln k|
    (i.e. the worst relative error of k) is minimized with Lawson's iteratively reweighted least squares,
    starting from the plain least-squares solution and keeping the best iterate.
    Returns (ln A, b, E/R, max relative error of the fitted k over the table).
    """
    X = np.column_stack([np.ones(len(T_arr)), np.log(T_arr), -1.0 / np.asarray(T_arr)])
    w = np.ones(len(T_arr))
    best_coef, best_max = None, np.inf
    for _ in range(n_iter):
        sw = np.sqrt(w)
        coef, *_ = np.linalg.lstsq(X * sw[:, None], ln_k * sw, rcond=None)
        res = np.abs(X @ coef - ln_k)
        if res.max() < best_max:
            best_coef, best_max = coef, res.max()
        if res.max() < 1e-12:
            break
        w = w * res
        w /= w.sum()
    err = float(np.max(np.abs(np.exp(X @ best_coef - ln_k) - 1.0)))
    return float(best_coef[0]), float(best_coef[1]), float(best_coef[2]), err


def fit_reverse_rates(rxns, mech: dict, T_min: float, T_max: float, n_pts: int) -> dict:
    """
    Arrhenius fit of 1/K_c(T) for every reversible reaction, sampled uniformly in 1/T:
        ln(1/K_c) ~ lnA_kc + b_kc ln T - E_R_kc/T,   so   k_r = k_f * exp(lnA_kc + b_kc ln T - E_R_kc/T).
    Returns {rid: {"lnA","b","E_R","err"}}. For elementary/three-body reactions the generator folds the
    fit into the forward Arrhenius parameters (k_f is exactly Arrhenius, so k_r is one Arrhenius
    expression); for falloff/PLOG it multiplies the forward k_f by exp(fit). "err" is the largest
    relative error of k_r over [T_min, T_max] (identical for all three uses).
    """
    T_arr = 1.0 / np.linspace(1.0 / T_max, 1.0 / T_min, n_pts)
    fits = {}
    for rx in rxns:
        if not rx.rev:
            continue
        lnA, b, E_R, err = fit_arrhenius(T_arr, -ln_kc_table(rx, T_arr, mech))
        fits[rx.rid] = {"lnA": lnA, "b": b, "E_R": E_R, "err": err}
    return fits


def yaml2nga(mech: dict, output_path, yaml_path=None, fit_reverse=None) -> None:
    """
    Generate the fcmech Fortran module and the reaction index file.
    fit_reverse: None (reverse rates from detailed balance at run time) or a dict
    {"T_min", "T_max", "n_pts"} to emit fitted reverse Arrhenius parameters instead.
    """
    species_names = mech["species_names"]
    nS = len(species_names)
    nA = mech["nA"]
    f90n = build_f90_species_names(species_names)
    rxns = parse_reactions(mech)

    # NASA7 table: thermo_coeffs(nS,14) = [low a1..a7, high a1..a7], T_mid(nS)
    thermo_coeffs = np.zeros((nS, 14))
    T_mid_sp = np.zeros(nS)
    for i, (low, high, _, T_mid, _) in enumerate(mech["thermo"]):
        thermo_coeffs[i, 0:7] = low
        thermo_coeffs[i, 7:14] = high
        T_mid_sp[i] = T_mid

    transport = build_transport_arrays(mech)
    if transport is None:
        missing = [nm for nm, t in zip(species_names, mech["transport"]) if t is None]
        warn(f"{len(missing)} species without transport data ({', '.join(missing[:6])}"
             f"{', ...' if len(missing) > 6 else ''}): transport arrays and routines are NOT generated")

    # ---- k / w layout: arr fwd | arr rev | tb fwd | tb rev | fo fwd | fo rev | plog fwd | plog rev ----
    arr = [r for r in rxns if r.kind == "arrhenius"]
    tb = [r for r in rxns if r.kind == "three-body"]
    fo = [r for r in rxns if r.kind == "falloff"]
    pl = [r for r in rxns if r.kind == "plog"]
    arr_rev = [r for r in arr if r.rev]
    tb_rev = [r for r in tb if r.rev]
    fo_rev = [r for r in fo if r.rev]
    pl_rev = [r for r in pl if r.rev]
    n_arr_fwd, n_arr_rev = len(arr), len(arr_rev)
    n_arr = n_arr_fwd + n_arr_rev
    n_tb, n_tb_rev, n_fo, n_fo_rev, n_plog, n_plog_rev = len(tb), len(tb_rev), len(fo), len(fo_rev), len(pl), len(pl_rev)
    n_m = n_tb + n_fo
    i_arr_end = n_arr
    i_tb_end = i_arr_end + n_tb
    i_tb_bwd_end = i_tb_end + n_tb_rev
    i_fo_end = i_tb_bwd_end + n_fo
    i_fo_bwd_end = i_fo_end + n_fo_rev
    i_plog_start = i_fo_bwd_end + 1
    i_plog_end = i_fo_bwd_end + n_plog
    i_plog_bwd_start = i_plog_end + 1
    i_plog_bwd_end = i_plog_end + n_plog_rev
    nR = i_plog_bwd_end

    entries = []  # (w_idx0, rxn, backward)

    def place(group, start, backward):
        for j, rx in enumerate(group):
            if backward:
                rx.w_rev = start + j
            else:
                rx.w_fwd = start + j
            entries.append((start + j, rx, backward))

    place(arr, 0, False)
    place(arr_rev, n_arr_fwd, True)
    place(tb, i_arr_end, False)
    place(tb_rev, i_tb_end, True)
    place(fo, i_tb_bwd_end, False)
    place(fo_rev, i_fo_end, True)
    place(pl, i_fo_bwd_end, False)
    place(pl_rev, i_plog_end, True)
    for j, rx in enumerate(tb):
        rx.m_idx = j + 1
    for j, rx in enumerate(fo):
        rx.m_idx = n_tb + j + 1
    entries.sort(key=lambda e: e[0])
    assert [e[0] for e in entries] == list(range(nR))

    # ---- Optional Arrhenius fits of the reverse rates ----
    fits = {}
    fit_note = "detailed balance at run time (k_r = k_f/K_c, K_c from the NASA7 table)"
    if fit_reverse is not None:
        T_min, T_max, n_pts = fit_reverse["T_min"], fit_reverse["T_max"], fit_reverse["n_pts"]
        fits = fit_reverse_rates(rxns, mech, T_min, T_max, n_pts)
        if fits:
            errs = np.array([f["err"] for f in fits.values()])
            worst = max(fits, key=lambda rid: fits[rid]["err"])
            fit_note = (f"Arrhenius fits of k_f/K_c over {T_min:.0f}-{T_max:.0f} K "
                        f"(max relative error {100 * errs.max():.2f}%)")
            print(f"  Reverse-rate fits over {T_min:.0f}-{T_max:.0f} K: median max-error {100 * np.median(errs):.2f}%, "
                  f"{int(np.sum(errs > 0.1))} of {len(errs)} reactions above 10%, worst {100 * errs.max():.1f}% "
                  f"(YAML reaction {worst})")
            if errs.max() > 0.25:
                warn(f"Reverse-rate Arrhenius fit error exceeds 25% for reaction {worst}; consider a narrower "
                     f"--fit-range or the default run-time detailed balance")

    # ---- Third-body groups: one M per distinct (default efficiency, explicit efficiencies) set ----
    # Reactions with identical efficiency sets share a group (gri30: 41 reactions -> 10 groups).
    eff_species_names = sorted({sp for rx in tb + fo for sp in rx.eff})
    n_eff = len(eff_species_names)
    eff_species_idx = [species_names.index(sp) + 1 for sp in eff_species_names]
    groups = {}
    for rx in tb + fo:
        key = (rx.default_eff, tuple(rx.eff.get(sp, rx.default_eff) for sp in eff_species_names))
        rx.m_idx = groups.setdefault(key, len(groups) + 1)
    n_m = len(groups)
    eff_matrix_reduced = np.array([list(key[1]) for key in groups], dtype=float).reshape(n_m, n_eff)
    default_eff_arr = np.array([key[0] for key in groups], dtype=float)
    for j, rx in enumerate(tb):
        rx.tb_pos = j + 1
    for j, rx in enumerate(fo):
        rx.fo_pos = j + 1
    rev_all = arr_rev + tb_rev + fo_rev + pl_rev  # reversible reactions in k-layout order
    n_rev = len(rev_all)
    n_kcinv = n_fo_rev + n_plog_rev  # reverse rates evaluated as k_f * kcinv(T)
    # ---- Reaction index file ----
    out_path = Path(output_path)
    stem = out_path.stem
    index_path = out_path.parent / ("chem_data_reactions.txt" if stem == "chem_data_fc" else stem + "_reactions.txt")
    yaml_labels = extract_yaml_labels(yaml_path) if yaml_path else {}
    eq_strs = []
    for w_idx, rx, backward in entries:
        lhs, rhs = (rx.right, rx.left) if backward else (rx.left, rx.right)
        text = f"{format_side(lhs, rx.collider)} => {format_side(rhs, rx.collider)}"
        if backward:
            text += f"  [reverse of w({rx.w_fwd + 1})"
            text += f"; fit err {100 * fits[rx.rid]['err']:.2f}%]" if fits else "]"
        elif rx.rev:
            text += f"  [reversible; reverse in w({rx.w_rev + 1})]"
        eq_strs.append(text)
    idx_width = max(4, len(str(nR)))
    eq_width = max([len(s) for s in eq_strs] + [8])
    with open(index_path, "w") as fidx:
        fidx.write("fcmech reaction index: w(i) = k(i) * product of reactant concentrations [* M]\n")
        fidx.write(f"reverse rates: {fit_note}\n")
        fidx.write("=" * (idx_width + eq_width + 30) + "\n")
        fidx.write(f"{'w':>{idx_width}}  {'Rate':<{eq_width}}  YAML\n")
        fidx.write("-" * (idx_width + eq_width + 30) + "\n")
        for (w_idx, rx, backward), text in zip(entries, eq_strs):
            fidx.write(f"{w_idx + 1:>{idx_width}}  {text:<{eq_width}}  {yaml_labels.get(rx.rid, f'# Reaction {rx.rid}')}\n")
    print(f"  Reaction index: {index_path}")

    # ---- Fortran source ----
    L = []
    yaml_name = Path(yaml_path).name if yaml_path else "(unknown)"
    L += [
        "!--------------------------------------------------------------------------------------------------",
        f"!  FILE {out_path.name}",
        "!  Module for chemical kinetics in NGA2 (finite chemistry)",
        "!  Generated by yaml2nga.py from a YAML mechanism. Do not edit manually.",
        f"!  Source: {yaml_name}",
        f"!  Reverse rates: {fit_note}",
        "!--------------------------------------------------------------------------------------------------",
        "",
        "module fcmech",
        "  use precision",
        "  use string",
        "  implicit none",
        "",
        "  !--------------------------------------------------------------------------------------------------",
        "  !  VARIABLE DEFINITIONS (SI units: m, s, mol, J, K, Pa)",
        "  !--------------------------------------------------------------------------------------------------",
        "  !  Rcst, P_ref    : universal gas constant (J/(mol K)), standard-state pressure of the thermo data (Pa)",
        "  !  nS, nR, nA     : number of species, of one-directional rates w(i), of atom types",
        "  !  sXXX           : species index (1..nS) for species XXX",
        "  !  W_sp           : molar mass (kg/mol), dimension(nS)",
        "  !  atom_names, atom_masses, comp_matrix(a,s) : atoms, their molar masses (kg/mol), count of atom a in species s",
        "  !  koveps         : 1/(eps/k) for Lennard-Jones (1/K), dimension(nS); T* = T*koveps",
        "  !  mucoeff        : viscosity prefactor, mu = mucoeff*sqrt(T)/Omega_mu(T*), dimension(nS)",
        "  !  Dcoeffs        : binary diffusion prefactor, D_ij = Dcoeffs*T^1.5/(P*Omega_D), dimension(nS,nS)",
        "  !  Ocoeffs        : sqrt(koveps_i*koveps_j) giving T* for the pair, dimension(nS,nS)",
        "  !  Cpsp, hsp      : species Cp (J/(mol K)) and enthalpy (J/mol), dimension(nS), set by fcmech_thermodata",
        "  !  T_mid          : NASA7 switch temperature (K), dimension(nS)",
        "  !  thermo_coeffs  : NASA7 coefficients, dimension(nS,14): cols 1-7 low T (a1..a7), cols 8-14 high T",
        "  !  k(1:nR), w(1:nR) layout (every reversible reaction has a forward and a reverse entry):",
        "  !     1..n_arr_fwd             Arrhenius forward          n_arr_fwd+1..i_arr_end   Arrhenius reverse",
        "  !     i_arr_end+1..i_tb_end    three-body forward         ..i_tb_bwd_end           three-body reverse",
        "  !     i_tb_bwd_end+1..i_fo_end falloff forward            ..i_fo_bwd_end           falloff reverse",
        "  !     i_plog_start..i_plog_end PLOG forward               i_plog_bwd_start..i_plog_bwd_end PLOG reverse",
        "  !  Reverse rates: see the header line above (run-time detailed balance, or Arrhenius fits of k_f/K_c",
        "  !     stored after the forward parameters: ln_A_arr(n_arr_fwd+1:n_arr), *_tb_b, *_b_fo, *_b_plog).",
        "  !  ln_A_*, b_*, E_R_*  : Arrhenius parameters ln(A) [SI], temperature exponent, Ea/R (K)",
        "  !  A_sign_*            : sign of A (present only if the mechanism uses negative-A reactions)",
        "  !  Troe_a,Troe_T3,Troe_T1,Troe_T2 : Troe parameters (Lindemann: a=1, T3=T1=1e30; T2=1e30 if absent)",
        "  !  n_m, default_eff, n_eff, eff_species_idx, eff_matrix_reduced : third-body groups",
        "  !     M(g) = default_eff(g)*(sum(c) - sum(c_eff)) + eff_matrix_reduced(g,:).c_eff; reactions with identical",
        "  !     efficiency sets share a group (fo_m_idx(i) = group of falloff reaction i)",
        "  !  P_plog, ln_P_plog, ln_A_plog, b_plog, E_R_plog, plog_reac_start, plog_reac_len : PLOG tables (P in Pa)",
        "  !  n_rev, rev_fwd_w, rev_bwd_w : reversible reactions and the w indices of their forward/reverse rates;",
        "  !     get_production_rates sums them as net rates w(fwd) - w(rev); get_destruction_rates sums the",
        "  !     one-directional rates consuming each species (fcmech_get_ydot_ddot returns both in mass units)",
        "  !  fo_rev_src, plog_rev_src, kcinv_cache : falloff/PLOG reverse rates are k_f(src) * kcinv(T)",
        "  !  *_cache : temperature-only data (k(T), falloff limits, F_cent, NASA7) kept from the previous call and",
        "  !     reused while T (and P for PLOG) is unchanged; exact comparison; threadprivate under OpenMP",
        "  !  fcmech_vexp : vectorizable exp used for all rate-coefficient exponentials (relative error < 1e-14)",
        "  !--------------------------------------------------------------------------------------------------",
        "",
        f"  real(WP), parameter :: Rcst = {f90_real(R_UNIV)}",
        f"  real(WP), parameter :: P_ref = {f90_real(P_REF)}",
        f"  integer, parameter :: nS = {nS}",
        f"  integer, parameter :: nR = {nR}",
        f"  integer, parameter :: nA = {nA}",
        "",
        "  ! --- Species indices ---",
    ]
    for i, nm in enumerate(species_names):
        L.append(f"  integer, parameter :: {f90n[nm]} = {i + 1}")
    L.append("")
    L.append("  ! --- Molar masses (kg/mol) ---")
    L += f90_array_parameter("W_sp", mech["W_sp"], comments=species_names, dim_expr="nS")
    L.append("")
    if nA > 0:
        L.append("  ! --- Atoms: names, molar masses (kg/mol), composition comp_matrix(a,s) ---")
        L += f90_array_parameter("atom_names", [f"'{a[:2].ljust(2)}'" for a in mech["atom_names"]],
                                 kind="character(len=2)", fmt=str, per_line=8, dim_expr="nA")
        L += f90_array_parameter("atom_masses", mech["atom_masses"], comments=mech["atom_names"], dim_expr="nA")
        L += f90_2d_array_parameter("comp_matrix", mech["comp_matrix"], nA, nS, kind="integer", fmt=f90_int,
                                    per_line=8, dim_expr="nA,nS")
        L.append("")
    if transport is not None:
        L.append("  ! --- Transport (Lennard-Jones / Chapman-Enskog, Neufeld collision integrals) ---")
        L += f90_array_parameter("koveps", transport["koveps"], comments=species_names, dim_expr="nS")
        L += f90_array_parameter("mucoeff", transport["mucoeff"], comments=species_names, dim_expr="nS")
        L.append("  ! nS x nS pair tables: DATA-initialized protected variables (a parameter array constructor is")
        L.append("  ! limited to 65535 elements, i.e. 255 species)")
        L += f90_2d_data_variable("Dcoeffs", transport["Dcoeffs"], nS, nS, dim_expr="nS,nS")
        L += f90_2d_data_variable("Ocoeffs", transport["Ocoeffs"], nS, nS, dim_expr="nS,nS")
        L.append("")
    L.append("  ! --- Thermodynamics: Cp and enthalpy (module variables, set by fcmech_thermodata) ---")
    L.append("  real(WP), dimension(nS) :: Cpsp, hsp")
    L.append("")
    L.append("  ! --- NASA7 thermo: switch temperature (K) and coefficients (low cols 1-7, high cols 8-14) ---")
    L += f90_array_parameter("T_mid", T_mid_sp, comments=species_names, dim_expr="nS")
    L += f90_2d_array_parameter("thermo_coeffs", thermo_coeffs, nS, 14, dim_expr="nS,14")
    L.append("")

    L.append("  ! --- Reaction counts and k/w index boundaries ---")
    L += [
        f"  integer, parameter :: n_arr_fwd = {n_arr_fwd}, n_arr_rev = {n_arr_rev}, n_arr = {n_arr}",
        f"  integer, parameter :: n_tb = {n_tb}, n_tb_rev = {n_tb_rev}",
        f"  integer, parameter :: n_fo = {n_fo}, n_fo_rev = {n_fo_rev}",
        f"  integer, parameter :: n_plog = {n_plog}, n_plog_rev = {n_plog_rev}, n_plog_pts = {sum(len(r.plog) for r in pl)}",
        f"  integer, parameter :: n_m = {n_m}  ! third-body/falloff M groups",
        f"  integer, parameter :: i_arr_end = {i_arr_end}, i_tb_end = {i_tb_end}, i_tb_bwd_end = {i_tb_bwd_end}",
        f"  integer, parameter :: i_fo_end = {i_fo_end}, i_fo_bwd_end = {i_fo_bwd_end}",
        f"  integer, parameter :: i_plog_start = {i_plog_start}, i_plog_end = {i_plog_end}",
        f"  integer, parameter :: i_plog_bwd_start = {i_plog_bwd_start}, i_plog_bwd_end = {i_plog_bwd_end}",
        "",
    ]
    # Arrhenius block: forward parameters, plus the fitted reverse ones (n_arr entries) with --fit-reverse.
    # Signs of A are applied after the exponential, so one sign array covers forward and reverse entries.
    def rev_params(r):
        """Fitted reverse Arrhenius parameters of an elementary/three-body reaction: forward + fit of 1/K_c."""
        f = fits[r.rid]
        return math.log(r.A) + f["lnA"], r.b + f["b"], r.Ea / R_UNIV + f["E_R"]

    arr_k = arr + arr_rev if fits else arr
    n_arr_k = "n_arr" if fits else "n_arr_fwd"
    any_negative_arr = any(r.sign < 0 for r in arr + arr_rev)
    any_negative_tb = any(r.sign < 0 for r in tb + tb_rev)
    if arr_k:
        L.append("  ! --- Arrhenius: ln(A), b, E/R" + (" (forward 1..n_arr_fwd, fitted reverse after) ---" if fits else " ---"))
        L += f90_array_parameter("ln_A_arr", [math.log(r.A) for r in arr] + [rev_params(r)[0] for r in arr_rev if fits], dim_expr=n_arr_k)
        L += f90_array_parameter("b_arr", [r.b for r in arr] + [rev_params(r)[1] for r in arr_rev if fits], dim_expr=n_arr_k)
        L += f90_array_parameter("E_R_arr", [r.Ea / R_UNIV for r in arr] + [rev_params(r)[2] for r in arr_rev if fits], dim_expr=n_arr_k)
        L.append("")
    if any_negative_arr:
        L.append("  ! --- Sign of A for every Arrhenius rate (forward, then reverse) ---")
        L += f90_array_parameter("A_sign_arr", [r.sign for r in arr] + [r.sign for r in arr_rev], dim_expr="n_arr")
        L.append("")
    if n_tb > 0:
        L.append("  ! --- Three-body: ln(A), b, E/R ---")
        L += f90_array_parameter("ln_A_tb", [math.log(r.A) for r in tb], dim_expr="n_tb")
        L += f90_array_parameter("b_tb", [r.b for r in tb], dim_expr="n_tb")
        L += f90_array_parameter("E_R_tb", [r.Ea / R_UNIV for r in tb], dim_expr="n_tb")
        L.append("")
    if fits and n_tb_rev > 0:
        L.append("  ! --- Three-body fitted reverse: ln(A), b, E/R ---")
        L += f90_array_parameter("ln_A_tb_b", [rev_params(r)[0] for r in tb_rev], dim_expr="n_tb_rev")
        L += f90_array_parameter("b_tb_b", [rev_params(r)[1] for r in tb_rev], dim_expr="n_tb_rev")
        L += f90_array_parameter("E_R_tb_b", [rev_params(r)[2] for r in tb_rev], dim_expr="n_tb_rev")
        L.append("")
    if any_negative_tb:
        L.append("  ! --- Sign of A for every three-body rate (forward, then reverse) ---")
        L += f90_array_parameter("A_sign_tb", [r.sign for r in tb] + [r.sign for r in tb_rev], dim_expr="n_tb+n_tb_rev")
        L.append("")
    if n_fo > 0:
        L.append("  ! --- Falloff: low-P (ln_A0,b0,E_R0), high-P (ln_Ainf,binf,E_Rinf), Troe; fo_m_idx = third-body group ---")
        L += f90_array_parameter("ln_A0_fo", [math.log(r.A0) for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("b0_fo", [r.b0 for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("E_R0_fo", [r.Ea0 / R_UNIV for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("ln_Ainf_fo", [math.log(r.Ainf) for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("binf_fo", [r.binf for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("E_Rinf_fo", [r.Eainf / R_UNIV for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("Troe_a", [r.troe[0] for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("Troe_T3", [r.troe[1] for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("Troe_T1", [r.troe[2] for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("Troe_T2", [r.troe[3] for r in fo], dim_expr="n_fo")
        L += f90_array_parameter("fo_m_idx", [r.m_idx for r in fo], kind="integer", fmt=f90_int, per_line=8, dim_expr="n_fo")
        L.append("")
    if n_fo_rev > 0:
        L.append("  ! --- Falloff reverse: fo_rev_src(j) = falloff reaction of reverse rate j (k_r = k_f * kcinv) ---")
        L += f90_array_parameter("fo_rev_src", [r.fo_pos for r in fo_rev], kind="integer", fmt=f90_int, per_line=8, dim_expr="n_fo_rev")
        if fits:
            L += f90_array_parameter("rev_lnA_fo", [fits[r.rid]["lnA"] for r in fo_rev], dim_expr="n_fo_rev")
            L += f90_array_parameter("rev_b_fo", [fits[r.rid]["b"] for r in fo_rev], dim_expr="n_fo_rev")
            L += f90_array_parameter("rev_E_R_fo", [fits[r.rid]["E_R"] for r in fo_rev], dim_expr="n_fo_rev")
        L.append("")
    if n_m > 0:
        L.append("  ! --- Third-body efficiencies: species with explicit efficiencies ---")
        L.append(f"  integer, parameter :: n_eff = {n_eff}")
        if n_eff > 0:
            L += f90_array_parameter("eff_species_idx", eff_species_idx, kind="integer", fmt=f90_int, per_line=8,
                                     comments=eff_species_names, dim_expr="n_eff")
            L += f90_2d_array_parameter("eff_matrix_reduced", eff_matrix_reduced, n_m, n_eff, dim_expr="n_m,n_eff")
        L += f90_array_parameter("default_eff", default_eff_arr, dim_expr="n_m")
        L.append("")
    if n_plog > 0:
        P_plog, lnA_plog, b_plog, ER_plog, starts, lens = [], [], [], [], [], []
        for rx in pl:
            starts.append(len(P_plog) + 1)
            lens.append(len(rx.plog))
            for P, A, b, Ea in rx.plog:
                P_plog.append(P)
                lnA_plog.append(math.log(A))
                b_plog.append(b)
                ER_plog.append(Ea / R_UNIV)
        L.append("  ! --- PLOG: tabulated pressures (Pa and ln Pa, ascending per reaction) and Arrhenius parameters ---")
        L += f90_array_parameter("P_plog", P_plog, dim_expr="n_plog_pts")
        L += f90_array_parameter("ln_P_plog", [math.log(P) for P in P_plog], dim_expr="n_plog_pts")
        L += f90_array_parameter("ln_A_plog", lnA_plog, dim_expr="n_plog_pts")
        L += f90_array_parameter("b_plog", b_plog, dim_expr="n_plog_pts")
        L += f90_array_parameter("E_R_plog", ER_plog, dim_expr="n_plog_pts")
        L += f90_array_parameter("plog_reac_start", starts, kind="integer", fmt=f90_int, per_line=8, dim_expr="n_plog")
        L += f90_array_parameter("plog_reac_len", lens, kind="integer", fmt=f90_int, per_line=8, dim_expr="n_plog")
        if n_plog_rev > 0:
            L.append("  ! --- PLOG reverse: plog_rev_src(j) = PLOG reaction of reverse rate j (k_r = k_f * kcinv) ---")
            L += f90_array_parameter("plog_rev_src", [pl.index(r) + 1 for r in pl_rev], kind="integer", fmt=f90_int, per_line=8,
                                     dim_expr="n_plog_rev")
            if fits:
                L += f90_array_parameter("rev_lnA_plog", [fits[r.rid]["lnA"] for r in pl_rev], dim_expr="n_plog_rev")
                L += f90_array_parameter("rev_b_plog", [fits[r.rid]["b"] for r in pl_rev], dim_expr="n_plog_rev")
                L += f90_array_parameter("rev_E_R_plog", [fits[r.rid]["E_R"] for r in pl_rev], dim_expr="n_plog_rev")
        L.append("")
    L.append(f"  integer, parameter :: n_rev = {n_rev}  ! reversible reactions")
    L.append(f"  integer, parameter :: n_kcinv = {n_kcinv}  ! falloff + PLOG reverse rates evaluated as k_f * kcinv(T)")
    if n_rev > 0:
        L.append("  ! --- w indices of the forward and reverse rate of each reversible reaction (net rate = w(fwd) - w(rev)) ---")
        L += f90_array_parameter("rev_fwd_w", [rx.w_fwd + 1 for rx in rev_all], kind="integer", fmt=f90_int, per_line=8, dim_expr="n_rev")
        L += f90_array_parameter("rev_bwd_w", [rx.w_rev + 1 for rx in rev_all], kind="integer", fmt=f90_int, per_line=8, dim_expr="n_rev")
    L.append("")
    L += [
        "  ! --- Caches of temperature-dependent data, reused by calls at the same T (Newton iterations, finite-",
        "  !     difference Jacobian columns that perturb the composition only); PLOG data are keyed on (T,P).",
        "  !     The exact comparisons on T and P are intentional. Serial use, or per-thread with OpenMP through",
        "  !     the threadprivate directive below. ---",
        "  real(WP), save :: T_nasa_cache = -1.0_WP",
        "  real(WP), dimension(nS), save :: cp_R_cache, h_RT_cache, s_R_cache",
        "  real(WP), save :: T_rate_cache = -1.0_WP",
        "  real(WP), save :: P_rate_cache = -1.0_WP",
        "  real(WP), dimension(nR), save :: k_T_cache  ! Arrhenius/three-body (forward and reverse) and PLOG forward k",
    ]
    threadprivate = ["T_nasa_cache", "cp_R_cache", "h_RT_cache", "s_R_cache", "T_rate_cache", "P_rate_cache", "k_T_cache"]
    if n_fo > 0:
        L.append("  real(WP), dimension(n_fo), save :: k0_cache, kinf_cache, fc_cache  ! falloff limits and F_cent")
        threadprivate += ["k0_cache", "kinf_cache", "fc_cache"]
    if n_kcinv > 0:
        L.append("  real(WP), dimension(n_kcinv), save :: kcinv_cache  ! k_r/k_f of the falloff and PLOG reverse rates")
        threadprivate.append("kcinv_cache")
    for c in range(0, len(threadprivate), 4):  # several directives keep every line under 132 columns
        L.append("  !$omp threadprivate(" + ", ".join(threadprivate[c:c + 4]) + ")")
    L.append("")
    L.append("  contains")
    L.append("")
    # ---- Collision integrals ----
    def neufeld_function(name, coeffs, what):
        body = [
            f"  ! Reduced collision integral {what} (Neufeld et al. 1972 rational fit)",
            f"  real(WP) function {name}(T_)",
            "    implicit none",
            "    real(WP), intent(in) :: T_  ! reduced temperature T* = T*k/eps",
            "    real(WP), parameter, dimension(9) :: mArray = [ &",
            f"         {f90_real(coeffs[0])}, {f90_real(coeffs[1])}, {f90_real(coeffs[2])}, &",
            f"         {f90_real(coeffs[3])}, {f90_real(coeffs[4])}, {f90_real(coeffs[5])}, &",
            f"         {f90_real(coeffs[6])}, {f90_real(coeffs[7])}, {f90_real(coeffs[8])} ]",
            "    integer :: i",
            "    real(WP) :: num, den",
            "    den = mArray(9)",
            "    do i = 1, 4",
            "      den = mArray(9-i) + T_*den",
            "    end do",
            "    num = mArray(4)",
            "    do i = 1, 3",
            "      num = mArray(4-i) + T_*num",
            "    end do",
            f"    {name} = num/den",
            f"  end function {name}",
            "",
        ]
        return body

    if transport is not None:
        L += neufeld_function("fcmech_omegamu", NEUFELD_OMEGA_MU, "Omega^(2,2)* for viscosity")
        L += neufeld_function("fcmech_omegaD", NEUFELD_OMEGA_D, "Omega^(1,1)* for diffusion")

    # ---- Falloff blending ----
    L += [
        "  !--------------------------------------------------------------------------------------------------",
        "  !  function getlindratecoeff: pressure-dependent rate coefficient (Lindemann / Troe blending)",
        "  !  k = kinf * F * Pr/(1+Pr), Pr = k0*M/kinf, log10 F = log10 Fcent / (1 + ((log10 Pr + c)/(n - d*(log10 Pr + c)))^2)",
        "  !--------------------------------------------------------------------------------------------------",
        "  real(WP) function getlindratecoeff(k0,kinf,fc,conc)",
        "    implicit none",
        "    ! k0, kinf : low-P and high-P rate coefficients (> 0)",
        "    ! fc       : Troe centering factor F_cent (1 for Lindemann); clamped to a tiny positive value like Cantera",
        "    ! conc     : third-body concentration M (mol/m^3); M <= 0 gives k = 0 (no fallback to P/(R*T))",
        "    real(WP), intent(in) :: k0,kinf,fc,conc",
        "    real(WP), parameter :: ln10 = 2.302585092994046_WP",
        "    real(WP) :: redP,lgfc,ntmp,ccoeff,lgPr,f",
        "    redP = k0 * MAX(conc, 0.0_WP) / MAX(kinf, TINY(1.0_WP))",
        "    if (redP <= 0.0_WP) then",
        "      getlindratecoeff = 0.0_WP",
        "      return",
        "    end if",
        "    lgfc = LOG10(MAX(fc, 1.0e-300_WP))  ! same clamp as Cantera (Fcent <= 0 is possible with odd Troe data)",
        "    ntmp = 0.75_WP - 1.27_WP*lgfc",
        "    ccoeff = -0.4_WP - 0.67_WP*lgfc",
        "    lgPr = LOG10(redP)",
        "    f = (lgPr + ccoeff)/(ntmp - 0.14_WP*(lgPr + ccoeff))",
        "    f = exp(ln10*lgfc/(f*f + 1.0_WP))",
        "    getlindratecoeff = kinf * f * redP/(1.0_WP + redP)",
        "  end function getlindratecoeff",
        "",
    ]

    # ---- Third bodies ----
    L += [
        "  ! Effective third-body concentrations M(r) = default_eff(r)*sum_{i not in eff}(c_i) + sum_{j in eff} eff(r,j)*c_j",
        "  subroutine get_thirdbodies(M,c)",
        "    implicit none",
        "    real(WP), dimension(nS), intent(in) :: c   ! species concentrations (mol/m^3)",
        "    real(WP), dimension(n_m), intent(out) :: M  ! third-body concentration per reaction (mol/m^3)",
    ]
    if n_m > 0 and n_eff > 0:
        L += [
            "    real(WP), dimension(n_eff) :: c_eff",
            "    real(WP) :: sum_c_rest",
            "    integer :: j",
            "    do j = 1, n_eff",
            "      c_eff(j) = c(eff_species_idx(j))",
            "    end do",
            "    sum_c_rest = sum(c) - sum(c_eff)",
            "    M(:) = default_eff(:) * sum_c_rest + matmul(eff_matrix_reduced, c_eff)",
        ]
    elif n_m > 0:
        L.append("    M(:) = default_eff(:) * sum(c)")
    else:
        L.append("    ! no third-body or falloff reactions in this mechanism (n_m = 0)")
    L += ["  end subroutine get_thirdbodies", ""]

    # ---- Vectorizable exponential ----
    L += [
        "  !--------------------------------------------------------------------------------------------------",
        "  !  subroutine fcmech_vexp: y = exp(x) for whole arrays, written so the compiler can SIMD-vectorize",
        "  !  it (libm exp is scalar): Cody-Waite range reduction, degree-11 Taylor polynomial, exponent",
        "  !  assembly. Relative error < 1e-14. Arguments are clamped to [-708, 709], so there is never an",
        "  !  Inf/NaN; results below ~1e-307 are inexact (rates that are zero for any practical purpose).",
        "  !--------------------------------------------------------------------------------------------------",
        "  pure subroutine fcmech_vexp(x, y)",
        "    implicit none",
        "    real(WP), dimension(:), intent(in) :: x",
        "    real(WP), dimension(:), intent(out) :: y",
        "    integer, parameter :: dp = selected_real_kind(15, 307)  ! the algorithm needs IEEE double internally",
        "    real(dp), parameter :: ln2_hi = 6.93147180369123816490e-01_dp, ln2_lo = 1.90821492927058770002e-10_dp",
        "    real(dp), parameter :: inv_ln2 = 1.44269504088896338700_dp",
        "    real(dp) :: xc, r, p",
        "    integer(8) :: n",
        "    integer :: i",
        "    do i = 1, size(x)",
        "      xc = min(max(real(x(i), dp), -708.0_dp), 709.0_dp)",
        "      n = nint(xc*inv_ln2, 8)",
        "      r = (xc - real(n, dp)*ln2_hi) - real(n, dp)*ln2_lo",
        "      p = 1.0_dp + r*(1.0_dp + r*(0.5_dp + r*(1.0_dp/6.0_dp + r*(1.0_dp/24.0_dp + r*(1.0_dp/120.0_dp &",
        "        + r*(1.0_dp/720.0_dp + r*(1.0_dp/5040.0_dp + r*(1.0_dp/40320.0_dp + r*(1.0_dp/362880.0_dp &",
        "        + r*(1.0_dp/3628800.0_dp + r*(1.0_dp/39916800.0_dp)))))))))))",
        "      y(i) = real(p*transfer(ishft(n + 1023_8, 52), 1.0_dp), WP)",
        "    end do",
        "  end subroutine fcmech_vexp",
        "",
    ]

    # ---- Rate coefficients ----
    runtime_kc = n_rev > 0 and not fits
    L += [
        "  !--------------------------------------------------------------------------------------------------",
        "  !  subroutine get_rate_coefficients: k(1:nR) at (T,P), forward and reverse. Everything that depends",
        "  !  on T only (and on P for PLOG) comes from the caches, refreshed when T (or P) changes; only the",
        "  !  falloff blending with the third-body concentration is evaluated on every call.",
        "  !--------------------------------------------------------------------------------------------------",
        "  subroutine get_rate_coefficients(k,M,Tloc,Ploc)",
        "    implicit none",
        "    real(WP), dimension(nR), intent(out) :: k   ! rate coefficients ((m^3/mol)^(n-1)/s)",
        "    real(WP), dimension(:), intent(in) :: M     ! third-body concentrations from get_thirdbodies",
        "    real(WP), intent(in) :: Tloc, Ploc          ! temperature (K), pressure (Pa)",
    ]
    if n_fo > 0 or n_kcinv > 0:
        L.append("    integer :: i")
    L.append("    if (Tloc /= T_rate_cache) then")
    L.append("      call update_rate_cache(Tloc)")
    if n_plog > 0:
        L.append("      P_rate_cache = -1.0_WP")
    L.append("    end if")
    if n_plog > 0:
        L.append("    if (Ploc /= P_rate_cache) call update_plog_cache(Tloc, Ploc)")
    if i_tb_bwd_end > 0:
        L.append("    k(1:i_tb_bwd_end) = k_T_cache(1:i_tb_bwd_end)")
    if n_fo > 0:
        L += [
            "    ! Falloff: blend the cached low/high-pressure limits with the current third-body concentration",
            "    do i = 1, n_fo",
            "      k(i_tb_bwd_end+i) = getlindratecoeff(k0_cache(i), kinf_cache(i), fc_cache(i), M(fo_m_idx(i)))",
            "    end do",
        ]
    if n_fo_rev > 0:
        L += [
            "    do i = 1, n_fo_rev",
            "      k(i_fo_end+i) = k(i_tb_bwd_end+fo_rev_src(i)) * kcinv_cache(i)",
            "    end do",
        ]
    if n_plog > 0:
        L.append("    k(i_plog_start:i_plog_end) = k_T_cache(i_plog_start:i_plog_end)")
    if n_plog_rev > 0:
        L += [
            "    do i = 1, n_plog_rev",
            "      k(i_plog_bwd_start+i-1) = k(i_plog_start+plog_rev_src(i)-1) * kcinv_cache(n_fo_rev+i)",
            "    end do",
        ]
    L += ["  end subroutine get_rate_coefficients", ""]

    # ---- Temperature-only cache refresh ----
    L += [
        "  ! Refresh every rate quantity that depends on the temperature only",
        "  subroutine update_rate_cache(Tloc)",
        "    implicit none",
        "    real(WP), intent(in) :: Tloc",
        "    real(WP) :: T_log, T_inv",
    ]
    if runtime_kc:
        L.append("    real(WP) :: ln_P0RT")
        L.append("    real(WP), dimension(nS) :: cp_R, h_RT, s_R, g_RT")
    if n_arr > 0:
        L.append("    real(WP), dimension(n_arr) :: ln_k_arr")
    if n_tb + n_tb_rev > 0:
        L.append("    real(WP), dimension(n_tb+n_tb_rev) :: ln_k_tb")
    if n_fo > 0:
        L.append("    real(WP), dimension(n_fo) :: ln_k0, ln_kinf, e1, e2, e3")
    if n_kcinv > 0:
        L.append("    real(WP), dimension(n_kcinv) :: ln_kcinv")
    L += ["    T_log = log(Tloc)", "    T_inv = 1.0_WP/Tloc"]
    if n_arr_fwd > 0:
        L.append("    ! Arrhenius" + (" (forward and fitted reverse)" if fits else " forward") + ": ln(k) = ln(A) + b*ln(T) - E/(R*T)")
        L.append(f"    ln_k_arr(1:{n_arr_k}) = ln_A_arr(:) + b_arr(:)*T_log - E_R_arr(:)*T_inv")
    if n_tb > 0:
        L.append("    ln_k_tb(1:n_tb) = ln_A_tb(:) + b_tb(:)*T_log - E_R_tb(:)*T_inv")
    if fits and n_tb_rev > 0:
        L.append("    ln_k_tb(n_tb+1:n_tb+n_tb_rev) = ln_A_tb_b(:) + b_tb_b(:)*T_log - E_R_tb_b(:)*T_inv")
    if n_fo > 0:
        L.append("    ln_k0(:) = ln_A0_fo(:) + b0_fo(:)*T_log - E_R0_fo(:)*T_inv")
        L.append("    ln_kinf(:) = ln_Ainf_fo(:) + binf_fo(:)*T_log - E_Rinf_fo(:)*T_inv")
    if fits and n_fo_rev > 0:
        L.append("    ln_kcinv(1:n_fo_rev) = rev_lnA_fo(:) + rev_b_fo(:)*T_log - rev_E_R_fo(:)*T_inv")
    if fits and n_plog_rev > 0:
        L.append("    ln_kcinv(n_fo_rev+1:n_kcinv) = rev_lnA_plog(:) + rev_b_plog(:)*T_log - rev_E_R_plog(:)*T_inv")
    if runtime_kc:
        L += [
            "    ! Reverse rates from detailed balance: ln k_r = ln k_f - ln K_c, with",
            "    ! ln K_c = -sum_i nu_i*g_i/(R*T) + dnu*ln(P_ref/(R*T)) and g/(RT) = h/(RT) - s/R from NASA7",
            "    call fcmech_nasa7(Tloc, cp_R, h_RT, s_R)",
            "    g_RT(:) = h_RT(:) - s_R(:)",
            "    ln_P0RT = log(P_ref*T_inv/Rcst)",
        ]

        def kc_expr_terms(rx, first=None):
            terms = _kc_terms(rx, f90n, species_names)
            if first is not None:
                return [first] + terms
            if not terms:
                return ["0.0_WP"]
            return [terms[0][2:] if terms[0].startswith("+ ") else terms[0]] + terms[1:]

        for j, rx in enumerate(arr_rev):
            L.append(f"    ! {rx.label()}")
            L += f90_wrap(f"    ln_k_arr({n_arr_fwd + j + 1}) = ", kc_expr_terms(rx, f"ln_k_arr({rx.w_fwd + 1})"), indent="        ")
        for j, rx in enumerate(tb_rev):
            L.append(f"    ! {rx.label()}")
            L += f90_wrap(f"    ln_k_tb({n_tb + j + 1}) = ", kc_expr_terms(rx, f"ln_k_tb({rx.tb_pos})"), indent="        ")
        for j, rx in enumerate(fo_rev + pl_rev):
            L.append(f"    ! {rx.label()}")
            L += f90_wrap(f"    ln_kcinv({j + 1}) = ", kc_expr_terms(rx), indent="        ")
    L.append("    ! Exponentials (vectorized); signs of negative pre-exponential factors are applied afterwards")
    if n_arr > 0:
        L.append("    call fcmech_vexp(ln_k_arr, k_T_cache(1:n_arr))")
        if any_negative_arr:
            L.append("    k_T_cache(1:n_arr) = A_sign_arr(:) * k_T_cache(1:n_arr)")
    if n_tb + n_tb_rev > 0:
        L.append("    call fcmech_vexp(ln_k_tb, k_T_cache(i_arr_end+1:i_tb_bwd_end))")
        if any_negative_tb:
            L.append("    k_T_cache(i_arr_end+1:i_tb_bwd_end) = A_sign_tb(:) * k_T_cache(i_arr_end+1:i_tb_bwd_end)")
    if n_fo > 0:
        L += [
            "    call fcmech_vexp(ln_k0, k0_cache)",
            "    call fcmech_vexp(ln_kinf, kinf_cache)",
            "    ! Troe centering F_cent = (1-a)*exp(-T/T3) + a*exp(-T/T1) + exp(-T2/T); T2 = 1e30 when absent (term vanishes)",
            "    call fcmech_vexp(-Tloc/Troe_T3(:), e1)",
            "    call fcmech_vexp(-Tloc/Troe_T1(:), e2)",
            "    call fcmech_vexp(-Troe_T2(:)*T_inv, e3)",
            "    fc_cache(:) = (1.0_WP - Troe_a(:))*e1(:) + Troe_a(:)*e2(:) + e3(:)",
        ]
    if n_kcinv > 0:
        L.append("    call fcmech_vexp(ln_kcinv, kcinv_cache)")
    L += ["    T_rate_cache = Tloc", "  end subroutine update_rate_cache", ""]

    # ---- PLOG (forward; depends on T and P) ----
    if n_plog > 0:
        L += [
            "  !--------------------------------------------------------------------------------------------------",
            "  !  subroutine update_plog_cache: PLOG forward rates at (T,P), log-log interpolation in pressure.",
            "  !  Entries tabulated at the same pressure are summed; outside the table the edge value is used.",
            "  !  The bracketing entries are found by index (i_inf/i_sup point at the last entry of each level).",
            "  !--------------------------------------------------------------------------------------------------",
            "  subroutine update_plog_cache(Tloc, Ploc)",
            "    implicit none",
            "    real(WP), intent(in) :: Tloc, Ploc",
            "    real(WP), parameter :: k_floor = 1.0e-300_WP  ! keeps log() finite when a level underflows (Cantera does the same)",
            "    real(WP) :: T_log, T_inv, P_log, k_inf, k_sup, lnk, w_lnP",
            "    integer :: j, is, ie, idx, i_inf, i_sup",
            "    T_log = log(Tloc)",
            "    T_inv = 1.0_WP/Tloc",
            "    P_log = log(Ploc)",
            "    do j = 1, n_plog",
            "      is = plog_reac_start(j)",
            "      ie = is + plog_reac_len(j) - 1",
            "      ! Last entry with P <= Ploc (or the first entry) and first entry with P >= Ploc (or the last entry)",
            "      i_inf = is",
            "      do idx = is, ie",
            "        if (ln_P_plog(idx) <= P_log) i_inf = idx",
            "      end do",
            "      i_sup = ie",
            "      do idx = ie, is, -1",
            "        if (ln_P_plog(idx) >= P_log) i_sup = idx",
            "      end do",
            "      ! Sum the entries tabulated at each bracketing pressure",
            "      k_inf = k_floor",
            "      do idx = i_inf, is, -1",
            "        if (ln_P_plog(idx) < ln_P_plog(i_inf)) exit",
            "        lnk = ln_A_plog(idx) + b_plog(idx)*T_log - E_R_plog(idx)*T_inv",
            "        k_inf = k_inf + exp(lnk)",
            "      end do",
            "      if (ln_P_plog(i_sup) > ln_P_plog(i_inf)) then",
            "        k_sup = k_floor",
            "        do idx = i_sup, ie",
            "          if (ln_P_plog(idx) > ln_P_plog(i_sup)) exit",
            "          lnk = ln_A_plog(idx) + b_plog(idx)*T_log - E_R_plog(idx)*T_inv",
            "          k_sup = k_sup + exp(lnk)",
            "        end do",
            "        w_lnP = (P_log - ln_P_plog(i_inf))/(ln_P_plog(i_sup) - ln_P_plog(i_inf))",
            "        k_T_cache(i_plog_start+j-1) = exp(log(k_inf) + (log(k_sup) - log(k_inf))*w_lnP)",
            "      else",
            "        ! exact hit on a tabulated pressure, or outside the table",
            "        k_T_cache(i_plog_start+j-1) = k_inf",
            "      end if",
            "    end do",
            "    P_rate_cache = Ploc",
            "  end subroutine update_plog_cache",
            "",
        ]
    # ---- Reaction rates w(i) ----
    L += [
        "  ! Reaction rates w(i) = k(i) * product of reactant concentrations [* M for three-body reactions]",
        "  subroutine get_reaction_rates(w,k,m,c)",
        "    implicit none",
        "    real(WP), dimension(nS), intent(in) :: c   ! species concentrations (mol/m^3)",
        "    real(WP), dimension(nR), intent(in) :: k   ! rate coefficients",
        "    real(WP), dimension(:), intent(in) :: m    ! third-body concentrations from get_thirdbodies",
        "    real(WP), dimension(nR), intent(out) :: w  ! reaction rates (mol/(m^3 s))",
    ]
    for (w_idx, rx, backward), text in zip(entries, eq_strs):
        side = rx.right if backward else rx.left
        terms = ["* " + t if k == 0 else t for k, t in enumerate(_reactant_terms(side, f90n))]
        if rx.kind == "three-body":
            terms.append(f"* m({rx.m_idx})")
        L.append(f"    ! {text}")
        L += f90_wrap(f"    w({w_idx + 1}) = k({w_idx + 1})", terms, indent="        ")
    L += ["  end subroutine get_reaction_rates", ""]

    # ---- Production rates ----
    L += [
        "  ! Species production rates cdot(i) = sum_r nu(i,r)*w(r); reversible reactions enter once, as net rates",
        "  subroutine get_production_rates(cdot,w)",
        "    implicit none",
        "    real(WP), dimension(nR), intent(in) :: w      ! reaction rates (mol/(m^3 s))",
        "    real(WP), dimension(nS), intent(out) :: cdot  ! production rates (mol/(m^3 s))",
    ]
    if n_rev > 0:
        L.append("    real(WP), dimension(n_rev) :: wn  ! net rates w(forward) - w(reverse) of the reversible reactions")
        L.append("    wn(:) = w(rev_fwd_w(:)) - w(rev_bwd_w(:))")
    rev_pos = {id(rx): j + 1 for j, rx in enumerate(rev_all)}
    cdot_terms = [[] for _ in range(nS)]
    for w_idx, rx, backward in entries:
        if backward:
            continue  # accounted for through the net rate
        ref = f"wn({rev_pos[id(rx)]})" if rx.rev else f"w({w_idx + 1})"
        for idx, nu in rx.nu.items():
            cdot_terms[idx].append((ref, nu))
    for i in range(nS):
        lhs = f"    cdot({f90n[species_names[i]]})"
        if not cdot_terms[i]:
            L.append(f"{lhs} = 0.0_WP")
            continue
        L += f90_sum_statements(lhs, [_signed_term(nu, ref) for ref, nu in cdot_terms[i]], indent="        ")
    L += ["  end subroutine get_production_rates", ""]

    # ---- Destruction rates ----
    L += [
        "  ! Species destruction rates cdes(i) = sum over the one-directional rates consuming species i of |nu|*w(r)",
        "  ! (>= 0). -cdes(i)/c(i) approximates the diagonal of the chemical Jacobian: the stiffness measure that",
        "  ! bounds the stable time step of an explicit integration of the chemical source terms.",
        "  subroutine get_destruction_rates(cdes,w)",
        "    implicit none",
        "    real(WP), dimension(nR), intent(in) :: w      ! reaction rates (mol/(m^3 s))",
        "    real(WP), dimension(nS), intent(out) :: cdes  ! destruction rates (mol/(m^3 s))",
    ]
    name_to_idx = {nm: i for i, nm in enumerate(species_names)}
    cdes_terms = [[] for _ in range(nS)]
    for w_idx, rx, backward in entries:
        side = rx.right if backward else rx.left
        for sp, c in side:
            cdes_terms[name_to_idx[sp]].append((f"w({w_idx + 1})", c))
    for i in range(nS):
        lhs = f"    cdes({f90n[species_names[i]]})"
        if not cdes_terms[i]:
            L.append(f"{lhs} = 0.0_WP")
            continue
        L += f90_sum_statements(lhs, [_signed_term(c, ref) for ref, c in cdes_terms[i]], indent="        ")
    L += ["  end subroutine get_destruction_rates", ""]
    # ---- Concentrations, thermo, names ----
    L += [
        "  ! Molar concentrations c = rho*y/W from mass fractions, ideal gas rho = P*W_mix/(R*T)",
        "  subroutine y2c(y, W_sp, P, T, c)",
        "    implicit none",
        "    real(WP), dimension(nS), intent(in) :: y     ! mass fractions",
        "    real(WP), dimension(nS), intent(in) :: W_sp  ! molar masses (kg/mol)",
        "    real(WP), intent(in) :: P, T                 ! pressure (Pa), temperature (K)",
        "    real(WP), dimension(nS), intent(out) :: c    ! concentrations (mol/m^3)",
        "    real(WP) :: W_mix, rho",
        "    W_mix = 1.0_WP / sum(y(1:nS) / W_sp(1:nS))",
        "    rho = P * W_mix / (Rcst * T)",
        "    c(1:nS) = rho * y(1:nS) / W_sp(1:nS)",
        "  end subroutine y2c",
        "",
        "  ! NASA7 polynomials for all species at T: cp/R, h/(RT), s/R (low-T coefficients for T <= T_mid).",
        "  ! Results are cached: repeated calls at the same T (thermo + kinetics in one RHS evaluation) cost a copy.",
        "  subroutine fcmech_nasa7(T, cp_R, h_RT, s_R)",
        "    implicit none",
        "    real(WP), intent(in) :: T",
        "    real(WP), dimension(nS), intent(out) :: cp_R, h_RT, s_R",
        "    real(WP) :: T2, T3, T4, lnT, Tinv",
        "    integer :: i, j",
        "    if (T /= T_nasa_cache) then",
        "      T2 = T*T",
        "      T3 = T2*T",
        "      T4 = T3*T",
        "      lnT = log(T)",
        "      Tinv = 1.0_WP/T",
        "      do i = 1, nS",
        "        if (T <= T_mid(i)) then",
        "          j = 0",
        "        else",
        "          j = 7",
        "        end if",
        "        cp_R_cache(i) = thermo_coeffs(i,j+1) + thermo_coeffs(i,j+2)*T + thermo_coeffs(i,j+3)*T2 &",
        "                      + thermo_coeffs(i,j+4)*T3 + thermo_coeffs(i,j+5)*T4",
        "        h_RT_cache(i) = thermo_coeffs(i,j+1) + 0.5_WP*thermo_coeffs(i,j+2)*T + thermo_coeffs(i,j+3)*T2/3.0_WP &",
        "                      + 0.25_WP*thermo_coeffs(i,j+4)*T3 + 0.2_WP*thermo_coeffs(i,j+5)*T4 + thermo_coeffs(i,j+6)*Tinv",
        "        s_R_cache(i) = thermo_coeffs(i,j+1)*lnT + thermo_coeffs(i,j+2)*T + 0.5_WP*thermo_coeffs(i,j+3)*T2 &",
        "                     + thermo_coeffs(i,j+4)*T3/3.0_WP + 0.25_WP*thermo_coeffs(i,j+5)*T4 + thermo_coeffs(i,j+7)",
        "      end do",
        "      T_nasa_cache = T",
        "    end if",
        "    cp_R(:) = cp_R_cache(:)",
        "    h_RT(:) = h_RT_cache(:)",
        "    s_R(:) = s_R_cache(:)",
        "  end subroutine fcmech_nasa7",
        "",
        "  ! Set module arrays Cpsp (J/(mol K)) and hsp (J/mol) at T",
        "  subroutine fcmech_thermodata(T)",
        "    implicit none",
        "    real(WP), intent(in) :: T",
        "    real(WP), dimension(nS) :: cp_R, h_RT, s_R",
        "    call fcmech_nasa7(T, cp_R, h_RT, s_R)",
        "    Cpsp(:) = Rcst * cp_R(:)",
        "    hsp(:) = Rcst * T * h_RT(:)",
        "  end subroutine fcmech_thermodata",
        "",
        "  subroutine fcmech_get_thermodata(h,cp,T)",
        "    implicit none",
        "    real(WP), dimension(nS), intent(out) :: h, cp  ! enthalpy (J/mol), heat capacity (J/(mol K))",
        "    real(WP), intent(in) :: T",
        "    call fcmech_thermodata(T)",
        "    h = hsp",
        "    cp = Cpsp",
        "  end subroutine fcmech_get_thermodata",
        "",
        "  subroutine fcmech_get_speciesnames(names)",
        "    implicit none",
        "    character(len=*), dimension(nS), intent(out) :: names",
    ]
    for i, nm in enumerate(species_names):
        L.append(f"    names({i + 1}) = '{nm.replace(chr(39), chr(39) * 2)}'")
    L += ["  end subroutine fcmech_get_speciesnames", ""]
    if nA > 0:
        L += [
            "  subroutine fcmech_get_atomnames(names)",
            "    implicit none",
            "    character(len=*), dimension(nA), intent(out) :: names",
            "    names(1:nA) = atom_names(1:nA)",
            "  end subroutine fcmech_get_atomnames",
            "",
            "  subroutine fcmech_get_atommasses(masses)",
            "    implicit none",
            "    real(WP), dimension(nA), intent(out) :: masses  ! atom molar masses (kg/mol)",
            "    masses(1:nA) = atom_masses(1:nA)",
            "  end subroutine fcmech_get_atommasses",
            "",
            "  subroutine fcmech_get_composition(comp)",
            "    implicit none",
            "    integer, dimension(nA,nS), intent(out) :: comp  ! comp(a,s) = count of atom a in species s",
            "    comp(1:nA,1:nS) = comp_matrix(1:nA,1:nS)",
            "  end subroutine fcmech_get_composition",
            "",
        ]

    # ---- Transport ----
    if transport is not None:
        L += [
            "  ! Pure-species viscosity (Pa s), Chapman-Enskog with Lennard-Jones potential",
            "  subroutine fcmech_get_viscosity(mu, T)",
            "    implicit none",
            "    real(WP), dimension(nS), intent(out) :: mu",
            "    real(WP), intent(in) :: T",
            "    integer :: i",
            "    do i = 1, nS",
            "      mu(i) = mucoeff(i)*sqrt(T)/fcmech_omegamu(T*koveps(i))",
            "    end do",
            "  end subroutine fcmech_get_viscosity",
            "",
            "  ! Pure-species thermal conductivity (W/(m K)), modified Eucken: lambda = mu*(cp + 1.25*R)/W",
            "  subroutine fcmech_get_conductivity(lambda, T, mu)",
            "    implicit none",
            "    real(WP), dimension(nS), intent(out) :: lambda",
            "    real(WP), intent(in) :: T",
            "    real(WP), dimension(nS), intent(in) :: mu  ! viscosity from fcmech_get_viscosity",
            "    real(WP), dimension(nS) :: cp_R, h_RT, s_R",
            "    call fcmech_nasa7(T, cp_R, h_RT, s_R)",
            "    lambda(:) = mu(:) * Rcst * (cp_R(:) + 1.25_WP) / W_sp(:)",
            "  end subroutine fcmech_get_conductivity",
            "",
            "  ! Inverse binary diffusion coefficients 1/D_ij (s/m^2), Chapman-Enskog",
            "  subroutine fcmech_get_invDij(invDij, T, P)",
            "    implicit none",
            "    real(WP), dimension(nS,nS), intent(out) :: invDij",
            "    real(WP), intent(in) :: T, P  ! temperature (K), pressure (Pa)",
            "    real(WP) :: TPterm",
            "    integer :: i, j",
            "    TPterm = P/(T*sqrt(T))",
            "    do i = 1, nS",
            "      do j = 1, i-1",
            "        invDij(i,j) = TPterm*fcmech_omegaD(T*Ocoeffs(i,j))/Dcoeffs(i,j)",
            "        invDij(j,i) = invDij(i,j)",
            "      end do",
            "      invDij(i,i) = 0.0_WP",
            "    end do",
            "  end subroutine fcmech_get_invDij",
            "",
        ]

    # ---- Main entry ----
    L += [
        "  ! Mass-based production rates dY_i/dt (1/s) at (P,T,Y): ydot_i = W_i*cdot_i/rho",
        "  subroutine fcmech_get_ydot(P, T, Y, ydot)",
        "    implicit none",
        "    real(WP), intent(in) :: P, T                 ! pressure (Pa), temperature (K)",
        "    real(WP), dimension(nS), intent(in) :: Y     ! mass fractions",
        "    real(WP), dimension(nS), intent(out) :: ydot",
        "    real(WP), dimension(nS) :: c, wdot",
        "    real(WP), dimension(nR) :: k, w",
        "    real(WP), dimension(n_m) :: M",
        "    real(WP) :: W_mix, rho",
        "    call y2c(Y, W_sp, P, T, c)",
        "    call get_thirdbodies(M, c)",
        "    call get_rate_coefficients(k, M, T, P)",
        "    call get_reaction_rates(w, k, M, c)",
        "    call get_production_rates(wdot, w)",
        "    W_mix = 1.0_WP / sum(Y(1:nS) / W_sp(1:nS))",
        "    rho = P * W_mix / (Rcst * T)",
        "    ydot(1:nS) = W_sp(1:nS) * wdot(1:nS) / rho",
        "  end subroutine fcmech_get_ydot",
        "",
        "  ! Mass-based production and destruction rates (1/s) at (P,T,Y): ydot_i = W_i*cdot_i/rho and",
        "  ! ddot_i = W_i*cdes_i/rho >= 0; ddot_i/Y_i is the destruction-rate coefficient of species i (explicit",
        "  ! time-step control). Same rate evaluation as fcmech_get_ydot, one extra pass over the reactions.",
        "  subroutine fcmech_get_ydot_ddot(P, T, Y, ydot, ddot)",
        "    implicit none",
        "    real(WP), intent(in) :: P, T                 ! pressure (Pa), temperature (K)",
        "    real(WP), dimension(nS), intent(in) :: Y     ! mass fractions",
        "    real(WP), dimension(nS), intent(out) :: ydot, ddot",
        "    real(WP), dimension(nS) :: c, wdot, wdes",
        "    real(WP), dimension(nR) :: k, w",
        "    real(WP), dimension(n_m) :: M",
        "    real(WP) :: W_mix, rho",
        "    call y2c(Y, W_sp, P, T, c)",
        "    call get_thirdbodies(M, c)",
        "    call get_rate_coefficients(k, M, T, P)",
        "    call get_reaction_rates(w, k, M, c)",
        "    call get_production_rates(wdot, w)",
        "    call get_destruction_rates(wdes, w)",
        "    W_mix = 1.0_WP / sum(Y(1:nS) / W_sp(1:nS))",
        "    rho = P * W_mix / (Rcst * T)",
        "    ydot(1:nS) = W_sp(1:nS) * wdot(1:nS) / rho",
        "    ddot(1:nS) = W_sp(1:nS) * wdes(1:nS) / rho",
        "  end subroutine fcmech_get_ydot_ddot",
        "",
        "end module fcmech",
    ]

    check_fortran_limits(L)
    with open(out_path, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"  Species: {nS}, atoms: {nA}, reactions: {len(rxns)} "
          f"(Arrhenius {n_arr_fwd}, three-body {n_tb}, falloff {n_fo}, PLOG {n_plog}), rates w(1:{nR})")


def main():
    parser = argparse.ArgumentParser(description="Generate the fcmech Fortran module (chem_data_fc.f90) from a YAML mechanism")
    parser.add_argument("input", help="Input YAML mechanism file")
    parser.add_argument("output", nargs="?", default="chem_data_fc.f90", help="Output Fortran file (default: chem_data_fc.f90)")
    parser.add_argument("--fit-reverse", action="store_true",
                        help="Fit k_f/K_c of every reversible reaction to an Arrhenius form and evaluate reverse rates "
                             "like forward ones (faster) instead of run-time detailed balance (exact)")
    parser.add_argument("--fit-range", nargs=2, type=float, metavar=("TMIN", "TMAX"), default=(300.0, 3000.0),
                        help="Temperature range (K) for --fit-reverse (default: 300 3000)")
    parser.add_argument("--fit-points", type=int, default=100,
                        help="Number of fit points, uniform in 1/T, for --fit-reverse (default: 100)")
    args = parser.parse_args()
    fit_reverse = None
    if args.fit_reverse:
        T_min, T_max = args.fit_range
        if not (0.0 < T_min < T_max) or args.fit_points < 4:
            print("Error: --fit-range needs 0 < TMIN < TMAX and --fit-points >= 4", file=sys.stderr)
            sys.exit(1)
        fit_reverse = {"T_min": T_min, "T_max": T_max, "n_pts": args.fit_points}

    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading mechanism from {input_path}")
    try:
        mech = load_mechanism(input_path)
        print(f"Writing {output_path}")
        yaml2nga(mech, output_path, yaml_path=input_path, fit_reverse=fit_reverse)
    except (MechanismError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    print("Done.")


if __name__ == "__main__":
    main()
