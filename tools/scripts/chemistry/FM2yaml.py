#!/usr/bin/env python3
"""
FM2yaml - FlameMaster to Cantera-YAML converter

Converts a FlameMaster mechanism (.mech, ScanMan syntax) with its CHEMKIN-format thermo (.thermo /
.chthermo) and transport (.trans / .chtrans) files into a Cantera-format YAML file for yaml2nga.

Reaction labels
  - number only (25):            irreversible
  - number + f (1f, i66f):       forward; reversible through K_eq unless a matching 'b' exists
  - number + b (i66b):           explicit backward rate (written as its own irreversible reaction)

Rate parameters in braces (FlameMaster units: A in cm, mol, s; E in kJ/mol unless 'Let units for E be [...]')
  - a, n, E                      Arrhenius; for pressure-dependent reactions this is the LOW-pressure limit k0
  - ai, ni, Ei                   HIGH-pressure limit k_inf (the 'i' stands for infinity)
  - fca, fcta, fcb, fctb, fcc, fctc
                                 F_cent = fca*exp(-T/fcta) + fcb*exp(-T/fctb) + fcc*exp(-fctc/T),
                                 mapped onto Cantera's Troe form. FlameMaster allows a coefficient on the
                                 exp(-fctc/T) term while Cantera fixes it at 1, so that case (and any other
                                 F_cent outside the Troe family) is least-squares fitted over --fit-range;
                                 the fit error is printed, written into the YAML and capped by --fit-tol.
Third bodies
  - 'M', "M'", "M''" tokens in equations, defined by  Let [M'] = 6.5 [CH4] + ... + 1.0 [OTHER]
    -> Cantera 'efficiencies' with [OTHER] as 'default-efficiency' (0 if [OTHER] is absent)
  - pressure-dependent reactions without an explicit M use the total concentration (all efficiencies 1)

CHEMKIN thermo cards list the HIGH-temperature coefficients first; they are re-ordered to Cantera's
low-T-first convention.
"""

from __future__ import annotations

import argparse
import io
import math
import re
import sys
from datetime import datetime
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))
from yaml2nga import MechanismError, canonical_element as _canonical_element  # noqa: E402

KJ_TO_CAL = 1000.0 / 4.184  # kJ/mol -> cal/mol
# 'Let units for E be [...]' -> factor to cal/mol (Cantera output unit); default FlameMaster unit is kJ/mole
_E_UNITS_TO_CAL = {"kJ/mole": KJ_TO_CAL, "kJ/mol": KJ_TO_CAL, "kcal/mole": 1000.0, "kcal/mol": 1000.0,
                   "cal/mole": 1.0, "cal/mol": 1.0, "J/mole": 1.0 / 4.184, "J/mol": 1.0 / 4.184}
_TROE_KEYS = ("fca", "fcta", "fcb", "fctb", "fcc", "fctc")
T_KILL = 1.0e-30  # Troe time scale that switches a term off (exp(-T/1e-30) = 0), as used by Cantera
T_INF = 1.0e30


class FMError(ValueError):
    pass


def fmt(val: float) -> str:
    """Compact, round-trip safe number for YAML."""
    return repr(float(val))


def upper_species(name: str) -> str:
    return name.strip().upper()


# ---------------------------------------------------------------------------------------------------
# Mechanism (.mech)
# ---------------------------------------------------------------------------------------------------
_KV_RE = re.compile(r"([A-Za-z]+)\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")
_KNOWN_KEYS = {"a", "n", "E", "ai", "ni", "Ei", "fca", "fcta", "fcb", "fctb", "fcc", "fctc"}
_M_TOKEN_RE = re.compile(r"^M'*$")


def parse_rate_params(text: str, ctx: str) -> dict:
    """'a = 2.e14 n = 0 E = 70.3 ...' -> {key: float}; unknown keys are an error."""
    text = text.replace("\t", " ")
    params = {}
    for m in _KV_RE.finditer(text):
        key, val = m.group(1), float(m.group(2))
        if key not in _KNOWN_KEYS:
            raise FMError(f"{ctx}: unknown rate parameter '{key}'")
        if key in params:
            raise FMError(f"{ctx}: rate parameter '{key}' given twice")
        params[key] = val
    leftover = _KV_RE.sub("", text).strip()
    if leftover:
        raise FMError(f"{ctx}: cannot parse rate parameters near '{leftover}'")
    for key in ("a", "n", "E"):
        if key not in params:
            raise FMError(f"{ctx}: missing rate parameter '{key}'")
    return params


def parse_side(text: str, ctx: str):
    """'O2 + H + M'' -> ([(coeff, SPECIES)], collider) with species uppercased."""
    species = []
    collider = None
    for term in text.replace("\t", " ").split("+"):
        term = term.strip()
        if not term:
            raise FMError(f"{ctx}: empty term in '{text}'")
        parts = term.split()
        if len(parts) == 2 and re.fullmatch(r"\d+(\.\d*)?", parts[0]):
            coeff, name = float(parts[0]), parts[1]
        elif len(parts) == 1:
            m = re.fullmatch(r"(\d+)([A-Za-z].*)", parts[0])  # '2OH'
            if m:
                coeff, name = float(m.group(1)), m.group(2)
            else:
                coeff, name = 1.0, parts[0]
        else:
            raise FMError(f"{ctx}: cannot parse term '{term}'")
        if _M_TOKEN_RE.match(name):
            if collider is not None:
                raise FMError(f"{ctx}: more than one third body on one side")
            collider = name
            continue
        species.append((coeff, upper_species(name)))
    if not species:
        raise FMError(f"{ctx}: no species on one side of '{text}'")
    return species, collider


def format_side(species, collider_out: str | None) -> str:
    text = " + ".join(f"{c:g} {nm}" if c != 1.0 else nm for c, nm in species)
    if collider_out is not None:
        text += (" " if collider_out.startswith("(") else " + ") + collider_out
    return text


def fm_fcent_terms(p: dict, ctx: str):
    """
    FlameMaster F_cent = fca*exp(-T/fcta) + fcb*exp(-T/fctb) + fcc*exp(-fctc/T), returned as the term
    lists (decay, rise): decay holds (c, tau) for c*exp(-T/tau), rise holds (c, T2) for c*exp(-T2/T).
    A missing fctc means exp(0) = 1, i.e. a constant term (written as a decay term with tau = inf).
    """
    decay, rise = [], []
    for c_key, t_key in (("fca", "fcta"), ("fcb", "fctb")):
        if c_key in p or t_key in p:
            if c_key not in p or t_key not in p:
                raise FMError(f"{ctx}: {c_key} and {t_key} must be given together")
            if p[t_key] <= 0.0:
                raise FMError(f"{ctx}: {t_key} must be positive")
            decay.append((p[c_key], p[t_key]))
    if "fcc" in p or "fctc" in p:
        c3, t3 = p.get("fcc", 1.0), p.get("fctc", 0.0)
        if t3 < 0.0:
            raise FMError(f"{ctx}: fctc must not be negative")
        if t3 == 0.0:
            decay.append((c3, T_INF))  # exp(-0/T) = 1: a constant term
        else:
            rise.append((c3, t3))
    return decay, rise


def fcent_fm(terms, T: float) -> float:
    decay, rise = terms
    return (sum(c * math.exp(-T / tau) for c, tau in decay)
            + sum(c * math.exp(-t2 / T) for c, t2 in rise))


def fcent_troe(A: float, T3: float, T1: float, T2, T: float) -> float:
    """Cantera/Troe F_cent = (1-A)*exp(-T/T3) + A*exp(-T/T1) [+ exp(-T2/T)]."""
    v = (1.0 - A) * math.exp(-T / T3) + A * math.exp(-T / T1)
    if T2 is not None:
        v += math.exp(-T2 / T)
    return v


def _nelder_mead(fun, x0, step, maxiter: int = 2000, tol: float = 1.0e-13):
    """Minimal Nelder-Mead simplex (the repo's Python tools are standard-library only)."""
    n = len(x0)
    pts = [list(x0)] + [list(x0) for _ in range(n)]
    for i in range(n):
        pts[i + 1][i] += step[i]
    val = [fun(x) for x in pts]
    for _ in range(maxiter):
        order = sorted(range(n + 1), key=lambda k: val[k])
        pts = [pts[k] for k in order]
        val = [val[k] for k in order]
        if abs(val[-1] - val[0]) <= tol * (abs(val[0]) + tol):
            break
        cen = [sum(pt[i] for pt in pts[:-1]) / n for i in range(n)]
        xr = [cen[i] + (cen[i] - pts[-1][i]) for i in range(n)]
        fr = fun(xr)
        if fr < val[0]:
            xe = [cen[i] + 2.0 * (cen[i] - pts[-1][i]) for i in range(n)]
            fe = fun(xe)
            pts[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            pts[-1], val[-1] = xr, fr
        else:
            xc = [cen[i] + 0.5 * (pts[-1][i] - cen[i]) for i in range(n)]
            fc = fun(xc)
            if fc < val[-1]:
                pts[-1], val[-1] = xc, fc
            else:
                for k in range(1, n + 1):
                    pts[k] = [pts[0][i] + 0.5 * (pts[k][i] - pts[0][i]) for i in range(n)]
                    val[k] = fun(pts[k])
    k = min(range(n + 1), key=lambda i: val[i])
    return pts[k], val[k]


def fit_troe(terms, ctx: str, t_min: float, t_max: float):
    """
    Least-squares Troe {A, T3, T1[, T2]} for an F_cent the Troe form cannot reproduce exactly
    (FlameMaster allows a coefficient on exp(-fctc/T); Cantera fixes it at 1). Fitted on log F_cent
    over [t_min, t_max]. Returns ((A, T3, T1, T2), largest relative F_cent error).
    """
    grid = [t_min * (t_max / t_min) ** (k / 24.0) for k in range(25)]
    ref = [fcent_fm(terms, T) for T in grid]
    if min(ref) <= 0.0:
        raise FMError(f"{ctx}: F_cent is not positive over {t_min:g}-{t_max:g} K; cannot fit a Troe form")
    lref = [math.log(v) for v in ref]

    def pow10(x: float) -> float:
        """10**x with the simplex kept inside a sane, overflow-free range."""
        return 10.0 ** min(max(x, -30.0), 30.0)

    def make_obj(use_T2):
        def obj(x):
            if not all(-1.0e6 < v < 1.0e6 for v in x):
                return 1.0e30
            A, T3, T1 = x[0], pow10(x[1]), pow10(x[2])
            T2 = pow10(x[3]) if use_T2 else None
            s = 0.0
            try:
                for T, lr in zip(grid, lref):
                    v = fcent_troe(A, T3, T1, T2, T)
                    if v <= 0.0:
                        return 1.0e30
                    d = math.log(v) - lr
                    s += d * d
            except (OverflowError, ValueError):
                return 1.0e30
            return s
        return obj

    best = None
    for use_T2 in (True, False):
        obj = make_obj(use_T2)
        for A0 in (0.3, 0.9, 1.2):
            for l3 in (0.0, 3.0):
                for l1 in (2.5, 6.0):
                    for l2 in (1.5, 3.5):
                        x, v = _nelder_mead(obj, [A0, l3, l1, l2], [0.2, 0.5, 0.5, 0.5])
                        if best is None or v < best[0]:
                            best = (v, x, use_T2)
    _, x, use_T2 = best
    A, T3, T1 = x[0], pow10(x[1]), pow10(x[2])
    T2 = pow10(x[3]) if use_T2 else None
    if T3 <= 0.0 or T1 <= 0.0:
        raise FMError(f"{ctx}: Troe fit produced a non-positive time scale")
    fine = [t_min * (t_max / t_min) ** (k / 200.0) for k in range(201)]
    err = max(abs(fcent_troe(A, T3, T1, T2, T) / fcent_fm(terms, T) - 1.0) for T in fine)
    return (A, T3, T1, T2), err


def troe_from_fm(p: dict, ctx: str, t_min: float = 300.0, t_max: float = 3000.0, tol: float = 0.10):
    """
    Map the FlameMaster F_cent onto Cantera's Troe form. Returns (troe_or_None, fit_error_or_None);
    None means Lindemann (no centering parameters at all). The mapping is exact when F_cent is a sum
    of at most two exp(-T/tau) terms whose coefficients sum to 1, plus at most one exp(-fctc/T) term
    with coefficient 1 -- that covers every F_cent a Troe form can represent. Anything else (in
    practice fcc != 1 with fctc != 0) is least-squares fitted and the error is reported to the caller.
    """
    if not any(k in p for k in _TROE_KEYS):
        return None, None
    decay, rise = fm_fcent_terms(p, ctx)
    T2 = None
    exact = True
    if len(rise) > 1:
        exact = False
    elif len(rise) == 1:
        if abs(rise[0][0] - 1.0) < 1.0e-9:
            T2 = rise[0][1]
        else:
            exact = False  # FlameMaster allows a coefficient here, Cantera does not
    if exact:
        if len(decay) == 0:
            A, T3, T1 = 1.0, T_INF, T_KILL
        elif len(decay) == 1:
            (c, tau), = decay
            A, T3, T1 = 1.0 - c, tau, T_KILL
        elif len(decay) == 2:
            (c1, tau1), (c2, tau2) = decay
            if abs(c1 + c2 - 1.0) > 1.0e-6:
                exact = False
            else:
                A, T3, T1 = c2, tau1, tau2
        else:
            exact = False
    if exact:
        return (A, T3, T1, T2), None
    troe, err = fit_troe((decay, rise), ctx, t_min, t_max)
    if err > tol:
        raise FMError(f"{ctx}: F_cent cannot be represented as a Troe form (best fit is off by "
                      f"{100.0 * err:.1f} % over {t_min:g}-{t_max:g} K, tolerance {100.0 * tol:.1f} %)")
    return troe, err


def parse_mech_file(mech_file: Path):
    """Returns (reactions, third_bodies). third_bodies: {"M'": ({SPECIES: eff}, default_or_None)}"""
    with open(mech_file) as fp:
        raw_lines = fp.readlines()

    # Strip comments and join continuation lines (braces and trailing '+')
    statements = []
    buf = ""
    for line in raw_lines:
        line = line.split("#", 1)[0].rstrip("\n")
        if not line.strip():
            if buf and (("{" in buf and "}" not in buf) or buf.rstrip().endswith("+")):
                continue  # blank line inside a brace block or a '+'-continued list
            if buf:
                statements.append(buf)
                buf = ""
            continue
        buf = (buf + " " + line.strip()) if buf else line.strip()
        if "{" in buf and "}" not in buf:
            continue
        if buf.rstrip().endswith("+"):
            continue
        statements.append(buf)
        buf = ""
    if buf:
        statements.append(buf)

    reactions = []
    third_bodies = {}
    additional = []
    E_to_cal = KJ_TO_CAL
    rxn_re = re.compile(r"^(\S+?)\s*:\s*(.*?)\s*->\s*(.*?)\s*\{(.*)\}\s*\.?\s*$")
    let_re = re.compile(r"^Let\s+\[(M'*)\]\s*=\s*(.*?)\s*\.?\s*$", re.IGNORECASE)
    add_re = re.compile(r"^Let\s+additional\s+species\s+be\s+(.*?)\s*\.?\s*$", re.IGNORECASE)
    units_re = re.compile(r"^Let\s+units\s+for\s+(\w+)\s+be\s*\[(.*)\]\s*\.?\s*$", re.IGNORECASE)
    for st in statements:
        m = add_re.match(st)
        if m:
            additional += [upper_species(s) for s in re.split(r"[,\s]+", m.group(1)) if s]
            continue
        m = units_re.match(st)
        if m:
            what, unit = m.group(1), re.sub(r"\s+", "", m.group(2))
            if what.upper() == "E":
                if unit not in _E_UNITS_TO_CAL:
                    raise FMError(f"Unsupported activation-energy unit '{unit}' (supported: {', '.join(_E_UNITS_TO_CAL)})")
                E_to_cal = _E_UNITS_TO_CAL[unit]
            elif what.upper() == "A":
                if "cm" not in unit or "mole" not in unit.lower():
                    raise FMError(f"Unsupported pre-exponential unit '{unit}' (expected cm^(3(n-1))/(s*mole^(n-1)*K^n))")
            continue
        m = let_re.match(st)
        if m:
            third_bodies[m.group(1)] = parse_third_body(m.group(2), st)
            continue
        if st.lower().startswith("let "):
            continue  # allowed atoms, temperature exponent, order of reaction, ...: no effect on the output
        m = rxn_re.match(st)
        if not m:
            raise FMError(f"Cannot parse statement (expected 'id: A + B -> C {{ a = ... }}' or 'Let ...'): '{st}'")
        rid, lhs, rhs, params = m.groups()
        ctx = f"Reaction {rid}"
        suffix = rid[-1] if rid[-1] in "fb" else ""
        base_id = rid[:-1] if suffix else rid
        left, col_l = parse_side(lhs, ctx)
        right, col_r = parse_side(rhs, ctx)
        if col_l != col_r:
            raise FMError(f"{ctx}: third body must appear on both sides ('{lhs} -> {rhs}')")
        p = parse_rate_params(params, ctx)
        reactions.append({
            "id": rid, "base_id": base_id, "suffix": suffix,
            "left": left, "right": right, "collider": col_l, "params": p,
        })
    if not reactions:
        raise FMError(f"No FlameMaster reactions ('id: A + B -> C {{ a = ... }}') found in {mech_file}")
    return reactions, third_bodies, additional, E_to_cal


def parse_third_body(body: str, st: str):
    """'6.5 [CH4] + [H2] + 1.0 [OTHER]' -> ({SPECIES: eff}, default_or_None); a missing coefficient means 1."""
    effs = {}
    default = None
    rest = body
    for cm in re.finditer(r"(?:([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*)?\[\s*(\S+?)\s*\]", body):
        val = float(cm.group(1)) if cm.group(1) is not None else 1.0
        sp = cm.group(2)
        if val < 0.0:
            raise FMError(f"Negative third-body efficiency for [{sp}] in '{st}'")
        if sp.upper() == "OTHER":
            if default is not None:
                raise FMError(f"[OTHER] given twice in '{st}'")
            default = val
        else:
            if upper_species(sp) in effs:
                raise FMError(f"Species [{sp}] listed twice in '{st}'")
            effs[upper_species(sp)] = val
        rest = rest.replace(cm.group(0), " ", 1)
    if rest.replace("+", "").strip():
        raise FMError(f"Cannot parse third-body definition near '{rest.strip()}' in '{st}'")
    if not effs and default is None:
        raise FMError(f"Empty third-body definition: '{st}'")
    return effs, default


def fcent_text(p: dict) -> str:
    """The FlameMaster centering function as written in the .mech, for the YAML comment."""
    bits = []
    for c_key, t_key in (("fca", "fcta"), ("fcb", "fctb")):
        if c_key in p:
            bits.append(f"{p[c_key]:g}*exp(-T/{p[t_key]:g})")
    if "fcc" in p or "fctc" in p:
        bits.append(f"{p.get('fcc', 1.0):g}*exp(-{p.get('fctc', 0.0):g}/T)")
    return " + ".join(bits)


def write_reactions(f, reactions, third_bodies, species_names, E_to_cal: float,
                    fit_range=(300.0, 3000.0), fit_tol: float = 0.10) -> int:
    fitted = []
    f.write("reactions:\n")
    has_b = {r["base_id"] for r in reactions if r["suffix"] == "b"}
    # Cantera requires 'duplicate: true' on reactions with identical reactants and products
    def dup_key(r):
        pressure_dep = any(k in r["params"] for k in ("ai", "ni", "Ei"))

        def side(terms):  # summed stoichiometry: 'CO + CO' and '2 CO' are the same side
            tot = {}
            for c, nm in terms:
                tot[nm] = tot.get(nm, 0.0) + c
            return tuple(sorted(tot.items()))
        return (side(r["left"]), side(r["right"]), r["collider"] is not None, pressure_dep)
    counts = {}
    for r in reactions:
        counts[dup_key(r)] = counts.get(dup_key(r), 0) + 1
    n_out = 0
    for r in reactions:
        rid, p, ctx = r["id"], r["params"], f"Reaction {r['id']}"
        reversible = r["suffix"] == "f" and r["base_id"] not in has_b
        pressure_dep = any(k in p for k in ("ai", "ni", "Ei"))
        collider = r["collider"]
        if pressure_dep:
            for k in ("ai", "ni", "Ei"):
                if k not in p:
                    raise FMError(f"{ctx}: pressure-dependent reaction needs ai, ni and Ei")
            # FlameMaster falloff reactions usually carry no explicit M: the third body is then the total
            # concentration (all efficiencies 1). An explicit M' selects the corresponding efficiency set.
            col_out = "(+M)"
        else:
            col_out = "M" if collider is not None else None
        arrow = "<=>" if reversible else "=>"
        n_out += 1
        eq = f"{format_side(r['left'], col_out)} {arrow} {format_side(r['right'], col_out)}"
        f.write(f"- equation: {eq}  # Reaction {n_out} (FlameMaster {rid})\n")
        if counts[dup_key(r)] > 1:
            f.write("  duplicate: true\n")
        if pressure_dep:
            f.write("  type: falloff\n")
            f.write(f"  low-P-rate-constant: {{A: {fmt(p['a'])}, b: {fmt(p['n'])}, Ea: {fmt(p['E'] * E_to_cal)}}}\n")
            f.write(f"  high-P-rate-constant: {{A: {fmt(p['ai'])}, b: {fmt(p['ni'])}, Ea: {fmt(p['Ei'] * E_to_cal)}}}\n")
            troe, fit_err = troe_from_fm(p, ctx, fit_range[0], fit_range[1], fit_tol)
            if troe is not None:
                A, T3, T1, T2 = troe
                t2 = f", T2: {fmt(T2)}" if T2 is not None else ""
                if fit_err is not None:
                    fitted.append((ctx, 100.0 * fit_err))
                    f.write(f"  # F_cent = {fcent_text(p)} is not a Troe form; least-squares fit over\n")
                    f.write(f"  # {fit_range[0]:g}-{fit_range[1]:g} K, largest F_cent error {100.0 * fit_err:.2f} %\n")
                f.write(f"  Troe: {{A: {fmt(A)}, T3: {fmt(T3)}, T1: {fmt(T1)}{t2}}}\n")
        else:
            if any(k in p for k in _TROE_KEYS):
                raise FMError(f"{ctx}: Troe parameters given for a reaction without a low/high pressure pair")
            if collider is not None:
                f.write("  type: three-body\n")
            f.write(f"  rate-constant: {{A: {fmt(p['a'])}, b: {fmt(p['n'])}, Ea: {fmt(p['E'] * E_to_cal)}}}\n")
        if collider is not None:
            if collider in third_bodies:
                effs, default = third_bodies[collider]
                known = {sp: v for sp, v in effs.items() if sp in species_names}
                dropped = sorted(set(effs) - set(known))
                if dropped:
                    print(f"  Warning: {ctx}: third-body efficiencies for species not in the thermo file dropped: {', '.join(dropped)}")
                if known:
                    f.write("  efficiencies: {" + ", ".join(f"{sp}: {fmt(v)}" for sp, v in known.items()) + "}\n")
                if default is None:
                    print(f"  Note: [{collider}] has no [OTHER] term; unlisted species get efficiency 0")
                    default = 0.0
                if default != 1.0:
                    f.write(f"  default-efficiency: {fmt(default)}\n")
            elif collider != "M":
                raise FMError(f"{ctx}: third body [{collider}] is used but never defined with 'Let [{collider}] = ...'")
    for ctx, err in fitted:
        print(f"  Warning: {ctx}: F_cent has a coefficient on exp(-fctc/T) that Cantera's Troe form"
              f" cannot hold; least-squares Troe fit over {fit_range[0]:g}-{fit_range[1]:g} K,"
              f" largest F_cent error {err:.2f} %")
    return n_out


# ---------------------------------------------------------------------------------------------------
# Transport (.trans)
# ---------------------------------------------------------------------------------------------------
def parse_transport_file(trans_file: Path) -> dict:
    """
    {SPECIES: {"geom", "eps", "sigma", "dipole", "polar", "rot"}} from a CHEMKIN-style transport table:
    name, geometry index (0 atom / 1 linear / 2 nonlinear; FlameMaster writes a dummy 0), eps/k [K],
    sigma [A], then optionally dipole [Debye], polarizability [A^3], rotational relaxation number.
    """
    data = {}
    if not trans_file.exists():
        return data
    with open(trans_file) as fp:
        for ln, line in enumerate(fp, 1):
            s = line.split("#", 1)[0].strip()
            if not s or "specname" in s.lower():
                continue
            parts = s.split()
            if not re.match(r"[A-Za-z]", parts[0]):
                continue
            try:
                nums = [float(x) for x in parts[1:]]
            except ValueError:
                raise FMError(f"{trans_file.name}:{ln}: cannot parse transport data '{s}'")
            if len(nums) < 2:
                raise FMError(f"{trans_file.name}:{ln}: transport line needs at least eps/k and sigma: '{s}'")
            if len(nums) == 2:
                geom, eps, sigma, extra = None, nums[0], nums[1], []
            else:
                geom, eps, sigma, extra = nums[0], nums[1], nums[2], nums[3:6]
            if eps <= 0.0 or sigma <= 0.0:
                raise FMError(f"{trans_file.name}:{ln}: eps/k and sigma must be positive: '{s}'")
            extra += [0.0] * (3 - len(extra))
            name = upper_species(parts[0])
            if name in data:
                raise FMError(f"{trans_file.name}:{ln}: duplicate transport entry for {name}")
            data[name] = {"geom": geom, "eps": eps, "sigma": sigma, "dipole": extra[0], "polar": extra[1], "rot": extra[2]}
    print(f"  Parsed {len(data)} transport entries")
    return data


# ---------------------------------------------------------------------------------------------------
# Thermo (.thermo, CHEMKIN cards)
# ---------------------------------------------------------------------------------------------------
def canonical_element(sym: str, ctx: str) -> str:
    """Canonical element symbol from the shared table (yaml2nga.ATOMIC_MASSES_KG); unknown -> error."""
    try:
        return _canonical_element(sym)
    except MechanismError as exc:
        raise FMError(f"{ctx}: {exc}") from None


def parse_card_composition(card: str, ctx: str) -> dict:
    """Element/count pairs in columns 25-44 (four 5-character fields) of CHEMKIN species card 1."""
    comp = {}
    card = card.ljust(80)
    for k in range(4):
        sym = card[24 + 5 * k : 26 + 5 * k].strip()
        cnt = card[26 + 5 * k : 29 + 5 * k].strip()
        if not sym or sym == "0" or (not cnt and not re.match(r"[A-Za-z]", sym)):
            continue
        try:
            n = int(float(cnt)) if cnt else 0
        except ValueError:
            raise FMError(f"{ctx}: cannot read the count '{cnt}' of element '{sym}' on card 1")
        if n < 0:
            raise FMError(f"{ctx}: negative element count on card 1")
        if n > 0:
            e = canonical_element(sym, ctx)
            comp[e] = comp.get(e, 0) + n
    return comp


_LINEAR_FORMULAS = {  # polyatomic molecules that CHEMKIN transport tables list as linear (geometry 1)
    frozenset({("C", 1), ("O", 2)}), frozenset({("N", 2), ("O", 1)}), frozenset({("C", 2), ("H", 1)}),
    frozenset({("C", 1), ("H", 2)}), frozenset({("H", 1), ("C", 1), ("N", 1)}), frozenset({("N", 1), ("C", 1), ("O", 1)}),
    frozenset({("C", 1), ("H", 3)}), frozenset({("C", 2), ("H", 2)}),
}


def geometry_for(comp: dict) -> str:
    """Geometry guessed from the formula when the transport table carries no geometry column."""
    n_atoms = sum(comp.values())
    if n_atoms == 1:
        return "atom"
    if n_atoms == 2 or frozenset(comp.items()) in _LINEAR_FORMULAS:
        return "linear"
    return "nonlinear"


def parse_thermo_file(thermo_file: Path, transport: dict):
    """Returns list of species dicts {name, comp, T_low, T_mid, T_high, low[7], high[7], trans}."""
    species = []
    with open(thermo_file) as fp:
        lines = [ln.rstrip("\n") for ln in fp]

    def is_card1(ln):
        return len(ln) >= 80 and ln[79:80] == "1"

    # Header: 'THERMO' then the global temperature line 'T_low T_mid T_high' used when a card leaves them blank
    global_T = None
    i = 0
    while i < len(lines) and not is_card1(lines[i]):
        vals = lines[i].split("!")[0].split()
        if len(vals) == 3:
            try:
                global_T = tuple(float(v) for v in vals)
            except ValueError:
                pass
        i += 1

    def temperature(field_text, which, name):
        s = field_text.strip()
        if s:
            try:
                return float(s)
            except ValueError:
                raise FMError(f"Species '{name}': cannot read {which} temperature '{s}' from card 1")
        if global_T is None:
            raise FMError(f"Species '{name}': {which} temperature is blank on card 1 and the THERMO header has no defaults")
        return global_T[{"low": 0, "common": 1, "high": 2}[which]]

    def field(line, k, name, card_no):
        s = line[15 * k : 15 * (k + 1)].strip()
        if not s:
            raise FMError(f"Species '{name}': blank NASA7 coefficient (field {k + 1} of card {card_no})")
        try:
            return float(s.replace("D", "E").replace("d", "e"))
        except ValueError:
            raise FMError(f"Species '{name}': cannot read NASA7 coefficient '{s}' (card {card_no})")

    while i < len(lines):
        card = lines[i]
        if not is_card1(card):
            if card.strip().upper() == "END":
                break
            i += 1
            continue
        name = card[:18].strip().split()[0].upper()
        ctx = f"Species '{name}'"
        if "&" in card[44:]:
            raise FMError(f"{ctx}: extended-composition continuation ('&') is not supported")
        comp = parse_card_composition(card, ctx)
        if not comp:
            raise FMError(f"{ctx}: no elemental composition on card 1 (columns 25-44)")
        T_low = temperature(card[45:55], "low", name)
        T_high = temperature(card[55:65], "high", name)
        T_mid = temperature(card[65:73], "common", name)
        if i + 3 >= len(lines):
            raise FMError(f"{ctx}: thermo file ends before card 4")
        c2, c3, c4 = (lines[i + 1].ljust(80), lines[i + 2].ljust(80), lines[i + 3].ljust(80))
        # CHEMKIN order: card 2 = high[0:5], card 3 = high[5:7] + low[0:3], card 4 = low[3:7]
        high = [field(c2, k, name, 2) for k in range(5)] + [field(c3, 0, name, 3), field(c3, 1, name, 3)]
        low = [field(c3, 2, name, 3), field(c3, 3, name, 3), field(c3, 4, name, 3)] + [field(c4, k, name, 4) for k in range(4)]
        if any(sp["name"] == name for sp in species):
            raise FMError(f"{ctx}: defined twice in the thermo file")
        species.append({"name": name, "comp": comp, "T_low": T_low, "T_mid": T_mid, "T_high": T_high,
                        "low": low, "high": high, "trans": transport.get(name)})
        i += 4
    if not species:
        raise FMError(f"No species cards found in {thermo_file}")
    print(f"  Parsed {len(species)} species")
    return species


def write_species(f, species, use_geometry_column: bool) -> None:
    f.write("species:\n")
    for sp in species:
        f.write(f"- name: {sp['name']}\n")
        f.write("  composition: {" + ", ".join(f"{e}: {n}" for e, n in sp["comp"].items()) + "}\n")
        f.write("  thermo:\n")
        f.write("    model: NASA7\n")
        f.write(f"    temperature-ranges: [{fmt(sp['T_low'])}, {fmt(sp['T_mid'])}, {fmt(sp['T_high'])}]\n")
        f.write("    data:\n")
        f.write("    - [" + ", ".join(fmt(c) for c in sp["low"]) + "]\n")
        f.write("    - [" + ", ".join(fmt(c) for c in sp["high"]) + "]\n")
        tr = sp["trans"]
        if tr is not None:
            if use_geometry_column and tr["geom"] in (0.0, 1.0, 2.0):
                geometry = {0.0: "atom", 1.0: "linear", 2.0: "nonlinear"}[tr["geom"]]
            else:
                geometry = geometry_for(sp["comp"])
            f.write("  transport:\n")
            f.write("    model: gas\n")
            f.write(f"    geometry: {geometry}\n")
            f.write(f"    diameter: {fmt(tr['sigma'])}\n")
            f.write(f"    well-depth: {fmt(tr['eps'])}\n")
            if tr["dipole"] != 0.0:
                f.write(f"    dipole: {fmt(tr['dipole'])}\n")
            if tr["polar"] != 0.0:
                f.write(f"    polarizability: {fmt(tr['polar'])}\n")
            if tr["rot"] != 0.0:
                f.write(f"    rotational-relaxation: {fmt(tr['rot'])}\n")


# ---------------------------------------------------------------------------------------------------
def convert_mechanism(mech_file: Path, thermo_file: Path, trans_file: Path, yaml_file: Path,
                      fit_range=(300.0, 3000.0), fit_tol: float = 0.10) -> None:
    transport = parse_transport_file(trans_file) if trans_file.exists() else {}
    # FlameMaster writes a dummy 0 in the CHEMKIN geometry column; only trust it if it is actually used
    use_geometry_column = any(t["geom"] not in (None, 0.0) for t in transport.values())
    thermo_species = parse_thermo_file(thermo_file, transport)
    reactions, third_bodies, additional, E_to_cal = parse_mech_file(mech_file)

    # The mechanism's species are those in its reactions plus 'Let additional species be ...' (inerts);
    # the thermo file may be a superset. Keep the thermo file order.
    used = set(additional)
    for r in reactions:
        for _, nm in r["left"] + r["right"]:
            used.add(nm)
    thermo_names = [sp["name"] for sp in thermo_species]
    unknown = sorted(used - set(thermo_names))
    if unknown:
        raise FMError(f"Species used in the mechanism but absent from the thermo file: {', '.join(unknown)}")
    species = [sp for sp in thermo_species if sp["name"] in used]
    species_names = [sp["name"] for sp in species]
    if len(species) < len(thermo_species):
        print(f"  Using {len(species)} of the {len(thermo_species)} thermo-file species (those in the mechanism)")
    missing_trans = [sp["name"] for sp in species if sp["trans"] is None]
    if transport and missing_trans:
        print(f"  Warning: no transport data for: {', '.join(missing_trans)}")
    elements = []
    for sp in species:
        for e in sp["comp"]:
            if e not in elements:
                elements.append(e)

    # Build the whole document first so a conversion error never leaves a truncated YAML behind
    f = io.StringIO()
    if True:
        f.write("description: |-\n")
        f.write("  Converted from FlameMaster format by FM2yaml.py\n")
        f.write(f"  Mechanism file: {mech_file.name}\n")
        f.write(f"  Thermo file: {thermo_file.name}\n")
        if trans_file.exists():
            f.write(f"  Transport file: {trans_file.name}\n")
        f.write("\n")
        f.write("generator: FM2yaml\n")
        f.write(f"date: {datetime.now().strftime('%a, %d %b %Y %H:%M:%S')}\n")
        f.write("input-files: [" + ", ".join(p.name for p in (mech_file, thermo_file, trans_file) if p.exists()) + "]\n")
        f.write("\n")
        f.write("units: {length: cm, quantity: mol, activation-energy: cal/mol}\n")
        f.write("\n")
        f.write("phases:\n")
        f.write("- name: gas\n")
        f.write("  thermo: ideal-gas\n")
        f.write("  elements: [" + ", ".join(elements) + "]\n")
        f.write("  species: [" + ", ".join(species_names) + "]\n")
        f.write("  kinetics: gas\n")
        if transport and not missing_trans:
            f.write("  transport: mixture-averaged\n")
        f.write("  reactions: all\n")
        f.write("  state:\n")
        f.write("    T: 300.0\n")
        f.write("    P: 1.01325e+05\n")
        f.write("\n")
        write_species(f, species, use_geometry_column)
        f.write("\n")
        n_out = write_reactions(f, reactions, third_bodies, set(species_names), E_to_cal, fit_range, fit_tol)
    with open(yaml_file, "w") as out:
        out.write(f.getvalue())
    print(f"  Parsed {len(reactions)} reactions -> {n_out} YAML reactions")


def main():
    parser = argparse.ArgumentParser(description="FM2yaml - FlameMaster to Cantera-YAML converter")
    parser.add_argument("mech_file", type=Path, help="Mechanism file (.mech)")
    parser.add_argument("thermo_file", type=Path, nargs="?", default=None, help="Thermo file. Default: <mech_base>.thermo or .chthermo")
    parser.add_argument("trans_file", type=Path, nargs="?", default=None, help="Transport file. Default: <mech_base>.trans or .chtrans")
    parser.add_argument("output", type=Path, nargs="?", default=None, help="Output YAML file. Default: <mech_base>.yaml")
    parser.add_argument("--fit-range", type=float, nargs=2, metavar=("TMIN", "TMAX"), default=[300.0, 3000.0],
                        help="Temperature range of the Troe fit used for centering functions that Cantera's"
                             " Troe form cannot hold exactly (default: 300 3000 K)")
    parser.add_argument("--fit-tol", type=float, default=10.0, metavar="PERCENT",
                        help="Largest accepted F_cent error of such a fit, in %% (default: 10)")
    args = parser.parse_args()
    if args.fit_range[0] <= 0.0 or args.fit_range[1] <= args.fit_range[0]:
        print("ERROR: --fit-range needs 0 < TMIN < TMAX", file=sys.stderr)
        sys.exit(1)

    mech_file = args.mech_file
    if not mech_file.exists():
        print(f"ERROR: Mechanism file not found: {mech_file}", file=sys.stderr)
        sys.exit(1)
    # strip only the real extension: 'CH4.Igni73.mech' -> 'CH4.Igni73' (Path.with_suffix would give 'CH4')
    base = str(mech_file)[: -len(mech_file.suffix)] if mech_file.suffix else str(mech_file)

    def default_file(explicit, *suffixes):
        if explicit is not None:
            return explicit
        for sfx in suffixes:
            if Path(base + sfx).exists():
                return Path(base + sfx)
        return Path(base + suffixes[0])

    thermo_file = default_file(args.thermo_file, ".thermo", ".chthermo")
    trans_file = default_file(args.trans_file, ".trans", ".chtrans")
    yaml_file = args.output or Path(base + ".yaml")

    print("FM2yaml - FlameMaster to Cantera-YAML converter")
    print(f"Mechanism file: {mech_file}")
    print(f"Thermo file:    {thermo_file}")
    print(f"Transport file: {trans_file if trans_file.exists() else '(not found, skipping)'}")
    print(f"Output file:    {yaml_file}")
    if not thermo_file.exists():
        print(f"ERROR: Thermo file not found: {thermo_file}", file=sys.stderr)
        sys.exit(1)
    try:
        convert_mechanism(mech_file, thermo_file, trans_file, yaml_file, tuple(args.fit_range), args.fit_tol / 100.0)
    except FMError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    print("Conversion complete!")


if __name__ == "__main__":
    main()
