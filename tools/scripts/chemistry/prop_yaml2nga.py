#!/usr/bin/env python3
"""
Generate Fortran fluid property file from a Peng-Robinson fluid YAML.

YAML lists one row per component.
Non-condensables (NCG) are flagged with ncg: true (vapor-only in flash).

Usage:
  python prop_yaml2nga.py <properties.yaml> [output.f90]
  python prop_yaml2nga.py kinetics/properties.yaml src/prop_data.f90
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
    check_fortran_limits,
    f90_1d_data_variable,
    f90_2d_array_parameter as _f90_2d_array_parameter,
    f90_array_parameter as _f90_array_parameter,
    f90_real,
    load_yaml,
    species_thermo,
)


def _f90_logical(v) -> str:
    return ".true." if v else ".false."

R_UNIV_SI = 8.314462618  # J/(mol K)
_OMEGA_MIN = -0.5
_OMEGA_MAX = 2.0


def f90_array_parameter(name, values, dim_name=None, max_per_line=4):
    return _f90_array_parameter(name, values, dim_expr=dim_name, per_line=max_per_line)


def _fort_dim_expr(fort_dim, nrows, ncols):
    if fort_dim:
        return ", ".join(p.strip() for p in fort_dim.split(","))
    return f"{nrows}, {ncols}"


def f90_real_2d_parameter(name, arr, nrows, ncols, max_per_line=4, fort_dim=None):
    return _f90_2d_array_parameter(name, arr, nrows, ncols, per_line=max_per_line,
                                   dim_expr=_fort_dim_expr(fort_dim, nrows, ncols))


def fortran_safe_symbol(sym: str) -> str:
    s = re.sub(r"[^0-9A-Za-z_]", "_", str(sym))
    if s and s[0].isdigit():
        s = "e_" + s
    return s or "X"


def canon_species_name(name: str) -> str:
    return str(name).strip()


def _require_float(val, field: str, species: str) -> float:
    if val is None:
        raise ValueError(f"species '{species}': missing '{field}'")
    try:
        x = float(val)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"species '{species}': '{field}' must be a number") from exc
    if not np.isfinite(x):
        raise ValueError(f"species '{species}': '{field}' must be finite")
    return x


def extract_nasa7_thermo(spec: dict) -> tuple[float, list[float], list[float]]:
    """
    (T_switch, low[7], high[7]) from yaml2nga's validated NASA7 reader (exactly 7 coefficients per row,
    increasing bounds, cp > 0 over the range, T_mid discontinuity warning; single range -> row repeated).
    """
    try:
        low, high, _, T_mid, _ = species_thermo(spec, str(spec.get("name")))
    except MechanismError as exc:
        raise ValueError(str(exc)) from None
    return T_mid, list(low), list(high)


def nasa7_h_mol(T: float, t_switch: float, lo: list[float], hi: list[float]) -> float:
    """NASA7 molar enthalpy: h/(RT) = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T (a7 is the entropy constant)."""
    a = lo if T <= t_switch else hi
    return R_UNIV_SI * T * (
        a[0]
        + a[1] * T / 2.0
        + a[2] * T**2 / 3.0
        + a[3] * T**3 / 4.0
        + a[4] * T**4 / 5.0
        + a[5] / T
    )


def load_pr_species(path: Path) -> dict:
    data = load_yaml(path)

    if not isinstance(data, dict):
        raise ValueError("YAML root must be a mapping")

    u_ig_Tref = float(data.get("ideal-gas-reference-temperature", 298.15))
    if not u_ig_Tref > 0.0:
        raise ValueError("'ideal-gas-reference-temperature' must be > 0 K")

    raw = data.get("species", [])
    if not raw:
        raise ValueError("YAML must contain a non-empty 'species' list")

    species_dicts = [s for s in raw if isinstance(s, dict)]
    if len(species_dicts) != len(raw):
        raise ValueError("every 'species' entry must be a mapping")

    names = [canon_species_name(s.get("name", "")) for s in species_dicts]
    if any(not nm for nm in names):
        raise ValueError("every species must have a non-empty 'name'")
    if len(names) != len(set(names)):
        raise ValueError("duplicate species names in YAML")

    nS = len(names)
    name_to_idx = {nm: i for i, nm in enumerate(names)}

    Tc = np.zeros(nS)
    Pc = np.zeros(nS)
    omega = np.zeros(nS)
    MM = np.zeros(nS)
    c_pen = np.zeros(nS)
    is_ncg = np.zeros(nS, dtype=bool)
    thermo_rows: list[tuple[float, list[float], list[float]]] = []
    thermo_data = np.zeros((nS, 15))
    n_cond = 0

    for i, s in enumerate(species_dicts):
        nm = names[i]
        crit = s.get("critical") or {}
        if not isinstance(crit, dict):
            raise ValueError(f"species '{nm}': 'critical' must be a mapping")

        Tc[i] = _require_float(crit.get("temperature"), "critical.temperature", nm)
        Pc[i] = _require_float(crit.get("pressure"), "critical.pressure", nm)
        omega[i] = _require_float(crit.get("acentric-factor"), "critical.acentric-factor", nm)
        MM[i] = _require_float(s.get("molecular-weight"), "molecular-weight", nm)

        pen = s.get("peneloux-volume-correction", 0.0)
        if pen is not None:
            c_pen[i] = float(pen)

        ts, lo, hi = extract_nasa7_thermo(s)
        thermo_rows.append((ts, lo, hi))
        thermo_data[i, 0] = ts
        thermo_data[i, 1:8] = lo
        thermo_data[i, 8:15] = hi

        if Tc[i] <= 0.0:
            raise ValueError(f"species '{nm}': critical.temperature must be > 0")
        if Pc[i] <= 0.0:
            raise ValueError(f"species '{nm}': critical.pressure must be > 0")
        if MM[i] <= 0.0:
            raise ValueError(f"species '{nm}': molecular-weight must be > 0")
        if omega[i] < _OMEGA_MIN or omega[i] > _OMEGA_MAX:
            print(
                f"Warning: species '{nm}': acentric-factor {omega[i]} outside [{_OMEGA_MIN}, {_OMEGA_MAX}]",
                file=sys.stderr,
            )

        ncg = s.get("ncg", False)
        if isinstance(ncg, str):
            ncg = ncg.strip().lower() in ("true", "t", "yes", "1")
        is_ncg[i] = bool(ncg)
        if not is_ncg[i]:
            n_cond += 1

    if n_cond < 1:
        raise ValueError("at least one species must have ncg: false (condensable)")

    kij = np.zeros((nS, nS))
    bi = data.get("binary-interaction") or {}
    kij_list = bi.get("kij", []) if isinstance(bi, dict) else []
    seen_pairs: set[tuple[int, int]] = set()

    for entry in kij_list:
        if not isinstance(entry, (list, tuple)) or len(entry) != 3:
            raise ValueError("binary-interaction.kij entries must be [name_i, name_j, value]")
        ni, nj, val = canon_species_name(entry[0]), canon_species_name(entry[1]), entry[2]
        if ni not in name_to_idx or nj not in name_to_idx:
            raise ValueError(f"kij pair ({ni}, {nj}): unknown species name")
        if ni == nj:
            raise ValueError(f"kij pair ({ni}, {nj}): diagonal must not appear in YAML")
        try:
            kval = float(val)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"kij ({ni}, {nj}): value must be a number") from exc
        ii, jj = name_to_idx[ni], name_to_idx[nj]
        pair = (min(ii, jj), max(ii, jj))
        if pair in seen_pairs:
            raise ValueError(f"duplicate kij entry for ({ni}, {nj})")
        seen_pairs.add(pair)
        kij[ii, jj] = kval
        kij[jj, ii] = kval

    nNCG = int(np.sum(is_ncg))
    i_first_cond = int(np.argmax(~is_ncg)) + 1  # 1-based Fortran index

    Cv_ig = np.zeros(nS)
    for i, nm in enumerate(names):
        ts, lo, hi = thermo_rows[i]
        Tcv = 300.0  # NOTE: the emitted Cv_ig is a single constant evaluated here (design of the vle consumers)
        a = lo if Tcv <= ts else hi
        cp_over_R = a[0] + a[1] * Tcv + a[2] * Tcv**2 + a[3] * Tcv**3 + a[4] * Tcv**4
        Cv_ig[i] = cp_over_R * R_UNIV_SI - R_UNIV_SI

    return {
        "names": names,
        "nS": nS,
        "nNCG": nNCG,
        "n_cond": n_cond,
        "i_first_cond": i_first_cond,
        "Tc": Tc,
        "Pc": Pc,
        "omega": omega,
        "MM": MM,
        "c_pen": c_pen,
        "is_ncg": is_ncg,
        "Cv_ig": Cv_ig,
        "kij": kij,
        "thermo_data": thermo_data,
        "u_ig_Tref": u_ig_Tref,
        "yaml_source": path.name,
    }


def emit_fortran(mech: dict, out_path: Path) -> None:
    names = mech["names"]
    nS = mech["nS"]
    nNCG = mech["nNCG"]

    lines = []
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append(f"!  FILE {out_path.name}")
    lines.append("!  Peng-Robinson fluid data")
    lines.append("!  Generated by prop_yaml2nga.py. Do not edit manually.")
    lines.append(f"!  Source: {mech['yaml_source']}")
    if mech.get("yaml_relpath"):
        lines.append(f"!  YAML:   {mech['yaml_relpath']}")
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append("!  Arrays are indexed 1..nS_pr in YAML species list order.")
    lines.append("!  is_ncg_pr(i) = .true. => non-condensable (vapor-only, x_i^L = 0).")
    lines.append("!  kij_pr(i,j) = kij_pr(j,i); diagonal zero; van der Waals one-fluid mixing.")
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append("")
    lines.append("module vle_pr_data")
    lines.append("   use precision, only: WP")
    lines.append("   implicit none")
    lines.append("   private")
    lines.append("   public :: nS_pr, nNCG_pr, n_cond_pr, i_first_cond_pr, is_ncg_pr")
    lines.append("   public :: Tc_pr, Pc_pr, omega_pr, MM, Cv_ig_pr, kij_pr")
    lines.append("   public :: Rcst, vle_pr_species_name, vle_pr_is_condensable")
    lines.append("   public :: vle_pr_z_ncg, vle_pr_z_condensable, vle_pr_is_pure_condensable")
    lines.append("   public :: vle_pr_is_cond_ncg_feed, vle_pr_i_dominant_cond")
    lines.append("   public :: vle_pr_set_pure_species, vle_pr_set_pure_ncg, vle_pr_set_feed_cond_ncg")
    lines.append("   public :: vle_pr_i_pure_species")
    for nm in names:
        sym = fortran_safe_symbol(nm)
        lines.append(f"   public :: s{sym}_pr")
    lines.append("")
    lines.append(f"   integer, parameter :: nS_pr = {nS}")
    lines.append(f"   integer, parameter :: nNCG_pr = {nNCG}")
    lines.append(f"   integer, parameter :: n_cond_pr = {mech['n_cond']}")
    lines.append(f"   integer, parameter :: i_first_cond_pr = {mech['i_first_cond']}")
    lines.append("")
    lines.append(f"   real(WP), parameter :: Rcst = {f90_real(R_UNIV_SI)}  ! J/(mol K)")
    lines.append("")
    lines.append("   ! --- Named species indices (YAML list order; use only in generated code) ---")
    for i, nm in enumerate(names):
        sym = fortran_safe_symbol(nm)
        ncg_note = "NCG" if mech["is_ncg"][i] else "condensable"
        lines.append(f"   integer, parameter :: s{sym}_pr = {i + 1}  ! {nm} ({ncg_note})")
    lines.append("")

    lines.extend(f90_array_parameter("Tc_pr", mech["Tc"], dim_name="nS_pr"))
    lines.append("")
    lines.extend(f90_array_parameter("Pc_pr", mech["Pc"], dim_name="nS_pr"))
    lines.append("")
    lines.extend(f90_array_parameter("omega_pr", mech["omega"], dim_name="nS_pr"))
    lines.append("")
    lines.extend(f90_array_parameter("MM", mech["MM"], dim_name="nS_pr"))
    lines.append("")
    lines.extend(f90_array_parameter("Cv_ig_pr", mech["Cv_ig"], dim_name="nS_pr"))
    lines.append("")

    lines.extend(_f90_array_parameter("is_ncg_pr", mech["is_ncg"], kind="logical", fmt=_f90_logical, per_line=8, dim_expr="nS_pr"))
    lines.append("")
    lines.extend(f90_real_2d_parameter("kij_pr", mech["kij"], nS, nS, fort_dim="nS_pr, nS_pr"))
    lines.append("")
    lines.append("contains")
    lines.append("")
    max_len = max(len(nm) for nm in names)
    lines.append("   pure function vle_pr_species_name(i) result(nm)")
    lines.append("      integer, intent(in) :: i")
    lines.append(f"      character(len={max_len}) :: nm")
    lines.append("      select case (i)")
    for i, nm in enumerate(names):
        q = nm.replace("'", "''")
        lines.append(f"      case ({i + 1}); nm = '{q}'")
    lines.append("      case default; nm = '?'")
    lines.append("      end select")
    lines.append("   end function vle_pr_species_name")
    lines.append("")
    lines.append("   pure function vle_pr_is_condensable(i) result(cond)")
    lines.append("      integer, intent(in) :: i")
    lines.append("      logical :: cond")
    lines.append("      if (i < 1 .or. i > nS_pr) then")
    lines.append("         cond = .false.")
    lines.append("      else")
    lines.append("         cond = .not. is_ncg_pr(i)")
    lines.append("      end if")
    lines.append("   end function vle_pr_is_condensable")
    lines.append("")
    lines.append("   pure function vle_pr_z_ncg(z) result(zn)")
    lines.append("      real(WP), intent(in) :: z(nS_pr)")
    lines.append("      real(WP) :: zn")
    lines.append("      integer :: i")
    lines.append("      zn = 0.0_WP")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (is_ncg_pr(i)) zn = zn + z(i)")
    lines.append("      end do")
    lines.append("   end function vle_pr_z_ncg")
    lines.append("")
    lines.append("   pure function vle_pr_z_condensable(z) result(zc)")
    lines.append("      real(WP), intent(in) :: z(nS_pr)")
    lines.append("      real(WP) :: zc")
    lines.append("      integer :: i")
    lines.append("      zc = 0.0_WP")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (.not. is_ncg_pr(i)) zc = zc + z(i)")
    lines.append("      end do")
    lines.append("   end function vle_pr_z_condensable")
    lines.append("")
    lines.append("   pure function vle_pr_is_pure_condensable(z) result(pure_c)")
    lines.append("      real(WP), intent(in) :: z(nS_pr)")
    lines.append("      logical :: pure_c")
    lines.append("      pure_c = vle_pr_z_condensable(z) >= 1.0_WP - 1.0e-10_WP")
    lines.append("   end function vle_pr_is_pure_condensable")
    lines.append("")
    lines.append("   pure function vle_pr_is_cond_ncg_feed(z) result(mix)")
    lines.append("      real(WP), intent(in) :: z(nS_pr)")
    lines.append("      logical :: mix")
    lines.append("      mix = (vle_pr_z_condensable(z) > 1.0e-15_WP) .and. (vle_pr_z_ncg(z) > 1.0e-15_WP)")
    lines.append("   end function vle_pr_is_cond_ncg_feed")
    lines.append("")
    lines.append("   pure function vle_pr_i_dominant_cond(z) result(ic)")
    lines.append("      real(WP), intent(in) :: z(nS_pr)")
    lines.append("      integer :: ic, i")
    lines.append("      real(WP) :: zmax")
    lines.append("      ic = i_first_cond_pr")
    lines.append("      zmax = 0.0_WP")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (.not. is_ncg_pr(i) .and. z(i) > zmax) then")
    lines.append("            zmax = z(i)")
    lines.append("            ic = i")
    lines.append("         end if")
    lines.append("      end do")
    lines.append("   end function vle_pr_i_dominant_cond")
    lines.append("")
    lines.append("   pure subroutine vle_pr_set_pure_species(z, is)")
    lines.append("      real(WP), intent(out) :: z(nS_pr)")
    lines.append("      integer, intent(in) :: is")
    lines.append("      integer :: i")
    lines.append("      z = 0.0_WP")
    lines.append("      if (is >= 1 .and. is <= nS_pr) z(is) = 1.0_WP")
    lines.append("   end subroutine vle_pr_set_pure_species")
    lines.append("")
    lines.append("   pure subroutine vle_pr_set_pure_ncg(z)")
    lines.append("      real(WP), intent(out) :: z(nS_pr)")
    lines.append("      integer :: i")
    lines.append("      z = 0.0_WP")
    lines.append("      if (nNCG_pr < 1) return")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (is_ncg_pr(i)) z(i) = 1.0_WP / real(nNCG_pr, WP)")
    lines.append("      end do")
    lines.append("   end subroutine vle_pr_set_pure_ncg")
    lines.append("")
    lines.append("   pure subroutine vle_pr_set_feed_cond_ncg(z, frac_ncg)")
    lines.append("      real(WP), intent(out) :: z(nS_pr)")
    lines.append("      real(WP), intent(in) :: frac_ncg")
    lines.append("      real(WP) :: fc, fn, zc, zn")
    lines.append("      integer :: i, nc, nn")
    lines.append("      z = 0.0_WP")
    lines.append("      fc = max(0.0_WP, min(1.0_WP, frac_ncg))")
    lines.append("      zc = 1.0_WP - fc")
    lines.append("      nc = 0")
    lines.append("      nn = 0")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (is_ncg_pr(i)) then")
    lines.append("            nn = nn + 1")
    lines.append("         else")
    lines.append("            nc = nc + 1")
    lines.append("         end if")
    lines.append("      end do")
    lines.append("      if (nc < 1 .or. nn < 1) return")
    lines.append("      zn = fc / real(nn, WP)")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (is_ncg_pr(i)) z(i) = zn")
    lines.append("      end do")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (.not. is_ncg_pr(i)) z(i) = zc / real(nc, WP)")
    lines.append("      end do")
    lines.append("   end subroutine vle_pr_set_feed_cond_ncg")
    lines.append("")
    lines.append("   pure function vle_pr_i_pure_species(x) result(ic)")
    lines.append("      real(WP), intent(in) :: x(nS_pr)")
    lines.append("      integer :: ic, i")
    lines.append("      ic = 0")
    lines.append("      do i = 1, nS_pr")
    lines.append("         if (x(i) >= 1.0_WP - 1.0e-10_WP) ic = i")
    lines.append("      end do")
    lines.append("   end function vle_pr_i_pure_species")
    lines.append("")
    lines.append("end module vle_pr_data")
    lines.append("")

    check_fortran_limits(lines)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Wrote {out_path} (nS_pr={nS}, nNCG_pr={nNCG})")


def _fortran_comp_names(names) -> tuple[int, list[str]]:
    """Quoted, blank-padded character literals for comp_name; width = longest name (>= 4), no truncation."""
    width = max(4, max(len(n.strip()) for n in names))
    lits = [f"'{n.strip():<{width}}'" for n in names]
    if len(set(lits)) != len(lits):
        raise ValueError("component names are not unique")
    return width, lits


def emit_vle_data_fortran(mech: dict, out_path: Path) -> None:
    names = mech["names"]
    nC = mech["nS"]
    nNCG = mech["nNCG"]
    # The vle_data index ranges (sGNCmin..sGNCmax, sLNCmin..sLNCmax) assume that all non-condensable
    # species come after all condensable ones in the YAML species list.
    is_ncg = list(mech["is_ncg"])
    first_ncg = is_ncg.index(True) if True in is_ncg else nC
    if not all(is_ncg[first_ncg:]):
        raise ValueError("vle-data target: species with 'ncg: true' must be listed after all condensable species")
    comp_width, comp_literals = _fortran_comp_names(names)

    lines = []
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append(f"!  FILE {out_path.name}")
    lines.append("!  Peng-Robinson fluid constants for standalone UV flash program.")
    lines.append("!  Generated by prop_yaml2nga.py. Do not edit manually.")
    lines.append(f"!  Source: {mech['yaml_source']}")
    if mech.get("yaml_relpath"):
        lines.append(f"!  YAML:   {mech['yaml_relpath']}")
    lines.append("!  Ideal-gas part: u_ig(T) = [h_NASA(T) - R T] - [h_NASA(T_ref) - R T_ref]; PR departure includes (Z-1)RT (ANN_VLE).")
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append("")
    lines.append("module vle_data")
    lines.append("   use precision, only: WP")
    lines.append("   implicit none")
    lines.append("   private")
    lines.append("")
    lines.append("   public :: nC, nS, nNCG, nThermoDataCols, sGmin, sGmax, sLmin, sLmax")
    lines.append("   public :: sGNCmin, sGNCmax, sLNCmin, sLNCmax")
    lines.append("   public :: Rcst, u_ig_Tref, comp_name, Tc, Pc, omega, MM, Cv_ig, c_pen, is_ncg, kij")
    lines.append("   public :: thermo_data, vle_nasa7_h_mol, vle_u_ig_mol")
    lines.append("   public :: vle_ic_from_is, vle_s_gas, vle_s_liq, vle_is_gas_species, vle_is_liq_species")
    lines.append("   public :: vle_component_name, vle_Y_pack, vle_Y_unpack, vle_species_name")
    lines.append("   public :: vle_is_condensable, vle_z_ncg, vle_z_condensable, vle_is_pure_condensable")
    lines.append("   public :: vle_is_cond_ncg_feed, vle_i_first_condensable, vle_i_dominant_cond")
    lines.append("   public :: vle_i_ref_condensable, vle_i_pure_species, vle_set_pure_species")
    lines.append("   public :: vle_set_pure_ncg, vle_set_feed_cond_ncg, vle_cmol_pen")
    lines.append("   public :: vle_c_pen, vle_set_peneloux, vle_use_peneloux")
    lines.append("")
    lines.append("   logical, save :: use_peneloux_flag = .true.")
    lines.append("")
    lines.append(f"   integer, parameter :: nC = {nC}")
    lines.append("   integer, parameter :: nS = 2 * nC")
    lines.append(f"   integer, parameter :: nNCG = {nNCG}")
    lines.append("   integer, parameter :: nThermoDataCols = 15")
    lines.append("")
    lines.append("   integer, parameter :: sGmin = 1")
    lines.append("   integer, parameter :: sGmax = nC")
    lines.append("   integer, parameter :: sLmin = nC + 1")
    lines.append("   integer, parameter :: sLmax = 2 * nC")
    lines.append("   integer, parameter :: sGNCmin = nC - nNCG + 1")
    lines.append("   integer, parameter :: sGNCmax = nC")
    lines.append("   integer, parameter :: sLNCmin = nC + nC - nNCG + 1")
    lines.append("   integer, parameter :: sLNCmax = 2 * nC")
    lines.append("")
    lines.append(f"   real(WP), parameter :: Rcst = {f90_real(R_UNIV_SI)}  ! J/(mol K)")
    lines.append(f"   real(WP), parameter :: u_ig_Tref = {f90_real(mech['u_ig_Tref'])}  ! K")
    lines.append("")
    lines.extend(_f90_array_parameter("comp_name", comp_literals, kind=f"character(len={comp_width})", fmt=str,
                                      per_line=6, dim_expr="nC"))
    lines.append("")
    lines.extend(f90_array_parameter("Tc", mech["Tc"], dim_name="nC"))
    lines.append("")
    lines.extend(f90_array_parameter("Pc", mech["Pc"], dim_name="nC"))
    lines.append("")
    lines.extend(f90_array_parameter("omega", mech["omega"], dim_name="nC"))
    lines.append("")
    lines.extend(f90_array_parameter("MM", mech["MM"], dim_name="nC"))
    lines.append("")
    lines.extend(f90_array_parameter("Cv_ig", mech["Cv_ig"], dim_name="nC"))
    lines.append("")
    lines.extend(f90_array_parameter("c_pen", mech["c_pen"], dim_name="nC"))
    lines.append("")
    lines.extend(f90_1d_data_variable("c_pen_eff", mech["c_pen"], dim_expr="nC", attrs="save"))
    lines.append("")
    lines.extend(_f90_array_parameter("is_ncg", mech["is_ncg"], kind="logical", fmt=_f90_logical, per_line=8, dim_expr="nC"))
    lines.append("")
    lines.extend(
        f90_real_2d_parameter(
            "thermo_data", mech["thermo_data"], nC, 15, fort_dim="nC, nThermoDataCols"
        )
    )
    lines.append("")
    lines.extend(f90_real_2d_parameter("kij", mech["kij"], nC, nC, fort_dim="nC, nC"))
    lines.append("")
    lines.append("contains")
    lines.append("")
    # --- helper functions (unchanged layout) ---
    lines.extend(_VLE_DATA_HELPERS_TEMPLATE.format(nC=nC).splitlines())
    lines.append("")
    lines.append("end module vle_data")
    lines.append("")

    check_fortran_limits(lines)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Wrote {out_path} (nC={nC}, nNCG={nNCG})")


_VLE_DATA_HELPERS_TEMPLATE = """
   pure subroutine vle_nasa7_coeffs(is, T, a)
      integer, intent(in) :: is
      real(WP), intent(in) :: T
      real(WP), intent(out) :: a(7)
      integer :: j
      if (is < 1 .or. is > nC) then
         a = 0.0_WP
         return
      end if
      if (T <= thermo_data(is, 1)) then
         do j = 1, 7
            a(j) = thermo_data(is, j + 1)
         end do
      else
         do j = 1, 7
            a(j) = thermo_data(is, j + 8)
         end do
      end if
   end subroutine vle_nasa7_coeffs

   pure function vle_nasa7_h_mol(T, is) result(h)
      real(WP), intent(in) :: T
      integer, intent(in) :: is
      real(WP) :: h, a(7)
      if (T <= 0.0_WP .or. is < 1 .or. is > nC) then
         h = 0.0_WP
         return
      end if
      call vle_nasa7_coeffs(is, T, a)
      ! NASA7: h/(RT) = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T  (a7 is the entropy constant)
      h = Rcst * T * (a(1) + a(2) * T / 2.0_WP + a(3) * T**2 / 3.0_WP &
         + a(4) * T**3 / 4.0_WP + a(5) * T**4 / 5.0_WP + a(6) / T)
   end function vle_nasa7_h_mol

   pure function vle_u_ig_mol(T, is) result(uig)
      real(WP), intent(in) :: T
      integer, intent(in) :: is
      real(WP) :: uig, h, href, Tref
      Tref = u_ig_Tref
      if (T <= 0.0_WP .or. is < 1 .or. is > nC) then
         uig = 0.0_WP
         return
      end if
      ! ideal-gas internal energy u = h - R*T, relative to T_ref
      h = vle_nasa7_h_mol(T, is)
      href = vle_nasa7_h_mol(Tref, is)
      uig = (h - Rcst * T) - (href - Rcst * Tref)
   end function vle_u_ig_mol

   pure function vle_ic_from_is(is) result(ic)
      integer, intent(in) :: is
      integer :: ic
      if (is >= sGmin .and. is <= sGmax) then
         ic = is
      else if (is >= sLmin .and. is <= sLmax) then
         ic = is - nC
      else
         ic = 0
      end if
   end function vle_ic_from_is

   pure function vle_s_gas(ic) result(is)
      integer, intent(in) :: ic
      integer :: is
      is = ic
   end function vle_s_gas

   pure function vle_s_liq(ic) result(is)
      integer, intent(in) :: ic
      integer :: is
      is = nC + ic
   end function vle_s_liq

   pure function vle_is_gas_species(is) result(gas)
      integer, intent(in) :: is
      logical :: gas
      gas = (is >= sGmin .and. is <= sGmax)
   end function vle_is_gas_species

   pure function vle_is_liq_species(is) result(liq)
      integer, intent(in) :: is
      logical :: liq
      liq = (is >= sLmin .and. is <= sLmax)
   end function vle_is_liq_species

   pure function vle_component_name(ic) result(nm)
      integer, intent(in) :: ic
      character(len=len(comp_name)) :: nm
      if (ic >= 1 .and. ic <= nC) then
         nm = comp_name(ic)
      else
         nm = '?'
      end if
   end function vle_component_name

   pure subroutine vle_Y_pack(Y, xL, yV)
      real(WP), intent(out) :: Y(nS)
      real(WP), intent(in) :: xL(nC), yV(nC)
      Y(sGmin:sGmax) = yV
      Y(sLmin:sLmax) = xL
   end subroutine vle_Y_pack

   pure subroutine vle_Y_unpack(Y, xL, yV)
      real(WP), intent(in) :: Y(nS)
      real(WP), intent(out) :: xL(nC), yV(nC)
      yV = Y(sGmin:sGmax)
      xL = Y(sLmin:sLmax)
   end subroutine vle_Y_unpack

   pure function vle_species_name(is) result(nm)
      integer, intent(in) :: is
      character(len=len(comp_name)+2) :: nm
      integer :: ic
      ic = vle_ic_from_is(is)
      if (ic < 1 .or. ic > nC) then
         nm = '?'
         return
      end if
      if (vle_is_gas_species(is)) then
         nm = trim(comp_name(ic)) // '-G'
      else
         nm = trim(comp_name(ic)) // '-L'
      end if
   end function vle_species_name

   pure function vle_is_condensable(ic) result(cond)
      integer, intent(in) :: ic
      logical :: cond
      if (ic < 1 .or. ic > nC) then
         cond = .false.
      else
         cond = .not. is_ncg(ic)
      end if
   end function vle_is_condensable

   pure function vle_z_ncg(z) result(zn)
      real(WP), intent(in) :: z(nC)
      real(WP) :: zn
      integer :: ic
      zn = 0.0_WP
      do ic = 1, nC
         if (is_ncg(ic)) zn = zn + z(ic)
      end do
   end function vle_z_ncg

   pure function vle_z_condensable(z) result(zc)
      real(WP), intent(in) :: z(nC)
      real(WP) :: zc
      integer :: ic
      zc = 0.0_WP
      do ic = 1, nC
         if (.not. is_ncg(ic)) zc = zc + z(ic)
      end do
   end function vle_z_condensable

   pure function vle_is_pure_condensable(z) result(pure_c)
      real(WP), intent(in) :: z(nC)
      logical :: pure_c
      pure_c = vle_z_condensable(z) >= 1.0_WP - 1.0e-10_WP
   end function vle_is_pure_condensable

   pure function vle_is_cond_ncg_feed(z) result(mix)
      real(WP), intent(in) :: z(nC)
      logical :: mix
      mix = (vle_z_condensable(z) > 1.0e-15_WP) .and. (vle_z_ncg(z) > 1.0e-15_WP)
   end function vle_is_cond_ncg_feed

   pure function vle_i_first_condensable() result(ic)
      integer :: ic, i
      ic = 0
      do i = 1, nC
         if (.not. is_ncg(i)) then
            ic = i
            return
         end if
      end do
   end function vle_i_first_condensable

   pure function vle_i_dominant_cond(z) result(ic)
      real(WP), intent(in) :: z(nC)
      integer :: ic, i
      real(WP) :: zmax
      ic = 0
      zmax = 0.0_WP
      do i = 1, nC
         if (.not. is_ncg(i)) then
            if (ic == 0) ic = i
            if (z(i) > zmax) then
               zmax = z(i)
               ic = i
            end if
         end if
      end do
      if (ic == 0) ic = vle_i_first_condensable()
   end function vle_i_dominant_cond

   pure function vle_i_ref_condensable(z) result(ic)
      real(WP), intent(in) :: z(nC)
      integer :: ic
      ic = vle_i_pure_species(z)
      if (ic > 0 .and. vle_is_condensable(ic)) return
      ic = vle_i_dominant_cond(z)
      if (ic > 0 .and. vle_is_condensable(ic)) return
      ic = vle_i_first_condensable()
   end function vle_i_ref_condensable

   pure function vle_i_pure_species(x) result(ic)
      real(WP), intent(in) :: x(nC)
      integer :: ic, i
      ic = 0
      do i = 1, nC
         if (x(i) >= 1.0_WP - 1.0e-10_WP) ic = i
      end do
   end function vle_i_pure_species

   pure subroutine vle_set_pure_species(z, ic)
      real(WP), intent(out) :: z(nC)
      integer, intent(in) :: ic
      z = 0.0_WP
      if (ic >= 1 .and. ic <= nC) z(ic) = 1.0_WP
   end subroutine vle_set_pure_species

   pure subroutine vle_set_pure_ncg(z)
      real(WP), intent(out) :: z(nC)
      integer :: ic
      z = 0.0_WP
      if (nNCG < 1) return
      do ic = 1, nC
         if (is_ncg(ic)) z(ic) = 1.0_WP / real(nNCG, WP)
      end do
   end subroutine vle_set_pure_ncg

   pure subroutine vle_set_feed_cond_ncg(z, frac_ncg)
      real(WP), intent(out) :: z(nC)
      real(WP), intent(in) :: frac_ncg
      real(WP) :: fc, zc, zn
      integer :: ic, nc, nn
      z = 0.0_WP
      fc = max(0.0_WP, min(1.0_WP, frac_ncg))
      zc = 1.0_WP - fc
      nc = 0
      nn = 0
      do ic = 1, nC
         if (is_ncg(ic)) then
            nn = nn + 1
         else
            nc = nc + 1
         end if
      end do
      if (nc < 1 .or. nn < 1) return
      zn = fc / real(nn, WP)
      do ic = 1, nC
         if (is_ncg(ic)) z(ic) = zn
      end do
      do ic = 1, nC
         if (.not. is_ncg(ic)) z(ic) = zc / real(nc, WP)
      end do
   end subroutine vle_set_feed_cond_ncg

   function vle_cmol_pen(x) result(cmix)
      real(WP), intent(in) :: x(nC)
      real(WP) :: cmix
      integer :: ic
      cmix = 0.0_WP
      do ic = 1, nC
         cmix = cmix + x(ic) * c_pen_eff(ic)
      end do
   end function vle_cmol_pen

   function vle_c_pen(is) result(cpen)
      integer, intent(in) :: is
      real(WP) :: cpen
      if (is >= 1 .and. is <= nC) then
         cpen = c_pen_eff(is)
      else
         cpen = 0.0_WP
      end if
   end function vle_c_pen

   subroutine vle_set_peneloux(flag)
      logical, intent(in) :: flag
      use_peneloux_flag = flag
      if (flag) then
         c_pen_eff = c_pen
      else
         c_pen_eff = 0.0_WP
      end if
   end subroutine vle_set_peneloux

   function vle_use_peneloux() result(flag)
      logical :: flag
      flag = use_peneloux_flag
   end function vle_use_peneloux
"""


def main():
    parser = argparse.ArgumentParser(
        description="Generate Fortran fluid data from Peng-Robinson YAML"
    )
    parser.add_argument("input", type=Path, help="Input YAML (species + binary-interaction)")
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=Path("vle_pr_data.f90"),
        help="Output Fortran file",
    )
    parser.add_argument(
        "--target",
        choices=("vle-pr-data", "vle-data"),
        default="vle-pr-data",
        help="Output module: vle_pr_data (default) or vle_data for vle_flash_uv",
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: input not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    try:
        mech = load_pr_species(args.input)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    try:
        mech["yaml_relpath"] = str(args.input.resolve().relative_to(args.output.resolve().parent))
    except ValueError:
        mech["yaml_relpath"] = str(args.input.resolve())

    print(f"Species: nS={mech['nS']} (nNCG={mech['nNCG']})")
    print("Order:", ", ".join(mech["names"]))
    try:
        if args.target == "vle-data":
            emit_vle_data_fortran(mech, args.output)
        else:
            emit_fortran(mech, args.output)
    except (ValueError, MechanismError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
