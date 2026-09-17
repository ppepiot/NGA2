"""
Noble-Abel Stiffened-Gas (NASG) and Stiffened-Gas (SG) EOS coefficient fitter + plotter.

Implements the methodology described in:
    Le Metayer, O. and Saurel, R., "The Noble-Abel Stiffened-Gas equation of state",
    Physics of Fluids 28, 046102 (2016).

NASG caloric EOS:   p = (gamma - 1) * (e - q) / (v - b)  -  gamma * p_inf
SG caloric EOS:     p = (gamma - 1) * (e - q) / v         -  gamma * p_inf  [b = 0]

This script:
  1. Pulls saturation data from IAPWS-IF97 (via the `iapws` package).
  2. Fits NASG parameters (gamma, Cv, p_inf, b, q, q_prime) for liquid and vapor.
  3. Fits SG parameters (same, but b forced to zero) for liquid and vapor.
  4. Plots a 2x3 saturation-curve figure comparing both EOS against IAPWS data.

Default: liquid water / steam in the range 300-500 K (paper Table V / Flatten-Lund 2011).

Usage:
    python fit_eos.py [--mode nasg|sg|both]

Dependencies: numpy, scipy, matplotlib, iapws  (pip install iapws)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from iapws import IAPWS97

plt.rcParams.update({
    "text.usetex":        True,
    "font.family":        "serif",
    "font.size":          13,
    "axes.labelsize":     14,
    "legend.fontsize":    12,
    "figure.dpi":         150,
    "axes.linewidth":     1.2,
    "xtick.major.size":   5,
    "ytick.major.size":   5,
    "xtick.major.width":  1.2,
    "ytick.major.width":  1.2,
    "xtick.minor.size":   3,
    "ytick.minor.size":   3,
})


# =============================================================================
# Reference data
# =============================================================================

WATER_REFERENCE_STATE = {
    "rho0": 999.9747,  # kg/m^3
    "P0":   1.0100e5,  # Pa
    "c0":   1421.53,   # m/s
}


def iapws_saturation_table(T_min, T_max, n_points=101):
    """
    Build a saturation-curve table from IAPWS-IF97.

    Returns
    -------
    ndarray, shape (N, 6) with columns (T, P_sat, h_l, h_g, v_l, v_g) in SI units.
    Rows where IAPWS97 raises (e.g. above the critical point) are skipped.
    """
    T = np.linspace(T_min, T_max, n_points)
    rows = []
    for Ti in T:
        try:
            liq = IAPWS97(T=Ti, x=0)
            vap = IAPWS97(T=Ti, x=1)
            rows.append([Ti, liq.P * 1e6, liq.h * 1e3, vap.h * 1e3,
                         liq.v, vap.v])
        except Exception:
            pass
    return np.array(rows)


# =============================================================================
# Vapor phase fit (shared by NASG and SG: b_g = 0, P_inf_g = 0)
# =============================================================================

def fit_vapor_phase(T, P, hg, vg):
    """
    Vapor-phase coefficients with b_g = 0 and P_inf_g = 0 (ideal-gas vapor).
    Implements Eqs. (50), (51), (55) of Le Metayer & Saurel (2016).
    """
    T_mean  = T.mean()
    hg_mean = hg.mean()

    # h_g(T) = CP_g * T + q_g  (least squares)
    CP_g = np.sum(T * (hg - hg_mean)) / np.sum(T * (T - T_mean))
    q_g  = hg_mean - CP_g * T_mean

    # v_g(T) = (CP_g - Cv_g) * T / P_sat(T)  (least squares)
    R_g  = np.sum(vg * T / P) / np.sum((T / P) ** 2)
    Cv_g = CP_g - R_g

    return {
        "CP":      CP_g,
        "Cv":      Cv_g,
        "gamma":   CP_g / Cv_g,
        "q":       q_g,
        "b":       0.0,
        "P_inf":   0.0,
    }


# =============================================================================
# NASG liquid phase fit (b_l solved jointly with CP_l - Cv_l from v_l data)
# =============================================================================

def _fit_liquid_nasg_for_given_Pinf(T, P, hl, vl, P_inf_l):
    """
    Given a candidate P_inf,l, solve Eqs. (60)-(61), (64)-(65) for
    (CP_l, Cv_l, q_l, b_l).  Both b_l and CP_l - Cv_l are free.
    """
    T_mean  = T.mean()
    P_mean  = P.mean()
    hl_mean = hl.mean()
    vl_mean = vl.mean()

    TP      = T / (P + P_inf_l)
    TP_mean = TP.mean()

    # Two-unknown least squares on v_l = (CP_l-Cv_l)*T/(P+P_inf) + b_l
    R_l = np.sum(TP * (vl - vl_mean)) / np.sum(TP * (TP - TP_mean))
    b_l = vl_mean - R_l * TP_mean

    # CP_l and q_l from h_l = CP_l*T + b_l*P + q_l
    num_CP = np.sum(T * (hl - hl_mean)) - b_l * np.sum(T * (P - P_mean))
    den_CP = np.sum(T * (T - T_mean))
    CP_l = num_CP / den_CP
    q_l  = hl_mean - CP_l * T_mean - b_l * P_mean

    return CP_l, CP_l - R_l, q_l, b_l


def fit_liquid_phase_nasg(T, P, hl, vl, ref_state, P_inf_bounds=(1.0e7, 5.0e9)):
    """
    Close the NASG liquid-phase system with the sound-speed reference state (Eq. 68).
    Uses Brent's method to find P_inf,l such that:
        P0 + P_inf,l - (Cv_l / CP_l) * rho0 * c0^2 * (1 - b_l * rho0) = 0
    """
    rho0 = ref_state["rho0"]
    P0   = ref_state["P0"]
    c0   = ref_state["c0"]

    def residual(P_inf_l):
        CP_l, Cv_l, _, b_l = _fit_liquid_nasg_for_given_Pinf(T, P, hl, vl, P_inf_l)
        return P0 + P_inf_l - (Cv_l / CP_l) * rho0 * c0**2 * (1.0 - b_l * rho0)

    a, b = P_inf_bounds
    fa, fb = residual(a), residual(b)
    if fa * fb > 0:
        for _ in range(8):
            b *= 5.0
            fb = residual(b)
            if fa * fb < 0:
                break
        else:
            raise RuntimeError(
                f"NASG: could not bracket root for P_inf_l. "
                f"f({a:.2e})={fa:.3e}, f({b:.2e})={fb:.3e}"
            )

    P_inf_l = brentq(residual, a, b, xtol=1.0, rtol=1e-10)
    CP_l, Cv_l, q_l, b_l = _fit_liquid_nasg_for_given_Pinf(T, P, hl, vl, P_inf_l)

    return {
        "CP":    CP_l,
        "Cv":    Cv_l,
        "gamma": CP_l / Cv_l,
        "q":     q_l,
        "b":     b_l,
        "P_inf": P_inf_l,
    }


# =============================================================================
# SG liquid phase fit (b_l = 0 forced throughout)
# =============================================================================

def _fit_liquid_sg_for_given_Pinf(T, P, hl, vl, P_inf_l):
    """
    Given a candidate P_inf,l, solve for (CP_l, Cv_l, q_l) with b_l = 0.

    For SG:
      h_l(T) = CP_l * T + q_l              [b=0, no P term]
      v_l(T) = (CP_l - Cv_l)*T/(P+P_inf)  [b=0, one-parameter fit]
    """
    T_mean  = T.mean()
    hl_mean = hl.mean()
    vl_mean = vl.mean()

    # CP_l and q_l from h_l = CP_l*T + q_l (independent of P_inf_l)
    CP_l = np.sum(T * (hl - hl_mean)) / np.sum(T * (T - T_mean))
    q_l  = hl_mean - CP_l * T_mean

    # CP_l - Cv_l from v_l = R_l * T/(P+P_inf_l), one-parameter least squares
    TP  = T / (P + P_inf_l)
    R_l = np.sum(vl * TP) / np.sum(TP ** 2)
    Cv_l = CP_l - R_l

    return CP_l, Cv_l, q_l, 0.0


def fit_liquid_phase_sg(T, P, hl, vl, ref_state, P_inf_bounds=(1.0e7, 5.0e9)):
    """
    Close the SG liquid-phase system with the sound-speed reference state.
    With b_l = 0 the closure (Eq. 68) simplifies to:
        P0 + P_inf,l - (Cv_l / CP_l) * rho0 * c0^2 = 0
    Uses Brent's method to find P_inf,l.
    """
    rho0 = ref_state["rho0"]
    P0   = ref_state["P0"]
    c0   = ref_state["c0"]

    def residual(P_inf_l):
        CP_l, Cv_l, _, _ = _fit_liquid_sg_for_given_Pinf(T, P, hl, vl, P_inf_l)
        return P0 + P_inf_l - (Cv_l / CP_l) * rho0 * c0**2

    a, b = P_inf_bounds
    fa, fb = residual(a), residual(b)
    if fa * fb > 0:
        for _ in range(8):
            b *= 5.0
            fb = residual(b)
            if fa * fb < 0:
                break
        else:
            raise RuntimeError(
                f"SG: could not bracket root for P_inf_l. "
                f"f({a:.2e})={fa:.3e}, f({b:.2e})={fb:.3e}"
            )

    P_inf_l = brentq(residual, a, b, xtol=1.0, rtol=1e-10)
    CP_l, Cv_l, q_l, b_l = _fit_liquid_sg_for_given_Pinf(T, P, hl, vl, P_inf_l)

    return {
        "CP":    CP_l,
        "Cv":    Cv_l,
        "gamma": CP_l / Cv_l,
        "q":     q_l,
        "b":     0.0,
        "P_inf": P_inf_l,
    }


# =============================================================================
# Entropy constants (shared by NASG and SG)
# Eq. (41)/(69-72).  E = (b_l - b_g)/Rg → 0 automatically for SG.
# =============================================================================

def fit_entropy_constants(T, P, liquid, vapor):
    """
    Determine q'_l and q'_g via the saturated-vapor-pressure relation.
    Convention: q'_l = 0.
    For SG, b_l = b_g = 0 so E = 0 and the formula reduces to the SG form.
    """
    CP_l, Cv_l, P_inf_l = liquid["CP"], liquid["Cv"], liquid["P_inf"]
    CP_g, Cv_g          = vapor["CP"],  vapor["Cv"]
    b_l, b_g            = liquid["b"],  vapor["b"]
    q_l, q_g            = liquid["q"],  vapor["q"]

    Rg = CP_g - Cv_g
    B  = (q_l - q_g) / Rg
    C  = (CP_g - CP_l) / Rg
    D  = (CP_l - Cv_l) / Rg
    E  = (b_l - b_g) / Rg          # zero for SG

    A = np.mean(
        np.log(P)
        - (B + E * P) / T
        - C * np.log(T)
        - D * np.log(P + P_inf_l)
    )

    q_prime_l = 0.0
    q_prime_g = A * Rg - (CP_l - CP_g) + q_prime_l
    return q_prime_l, q_prime_g


# =============================================================================
# Top-level fitters
# =============================================================================

def compute_nasg_coefficients(sat_data, ref_state, T_range):
    """Full NASG fit: vapor (ideal-gas) + liquid (with covolume b_l)."""
    Tmin, Tmax = T_range
    mask = (sat_data[:, 0] >= Tmin) & (sat_data[:, 0] <= Tmax)
    d = sat_data[mask]
    if len(d) < 4:
        raise ValueError(f"Need >=4 points in [{Tmin},{Tmax}] K, got {len(d)}.")

    T, P, hl, hg, vl, vg = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]

    vapor  = fit_vapor_phase(T, P, hg, vg)
    liquid = fit_liquid_phase_nasg(T, P, hl, vl, ref_state)
    qp_l, qp_g = fit_entropy_constants(T, P, liquid, vapor)
    liquid["q_prime"] = qp_l
    vapor["q_prime"]  = qp_g
    return liquid, vapor


def compute_sg_coefficients(sat_data, ref_state, T_range):
    """Full SG fit: vapor (ideal-gas) + liquid (b_l = 0 forced)."""
    Tmin, Tmax = T_range
    mask = (sat_data[:, 0] >= Tmin) & (sat_data[:, 0] <= Tmax)
    d = sat_data[mask]
    if len(d) < 4:
        raise ValueError(f"Need >=4 points in [{Tmin},{Tmax}] K, got {len(d)}.")

    T, P, hl, hg, vl, vg = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]

    vapor  = fit_vapor_phase(T, P, hg, vg)
    liquid = fit_liquid_phase_sg(T, P, hl, vl, ref_state)
    qp_l, qp_g = fit_entropy_constants(T, P, liquid, vapor)
    liquid["q_prime"] = qp_l
    vapor["q_prime"]  = qp_g
    return liquid, vapor


# =============================================================================
# Analytical saturation curves (work for both NASG and SG; b=0 ⟹ SG formulae)
# =============================================================================

def eos_psat(T, liquid, vapor):
    """
    Solve Eq. (41) for the saturated vapor pressure at each T.
    For SG, E = 0 (b_l = b_g = 0) so the equation reduces to:
        ln(P + P_inf,g) = A + B/T + C*ln(T) + D*ln(P + P_inf,l)
    """
    CP_l, Cv_l, P_inf_l = liquid["CP"], liquid["Cv"], liquid["P_inf"]
    CP_g, Cv_g, P_inf_g = vapor["CP"],  vapor["Cv"],  vapor["P_inf"]
    b_l, b_g            = liquid["b"],  vapor["b"]
    q_l, q_g            = liquid["q"],  vapor["q"]
    qp_l, qp_g          = liquid["q_prime"], vapor["q_prime"]

    Rg = CP_g - Cv_g
    A  = (CP_l - CP_g + qp_g - qp_l) / Rg
    B  = (q_l - q_g) / Rg
    C  = (CP_g - CP_l) / Rg
    D  = (CP_l - Cv_l) / Rg
    E  = (b_l - b_g) / Rg

    def residual(P, Tval):
        return (np.log(P + P_inf_g) - A - (B + E * P) / Tval
                - C * np.log(Tval) - D * np.log(P + P_inf_l))

    T = np.atleast_1d(T)
    out = np.empty_like(T, dtype=float)
    for i, Ti in enumerate(T):
        try:
            out[i] = brentq(residual, 1.0, 1.0e9, args=(Ti,),
                            xtol=1e-3, rtol=1e-10)
        except ValueError:
            out[i] = np.nan
    return out


def eos_vl(T, P, liquid):
    """Eq. (46): saturated liquid specific volume. b=0 for SG."""
    return (liquid["CP"] - liquid["Cv"]) * T / (P + liquid["P_inf"]) + liquid["b"]


def eos_vg(T, P, vapor):
    """Eq. (45): saturated vapor specific volume. b=0 for SG."""
    return (vapor["CP"] - vapor["Cv"]) * T / (P + vapor["P_inf"]) + vapor["b"]


def eos_hl(T, P, liquid):
    """Eq. (44): saturated liquid specific enthalpy. b=0 for SG."""
    return liquid["CP"] * T + liquid["b"] * P + liquid["q"]


def eos_hg(T, P, vapor):
    """Eq. (43): saturated vapor specific enthalpy. b=0 for SG."""
    return vapor["CP"] * T + vapor["b"] * P + vapor["q"]


# =============================================================================
# Output: text table and figure
# =============================================================================

def print_coefficients(liquid, vapor, fluid_name, T_range, eos_name="NASG"):
    Tmin, Tmax = T_range
    header = f"{eos_name} coefficients for {fluid_name} in [{Tmin:.0f}-{Tmax:.0f}] K"
    print(header)
    print("=" * len(header))
    rows = [
        ("C_P    [J/(kg.K)]", liquid["CP"],      vapor["CP"]),
        ("C_v    [J/(kg.K)]", liquid["Cv"],      vapor["Cv"]),
        ("gamma  [-]",        liquid["gamma"],   vapor["gamma"]),
        ("P_inf  [Pa]",       liquid["P_inf"],   vapor["P_inf"]),
        ("b      [m^3/kg]",   liquid["b"],       vapor["b"]),
        ("q      [J/kg]",     liquid["q"],       vapor["q"]),
        ("q'     [J/(kg.K)]", liquid["q_prime"], vapor["q_prime"]),
    ]
    print(f"{'Coefficient':<20}{'Liquid':>18}{'Vapor':>18}")
    print("-" * 56)
    for label, lv, gv in rows:
        print(f"{label:<20}{lv:>18.4e}{gv:>18.4e}")
    print()


def plot_saturation(results, sat_data,
                    T_plot_range=None, n_curve=300,
                    savepath=None):
    """
    2x3 figure: p, L_v, h_l, h_g, v_l, v_g vs T.
    IAPWS-IF97 as solid black line; EOS curves as colored markers (no line).

    Parameters
    ----------
    results : list of (label, liquid_dict, vapor_dict) tuples
    """
    T_ref  = sat_data[:, 0]
    p_ref  = sat_data[:, 1]
    hl_ref = sat_data[:, 2]
    hg_ref = sat_data[:, 3]
    vl_ref = sat_data[:, 4]
    vg_ref = sat_data[:, 5]
    Lv_ref = hg_ref - hl_ref

    if T_plot_range is None:
        T_plot_range = (T_ref.min(), T_ref.max())

    T_curve   = np.linspace(*T_plot_range, n_curve)
    markers   = ["o", "s", "^", "D"]
    colors    = ["b", "r", "c", "m"]
    filled    = [True, False, True, False]   # NASG filled, SG hollow
    markevery = max(1, n_curve // 20)

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    panel_defs = [
        (axes[0, 0], r"$p\;[\mathrm{Pa}]$",                    True ),
        (axes[0, 1], r"$L_v\;[\mathrm{J\,kg^{-1}}]$",          False),
        (axes[0, 2], r"$h_l\;[\mathrm{J\,kg^{-1}}]$",          False),
        (axes[1, 0], r"$h_v\;[\mathrm{J\,kg^{-1}}]$",          False),
        (axes[1, 1], r"$v_l\;[\mathrm{m^3\,kg^{-1}}]$",        False),
        (axes[1, 2], r"$v_v\;[\mathrm{m^3\,kg^{-1}}]$",        True ),
    ]

    ref_curves = [p_ref, Lv_ref, hl_ref, hg_ref, vl_ref, vg_ref]

    for (ax, ylabel, log_y), y_ref in zip(panel_defs, ref_curves):
        ax.plot(T_ref, y_ref, "k-", lw=2.5, label="IAPWS-IF97", zorder=3)
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel(r"$T\;[\mathrm{K}]$")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.3, lw=0.7)

    for idx, (label, liquid, vapor) in enumerate(results):
        p_c  = eos_psat(T_curve, liquid, vapor)
        hl_c = eos_hl(T_curve, p_c, liquid)
        hg_c = eos_hg(T_curve, p_c, vapor)
        vl_c = eos_vl(T_curve, p_c, liquid)
        vg_c = eos_vg(T_curve, p_c, vapor)
        Lv_c = hg_c - hl_c

        eos_curves = [p_c, Lv_c, hl_c, hg_c, vl_c, vg_c]
        mk    = markers[idx % len(markers)]
        color = colors[idx % len(colors)]
        mfc   = color if filled[idx % len(filled)] else "none"
        for (ax, _, _), y_c in zip(panel_defs, eos_curves):
            ax.plot(T_curve, y_c, lw=0, marker=mk, ms=5,
                    markerfacecolor=mfc, markeredgecolor=color, markeredgewidth=1.5,
                    markevery=markevery, label=label, zorder=4)

    for ax, _, _ in panel_defs:
        ax.legend(markerscale=1.5)

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, bbox_inches="tight")
        print(f"Saved {savepath}")
    return fig


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit NASG and/or SG EOS coefficients to IAPWS saturation data."
    )
    parser.add_argument(
        "--mode", choices=["nasg", "sg", "both"], default="both",
        help="Which EOS(es) to fit and plot (default: both).",
    )
    args = parser.parse_args()

    T_FIT     = (300.0, 500.0)   # fitting window
    T_PLOT    = (275.0, 640.0)   # wider window for the plot
    REF_STATE = WATER_REFERENCE_STATE

    sat_fit  = iapws_saturation_table(T_FIT[0],  T_FIT[1],  n_points=101)
    sat_wide = iapws_saturation_table(T_PLOT[0], T_PLOT[1], n_points=300)

    plot_results = []

    if args.mode in ("nasg", "both"):
        liq_nasg, vap_nasg = compute_nasg_coefficients(sat_fit, REF_STATE, T_FIT)
        print_coefficients(liq_nasg, vap_nasg, "water/steam", T_FIT, eos_name="NASG")
        print("Reference (Le Metayer & Saurel 2016, Table V):")
        print("  Liquid: CP=4285  Cv=3610  gamma=1.19  P_inf=7.028e8 Pa  b=6.61e-4 m^3/kg")
        print("          q=-1.177788e6 J/kg   q'=0")
        print("  Vapor : CP=1401  Cv=955   gamma=1.47  P_inf=0  b=0")
        print("          q=2.077616e6 J/kg   q'=14317 J/(kg.K)\n")
        plot_results.append(("NASG", liq_nasg, vap_nasg))

    if args.mode in ("sg", "both"):
        liq_sg, vap_sg = compute_sg_coefficients(sat_fit, REF_STATE, T_FIT)
        print_coefficients(liq_sg, vap_sg, "water/steam", T_FIT, eos_name="SG")
        print("Reference (Flatten & Lund 2011 / CLAUDE.md):")
        print("  Liquid: CP=4267  Cv=1816  gamma=2.35  P_inf=1.0e9 Pa  b=0")
        print("          q=-1.167e6 J/kg   q'=0")
        print("  Vapor : CP=1487  Cv=1040  gamma=1.43  P_inf=0  b=0")
        print("          q=2.030e6 J/kg   q'=-23400 J/(kg.K)\n")
        plot_results.append(("SG", liq_sg, vap_sg))

    savepath = f"fit_eos_{args.mode}_saturation.pdf"
    plot_saturation(
        plot_results, sat_wide,
        T_plot_range=T_PLOT,
        savepath=savepath,
    )