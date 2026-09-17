#!/usr/bin/env python3
"""Cantera reference for the pseudo-1D laminar premixed flame of examples/amrcomp_flame.

usage: cantera_freeflame.py MECH.yaml FUEL PHI T_u P [out.csv] [--transport unity-Lewis-number|mixture-averaged]

Solves a freely propagating flame, prints the flame speed, adiabatic flame temperature, thermal thickness and
the `Burnt Y <species>` / `Burnt temperature` / `Inflow velocity` lines to paste into the NGA2 input, and
writes the profiles (x, T, u, rho, HRR, Y_k) to out.csv (default flame_<fuel>.csv).
"""
import sys
import numpy as np
import cantera as ct


def main(argv):
    mech, fuel, phi, Tu, P = argv[0], argv[1], float(argv[2]), float(argv[3]), float(argv[4])
    out = argv[5] if len(argv) > 5 and not argv[5].startswith('--') else f'flame_{fuel}.csv'
    transport = 'unity-Lewis-number'
    if '--transport' in argv:
        transport = argv[argv.index('--transport') + 1]
    gas = ct.Solution(mech, transport_model='mixture-averaged')   # YAML files without a default transport model
    gas.TP = Tu, P
    gas.set_equivalence_ratio(phi, fuel, 'O2:1.0, N2:3.76')
    Yu = gas.Y.copy(); rho_u = gas.density
    flame = ct.FreeFlame(gas, width=0.03)
    flame.transport_model = transport
    flame.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.1)
    flame.solve(loglevel=0, auto=True)
    x, T, u = flame.grid, flame.T, flame.velocity
    sL = u[0]
    Tb = T[-1]
    delta = (T[-1] - T[0]) / np.max(np.gradient(T, x))
    Yb = flame.Y[:, -1]
    gas.TPY = Tu, P, Yu; gas.equilibrate('HP'); Tad = gas.T
    print(f'# mechanism {mech}, {fuel}/air phi={phi}, T_u={Tu} K, P={P} Pa, transport={transport}')
    print(f'# S_L = {sL:.5f} m/s, T_b = {Tb:.2f} K (T_ad equilibrium {Tad:.2f} K), delta_L = {delta*1e3:.4f} mm, rho_u = {rho_u:.5f} kg/m3')
    print(f'Inflow velocity : {sL:.6f}')
    print(f'Burnt temperature : {Tb:.3f}')
    for k, name in enumerate(gas.species_names):
        if Yb[k] > 1e-8:
            print(f'Burnt Y {name} : {Yb[k]:.8e}')
    with open(out, 'w') as f:
        f.write('x T u rho hrr ' + ' '.join(gas.species_names) + '\n')
        for i in range(len(x)):
            f.write(f'{float(x[i])!r} {float(T[i])!r} {float(u[i])!r} {float(flame.density[i])!r} {float(flame.heat_release_rate[i])!r} '
                    + ' '.join(repr(float(v)) for v in flame.Y[:, i]) + '\n')
    print(f'# wrote {out} ({len(x)} points)')


if __name__ == '__main__':
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    main(sys.argv[1:])
