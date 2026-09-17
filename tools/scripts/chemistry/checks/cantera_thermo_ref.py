#!/usr/bin/env python3
"""Reference states and properties for examples/chem_thermo_check (Cantera 3.x).

usage:
  cantera_thermo_ref.py gen     MECH.yaml states.dat
  cantera_thermo_ref.py compare MECH.yaml states.dat nga_props.dat [nga_transport.dat]

`gen` writes a table of (T, P, Y_1..Y_ns) states in the mechanism's species order.
`compare` reads the properties written by the Fortran program for those states and reports the maximum
relative deviation from Cantera for each quantity (thermo: W, rho, e, h, cp, cv, s, c, gamma; transport,
if the transport file is given: mixture-averaged mu and lambda against Cantera's mixture-averaged model,
and the unity-Lewis species diffusivity lambda/(rho cp)).
"""
import sys
import numpy as np
import cantera as ct

TEMPERATURES = [300.0, 600.0, 999.0, 1001.0, 1500.0, 2500.0]
THERMO_TOL = 1.0e-9      # same NASA7 polynomials; only R (8.314462618 vs 8.31446261815324) and round-off differ


def compositions(gas):
    """Named compositions as mass-fraction arrays in the mechanism's species order."""
    comps = {}
    for sp in ('H2', 'N2'):
        if sp in gas.species_names:
            y = np.zeros(gas.n_species); y[gas.species_index(sp)] = 1.0
            comps['pure_' + sp] = y
    fuel = 'CH4' if 'CH4' in gas.species_names else 'H2'
    gas.TP = 300.0, ct.one_atm
    gas.set_equivalence_ratio(1.0, fuel, 'O2:1.0, N2:3.76')
    comps['stoich_' + fuel + '_air'] = gas.Y.copy()
    gas.equilibrate('HP')
    comps['equilibrium_' + fuel + '_air'] = gas.Y.copy()
    return comps


def gen(mech, states_file):
    gas = ct.Solution(mech)
    rows = []
    for name, y in compositions(gas).items():
        pressures = [ct.one_atm, 5.0e5] if name.startswith('stoich') else [ct.one_atm]
        for P in pressures:
            for T in TEMPERATURES:
                rows.append((T, P, y))
    with open(states_file, 'w') as f:
        f.write(f'{len(rows)}\n')
        for T, P, y in rows:
            f.write(f'{T!r} {P!r} ' + ' '.join(repr(float(v)) for v in y) + '\n')
    print(f'wrote {len(rows)} states for {gas.n_species} species to {states_file}')


def read_table(path):
    header = None
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if header is None:
                    header = line[1:].split()
                elif line.startswith('# species:'):
                    header_species = line.split(':', 1)[1].split()
                    header = (header, header_species)
                continue
            if line.strip():
                data.append([float(v) for v in line.split()])
    if isinstance(header, tuple):
        cols, species = header
    else:
        cols, species = header, None
    return cols, species, np.array(data)


def compare(mech, states_file, props_file, transport_file=None):
    gas = ct.Solution(mech)
    with open(states_file) as f:
        n = int(f.readline())
        states = [np.array([float(v) for v in f.readline().split()]) for _ in range(n)]
    cols, species, nga = read_table(props_file)
    if species is not None and species != gas.species_names:
        sys.exit('species order in the Fortran output differs from the mechanism order')
    col = {c: i for i, c in enumerate(cols)}
    if nga.shape[0] != n:
        sys.exit(f'{props_file} has {nga.shape[0]} rows, expected {n}')
    names = ['W_mix', 'rho', 'e', 'h', 'cp', 'cv', 's', 'c', 'gamma']
    worst = {k: 0.0 for k in names}
    for row, st in zip(nga, states):
        T, P, y = st[0], st[1], st[2:]
        gas.TPY = T, P, y
        ref = {'W_mix': gas.mean_molecular_weight * 1.0e-3, 'rho': gas.density, 'e': gas.int_energy_mass,
               'h': gas.enthalpy_mass, 'cp': gas.cp_mass, 'cv': gas.cv_mass, 's': gas.entropy_mass,
               'c': gas.sound_speed, 'gamma': gas.cp_mass / gas.cv_mass}
        for k in names:
            worst[k] = max(worst[k], abs(row[col[k]] - ref[k]) / max(abs(ref[k]), 1.0e-300))
    ok = True
    print('thermo: max relative deviation from Cantera (tolerance %.1e)' % THERMO_TOL)
    for k in names:
        flag = 'ok' if worst[k] <= THERMO_TOL else 'FAIL'
        ok &= worst[k] <= THERMO_TOL
        print(f'  {k:6s} {worst[k]:10.3e} {flag}')
    # Solver diagnostics written by the Fortran program
    T = nga[:, col['T']]
    for k, lab in (('T_rt_cold', 'cold'), ('T_rt_warm', 'warm')):
        err = np.max(np.abs(nga[:, col[k]] - T) / T)
        its = nga[:, col['it_cold' if lab == 'cold' else 'it_warm']]
        flag = 'ok' if err <= 1.0e-9 else 'FAIL'; ok &= err <= 1.0e-9
        print(f'T(rho,e) round trip, {lab} start: max rel error {err:.3e}, iterations {int(its.min())}-{int(its.max())} {flag}')
    # accessor vs separate calls: both converge to tol_T=1e-10 relative in T, so agree to ~1e-9 in p and c
    for k, tol in (('dp_acc', 1.0e-9), ('dc_acc', 1.0e-9), ('dhk_max', 1.0e-12)):
        v = np.max(nga[:, col[k]])
        flag = 'ok' if v <= tol else 'FAIL'; ok &= v <= tol
        print(f'{k}: max {v:.3e} (tol {tol:.0e}) {flag}')
    if transport_file is not None:
        cols_t, _, tr = read_table(transport_file)
        ct_ = {c: i for i, c in enumerate(cols_t)}
        gas.transport_model = 'mixture-averaged'
        w_mu = w_lam = w_D = 0.0
        w_mu_polar = 0.0
        for row, st in zip(tr, states):
            T, P, y = st[0], st[1], st[2:]
            gas.TPY = T, P, y
            dmu = abs(row[ct_['mu_mix']] - gas.viscosity) / gas.viscosity
            dlam = abs(row[ct_['lambda_mix']] - gas.thermal_conductivity) / gas.thermal_conductivity
            D_unity = gas.thermal_conductivity / (gas.density * gas.cp_mass)
            dD = abs(row[ct_['D_unityLe']] - D_unity) / D_unity
            polar = y[gas.species_index('H2O')] > 0.01 if 'H2O' in gas.species_names else False
            if polar:
                w_mu_polar = max(w_mu_polar, dmu)
            else:
                w_mu = max(w_mu, dmu)
            w_lam = max(w_lam, dlam); w_D = max(w_D, dD)
        print('transport (mixavg mode vs Cantera mixture-averaged):')
        print(f'  mu     non-polar mixtures {w_mu:9.3e} (tol 2e-2) {"ok" if w_mu <= 2e-2 else "FAIL"}')
        print(f'  mu     H2O-rich mixtures  {w_mu_polar:9.3e} (informational: no Stockmayer correction)')
        print(f'  lambda {w_lam:9.3e} (tol 1.5e-1: modified Eucken vs Cantera species conductivities) {"ok" if w_lam <= 1.5e-1 else "FAIL"}')
        print(f'  D=lambda/(rho cp) {w_D:9.3e} (inherits the lambda tolerance) {"ok" if w_D <= 1.5e-1 else "FAIL"}')
        ok &= (w_mu <= 2e-2) and (w_lam <= 1.5e-1)
    print('RESULT:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == 'gen':
        gen(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 5 and sys.argv[1] == 'compare':
        sys.exit(compare(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else None))
    else:
        sys.exit(__doc__)
