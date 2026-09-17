#!/usr/bin/env python3
"""Post-process a slot-burner run the way Selle, Poinsot & Ferret (Combust. Flame 158 (2011) 146-154) do.

usage: compare_slotburner.py <run_dir> <input_file> [mechanism.yaml]

Reads monitor/flame and monitor/simulation and reports the flame length L_f, the mean consumption speed
s_L = U h / L_f (their Eq. 1) and the ratio eta_1 = s_L^0 / s_L_bar (their Eq. 2), which the paper finds to
be within 1.5% of unity for methane. With a mechanism file the unstretched reference s_L^0 is recomputed
with Cantera at the run's operating point.
"""
import sys, os, re
import numpy as np


def read_monitor(path):
    head, rows = [], []
    with open(path) as f:
        for line in f:
            try:
                rows.append([float(v) for v in line.split()])
            except ValueError:
                head.append(line.rstrip('\n'))
    return head, np.array(rows)


def read_input(path):
    """Operating point from the input file the run actually used."""
    par = {}
    for line in open(path):
        line = line.split('#')[0]
        if ':' in line:
            k, v = line.split(':', 1)
            par[k.strip()] = v.strip()
    return par


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__.strip().split('\n')[2])
    run, inp = sys.argv[1], sys.argv[2]
    mech = sys.argv[3] if len(sys.argv) > 3 else None
    _, fl = read_monitor(os.path.join(run, 'monitor', 'flame'))
    if fl.shape[0] < 3:
        sys.exit('not enough rows in monitor/flame')
    par = read_input(inp if os.path.exists(inp) else os.path.join(run, inp))
    U = float(par.get('Bulk velocity', 'nan'))
    h = float(par.get('Slot width', 'nan'))
    fuel = par.get('Fuel', '?')
    phi = float(par.get('Equivalence ratio', 'nan'))

    t, Hf, Lf, sL, hrr, ratio, Tmax = (fl[:, i] for i in range(1, 8))
    # Average over the last third, which is the closest thing to a steady state
    m = t >= t[0] + 2.0 * (t[-1] - t[0]) / 3.0
    print(f'run: {run}')
    print(f'  {fuel}/air, phi = {phi}, U = {U} m/s, h = {h*1e3:.1f} mm, t = {t[0]:.3e} .. {t[-1]:.3e} s '
          f'({len(t)} samples)')
    print(f'  flame tip height H_f : {Hf[m].mean()*1e3:8.3f} mm   (last value {Hf[-1]*1e3:.3f})')
    print(f'  flame length   L_f   : {Lf[m].mean()*1e3:8.3f} mm   TRIANG estimate, 2*sqrt((h/2)^2+H_f^2)')
    print(f'  consumption speed    : {sL[m].mean():8.4f} m/s  = U h / L_f  (paper Eq. 1)')
    print(f'  max temperature      : {Tmax[m].mean():8.2f} K')
    print(f'  heat release / fed   : {ratio[m].mean():8.4f}     (1 at steady state; a smaller value means '
          f'the flame is still building)')
    drift = np.polyfit(t[m], Hf[m], 1)[0] if m.sum() > 2 else float('nan')
    print(f'  tip drift            : {drift:8.3f} m/s   (should approach 0)')

    if mech:
        import cantera as ct
        g = ct.Solution(mech, transport_model='mixture-averaged')
        g.set_equivalence_ratio(phi, fuel, 'O2:0.21,N2:0.79')
        g.TP = float(par.get('Temperature', 300.0)), float(par.get('Pressure', 101325.0))
        f = ct.FreeFlame(g, width=0.05)
        f.transport_model = 'mixture-averaged'
        f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
        f.solve(loglevel=0, auto=True)
        sL0 = f.velocity[0]
        dL = (f.T[-1] - f.T[0]) / np.max(np.abs(np.gradient(f.T, f.grid)))
        print(f'  Cantera reference    : s_L^0 = {sL0:.4f} m/s, delta_L = {dL*1e3:.4f} mm, T_b = {f.T[-1]:.1f} K')
        print(f'  eta_1 = s_L^0/s_L    : {sL0/sL[m].mean():8.4f}     (paper Table 3: 1.01 for methane)')
        print(f'  T_max / T_b          : {Tmax[m].mean()/f.T[-1]:8.4f}')


if __name__ == '__main__':
    main()
