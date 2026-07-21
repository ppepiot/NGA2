# 0D constant-volume heat addition — `thermo0D`

Closed reactor with fixed total volume. Supports pure liquid water (NASG),
pure nitrogen (ideal gas), or a mixture (default: half / half by volume).
Specific internal energy of each present phase rises at a prescribed rate.
Mixture cases can apply mechanical (`p`) or mechanical+thermal (`pT`)
relaxation after each heat step.

## Build and run

```bash
cd examples/thermo0D
make
mpirun -np 1 ./thermo0D.dp.gnu.opt.mpi.exe -i input
python3 plot_thermo_interactive.py data/thermo.csv
```

## Physics

Fixed total volume \(V\). Liquid volume fraction \(\alpha\):
\[
V_L=\alpha V,\quad V_G=(1-\alpha)V,\qquad
\frac{de_k}{dt}=\dot{q}\quad(k=L,G\text{ if present}).
\]
For \(0<\alpha<1\), optional `relax_ig_nasg` equalization of \(P\) (and \(T\)).

Default NASG water: Le Métayer & Saurel (2016) Table V.  
Default N₂: \(\gamma=1.4\), \(c_v=742\) J/(kg·K) (\(R=296.8\)).

## Critical flags

Applied to the **liquid** when present (else gas). Thresholds are IAPWS water
by default:

| Column | Meaning |
|--------|---------|
| `flag_Tcrit` | latched 1 once \(T \ge T_\mathrm{crit}\) |
| `flag_Pcrit` | latched 1 once \(P \ge P_\mathrm{crit}\) |
| `flag_critical` | 1 when both have been crossed |
| `flag_packing` | 1 if liquid \(b\rho\) reaches the packing warn fraction |

## Input highlights

| Key | Unit | Default |
|-----|------|---------|
| `Liquid volume fraction` | — | 0.5 |
| `Volume` | m³ | 1 |
| `Pressure` / `Temperature` | Pa / K | 101325 / 298.15 |
| `Gas Gamma` / `Gas Cv` | — / J/(kg·K) | 1.4 / 742 |
| `Relaxation type` | none / p / pT | pT |
| `Heating rate` | J/(kg·s) | 1e4 |

Optional per-phase ICs: `Liquid pressure`, `Liquid temperature`, `Gas pressure`, `Gas temperature`.

## CSV columns

Phasic and mixture fields: `t, alpha, eL, eG, e, PL, PG, P, TL, TG, T, …`,
plus the four flag columns.

## Interactive plot

`plot_thermo_interactive.py` — choose the x-axis (`t` or mixture specific energy `e`)
with the radio buttons, and click property checkboxes for the y-axis.
