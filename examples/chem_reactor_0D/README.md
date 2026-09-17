# chem_reactor_0D — 0D Homogeneous Reactor

0D adiabatic constant-pressure chemistry solver. Integrates mass fraction and temperature evolution using CVODE (SUNDIALS) with adaptive stepping. Chemical source terms come from the `fcmech` kinetics module.

## Naming convention (for future variants)

Internal prefixes for homogeneous reactor variants:

| Prefix | Variant | Description |
|--------|---------|-------------|
| `hr_ib_` | Isobaric | Constant pressure (`Reactor : isobar`, default) |
| `hr_ic_` | Isochoric | Constant volume (`Reactor : isochor`) |

## Features

- **Adiabatic, constant pressure**: ideal gas, NASA7 thermo
- **Initial composition**: Fuel + equivalence ratio (standard air) or per-species mass fractions
- **Kinetics**: Arrhenius, three-body, falloff (Troe), plog
- **Output**: Time series of species mass fractions, T, P, ρ

## Requirements

- **Fortran compiler**: gfortran (or Intel, etc.)
- **SUNDIALS/CVODE**: For stiff ODE integration
  - macOS: `brew install sundials`
  - Linux: install via package manager or build from source

## Directory Structure

```
chem_reactor_0D/
├── GNUmakefile
├── README.md
├── input              # Input parameters
├── plot_results.py    # Plot results
└── src/
    ├── Make.package
    ├── simulation.f90      # Main program
    ├── chem_data_fc.f90      # Kinetics (generated)
    └── chem_data_reactions.txt # Reaction index (generated)
```

## Quick Start

### 1. Generate the kinetics module

`chem_data_fc.f90` must be generated from a Cantera YAML mechanism:

```bash
cd examples/chem_reactor_0D

# Generate from a mechanism (e.g. gri30)
python ../../tools/scripts/chemistry/yaml2nga.py \
    ../../tools/scripts/chemistry/kinetics/gri30.yaml \
    src/chem_data_fc.f90
```

This creates `src/chem_data_fc.f90` and `src/chem_data_reactions.txt`. The generator
rejects anything it cannot convert faithfully (unknown units, unsupported reaction
types, unbalanced reactions, ...) instead of guessing; fix the YAML if it complains.
Reverse rates of reversible reactions are computed at run time from detailed balance
(k_r = k_f/K_c), so equilibrium is exactly consistent with the NASA7 thermodynamics.
`--fit-reverse [--fit-range TMIN TMAX]` instead fits 1/K_c of every reversible reaction to
an Arrhenius form over the given range (default 300–3000 K) and evaluates reverse rates
exactly like forward ones; the generator prints the fit-error statistics and writes the
per-reaction error to `chem_data_reactions.txt` (gri30: worst 15% over 300–3000 K, 4%
over 600–2800 K; equilibrium is then only approximately reproduced). Since the generated
module caches all temperature-only rate data, the fit saves only the NASA7 evaluation
(a few percent); it is kept for cases that require an Arrhenius reverse form.

### Performance notes

The generated module evaluates all rate-coefficient exponentials with a vectorizable
`exp` and caches everything that depends on temperature only (Arrhenius and three-body
rates, falloff limits and Troe centering, reverse-rate factors, NASA7 polynomials), keyed
on the exact T (and P for PLOG). Calls at a repeated temperature — Newton iterations and
the composition columns of a finite-difference Jacobian — cost about a third of a full
evaluation (gri30: ~0.8 µs vs ~1.7 µs per `fcmech_get_ydot`, 2.7 µs before). The reactor
supplies its own finite-difference Jacobian (`hr_fd_jacobian`) that exploits this and
uses a central difference for the temperature column, and prints CVODE statistics and
the integration-loop time at the end of a run.

### 2. Build

```bash
make CVODE_DIR=$(brew --prefix sundials)
```

If SUNDIALS is elsewhere, set `CVODE_DIR` to the install root (with `include/`, `lib/`, `fortran/`).

### 3. Run

```bash
./chem_reactor_0D.dp.gnu.opt.exe -i input
```

With MPI (single rank):

```bash
mpirun -np 1 ./chem_reactor_0D.dp.gnu.opt.mpi.exe -i input
```

### 4. Plot results

```bash
python plot_results.py results_hr.out
```

## Input File

The `input` file uses NGA2 param format (`key : value`). Lines starting with `#` are comments.

### Required parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `Temperature` | Initial temperature (K) | `1000.0` |
| `Pressure` | Pressure (Pa) | `100000.0` |
| `End time` | Integration end time (s) | `0.1` |
| `Time step` | Output interval (s) | `1e-6` |
| `Output file` | Output filename | `results_hr.out` |

### Reactor type

| Parameter | Description | Example |
|-----------|-------------|---------|
| `Reactor` | `isobar` (constant P) or `isochor` (constant V) | `isobar` |

Default: `isobar`. Case independent.

### Initial composition (choose one)

**Option A — Fuel + equivalence ratio**

```text
Fuel : NH3
Equivalence ratio : 1.0
```

Uses standard air (N2=0.79, O2=0.21 mole frac). Works for fuels without carbon (e.g. NH3) or with carbon (e.g. CH4).

**Option B — Per-species mass fractions**

```text
Initial Y CH4 : 0.05
Initial Y O2 : 0.2
Initial Y N2 : 0.75
```

`Initial Y <species>` for each non-zero species; the values are normalized to sum to 1.

### Command-line input

```bash
./chem_reactor_0D.dp.gnu.opt.exe -i input
```

## Output Format

`results_hr.out` is a space-separated file:

- **Header**: `time`, `species1-Y1`, `species2-Y2`, ..., `T`, `P`, `rho`
- **Columns**: Fixed-width columns (18 chars for time, species, T, P, rho)

## Plotting

`plot_results.py` reads the output and plots temperature and selected species:

```bash
python plot_results.py [results_hr.out] [-o output_prefix]
```

Requires: `numpy`, `matplotlib`

## Available Mechanisms

YAML mechanisms are in `tools/scripts/chemistry/kinetics/`:

- `gri30.yaml` — GRI-Mech 3.0 (CH4/hydrocarbon), detailed and thermodynamically closed. What
  `src/chem_data_fc.f90` ships with, and what any case where **auto-ignition drives the dynamics**
  should use (this reactor, and `examples/amrcomp_flame` `input_ignition`).
- `FM/CH4.Igni73.yaml` — skeletal CH4 (28 species, 97 reactions), converted from FlameMaster with
  `FM2yaml.py`. Reduced for **flame propagation**: laminar flame speed within 8% of gri30, ~5x cheaper
  per rate evaluation and ~3.6x less stiff, so it is the mechanism for the CH4 flame and jet cases.
  Do **not** use it where ignition timing matters: its constant-volume ignition delay is 1.3x gri30 at
  1200 K, 4.3x at 1800 K and 5.1x at 2000 K. Also note that 95 of its 97 reactions are irreversible
  (FlameMaster supplies explicit forward/backward Arrhenius pairs), so it does not relax to
  thermodynamic equilibrium and is not meant for long-time post-flame states.
- `h2o2.yaml` — H2/O2
- `FM/*.yaml` — other FlameMaster mechanisms converted with `FM2yaml.py`

Any Cantera-format YAML mechanism can be used.

## Build Options

| Variable | Default | Description |
|----------|---------|-------------|
| `PRECISION` | `DOUBLE` | `DOUBLE` or `SINGLE` |
| `USE_MPI` | `TRUE` | MPI support |
| `COMP` | `gnu` | Compiler: `gnu`, `intel`, etc. |
| `DEBUG` | `FALSE` | Debug build |

Example:

```bash
make PRECISION=DOUBLE COMP=intel CVODE_DIR=/path/to/sundials
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `CVODE_DIR not set` | Install SUNDIALS and pass `CVODE_DIR=$(brew --prefix sundials)` |
| `Fuel, O2, or N2 not found` | Mechanism must define O2 and N2; fuel name must match species name |
| `No Initial Y specified` | Provide either Fuel+Equivalence ratio or Initial Y per species |
| `chem_data_fc.f90` missing | Run `yaml2nga.py` to generate from a YAML mechanism |

## See Also

- `tools/scripts/chemistry/yaml2nga.py` — YAML → chem_data_fc converter
- `tools/scripts/chemistry/check_diffusion.py`, `cantera_thermo_ref.py`, `cantera_freeflame.py` — Cantera comparison utilities
- `examples/amrcomp_flame` — the same `fcmech` module coupled to the AMR compressible solver (explicit chemistry)
