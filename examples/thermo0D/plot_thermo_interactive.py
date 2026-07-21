#!/usr/bin/env python3
"""Interactive thermo plotter for thermo0D CSV output.

Loads a single CSV (default: data/thermo.csv) and opens a matplotlib window
with:
  * radio buttons to choose the x-axis: time (t) or mixture specific energy (e)
  * checkboxes to pick which properties to plot on the y-axis

Usage:
  python3 plot_thermo_interactive.py [path/to/thermo.csv]

Requires: numpy, matplotlib
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import CheckButtons, RadioButtons


# Preferred display order; any extra CSV columns are appended.
PREFERRED = [
    "alpha",
    "eL", "eG", "e",
    "PL", "PG", "P",
    "TL", "TG", "T",
    "rhoL", "rhoG", "rho",
    "vL", "vG",
    "hL", "hG", "sL", "sG", "gL", "gG",
    "cL", "cG", "cvL", "cvG",
    "massL", "massG", "mass", "V", "Q_added", "b_rho",
    "flag_Tcrit", "flag_Pcrit", "flag_critical", "flag_packing",
]

FLAG_COLS = {"flag_Tcrit", "flag_Pcrit", "flag_critical", "flag_packing"}

# x-axis choices: (CSV column, axis label)
X_CHOICES = [
    ("t", "t [s]"),
    ("e", "e [J/kg]  (mixture specific energy)"),
]


def axis_limits(vals, rel_floor: float = 1e-4, rel_pad: float = 0.05):
    lo = float(np.nanmin(vals))
    hi = float(np.nanmax(vals))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return 0.0, 1.0
    span = hi - lo
    mid = 0.5 * (lo + hi)
    scale = max(abs(mid), abs(lo), abs(hi), 1.0)
    if span <= rel_floor * scale:
        half = max(0.01 * scale, 0.5 * rel_floor * scale, 1.0)
        return mid - half, mid + half
    pad = rel_pad * span
    return lo - pad, hi + pad


def load_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as f:
        rows = []
        header = None
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if header is None:
                header = next(csv.reader([s]))
                continue
            rows.append(next(csv.reader([s])))
    if header is None or not rows:
        raise SystemExit(f"No data rows in {path}")
    data = {h: [] for h in header}
    for row in rows:
        for h, v in zip(header, row):
            data[h].append(float(v))
    return {h: np.asarray(v, dtype=float) for h, v in data.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "csv",
        nargs="?",
        default="data/thermo.csv",
        help="Path to thermo CSV (default: data/thermo.csv)",
    )
    args = ap.parse_args()
    path = Path(args.csv)
    if not path.is_file():
        raise SystemExit(f"File not found: {path}")

    data = load_csv(path)
    for col, _ in X_CHOICES:
        if col not in data:
            raise SystemExit(f"CSV must contain a '{col}' column")

    x_keys = [c for c, _ in X_CHOICES]
    x_labels = [lab for _, lab in X_CHOICES]
    x_key = "t"

    props = [c for c in PREFERRED if c in data and c not in x_keys]
    props += [c for c in data if c not in props and c not in x_keys]
    if not props:
        raise SystemExit("No property columns to plot")

    # Default: show P and T
    active = {c: c in ("P", "T") for c in props}

    fig, ax = plt.subplots(figsize=(10, 6.0))
    plt.subplots_adjust(left=0.30, right=0.98, top=0.92, bottom=0.22)
    ax.set_ylabel("value")
    ax.set_title(f"thermo0D — {path.name}")
    ax.grid(True, alpha=0.3)

    lines: dict[str, object] = {}
    for col in props:
        (ln,) = ax.plot(
            data[x_key],
            data[col],
            label=col,
            lw=1.8 if col not in FLAG_COLS else 1.2,
            ls="-" if col not in FLAG_COLS else "--",
            visible=active[col],
        )
        lines[col] = ln

    # Critical-crossing markers (x follows the selected abscissa)
    vlines = []
    for flag, color, name in (
        ("flag_Tcrit", "C3", "Tcrit"),
        ("flag_Pcrit", "C1", "Pcrit"),
    ):
        if flag not in data:
            continue
        idx = np.where(data[flag] >= 0.5)[0]
        if not len(idx):
            continue
        i0 = int(idx[0])
        vl = ax.axvline(
            data[x_key][i0],
            color=color,
            ls=":",
            lw=1.0,
            alpha=0.8,
            label=f"{name} @ {x_key}={data[x_key][i0]:.3g}",
        )
        vlines.append((vl, i0, name, color))

    def refresh_legend():
        handles = [lines[c] for c in props if lines[c].get_visible()]
        labels = [c for c in props if lines[c].get_visible()]
        for vl, i0, name, _ in vlines:
            handles.append(vl)
            labels.append(f"{name} @ {x_key}={data[x_key][i0]:.3g}")
        if handles:
            ax.legend(handles, labels, loc="best", fontsize=9)
        elif ax.get_legend() is not None:
            ax.get_legend().remove()

    def rescale():
        x = data[x_key]
        ax.set_xlim(*axis_limits(x))
        ys = [data[c] for c in props if lines[c].get_visible()]
        if ys:
            ax.set_ylim(*axis_limits(np.concatenate(ys)))

    def set_x(key: str):
        nonlocal x_key
        x_key = key
        x = data[x_key]
        xlab = next(lab for k, lab in X_CHOICES if k == x_key)
        ax.set_xlabel(xlab)
        for col in props:
            lines[col].set_xdata(x)
        for vl, i0, name, _ in vlines:
            vl.set_xdata([x[i0], x[i0]])
            vl.set_label(f"{name} @ {x_key}={x[i0]:.3g}")
        refresh_legend()
        rescale()
        fig.canvas.draw_idle()

    set_x("t")

    # Y-property checkboxes
    yax = fig.add_axes([0.02, 0.22, 0.24, 0.68])
    yax.set_title("Y properties", fontsize=10)
    checks = CheckButtons(yax, props, [active[c] for c in props])

    def on_y_clicked(label: str):
        ln = lines[label]
        ln.set_visible(not ln.get_visible())
        refresh_legend()
        rescale()
        fig.canvas.draw_idle()

    checks.on_clicked(on_y_clicked)

    # X-axis radio buttons
    xax = fig.add_axes([0.30, 0.02, 0.55, 0.12])
    xax.set_title("X axis", fontsize=10)
    radios = RadioButtons(xax, x_labels, active=0)

    def on_x_clicked(label: str):
        key = next(k for k, lab in X_CHOICES if lab == label)
        set_x(key)

    radios.on_clicked(on_x_clicked)

    plt.show()


if __name__ == "__main__":
    main()
