#!/usr/bin/env python3
"""Interactive viewer for relax1D data/profiles_NNNNNN.csv dumps.

Loads every profiles_*.csv in a directory (default: ./data), reads the
'# time = ...' header comment, and opens a matplotlib window with:

  * radio buttons to pick the x-axis column
  * checkboxes to pick one or more y-axis columns
  * an animated line plot that cycles through the dumps in time order
  * play/pause button and a time slider for manual scrubbing

Usage:
  python3 plot_profiles_interactive.py [data_dir]

Example:
  cd examples/relax1D/run_1mm
  python3 ../plot_profiles_interactive.py data
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import CheckButtons, RadioButtons, Button, Slider
import numpy as np


TIME_RE = re.compile(r"#\s*time\s*=\s*([-+0-9.eE]+)")

# Pure-phase thresholds (match amrvof VFlo ~ 1e-12 in spirit; slightly looser for CSV)
_VF_LO = 1.0e-12
_VF_HI = 1.0 - _VF_LO


def axis_limits(vals, rel_floor: float = 1e-4, rel_pad: float = 0.05):
    """Stable (lo, hi) for a data array, without zooming into roundoff.

    Global min/max over a run is what we want when a field really varies. But
    for a field that is physically constant (e.g. composite P stuck at 1e5 to
    machine noise), nanmin/nanmax differ by ~1e-6 and autoscale turns that into
    a wild-looking plot. If the relative span is below ``rel_floor``, expand
    to a readable window around the mean instead.
    """
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


def join_phasic_fields(frame: dict) -> None:
    """Add joined P/T and mask inactive-phase dumps, in place.

    The solver stores phasic fields in every cell, but the inactive phase is
    undefined there and the CSV dump writes zeros (PL=0 in pure gas, PG=0 in
    pure liquid, etc.). Plotting those raw columns therefore jumps between the
    physical value and 0. Mask inactive cells with NaN so line/colormap plots
    only show the phase that is present.

    Also builds single-field ``P`` and ``T``: liquid value where VF=1, gas
    where VF=0, arithmetic average in mixed cells.
    """
    required = ("VF", "PL", "PG", "TL", "TG")
    if not all(k in frame for k in required):
        return
    vf = frame["VF"]
    pure_liq = vf >= _VF_HI
    pure_gas = vf <= _VF_LO
    mixed = ~(pure_liq | pure_gas)
    has_liq = vf >= _VF_LO  # liquid present (pure or mixed)
    has_gas = vf <= _VF_HI  # gas present (pure or mixed)

    # Mask inactive-phase raw columns (copy so we do not alias CSV arrays)
    for col, active in (("PL", has_liq), ("TL", has_liq), ("PG", has_gas), ("TG", has_gas)):
        vals = np.array(frame[col], dtype=float, copy=True)
        vals[~active] = np.nan
        frame[col] = vals
    if "Yv" in frame:
        yv = np.array(frame["Yv"], dtype=float, copy=True)
        yv[~has_gas] = np.nan
        frame["Yv"] = yv

    P = np.full_like(vf, np.nan, dtype=float)
    T = np.full_like(vf, np.nan, dtype=float)
    P[pure_liq] = frame["PL"][pure_liq]
    T[pure_liq] = frame["TL"][pure_liq]
    P[pure_gas] = frame["PG"][pure_gas]
    T[pure_gas] = frame["TG"][pure_gas]
    P[mixed] = 0.5 * (frame["PL"][mixed] + frame["PG"][mixed])
    T[mixed] = 0.5 * (frame["TL"][mixed] + frame["TG"][mixed])
    frame["P"] = P
    frame["T"] = T


def load_profiles(data_dir: Path):
    """Return (times, frames) where each frame is a dict of column -> ndarray.

    Besides the CSV columns, each frame also carries derived fields
    ``P`` and ``T`` (joined liquid/gas pressure and temperature).
    """
    files = sorted(data_dir.glob("profiles_*.csv"))
    if not files:
        raise SystemExit(f"No profiles_*.csv found in {data_dir}")

    times = []
    frames = []
    columns = None
    for path in files:
        text = path.read_text().splitlines()
        t = None
        header = None
        rows = []
        for line in text:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                m = TIME_RE.search(line)
                if m:
                    t = float(m.group(1))
                continue
            if header is None:
                header = [c.strip() for c in line.split(",")]
                continue
            rows.append([float(v) for v in line.split(",")])
        if header is None or not rows:
            continue
        if t is None:
            # Fall back to file index if the time comment is missing
            t = float(len(times))
        if columns is None:
            columns = list(header)
        elif list(header) != columns:
            raise SystemExit(f"Column mismatch in {path.name}: {header} vs {columns}")
        arr = np.asarray(rows, dtype=float)
        frame = {name: arr[:, i] for i, name in enumerate(header)}
        join_phasic_fields(frame)
        frames.append(frame)
        times.append(t)

    # Drop exact duplicate final dump if present (same time as previous)
    keep = [0]
    for i in range(1, len(times)):
        if abs(times[i] - times[keep[-1]]) > 1e-18:
            keep.append(i)
    times = [times[i] for i in keep]
    frames = [frames[i] for i in keep]

    # Expose derived fields in the UI column list (after the raw CSV ones)
    if frames and "P" in frames[0]:
        for name in ("P", "T"):
            if name not in columns:
                columns.append(name)
    return np.asarray(times), frames, columns


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="data",
        help="Directory containing profiles_*.csv (default: data)",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=200.0,
        help="Animation frame interval in milliseconds (default: 200)",
    )
    args = parser.parse_args()
    data_dir = Path(args.data_dir)
    if not data_dir.is_dir():
        raise SystemExit(f"Not a directory: {data_dir}")

    times, frames, columns = load_profiles(data_dir)
    nframes = len(frames)
    print(f"Loaded {nframes} profiles from {data_dir.resolve()}")
    print(f"  t in [{times[0]:.6e}, {times[-1]:.6e}]")
    print(f"  columns: {columns}")

    # Defaults: x = first column (usually 'x'), y = a useful scalar if present
    x_col = columns[0]
    preferred_y = ["Yv", "VF", "T", "P", "TG", "PG", "U", "RHOmix"]
    y_active = {c: False for c in columns}
    for name in preferred_y:
        if name in y_active and name != x_col:
            y_active[name] = True
            break
    if not any(y_active.values()):
        # Fall back to second column
        y_active[columns[min(1, len(columns) - 1)]] = True

    # Global axis limits over the whole run (so the scale does not jump frame-to-frame).
    # Use nanmin/nanmax: phasic columns mask inactive cells with NaN.
    # axis_limits() avoids zooming into roundoff when a field is flat.
    col_limits = {}
    for col in columns:
        vals = np.concatenate([fr[col] for fr in frames])
        col_limits[col] = axis_limits(vals)

    # ---- layout ----
    fig = plt.figure(figsize=(11, 6.5))
    fig.canvas.manager.set_window_title(f"relax1D profiles — {data_dir}")
    # Main axes leave room on the right for controls
    ax = fig.add_axes([0.08, 0.28, 0.58, 0.64])
    ax_xradio = fig.add_axes([0.72, 0.55, 0.24, 0.35])
    ax_ycheck = fig.add_axes([0.72, 0.18, 0.24, 0.32])
    ax_slider = fig.add_axes([0.08, 0.12, 0.58, 0.04])
    ax_play = fig.add_axes([0.08, 0.03, 0.12, 0.05])
    ax_info = fig.add_axes([0.72, 0.03, 0.24, 0.10])
    ax_info.axis("off")

    lines = {}  # col -> Line2D
    for col in columns:
        (ln,) = ax.plot([], [], "-o", ms=3, lw=1.2, label=col, visible=y_active[col])
        lines[col] = ln
    ax.set_xlabel(x_col)
    ax.set_ylabel("value")
    ax.grid(True, alpha=0.3)
    legend = ax.legend(loc="best", fontsize=8)

    info_text = ax_info.text(
        0.0,
        0.5,
        "",
        transform=ax_info.transAxes,
        va="center",
        fontsize=9,
        family="monospace",
    )

    state = {"iframe": 0, "playing": True, "x_col": x_col}

    def active_y_cols():
        return [c for c, on in y_active.items() if on and c != state["x_col"]]

    def redraw(iframe: int):
        iframe = int(np.clip(iframe, 0, nframes - 1))
        state["iframe"] = iframe
        frame = frames[iframe]
        xc = state["x_col"]
        x = frame[xc]
        any_y = False
        for col, ln in lines.items():
            on = y_active[col] and col != xc
            ln.set_visible(on)
            if on:
                ln.set_data(x, frame[col])
                any_y = True
            else:
                ln.set_data([], [])
        ax.set_xlabel(xc)
        # X limits: full spatial range of the chosen x column over the run
        xlo, xhi = col_limits[xc]
        ax.set_xlim(xlo, xhi)
        # Y limits: span of all active y fields over the entire simulation
        # (axis_limits already includes padding / roundoff floor)
        if any_y:
            ymin = min(col_limits[c][0] for c in active_y_cols())
            ymax = max(col_limits[c][1] for c in active_y_cols())
            ax.set_ylim(ymin, ymax)
        ax.set_title(f"frame {iframe + 1}/{nframes}")
        info_text.set_text(f"t = {times[iframe]:.6e} s\nx = {xc}\ny = {', '.join(active_y_cols()) or '(none)'}")
        # Refresh legend for visible lines only
        handles = [lines[c] for c in active_y_cols()]
        labels = active_y_cols()
        if handles:
            ax.legend(handles, labels, loc="best", fontsize=8)
        elif legend is not None:
            legend.remove()
        # Keep slider in sync without retriggering callbacks awkwardly
        if abs(time_slider.val - iframe) > 1e-9:
            time_slider.set_val(iframe)
        fig.canvas.draw_idle()

    # ---- widgets ----
    ax_xradio.set_title("X axis", fontsize=10)
    x_radio = RadioButtons(ax_xradio, columns, active=columns.index(x_col))
    for lab in x_radio.labels:
        lab.set_fontsize(9)

    ax_ycheck.set_title("Y axis", fontsize=10)
    y_check = CheckButtons(ax_ycheck, columns, [y_active[c] for c in columns])
    for lab in y_check.labels:
        lab.set_fontsize(9)

    time_slider = Slider(
        ax_slider,
        "frame",
        0,
        nframes - 1,
        valinit=0,
        valstep=1,
        valfmt="%d",
    )
    play_button = Button(ax_play, "Pause")

    def on_x(label):
        state["x_col"] = label
        # A column used as x should not also plot as y
        redraw(state["iframe"])

    def on_y(label):
        y_active[label] = not y_active[label]
        # Keep checkbox visual state in sync if user re-clicks x-as-y
        redraw(state["iframe"])

    def on_slider(val):
        iframe = int(round(val))
        if iframe != state["iframe"]:
            redraw(iframe)

    def on_play(_event):
        state["playing"] = not state["playing"]
        play_button.label.set_text("Pause" if state["playing"] else "Play")

    x_radio.on_clicked(on_x)
    y_check.on_clicked(on_y)
    time_slider.on_changed(on_slider)
    play_button.on_clicked(on_play)

    def animate(_frame_num):
        if state["playing"]:
            nxt = (state["iframe"] + 1) % nframes
            redraw(nxt)
        return tuple(lines.values())

    anim = animation.FuncAnimation(
        fig,
        animate,
        interval=args.interval,
        blit=False,
        cache_frame_data=False,
    )
    # Keep a reference so the animation is not garbage-collected
    fig._relax1d_anim = anim  # noqa: SLF001

    redraw(0)
    plt.show()


if __name__ == "__main__":
    main()
