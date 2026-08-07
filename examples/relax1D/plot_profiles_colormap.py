#!/usr/bin/env python3
"""Colormap / filled-rectangle viewer for relax1D profiles_*.csv dumps.

Same data and time controls as plot_profiles_interactive.py, but each 1D
profile is extruded into a 2D rectangle and shown with a colormap:

  * horizontal axis  = chosen spatial column (usually x)
  * vertical axis    = artificial transverse coordinate (fills a band)
  * color            = chosen field value at that x

Also includes a space–time panel (x vs t, color = field) so you can see
the whole run at a glance while the top panel animates.

Usage:
  python3 plot_profiles_colormap.py [data_dir]

Example:
  cd examples/relax1D/run_1mm
  python3 ../plot_profiles_colormap.py data
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import RadioButtons, Button, Slider
from matplotlib.colors import Normalize
import numpy as np

# Reuse the CSV loader from the line-plot viewer
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_profiles_interactive import load_profiles, axis_limits  # noqa: E402


def build_strip(x: np.ndarray, field: np.ndarray, n_transverse: int = 32):
    """Extrude a 1D profile into a 2D array for pcolormesh.

    Returns (X, Y, Z) suitable for ax.pcolormesh(X, Y, Z, shading='flat')
    where Y spans [0, 1] (dimensionless transverse fill).
    """
    # Cell edges from cell centers (or pass-through if already edges)
    if len(x) < 2:
        xe = np.array([x[0] - 0.5, x[0] + 0.5])
    else:
        dx = np.diff(x)
        xe = np.empty(len(x) + 1)
        xe[1:-1] = 0.5 * (x[:-1] + x[1:])
        xe[0] = x[0] - 0.5 * dx[0]
        xe[-1] = x[-1] + 0.5 * dx[-1]
    ye = np.linspace(0.0, 1.0, n_transverse + 1)
    X, Y = np.meshgrid(xe, ye)
    # Z is (n_transverse, nx) — constant in the transverse direction
    Z = np.tile(field[np.newaxis, :], (n_transverse, 1))
    return X, Y, Z


def stack_spacetime(times, frames, x_col: str, f_col: str):
    """Build (x_edges, t_edges, Z) for a space–time pcolormesh."""
    x = frames[0][x_col]
    if len(x) < 2:
        xe = np.array([x[0] - 0.5, x[0] + 0.5])
    else:
        dx = np.diff(x)
        xe = np.empty(len(x) + 1)
        xe[1:-1] = 0.5 * (x[:-1] + x[1:])
        xe[0] = x[0] - 0.5 * dx[0]
        xe[-1] = x[-1] + 0.5 * dx[-1]

    # Time edges mid-way between dumps; extend ends by half the neighbor dt
    t = np.asarray(times, dtype=float)
    if len(t) == 1:
        te = np.array([t[0] - 0.5, t[0] + 0.5])
    else:
        dt = np.diff(t)
        te = np.empty(len(t) + 1)
        te[1:-1] = 0.5 * (t[:-1] + t[1:])
        te[0] = t[0] - 0.5 * dt[0]
        te[-1] = t[-1] + 0.5 * dt[-1]

    Z = np.vstack([fr[f_col] for fr in frames])  # (nt, nx)
    return xe, te, Z


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
        help="Animation frame interval in ms (default: 200)",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap name (default: viridis)",
    )
    parser.add_argument(
        "--ny",
        type=int,
        default=48,
        help="Transverse resolution of the filled strip (default: 48)",
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

    x_col = columns[0]
    preferred = ["Yv", "VF", "T", "P", "TG", "PG", "U", "RHOmix", "TL", "PL"]
    f_col = next((c for c in preferred if c in columns and c != x_col), columns[min(1, len(columns) - 1)])

    # Global color limits over all times for a stable scale while animating.
    # axis_limits avoids a roundoff-zoomed colorbar on flat fields (e.g. P≡1e5).
    def field_limits(col: str):
        vals = np.concatenate([fr[col] for fr in frames])
        return axis_limits(vals)

    # ---- layout ----
    fig = plt.figure(figsize=(12, 7.5))
    fig.canvas.manager.set_window_title(f"relax1D colormap — {data_dir}")
    # Top: animated filled strip; bottom: space–time
    ax_strip = fig.add_axes([0.08, 0.52, 0.58, 0.40])
    ax_st = fig.add_axes([0.08, 0.18, 0.58, 0.28])
    ax_cbar = fig.add_axes([0.67, 0.18, 0.015, 0.74])
    ax_xradio = fig.add_axes([0.74, 0.55, 0.22, 0.35])
    ax_fradio = fig.add_axes([0.74, 0.18, 0.22, 0.32])
    ax_slider = fig.add_axes([0.08, 0.08, 0.58, 0.035])
    ax_play = fig.add_axes([0.08, 0.02, 0.12, 0.04])
    ax_info = fig.add_axes([0.74, 0.02, 0.22, 0.12])
    ax_info.axis("off")

    lo, hi = field_limits(f_col)
    norm = Normalize(vmin=lo, vmax=hi)

    # Initial strip
    X, Y, Z = build_strip(frames[0][x_col], frames[0][f_col], n_transverse=args.ny)
    mesh_strip = ax_strip.pcolormesh(X, Y, Z, cmap=args.cmap, norm=norm, shading="flat")
    ax_strip.set_ylabel("transverse (fill)")
    ax_strip.set_yticks([0, 1])
    ax_strip.set_yticklabels(["0", "1"])
    ax_strip.set_title("filled 1D profile")

    # Space–time
    xe, te, Zst = stack_spacetime(times, frames, x_col, f_col)
    mesh_st = ax_st.pcolormesh(xe, te, Zst, cmap=args.cmap, norm=norm, shading="flat")
    # Time cursor
    t_cursor = ax_st.axhline(times[0], color="w", lw=1.2, ls="--", alpha=0.9)
    ax_st.set_xlabel(x_col)
    ax_st.set_ylabel("t [s]")
    ax_st.set_title("space–time")

    cbar = fig.colorbar(mesh_strip, cax=ax_cbar)
    cbar.set_label(f_col)

    info_text = ax_info.text(
        0.0,
        0.5,
        "",
        transform=ax_info.transAxes,
        va="center",
        fontsize=9,
        family="monospace",
    )

    state = {"iframe": 0, "playing": True, "x_col": x_col, "f_col": f_col}

    def rebuild_spacetime():
        nonlocal mesh_st
        xe, te, Zst = stack_spacetime(times, frames, state["x_col"], state["f_col"])
        mesh_st.remove()
        mesh_st = ax_st.pcolormesh(xe, te, Zst, cmap=args.cmap, norm=norm, shading="flat")
        ax_st.set_xlabel(state["x_col"])

    def redraw(iframe: int):
        nonlocal mesh_strip
        iframe = int(np.clip(iframe, 0, nframes - 1))
        state["iframe"] = iframe
        fr = frames[iframe]
        xc, fc = state["x_col"], state["f_col"]
        X, Y, Z = build_strip(fr[xc], fr[fc], n_transverse=args.ny)
        mesh_strip.remove()
        mesh_strip = ax_strip.pcolormesh(X, Y, Z, cmap=args.cmap, norm=norm, shading="flat")
        ax_strip.set_xlim(float(X.min()), float(X.max()))
        ax_strip.set_ylim(0.0, 1.0)
        ax_strip.set_title(f"filled 1D profile — frame {iframe + 1}/{nframes}")
        t_cursor.set_ydata([times[iframe], times[iframe]])
        cbar.set_label(fc)
        info_text.set_text(
            f"t = {times[iframe]:.6e} s\n"
            f"x = {xc}\n"
            f"color = {fc}\n"
            f"range = [{norm.vmin:.3g}, {norm.vmax:.3g}]"
        )
        if abs(time_slider.val - iframe) > 1e-9:
            time_slider.set_val(iframe)
        fig.canvas.draw_idle()

    def on_x(label):
        state["x_col"] = label
        if state["f_col"] == label:
            # Pick another field if x and color collide
            for c in columns:
                if c != label:
                    state["f_col"] = c
                    f_radio.set_active(columns.index(c))
                    break
        lo, hi = field_limits(state["f_col"])
        norm.vmin, norm.vmax = lo, hi
        mesh_strip.set_norm(norm)
        mesh_st.set_norm(norm)
        rebuild_spacetime()
        redraw(state["iframe"])

    def on_f(label):
        if label == state["x_col"]:
            # Ignore coloring by the spatial coordinate itself
            f_radio.set_active(columns.index(state["f_col"]))
            return
        state["f_col"] = label
        lo, hi = field_limits(label)
        norm.vmin, norm.vmax = lo, hi
        mesh_strip.set_clim(lo, hi)
        mesh_st.set_clim(lo, hi)
        rebuild_spacetime()
        redraw(state["iframe"])

    def on_slider(val):
        iframe = int(round(val))
        if iframe != state["iframe"]:
            redraw(iframe)

    def on_play(_event):
        state["playing"] = not state["playing"]
        play_button.label.set_text("Pause" if state["playing"] else "Play")

    ax_xradio.set_title("X axis", fontsize=10)
    x_radio = RadioButtons(ax_xradio, columns, active=columns.index(x_col))
    for lab in x_radio.labels:
        lab.set_fontsize(9)

    ax_fradio.set_title("Color field", fontsize=10)
    f_radio = RadioButtons(ax_fradio, columns, active=columns.index(f_col))
    for lab in f_radio.labels:
        lab.set_fontsize(9)

    time_slider = Slider(ax_slider, "frame", 0, nframes - 1, valinit=0, valstep=1, valfmt="%d")
    play_button = Button(ax_play, "Pause")

    x_radio.on_clicked(on_x)
    f_radio.on_clicked(on_f)
    time_slider.on_changed(on_slider)
    play_button.on_clicked(on_play)

    def animate(_frame_num):
        if state["playing"]:
            redraw((state["iframe"] + 1) % nframes)
        return (mesh_strip,)

    anim = animation.FuncAnimation(
        fig, animate, interval=args.interval, blit=False, cache_frame_data=False
    )
    fig._relax1d_anim = anim  # noqa: SLF001

    redraw(0)
    plt.show()


if __name__ == "__main__":
    main()
