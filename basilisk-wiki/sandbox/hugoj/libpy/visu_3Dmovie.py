#!/usr/bin/env python3
"""
visu_3Dmovie — render a 3D PyVista movie from a Basilisk/xarray NetCDF output.

Usage:
    visu_3Dmovie file.nc wave
    visu_3Dmovie file.nc wave --fps 24 --speed-factor 0.5 --method interp
    visu_3Dmovie file.nc wave --var-top u.x --clim-top -2.5 2.5 --cmap-top bwr

Run `visu_3Dmovie --help` for the full list of options.

You might need to
chmod +x visu_3Dmovie.py
and add it to your PATH
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import xarray as xr

from tools import build_z_ds
from tools_3Dplots import (
    build_3Dgrid,
    add_mesh_domain_side,
    add_ticks,
    clip_physical_to_index,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="visu_3Dmovie",
        description="Render a 3D PyVista movie from a NetCDF simulation output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional args
    p.add_argument("input", type=str, help="Path to the input NetCDF file (.nc)")
    p.add_argument(
        "output",
        type=str,
        help="Output movie filename, with or without extension, e.g. 'wave' or 'wave.mp4'",
    )

    p.add_argument("--skip", type=int, default=0, help="Skip n first steps")

    # Movie timing
    p.add_argument(
        "--fps", type=int, default=30, help="Frames per second of output video"
    )
    p.add_argument(
        "--speed-factor",
        type=float,
        default=1.0,
        help="Playback speed: 1 sim-second = 1/speed_factor video-seconds",
    )
    p.add_argument(
        "--method",
        choices=["nearest", "interp"],
        default="nearest",
        help="How to sample simulation time onto the video time axis",
    )
    p.add_argument(
        "--t-start",
        type=float,
        default=None,
        help="Start time (default: first available)",
    )
    p.add_argument(
        "--t-end", type=float, default=None, help="End time (default: last available)"
    )

    # Domain geometry
    p.add_argument("--L0", type=float, default=200.0, help="Horizontal domain size (m)")
    p.add_argument("--H0", type=float, default=50.0, help="Vertical domain size (m)")

    # Fields to plot
    p.add_argument(
        "--var-top", type=str, default="u.x", help="Variable plotted on top surface"
    )
    p.add_argument(
        "--clim-top",
        type=float,
        nargs=2,
        default=[-2.5, 2.5],
        metavar=("MIN", "MAX"),
        help="Color limits for top-surface variable",
    )
    p.add_argument(
        "--cmap-top",
        type=str,
        default="seismic",
        help="Colormap for top-surface variable",
    )

    p.add_argument(
        "--var-side", type=str, default="u.x", help="Variable plotted on side surfaces"
    )
    p.add_argument(
        "--clim-side",
        type=float,
        nargs=2,
        default=[-2.5, 2.5],
        metavar=("MIN", "MAX"),
        help="Color limits for side-surface variable",
    )
    p.add_argument(
        "--cmap-side",
        type=str,
        default="seismic",
        help="Colormap for side-surface variable",
    )

    # Ticks
    p.add_argument(
        "--tick-len-min",
        type=float,
        default=3.0,
        help="Minor tick length (world units)",
    )
    p.add_argument(
        "--tick-len-maj",
        type=float,
        default=6.0,
        help="Major tick length (world units)",
    )
    p.add_argument(
        "--label-offset",
        type=float,
        default=10.0,
        help="Tick label offset (world units)",
    )

    # Camera / rendering
    p.add_argument(
        "--camera-position",
        type=str,
        default="yz",
        help="PyVista camera preset, e.g. 'yz', 'iso'",
    )
    p.add_argument(
        "--elevation", type=float, default=25.0, help="Camera elevation (degrees)"
    )
    p.add_argument(
        "--azimuth", type=float, default=-120.0, help="Camera azimuth (degrees)"
    )
    p.add_argument(
        "--background", type=str, default="lightgrey", help="Background color"
    )
    p.add_argument(
        "--window-size",
        type=int,
        nargs=2,
        default=[1024, 1024],
        metavar=("W", "H"),
        help="Render window size in pixels",
    )
    p.add_argument(
        "--off-screen", action="store_true", help="Render off-screen (no window popup)"
    )

    p.add_argument(
        "--xclip",
        type=float,
        nargs=2,
        default=None,
        metavar=("xmin", "xmax"),
        help="X direction clip. Default is full domain",
    )
    p.add_argument(
        "--yclip",
        type=float,
        nargs=2,
        default=None,
        metavar=("ymin", "ymax"),
        help="Y direction clip. Default is full domain",
    )
    p.add_argument(
        "--zclip",
        type=float,
        nargs=2,
        default=None,
        metavar=("zmin", "zmax"),
        help="Z direction clip. Default is full domain. The clip is made by guessing an index from the averaged height z",
    )

    return p.parse_args(argv)


def render_movie(
    input,
    output,
    skip=1,
    fps=30,
    speed_factor=1.0,
    method="nearest",
    t_start=None,
    t_end=None,
    L0=200.0,
    H0=50.0,
    var_top="u.x",
    clim_top=(-2.5, 2.5),
    cmap_top="seismic",
    var_side="u.x",
    clim_side=(-2.5, 2.5),
    cmap_side="seismic",
    tick_len_min=3.0,
    tick_len_maj=6.0,
    label_offset=10.0,
    camera_position="yz",
    elevation=25.0,
    azimuth=-120.0,
    background="lightgrey",
    window_size=(1024, 1024),
    off_screen=False,
    verbose=True,
    xclip=None,
    yclip=None,
    zclip=None,
):
    """
    Render a 3D PyVista movie from a NetCDF simulation output.

    Same logic as the `visu3D_movie` command-line tool, exposed as a
    plain function so it can be called from other Python code, e.g.:

        from visu3D_movie import render_movie
        render_movie("file.nc", "wave", fps=24, speed_factor=0.5)

    Parameters mirror the CLI flags (see `visu3D_movie --help`).
    Returns the path to the written movie file (str).
    """
    input_path = Path(input).expanduser()
    if not input_path.exists():
        raise FileNotFoundError(f"input file not found: {input_path}")

    output_path = Path(output)
    filemovie = (
        str(output_path)
        if output_path.suffix.lower() in (".mp4", ".avi", ".mov")
        else f"{output_path}.mp4"
    )

    vartop = {"name": var_top, "clim": list(clim_top), "cmap": cmap_top}
    varside = {"name": var_side, "clim": list(clim_side), "cmap": cmap_side}

    # --- Open dataset ---
    dst = xr.open_dataset(input_path)
    dst = dst.isel(time=slice(skip, len(dst.time)))
    H0 = -dst.zb.values.flatten()[0]

    # Clipping
    if xclip is None:
        xmin, xmax = -L0 / 2, L0 / 2
    else:
        xmin, xmax = xclip
    if yclip is None:
        ymin, ymax = -L0 / 2, L0 / 2
    else:
        ymin, ymax = yclip
    if zclip is None:
        zmin, zmax = -H0, 0.0
    else:
        zmin, zmax = zclip
    physical_clip = xmin, xmax, ymin, ymax, zmin, zmax

    t0 = t_start if t_start is not None else float(dst.time.values[0])
    t1 = t_end if t_end is not None else float(dst.time.values[-1])
    if t1 <= t0:
        raise ValueError(f"t_end ({t1}) must be greater than t_start ({t0})")

    duration_video = (t1 - t0) / speed_factor
    nframes = int(np.round(duration_video * fps))
    if nframes < 1:
        raise ValueError("computed 0 frames — check fps / speed_factor / time range")
    target_times = t0 + np.arange(nframes) * (speed_factor / fps)

    nx, ny, nz = dst.x.size, dst.y.size, dst.level.size

    if verbose:
        print(f"Input:        {input_path}")
        print(f"Output:       {filemovie}")
        print(f"Time range:   [{t0}, {t1}]  ({method})")
        print(f"fps:          {fps}")
        print(f"Frames:       {nframes}  (video duration: {nframes / fps:.2f} s)")
        print(f"Grid size:    {nx} x {ny} x {nz}")
        print("Domain size:")
        print(f"    X: [{xmin}, {xmax}] m")
        print(f"    Y: [{ymin}, {ymax}] m")
        print(f"    Z: [{zmin}, {zmax}] m")

    if nframes > len(dst.time):
        print(
            f"Warning: you asked for {nframes} frames but the file has only {len(dst.time)} points."
        )
        print("         results can behave strangely !")

    # --- Plot setup ---
    plotter = pv.Plotter(window_size=tuple(window_size), off_screen=off_screen)
    plotter.open_movie(filemovie, framerate=fps)

    # nice box around domain, default is full domain
    box = pv.Box(bounds=(xmin, xmax, ymin, ymax, zmin, zmax))
    plotter.add_mesh(box, style="wireframe", color="black", line_width=1)

    add_ticks(
        plotter, L0, H0, tick_len_min, tick_len_maj, label_offset, clip=physical_clip
    )

    # Camera position
    plotter.camera_position = camera_position
    plotter.camera.elevation = elevation
    plotter.camera.azimuth = azimuth
    plotter.enable_3_lights()
    plotter.set_background(background)

    # better view of the simu, with space for colorbar
    # increase shiftV : cube goes up
    # increase shiftH : cube goes left
    shiftV = 1 / 3 * H0
    shiftH = 1 / 5 * H0
    pos = plotter.camera.position
    fp = plotter.camera.focal_point
    plotter.camera.focal_point = (fp[0] + shiftH, fp[1], fp[2] - shiftV)
    plotter.camera.position = (pos[0] + shiftH, pos[1], pos[2] - shiftV)
    plotter.camera.view_angle = 35

    # --- Frame loop ---
    for it_video, t_target in enumerate(target_times):
        if method == "interp":
            ds = dst.interp(time=t_target)
        else:
            ds = dst.sel(time=t_target, method="nearest")

        ds = build_z_ds(ds, "float32")
        grid, xv, yv, zv = build_3Dgrid(ds)

        T = ds[var_side].values.transpose(2, 1, 0)

        grid.cell_data[var_side] = T.ravel(order="F")
        UX = ds[var_top].values.transpose(2, 1, 0)
        grid.cell_data[var_top] = UX.ravel(order="F")

        # Clipping
        if physical_clip is None:
            iclip = None
        else:
            zvm = np.mean(zv, axis=(1, 2))
            iclip = clip_physical_to_index(physical_clip, xv, yv, zvm)

        # Disable surface from previous timestep
        plotter.remove_actor("top_surface", render=False)
        plotter.remove_actor("side_surface", render=False)

        # Add mesh on the sides and the top
        add_mesh_domain_side(grid, plotter, vartop, varside, index_clip=iclip)

        plotter.write_frame()

        if verbose and (
            (it_video + 1) % max(1, nframes // 20) == 0 or it_video == nframes - 1
        ):
            pct = 100 * (it_video + 1) / nframes
            print(
                f"  frame {it_video + 1}/{nframes} ({pct:.0f}%)", end="\r", flush=True
            )

    if verbose:
        print()  # newline after progress
    plotter.close()
    if verbose:
        print(f"Done. Movie written to {filemovie}")
    return filemovie


def main(argv=None):
    args = parse_args(argv)
    try:
        render_movie(
            input=args.input,
            output=args.output,
            fps=args.fps,
            speed_factor=args.speed_factor,
            method=args.method,
            t_start=args.t_start,
            t_end=args.t_end,
            L0=args.L0,
            H0=args.H0,
            var_top=args.var_top,
            clim_top=args.clim_top,
            cmap_top=args.cmap_top,
            var_side=args.var_side,
            clim_side=args.clim_side,
            cmap_side=args.cmap_side,
            tick_len_min=args.tick_len_min,
            tick_len_maj=args.tick_len_maj,
            label_offset=args.label_offset,
            camera_position=args.camera_position,
            elevation=args.elevation,
            azimuth=args.azimuth,
            background=args.background,
            window_size=args.window_size,
            off_screen=args.off_screen,
            xclip=args.xclip,
            yclip=args.yclip,
            zclip=args.zclip,
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
