#!/usr/bin/env python3
"""
visu_3Dsnap — render a 3D PyVista snapshot from a Basilisk/xarray NetCDF output.

Usage:
    visu_3Dsnap file.nc wave
    visu_3Dsnap file.nc wave --method interp
    visu_3Dsnap file.nc wave --var-top u.x --clim-top -2.5 2.5 --cmap-top bwr

Run `visu_3Dsnap --help` for the full list of options.

You might need to
chmod +x visu_3Dsnap.py
and add it to your PATH
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import xarray as xr

from tools import build_z_ds
from tools_3Dplots import build_3Dgrid, add_mesh_domain_side, add_ticks


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="visu_3Dsnap",
        description="Render a 3D PyVista snapshot from a NetCDF simulation output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional args
    p.add_argument("input", type=str, help="Path to the input NetCDF file (.nc)")
    p.add_argument(
        "output",
        type=str,
        help="Output snapshot filename. Added is the time and format (pdf)",
    )

    p.add_argument("ttime", type=float, help="specific time to plot the snapshot")

    p.add_argument(
        "--method",
        choices=["nearest", "interp"],
        default="nearest",
        help="How to sample simulation time onto the video time axis",
    )

    # Domain geometry
    p.add_argument("--L0", type=float, default=200.0, help="Horizontal domain size (m)")
    p.add_argument("--H0", type=float, default=50.0, help="Vertical domain size (m)")
    p.add_argument(
        "--surf-level", type=int, default=-1, help="Level index of the free surface"
    )

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
        "--cmap-top", type=str, default="bwr", help="Colormap for top-surface variable"
    )

    p.add_argument(
        "--var-side", type=str, default="T", help="Variable plotted on side surfaces"
    )
    p.add_argument(
        "--clim-side",
        type=float,
        nargs=2,
        default=[19.95, 20.0],
        metavar=("MIN", "MAX"),
        help="Color limits for side-surface variable",
    )
    p.add_argument(
        "--cmap-side",
        type=str,
        default="plasma",
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

    return p.parse_args(argv)


def render_snapshot(
    input,
    output,
    ttime=0.0,
    method="nearest",
    L0=200.0,
    H0=50.0,
    var_top="u.x",
    clim_top=(-2.5, 2.5),
    cmap_top="seismic",
    var_side="T",
    clim_side=(19.95, 20.0),
    cmap_side="plasma",
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
):

    input_path = Path(input).expanduser()
    if not input_path.exists():
        raise FileNotFoundError(f"input file not found: {input_path}")
    output_path = Path(output)
    output_path = f"{Path(output)}_t{ttime}.pdf"

    vartop = {"name": var_top, "clim": list(clim_top), "cmap": cmap_top}
    varside = {"name": var_side, "clim": list(clim_side), "cmap": cmap_side}
    # --- Open dataset ---
    dst = xr.open_dataset(input_path, chunks="auto")
    nx, ny, nz = dst.x.size, dst.y.size, dst.level.size
    t0 = dst.time.values[0]
    if "time" not in dst.keys():
        ds = dst
    else:
        ds = dst.sel(time=ttime, method=method)

    if verbose:
        print(f"Input:        {input_path}")
        print(f"Output:       {output_path}")
        print(f"Time:         {t0}")
        print(f"Method:       {method}")
        print(f"Grid size:    {nx} x {ny} x {nz}")

    # --- Plot setup ---
    plotter = pv.Plotter(window_size=tuple(window_size), off_screen=off_screen)

    box = pv.Box(bounds=(-L0 / 2, L0 / 2, -L0 / 2, L0 / 2, -H0, 0))
    plotter.add_mesh(box, style="wireframe", color="black", line_width=1)

    add_ticks(plotter, L0, H0, tick_len_min, tick_len_maj, label_offset)

    plotter.camera_position = camera_position
    plotter.camera.elevation = elevation
    plotter.camera.azimuth = azimuth
    plotter.enable_3_lights()
    plotter.set_background(background)

    ds = build_z_ds(ds, "float32")
    grid = build_3Dgrid(ds)

    T = ds[var_side].values.transpose(2, 1, 0)

    grid.cell_data[var_side] = T.ravel(order="F")
    UX = ds[var_top].values.transpose(2, 1, 0)
    grid.cell_data[var_top] = UX.ravel(order="F")
    add_mesh_domain_side(grid, (nx, ny, nz), plotter, vartop, varside)

    plotter.save_graphic(output_path)
    plotter.close()
    if verbose:
        print(f"Done. Image written to {output_path}")
    return output_path


def main(argv=None):
    args = parse_args(argv)
    try:
        render_snapshot(
            input=args.input,
            output=args.output,
            ttime=args.ttime,
            method=args.method,
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
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
