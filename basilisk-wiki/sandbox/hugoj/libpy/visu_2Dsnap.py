#!/usr/bin/env python3
"""
visu_2Dsnap — render a 2D snapshot from a Basilisk/xarray NetCDF output.

Usage:
    visu_2Dsnap file.nc wave
    visu_2Dsnap file.nc wave --method interp
    visu_2Dsnap file.nc wave --var u.x --clim -2.5 2.5 --cmap bwr

Run `visu_2Dsnap --help` for the full list of options.

You might need to
chmod +x visu_2Dsnap.py
and add it to your PATH
"""

import argparse
import sys
from pathlib import Path
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="visu_2Dsnap",
        description="Render a 2D snapshot from a NetCDF simulation output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Positional args
    p.add_argument("input", type=str, help="Path to the input NetCDF file (.nc)")
    p.add_argument(
        "output",
        type=str,
        help="Output filename, with or without extension, e.g. 'wave' or 'wave.png'",
    )
    p.add_argument("ttime", type=float, help="specific time to plot the snapshot")
    p.add_argument(
        "--method",
        choices=["nearest", "interp"],
        default="nearest",
        help="How to sample simulation time onto the video time axis",
    )
    p.add_argument(
        "--shading",
        choices=["nearest", "gouraud"],
        default="gouraud",
        help="Matplotlib shading to use, nearest=faster, gouraud=cleanest",
    )

    # Domain geometry
    p.add_argument("--L0", type=float, default=200.0, help="Horizontal domain size (m)")
    p.add_argument("--H0", type=float, default=50.0, help="Vertical domain size (m)")
    p.add_argument("--level", type=int, default=-1, help="Layer")

    # Fields to plot
    p.add_argument("--var", type=str, default="u.x", help="Variable plotted")
    p.add_argument("--unit", type=str, default="m/s", help="Unit")
    p.add_argument("--detrend", type=bool, default=False, help="Removes XY avg")
    p.add_argument(
        "--clim",
        type=float,
        nargs=2,
        default=[-2.5, 2.5],
        metavar=("MIN", "MAX"),
        help="Color limits for the variable",
    )
    p.add_argument(
        "--cmap",
        type=str,
        default="Greys_r",
        help="Colormap for the variable",
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


def render_movie2D(
    input,
    output,
    ttime,
    method="nearest",
    shading="gouraud",
    L0=200.0,
    H0=50.0,
    var="u.x",
    unit="m/s",
    detrend=False,
    level=-1,
    clim=(-2.5, 2.5),
    cmap="Greys_r",
    xclip=None,
    yclip=None,
):
    ds = xr.open_dataset(input)
    ds = ds.isel(time=~ds.indexes["time"].duplicated(keep="first"))

    if (len(ds.y) == 1) or ("y" not in ds.keys()):
        raise Exception("You want to use a 2D plotter for a 1D simulation ... Aborting")

    datavar = ds[var].sel(time=ttime, method="nearest").isel(level=level)

    if detrend:
        datavarm = datavar.mean(dim=["x", "y"])
    else:
        datavarm = 0.0

    datavar -= datavarm

    if xclip is None:
        xclip = (-L0 / 2, L0 / 2)
    if yclip is None:
        yclip = (-L0 / 2, L0 / 2)

    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True, dpi=200)
    s = ax.pcolormesh(
        ds.x,
        ds.y,
        datavar,
        cmap=cmap,
        vmin=clim[0],
        vmax=clim[1],
        shading=shading,
    )
    plt.colorbar(s, ax=ax, label=f"{var} ({unit})")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_ylim(xclip)
    ax.set_xlim(yclip)
    ax.set_aspect(1)
    if level == -1:
        nlevel = "sfx"
    else:
        nlevel = level
    ax.set_title(f"t = %.1f (level={nlevel})" % (ttime))
    plt.savefig(f"{output}_{var}_t%.1d_l{nlevel}.png" % ttime)
    plt.close(fig)


def main(argv=None):
    args = parse_args(argv)
    try:
        render_movie2D(
            input=args.input,
            output=args.output,
            ttime=args.ttime,
            method=args.method,
            shading=args.shading,
            L0=args.L0,
            H0=args.H0,
            var=args.var,
            unit=args.unit,
            detrend=args.detrend,
            level=args.level,
            clim=args.clim,
            cmap=args.cmap,
            xclip=args.xclip,
            yclip=args.yclip,
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
