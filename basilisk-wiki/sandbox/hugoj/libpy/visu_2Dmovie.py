#!/usr/bin/env python3
"""
visu_2Dmovie — render a 2D snapshot from a Basilisk/xarray NetCDF output.

Usage:
    visu_2Dmovie file.nc wave
    visu_2Dmovie file.nc wave --method interp
    visu_2Dmovie file.nc wave --var u.x --clim -2.5 2.5 --cmap bwr

Run `visu_2Dmovie --help` for the full list of options.

You might need to
chmod +x visu_2Dmovie.py
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
        prog="visu_3Dsnap",
        description="Render a 3D PyVista snapshot from a NetCDF simulation output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Positional args
    p.add_argument("input", type=str, help="Path to the input NetCDF file (.nc)")
    p.add_argument(
        "output",
        type=str,
        help="Output movie filename, with or without extension, e.g. 'wave' or 'wave.mp4'",
    )

    # p.add_argument("--skip", type=int, default=0, help="Skip n first steps")

    # Movie timing
    p.add_argument(
        "--fps", type=int, default=30, help="Frames per second of output video"
    )
    p.add_argument(
        "--speed",
        type=float,
        default=1.0,
        help="Playback speed: 1 sim-second = 1/speed video-seconds",
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
    skip=1,
    fps=30,
    speed=1.0,
    method="nearest",
    t_start=None,
    t_end=None,
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
    zclip=None,
):
    nameout = output + ".mp4"
    ds = xr.open_dataset(input)
    # clean up duplicates
    ds = ds.isel(time=~ds.indexes["time"].duplicated(keep="first"))

    if len(ds.y) == 1:
        raise Exception("You want to use a 2D plotter for a 1D simulation ... Aborting")

    t0 = t_start if t_start is not None else float(ds.time.values[0])
    t1 = t_end if t_end is not None else float(ds.time.values[-1])
    if t1 <= t0:
        raise ValueError(f"t_end ({t1}) must be greater than t_start ({t0})")
    duration_video = (t1 - t0) / speed
    nframes = int(np.round(duration_video * fps))
    if nframes < 1:
        raise ValueError("computed 0 frames — check fps / speed_factor / time range")
    target_times = t0 + np.arange(nframes) * (speed / fps)

    os.system("mkdir -p tmp")
    os.system("rm tmp/*")
    for it in range(len(target_times)):
        # print("t=%f" % target_times[it])
        pct = 100 * (it) / nframes

        datavar = ds[var].sel(time=target_times[it], method="nearest").isel(level=level)

        if detrend:
            datavarm = datavar.mean(dim=["x", "y"])
        else:
            datavarm = 0.0

        datavar -= datavarm

        print(f"  frame {it}/{nframes} ({pct:.0f}%)", end="\r", flush=True)
        fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True, dpi=200)
        s = ax.pcolormesh(
            ds.x,
            ds.y,
            datavar,
            cmap=cmap,
            vmin=clim[0],
            vmax=clim[1],
        )
        plt.colorbar(s, ax=ax, label=f"{var} ({unit})")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        # ax.set_ylim([-0.1, 0.1])
        # ax.set_xlim([-L0 / 2, L0 / 2])
        ax.set_aspect(1)
        ax.set_title("t = %.1f" % (target_times[it]))
        plt.savefig(f"tmp/{var}_t%05d.png" % it)
        plt.close(fig)

    os.system(f"rm {nameout}")
    os.system(
        f"ffmpeg -framerate 30 -i ./tmp/{var}_t%05d.png -c:v libx264 -pix_fmt yuv420p -r 30 {nameout}"
    )


def main(argv=None):
    args = parse_args(argv)
    try:
        render_movie2D(
            input=args.input,
            output=args.output,
            fps=args.fps,
            speed=args.speed,
            method=args.method,
            t_start=args.t_start,
            t_end=args.t_end,
            L0=args.L0,
            H0=args.H0,
            var=args.var,
            unit=args.unit,
            detrend=args.detrend,
            clim=args.clim,
            cmap=args.cmap,
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
