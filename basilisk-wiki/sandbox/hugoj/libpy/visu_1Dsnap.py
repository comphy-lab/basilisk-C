#!/usr/bin/env python3
"""
visu_profile_snap — render a 1D profile snapshot (matplotlib) from an xarray
NetCDF output.

Usage:
    visu_profile_snap stokes/out.nc snap 2.5
    visu_profile_snap stokes/out.nc snap 2.5 --method interp --out svg
    visu_profile_snap stokes/out.nc snap 2.5 --var eta --L0 2 --T0 1.5

Run `visu_profile_snap --help` for the full list of options.

You might need to
    chmod +x visu_profile_snap.py
and add it to your PATH.
"""

import argparse
import sys
from pathlib import Path

from tools_1Dplots import (
    ProfilePlot,
    add_plot_arguments,
    make_title,
    open_profile_dataset,
    resolve_limits,
    select_time,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="visu_profile_snap",
        description="Render a 1D profile snapshot from a NetCDF simulation output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional args
    p.add_argument("input", type=str, help="Path to the input NetCDF file (.nc)")
    p.add_argument(
        "output",
        type=str,
        help="Output snapshot filename. Added is the time and format",
    )
    p.add_argument("ttime", type=float, help="Time at which to plot the snapshot")
    p.add_argument(
        "--out",
        choices=["pdf", "svg", "png"],
        default="pdf",
        help="Output format",
    )

    # Variable, limits, style (shared with visu_profile_movie)
    add_plot_arguments(p)
    return p.parse_args(argv)


def render_snapshot1D(
    input,
    output,
    ttime=0.0,
    method="nearest",
    var="eta",
    y_index=0,
    L0=None,
    T0=None,
    xlim=None,
    ylim=None,
    xlabel="x (m)",
    ylabel="z (m)",
    color="C0",
    linewidth=1.5,
    equal_aspect=False,
    figsize=(6, 3),
    dpi=100,
    out="pdf",
    verbose=True,
):
    """
    Render a 1D profile snapshot from a NetCDF simulation output.

    Same logic as the `visu_profile_snap` command-line tool, exposed as a
    plain function. Returns the path to the written file (str).
    """
    output_path = f"{Path(output)}_t{ttime}.{out}"

    # --- Open dataset ---
    ds = open_profile_dataset(input, y_index=y_index)
    xlim, ylim = resolve_limits(ds, var, xlim=xlim, ylim=ylim, L0=L0)

    # --- Select a timestamp ---
    da = select_time(ds[var], ttime, method)
    t_used = float(da.time.values) if "time" in da.coords else ttime

    if verbose:
        print(f"Input: {input}")
        print(f"Output: {output_path}")
        print(f"Variable: {var}")
        print(f"Time: {t_used}")
        print(f"Method: {method}")
        print(f"x-limits: {xlim}  y-limits: {ylim}")

    # --- Plot ---
    plot = ProfilePlot(
        ds.x.values,
        xlim,
        ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        figsize=figsize,
        dpi=dpi,
        color=color,
        linewidth=linewidth,
        equal_aspect=equal_aspect,
    )
    plot.update(da.values, make_title(t_used, T0))
    plot.save(output_path)
    plot.close()

    if verbose:
        print(f"Done. Image written to {output_path}")
    return output_path


def main(argv=None):
    args = parse_args(argv)
    try:
        render_snapshot1D(
            input=args.input,
            output=args.output,
            ttime=args.ttime,
            method=args.method,
            var=args.var,
            y_index=args.y_index,
            L0=args.L0,
            T0=args.T0,
            xlim=args.xlim,
            ylim=args.ylim,
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            color=args.color,
            linewidth=args.linewidth,
            equal_aspect=args.equal_aspect,
            figsize=args.figsize,
            dpi=args.dpi,
            out=args.out,
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
