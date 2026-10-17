#!/usr/bin/env python3
"""
visu_1Dmovie — render a 1D profile movie (matplotlib + ffmpeg) from an
xarray NetCDF output.

Usage:
    visu_1Dmovie stokes/out.nc movie
    visu_1Dmovie stokes/out.nc movie --fps 24 --speed-factor 0.5 --method interp
    visu_1Dmovie stokes/out.nc movie --var eta --L0 2 --T0 1.5 --ylim -0.1 0.1

Run `visu_profile_movie --help` for the full list of options.

You might need to
    chmod +x visu_profile_movie.py
and add it to your PATH. Requires `ffmpeg` to be installed.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

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
        prog="visu_profile_movie",
        description="Render a 1D profile movie from a NetCDF simulation output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional args
    p.add_argument("input", type=str, help="Path to the input NetCDF file (.nc)")
    p.add_argument(
        "output",
        type=str,
        help="Output movie filename, with or without extension, e.g. 'movie' or 'movie.mp4'",
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
        "--t-start",
        type=float,
        default=None,
        help="Start time (default: first available)",
    )
    p.add_argument(
        "--t-end", type=float, default=None, help="End time (default: last available)"
    )
    p.add_argument(
        "--keep-frames",
        type=str,
        default=None,
        metavar="DIR",
        help="Keep the PNG frames in DIR (default: temporary dir, deleted at the end)",
    )

    # Variable, limits, style (shared with visu_profile_snap)
    add_plot_arguments(p)
    return p.parse_args(argv)


def render_movie1D(
    input,
    output,
    skip=0,
    fps=30,
    speed_factor=1.0,
    method="nearest",
    t_start=None,
    t_end=None,
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
    equal_aspect=True,
    figsize=(6, 3),
    dpi=100,
    keep_frames=None,
    verbose=True,
):
    """
    Render a 1D profile movie from a NetCDF simulation output.

    Same logic as the `visu_profile_movie` command-line tool, exposed as a
    plain function so it can be called from other Python code, e.g.:

        from visu_profile_movie import render_movie
        render_movie("file.nc", "movie", fps=24, speed_factor=0.5)

    Parameters mirror the CLI flags (see `visu_profile_movie --help`).
    Returns the path to the written movie file (str).
    """
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg not found in PATH")

    output_path = Path(output)
    filemovie = (
        str(output_path)
        if output_path.suffix.lower() in (".mp4", ".avi", ".mov")
        else f"{output_path}.mp4"
    )

    # --- Open dataset ---
    ds = open_profile_dataset(input, skip=skip, y_index=y_index)
    xlim, ylim = resolve_limits(ds, var, xlim=xlim, ylim=ylim, L0=L0)
    da = ds[var]

    t0 = t_start if t_start is not None else float(ds.time.values[0])
    t1 = t_end if t_end is not None else float(ds.time.values[-1])
    if t1 <= t0:
        raise ValueError(f"t_end ({t1}) must be greater than t_start ({t0})")

    duration_video = (t1 - t0) / speed_factor
    nframes = int(np.round(duration_video * fps))
    if nframes < 1:
        raise ValueError("computed 0 frames — check fps / speed_factor / time range")
    target_times = t0 + np.arange(nframes) * (speed_factor / fps)

    if verbose:
        print(f"Input: {input}")
        print(f"Output: {filemovie}")
        print(f"Variable: {var}")
        print(f"Time range: [{t0}, {t1}] ({method})")
        print(f"fps: {fps}")
        print(f"Frames: {nframes} (video duration: {nframes / fps:.2f} s)")
        print(f"x-limits: {xlim}  y-limits: {ylim}")
        if nframes > ds.time.size:
            print(
                f"Warning: you asked for {nframes} frames but the file has only {ds.time.size} points."
            )
            print("         results can behave strangely !")

    # --- Plot setup (figure is created once, only the line is updated) ---
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

    # --- Frame loop ---
    if keep_frames is not None:
        tmpdir = Path(keep_frames)
        tmpdir.mkdir(parents=True, exist_ok=True)
        cleanup = None
    else:
        cleanup = tempfile.TemporaryDirectory(prefix="visu_profile_")
        tmpdir = Path(cleanup.name)

    try:
        for it_video, t_target in enumerate(target_times):
            values = select_time(da, t_target, method).values
            plot.update(values, make_title(t_target, T0))
            plot.save(tmpdir / f"frame_{it_video:05d}.png")

            if verbose and (
                (it_video + 1) % max(1, nframes // 20) == 0 or it_video == nframes - 1
            ):
                pct = 100 * (it_video + 1) / nframes
                print(
                    f"  frame {it_video + 1}/{nframes} ({pct:.0f}%)",
                    end="\r",
                    flush=True,
                )
        if verbose:
            print()  # newline after progress
        plot.close()

        # --- Encode ---
        cmd = [
            "ffmpeg", "-y", "-loglevel", "error",
            "-framerate", str(fps),
            "-i", str(tmpdir / "frame_%05d.png"),
            # libx264 + yuv420p need even width/height
            "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
            "-c:v", "libx264", "-pix_fmt", "yuv420p",
            "-r", str(fps),
            filemovie,
        ]  # fmt: skip
        subprocess.run(cmd, check=True)
    finally:
        if cleanup is not None:
            cleanup.cleanup()

    if verbose:
        print(f"Done. Movie written to {filemovie}")
    return filemovie


def main(argv=None):
    args = parse_args(argv)
    try:
        render_movie1D(
            input=args.input,
            output=args.output,
            skip=args.skip,
            fps=args.fps,
            speed_factor=args.speed_factor,
            method=args.method,
            t_start=args.t_start,
            t_end=args.t_end,
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
            keep_frames=args.keep_frames,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as e:
        print(f"Error: ffmpeg failed (exit code {e.returncode})", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
