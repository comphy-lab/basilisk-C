"""
tools_profile — helpers shared by visu_profile_snap.py and visu_profile_movie.py.

Matplotlib only. Plots a 1D profile (e.g. eta(x)) from an
xarray NetCDF output.
"""

from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")  # never opens a window: safe on clusters / over ssh
import matplotlib.pyplot as plt  # noqa: E402


# --------------------------------------------------------------------------
# Command-line arguments common to the snapshot and the movie scripts
# --------------------------------------------------------------------------
def add_plot_arguments(p):
    p.add_argument(
        "--method",
        choices=["nearest", "interp"],
        default="nearest",
        help="How to sample simulation time onto the requested time(s)",
    )
    p.add_argument("--var", type=str, default="eta", help="Variable to plot")
    p.add_argument(
        "--y-index",
        type=int,
        default=0,
        help="Index along y used to extract the 1D profile (ignored if no y dim)",
    )
    p.add_argument(
        "--L0",
        type=float,
        default=None,
        help="Horizontal domain size (m). x-limits become [-L0/2, L0/2]",
    )
    p.add_argument(
        "--T0",
        type=float,
        default=None,
        help="Reference period. If given, title shows t/T0 instead of t",
    )
    p.add_argument(
        "--xlim",
        type=float,
        nargs=2,
        default=None,
        metavar=("XMIN", "XMAX"),
        help="x-limits. Overrides --L0. Default: --L0 if given, else data range",
    )
    p.add_argument(
        "--ylim",
        type=float,
        nargs=2,
        default=None,
        metavar=("YMIN", "YMAX"),
        help="y-limits. Default: min/max of the variable over the whole file",
    )
    p.add_argument("--xlabel", type=str, default="x (m)", help="x-axis label")
    p.add_argument("--ylabel", type=str, default="z (m)", help="y-axis label")
    p.add_argument("--color", type=str, default="C0", help="Line colour")
    p.add_argument("--linewidth", type=float, default=1.5, help="Line width")
    p.add_argument("--equal-aspect", default=True, help="Use the same scale on x and y")
    p.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        default=[6.0, 3.0],
        metavar=("W", "H"),
        help="Figure size in inches",
    )
    p.add_argument("--dpi", type=int, default=100, help="Figure resolution")
    return p


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
def open_profile_dataset(path, skip=0, y_index=0):
    """Open the file, drop duplicated times, skip first steps, pick one y."""
    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"input file not found: {path}")
    ds = xr.open_dataset(path)
    if "time" in ds.dims:
        ds = ds.isel(time=~ds.indexes["time"].duplicated(keep="first"))
        if skip:
            ds = ds.isel(time=slice(skip, None))
    if y_index is not None and "y" in ds.dims:
        ds = ds.isel(y=y_index)
    return ds


def select_time(da, t, method="nearest"):
    """Sample a DataArray at time t (nearest neighbour or linear interpolation)."""
    if "time" not in da.dims:
        return da
    if method == "interp":
        return da.interp(time=t)
    return da.sel(time=t, method="nearest")


def auto_ylim(da, margin=0.05):
    """Fixed y-limits from the whole record, so the axis doesn't jump in a movie."""
    lo, hi = float(da.min()), float(da.max())
    if lo == hi:
        lo, hi = lo - 1.0, hi + 1.0
    pad = margin * (hi - lo)
    return lo - pad, hi + pad


def resolve_limits(ds, var, xlim=None, ylim=None, L0=None):
    if var not in ds:
        raise ValueError(f"variable '{var}' not found. Available: {list(ds.data_vars)}")
    if xlim is None:
        xlim = (
            (-L0 / 2, L0 / 2)
            if L0 is not None
            else (
                float(ds.x.min()),
                float(ds.x.max()),
            )
        )
    if ylim is None:
        ylim = auto_ylim(ds[var])
    return tuple(xlim), tuple(ylim)


def make_title(t, T0=None):
    return f"t/T0 = {t / T0:.1f}" if T0 else f"t = {t:.1f} s"


# --------------------------------------------------------------------------
# Plot: one figure created once, then only the line + title are updated
# --------------------------------------------------------------------------
class ProfilePlot:
    def __init__(
        self,
        x,
        xlim,
        ylim,
        xlabel="x (m)",
        ylabel="z (m)",
        figsize=(6, 3),
        dpi=100,
        color="C0",
        linewidth=1.5,
        equal_aspect=True,
    ):
        self.fig, self.ax = plt.subplots(
            figsize=tuple(figsize), constrained_layout=True, dpi=dpi
        )
        self.x = np.asarray(x)
        (self.line,) = self.ax.plot(
            self.x, np.full_like(self.x, np.nan, dtype=float), color=color, lw=linewidth
        )
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        if equal_aspect:
            self.ax.set_aspect(1)
        self.title = self.ax.set_title("")

    def update(self, values, title):
        self.line.set_ydata(np.asarray(values))
        self.title.set_text(title)

    def save(self, path, **kwargs):
        self.fig.savefig(path, **kwargs)

    def close(self):
        plt.close(self.fig)
