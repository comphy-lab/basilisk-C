"""
Plot a 3D volume (time, level, y, x) with its free surface (from eta) and
a colored interior, using xarray (to read the NetCDF) + PyVista (to render).

Assumes your NetCDF has, e.g.:
    T(time, level, y, x)   -- the variable to color by
    Z(time, level, y, x)   -- vertical coordinate of each level (depth/height,
                               terrain- or sigma-following -> varies with
                               level, y AND x)
    eta(time, y, x)         -- free-surface elevation (already in the file!)
    x(x), y(y)               -- horizontal coordinates (1D, regular or not)

Adjust variable/dimension names below to match the file.

conda env is visu3D
"""

import numpy as np
import xarray as xr
import pyvista as pv

# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)

from visu_3Dmovie import render_movie
from visu_3Dsnap import render_snapshot


# files
# myfile = "N512_nl30_0.000002_Tinizl/out.nc"
myfile = "~/basilisk/wiki/sandbox/hugoj/breaking_strat/ml_breaking_strat.gpu/out.nc"
# snapfile = "~/work/BasiLagrangian/models/N1024/out954.nc"
filemovie = "wave.mp4"
ttime = 50.0
snap = "wave"


if True:
    render_snapshot(
        myfile,
        "snapshot",
        50.0,
        var_side="T",
        clim_side=(19.98, 20.0),
        cmap_side="plasma",
    )
if False:
    render_movie(
        myfile,
        "wave",
        skip=1,
        speed_factor=2.0,
        var_side="T",
        clim_side=(19.98, 20.0),
        cmap_side="plasma",
        off_screen=True,
    )

# TODO:
# slice widget
