import xarray as xr
import numpy as np
from xgcm import Grid
from tools import build_z_ds


# read netcdf, set up for xgcm
def read_bas_data(filename, chunks="auto", dtype="float32"):
    """
    This function reads the output of Basilisk (using bderembl/libs/netcdf_pas.h)
        and return a dataset with xgcm coordinates to allow for easy differentiation

    INPUT:
        filename: a string where the .nc is
    OUPUT:
        dataset: a xr.Dataset
    """
    ds = xr.open_dataset(filename, chunks=chunks)

    # building z and left side of z
    ds = build_z_ds(ds, dtype)

    # adding left side coordinates
    #   and z and z_l
    dx = (ds.x[1] - ds.x[0]).values
    dy = (ds.y[1] - ds.y[0]).values
    coords = {
        "x_l": (["x_l"], ds.x.data - dx / 2),
        "y_l": (["y_l"], ds.y.data - dy / 2),
        # "zl":(["zl"], ds.level.data), # doing nothing but its cleaner
        "level_l": (["level_l"], ds.level.data - 0.5),
        # "z":(['z'], ds.z.data),
        # "z_l":(['z_l'], ds.z_l.data),
    }
    ds = ds.assign_coords(coords)
    # ds = ds.rename({'u.z':'w'}) # this could be set in bderemble/netcdf_bas.h
    ds = ds.rename({"w": "u.z"})

    # xgcm grid
    grid = build_grid(ds)

    # fix: sometimes Basilisk simulation has 2 times the initial timestamp.
    #       This removes the first one so that .sel() works properly
    if ds.time[0] == ds.time[1]:
        ds = ds.isel(time=slice(1, len(ds.time)))

    return ds, grid


def read_data(filename, chunks="auto", dtype="float32"):
    return read_bas_data(filename, chunks, dtype)


def build_grid(ds, zvar="level"):
    # xgcm grid
    grid = Grid(
        ds,
        coords={
            "X": {"center": "x", "left": "x_l"},
            "Y": {"center": "y", "left": "y_l"},
            "Z": {"center": zvar, "left": zvar + "_l"},
        },
        autoparse_metadata=False,
        # periodic={'X':'True','Y':'True','Z':'False'},
        padding={"X": "periodic", "Y": "periodic", "Z": "fill"},
        fill_value={"Z": 0},
    )
    return grid
