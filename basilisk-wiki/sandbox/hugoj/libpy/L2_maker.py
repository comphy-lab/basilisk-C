import numpy as np
import xarray as xr
import xgcm
import gc
from pathlib import Path

from diags import grad_velocities, interpz
from data_reader import build_grid


def check_var_in(ds, list_var):
    # fast check
    for var in list_var:
        if var not in ds.keys():
            raise Exception(f"{var} is not in the dataset you provided ! Aborting ...")


def common_diags(
    ds: xr.Dataset, var_vec: list, var_vec2: list, u_vec: list, u_vec2: list
):

    # mean
    for i, var in enumerate(var_vec):
        ds[var_vec2[i] + "_m"] = ds[var].mean(dim=["x", "y"])

    # fluctuations
    for i, var in enumerate(var_vec):
        ds[var_vec2[i] + "_f"] = ds[var] - ds[var_vec2[i] + "_m"]

    # mke
    mke = ds[var_vec2[0] + "_m"] * 0.0
    for i, var in enumerate(u_vec):
        mke += ds[u_vec2[i] + "_m"] ** 2
    ds["mke"] = mke

    # flux
    for uj in var_vec2:
        for ui in var_vec2:
            ds[uj + ui] = ds[uj + "_f"] * ds[ui + "_f"]

    # tke
    tke = ds[u_vec[0]] * 0.0
    for var in u_vec2:
        tke += ds[var + "_f"] ** 2
    ds["tke"] = 0.5 * tke

    del (mke, tke)
    gc.collect()
    return ds


def make_L2_layer(
    ds: xr.Dataset,
    grid: xgcm.Grid,
    outpath: str,
    outfile: str,
    save=True,
    verb=False,
):
    """
    Builds a xarray dataset and saves on disk. The average operator is the layer average
    (z = z(t,x,y)). This dataset contains mean profiles, correlations of fluctuations,
    tke, mke and gradients.

    ds and grid can be obtained using 'read_bas_data' from 'data_reader.py'

    INPUTS:
        - ds: the original dataset
        - grid: a xgcm grid
        - outpath: where to save the file
        - outfile: how to name the saved file
        - save: whether to trigger the computation by saving or not
    OUPUT:
        - dsL2: the dataset with the new diagnostics
    """
    input_path = Path(outpath + outfile).expanduser()
    if not input_path.exists():
        if verb:
            print("=====================================================")
            print("Original Dataset: \n")
            print(ds)
            print("\n =====================================================")
            print("-> Building L2 file !\n")

        check_var_in(ds, ["u.x", "u.y", "u.z", "T"])

        # let's gather names of variables in the initial dataset
        list_var_ini = list(ds.keys())

        # _____________________
        # compute diags

        u_vec = ["u.x", "u.y", "u.z"]
        u_vec2 = ["u", "v", "w"]
        var_vec = u_vec + ["T"]
        var_vec2 = u_vec2 + ["T"]

        ds = common_diags(ds, var_vec, var_vec2, u_vec, u_vec2)

        # Divergence at surface
        # ds, _ = grad_velocities(ds, grid, compute=False)

        # make the dataset lighter by removing already saved on disk variables
        dsL2 = ds.drop_vars(list_var_ini).drop_vars("time")

        # save also z
        dsL2["z"] = ds.z

        # ______________________
        # save file
        if save:
            print("saving ...")
            print("dataset saved:\n")
            print(dsL2)
            dsL2.to_netcdf(outpath + outfile)
            print("\ndone !")
            dsL2.close()
            dsL2 = xr.open_dataset(outpath + outfile)
        else:
            print("saving skipped !")
    else:
        dsL2 = xr.open_dataset(outpath + outfile)

    return dsL2


def make_L2_eulerian(
    ds: xr.Dataset,
    grid: xgcm.Grid,
    znew: np.ndarray,
    outpath: str,
    outfile: str,
    fill_value=0.0,
    save=True,
    verb=False,
):
    """
    Builds a xarray dataset and saves on disk. The average operator is the horizontal average
    (z independant of x,y). This dataset contains mean profiles, correlations of fluctuations,
    tke, mke and gradients.

    ds and grid can be obtained using 'read_bas_data' from 'data_reader.py'

    INPUTS:
        - ds: the original dataset
        - grid: a xgcm grid
        - outpath: where to save the file
        - outfile: how to name the saved file
        - save: whether to trigger the computation by saving or not
    OUPUT:
        - dsL2: dataset with newdiags
    """
    input_path = Path(outpath + outfile).expanduser()
    if not input_path.exists():
        if verb:
            print("=====================================================")
            print("Original Dataset: \n")
            print(ds)
            print("\n =====================================================")
            print("-> Building L2 file !\n")
        check_var_in(ds, ["u.x", "u.y", "u.z", "T"])

        # _____________________________________
        # step 1: interpolate on cartesian grid
        l_to_interp = ["u.x", "u.y", "u.z", "T"]

        for var in l_to_interp:
            ds[var + "i"] = interpz(ds.z, ds[var], znew, fill_value=fill_value)

        ds = ds.assign_coords(
            {
                "znew": (["znew"], znew),
                "znew_l": (["znew_l"], znew - 0.5 * (znew[1] - znew[0])),
            }
        )

        # vertical coordinate is now znew so rebuilding grid
        # gridnew = build_grid(ds, zvar="znew")

        # _____________________
        # step 2: compute diags

        u_vec = ["u.xi", "u.yi", "u.zi"]
        u_vec2 = ["u", "v", "w"]
        var_vec = u_vec + ["Ti"]
        var_vec2 = u_vec2 + ["T"]

        ds = common_diags(ds, var_vec, var_vec2, u_vec, u_vec2)

        # Divergence at surface
        # ds, _ = grad_velocities(ds, gridnew, uname=u_vec, zvar="znew", compute=False)

        # _________________________
        # step 3: cleanup, renaming
        l_to_drop = ["u.x", "u.y", "u.z", "z", "z_l", "h", "T", "level", "level_l"]
        dsL2 = ds.drop_vars(l_to_drop)  # remove var from initial dataset
        dsL2 = dsL2.transpose(
            "znew", "y", "x", "znew_l", "y_l", "x_l", ...
        )  # reorder as (z,y,x)
        dsL2.attrs["Note"] = (
            f"z coord: interpolated from semi-Lagrangian to cartesian z, with fill_value={fill_value}"
        )

        dsL2 = dsL2.rename({"u.xi": "u.x", "u.yi": "u.y", "u.zi": "u.z", "Ti": "T"})
        # ______________________
        # save file
        if save:
            print("saving ...")
            print("dataset saved:\n")
            print(dsL2)
            dsL2.to_netcdf(outpath + outfile)
            print("\ndone !")
            dsL2.close()
            dsL2 = xr.open_dataset(
                outpath + outfile,
                chunks={"znew": -1, "x": ds.chunks["x"], "y": ds.chunks["y"]},
            )
        else:
            print("saving skipped !")
    else:
        print(f"I found a L2 file here: {outpath + outfile}\n")
        dsL2 = xr.open_dataset(
            outpath + outfile,
            chunks={"znew": -1, "x": ds.chunks["x"], "y": ds.chunks["y"]},
        )

    return dsL2
