"""
breaking_structure.py

This script aims to detect turbulent structures associated with breaking events

METHOD:
1) detection of the breaking fronts using Wu et al. (2023) python function
2) group point detected as breaking into objects (both space and time)
3) diagnose breaking characteristics: dimensions (both space and time)
    -> do a composite of detected breakings
4) Use previous point to analyse a volume of water under the breaking
5) Expand diagnostics to higher ordere quantities:
    Q criterion, T structures ?
6) Gather results
    -> do a composite of breaking + underneath turbulence

HYPOTHESIS:
- breaking celerity = u.x of the surface layer
- detection of a neighbour is in a 3x3 area around a point


NOTES:


TODO:
- how to decide neighbour detection range ?
- this whole script needs to be efficient: to be used on many snapshots
- output frequency of the model ?

REFERENCES:

Wu, J., Popinet, S., & Deike, L. (2023). Breaking wave field statistics with
a multi-layer model. Journal of Fluid Mechanics, 968, A12.
https://doi.org/10.1017/jfm.2023.522


"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import ndimage as ndi


# pytest test_label_periodic.py -v add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)

from breaking import (
    get_bins,
    simple_mapping,
    label_periodic,
    btrack,
    simple_mapping_parallel,
)


path = "/home/jacqhugo/basilisk/wiki/sandbox/hugoj/breaking_strat/N512_40_RE40000_TiniZl_remap1.0_cuda_smooth/"
myfile = "out.nc"
seltime = slice(-4, -1, 1)
L0 = 200
kp = 5 * 2 * np.pi / L0
lvl_surf = -1
threshold = -2 * kp
hchunks = 256

# =========================================
print("Step 0: get the data")
# chunking along 'time' but no 'x' or 'y' as we compute XY diags and we want to // on time
ds = xr.open_dataset(path + myfile, chunks={"time": 1, "x": -1, "y": -1})

# ds = ds.isel(time=slice(4, 745))  # t/tp = [2,5] -> ~115s to be completed

# ds = ds.isel(time=seltime)  # only a subset for now to test the method
ds = ds.isel(time=slice(742, 745))
# ds = ds.isel(time=slice(150, 200))

delta = (ds.x[1] - ds.x[0]).values
nt = len(ds.time)
bins = get_bins(kp)
bins_center = bins[1:] - (bins[2] - bins[1]) / 2

DOPLOTS = False
if nt <= 3:
    DOPLOTS = True

# =========================================
print("Step 1: detection of breaking fronts")

dmask = simple_mapping_parallel(
    ds.eta,
    ds["u.x"],
    ds["u.y"],
    lvl_surf=lvl_surf,
    delta=delta,
    bins=bins,
    threshold=threshold,
)

# =========================================
if DOPLOTS:
    print("a first Step 2: label breakings (no time coherence)")
    structure = ndi.generate_binary_structure(2, 2)

    fig, ax = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True, dpi=200)
    for it in range(nt):
        # load only the slice you need (lazy with xarray/dask)
        u = ds["u.x"].isel(time=it, level=lvl_surf)
        s = ax[it].pcolormesh(
            u.x, u.y, u, cmap="Greys_r", vmin=-4, vmax=4, rasterized=True
        )

        # labels
        points = np.logical_not(dmask.mask[it])
        labels, num_features = label_periodic(points, structure)

        ax[it].pcolormesh(
            u.x,
            u.y,
            np.ma.masked_equal(labels, 0),
            cmap="tab20",
            vmin=0,
            vmax=max(num_features, 1),
            rasterized=True,
        )
        ax[it].set_aspect(1)
        ax[it].set_xlabel("x (m)")
        ax[it].set_title("t=%.02f" % u.time.values)
        del points, labels  # free before the next iteration

    plt.colorbar(s, ax=ax[-1])
    ax[0].set_ylabel("y (m)")
    fig.suptitle("Label: no time coherence")

# =========================================
print("Step 2: track in time, birth, merge, death (time coherence)")
structure = ndi.generate_binary_structure(2, 2)

# grid spacing (the tracker works on cell indices, so only the spacing matters)
dx = float(ds.x[1] - ds.x[0])
dy = float(ds.y[1] - ds.y[0])

# r_max: gap allowed between an object's edges in consecutive frames (metres)
r_max = 3 * max(dx, dy)  # starting guess, I need to tune it

# run the tracker; frames are produced lazily, one at a time
events = []
ids = btrack(
    dmask.mask,
    dx,
    dy,
    structure,
    out_file="crest_mask_labelled.nc",
    r_max=r_max,
    events=events,
)


# =========================================
if DOPLOTS:
    print("Step 3: plot")
    fig, ax = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True, dpi=200)
    for a, it in zip(ax, range(nt)):
        u = ds["u.x"].isel(time=it, level=lvl_surf)
        s = a.pcolormesh(u.x, u.y, u, cmap="Greys_r", vmin=-4, vmax=4, rasterized=True)

        tid = np.ma.masked_equal(ids[it], 0)
        a.pcolormesh(
            ds.x, ds.y, (tid - 1) % 20, cmap="tab20", vmin=0, vmax=19, rasterized=True
        )  # same id -> same colour across frames
        a.set_aspect(1)
        a.set_xlabel("x (m)")
        a.set_xlim([0, 20])
        a.set_ylim([35, 55])
        a.set_title("t=%.02f" % u.time.values)
    plt.colorbar(s, ax=ax[-1])
    ax[0].set_ylabel("y (m)")
    fig.suptitle("Label: with time coherence")


# fig, ax = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True, dpi=200)
# for it in range(nt):
#     # load only the slice you need (lazy with xarray/dask)
#     u = ds["u.x"].isel(time=it, level=lvl_surf)
#     s = ax[it].pcolormesh(u.x, u.y, u, cmap="Greys_r", vmin=-1, vmax=2, rasterized=True)
#
#     # labels
#     points = ~l_mask[it]
#     labels, num_features = label_periodic(points, structure)
#
#     ax[it].pcolormesh(
#         u.x,
#         u.y,
#         np.ma.masked_equal(labels, 0),
#         cmap="tab20",
#         vmin=0,
#         vmax=max(num_features, 1),
#         rasterized=True,
#     )
#     ax[it].set_aspect(1)
#     ax[it].set_xlabel("x (m)")
#     ax[it].set_title("t=%.02f" % u.time.values)
#     del points, labels  # free before the next iteration
# plt.colorbar(s, ax=ax[-1])
# ax[0].set_ylabel("y (m)")

plt.show()
