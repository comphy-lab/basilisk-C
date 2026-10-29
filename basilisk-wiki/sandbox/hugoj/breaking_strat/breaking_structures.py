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
 - see also hugoj/test_breaking_detection/test_breaking_shape.py

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


import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)

from breaking import (
    get_bins,
    label_periodic,
    btrack,
    simple_mapping_parallel,
    track_table,
    length1,
    bbarycenter,
    speed_direction,
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
# ds = ds.isel(time=slice(742, 745))
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
print("Step 2: track in time, birth, merge, death (time coherence)")

# we define here a structure for spatial coherence
structure = ndi.generate_binary_structure(2, 2)

# grid spacing (the tracker works on cell indices, so only the spacing matters)
dx = float(ds.x[1] - ds.x[0])
dy = float(ds.y[1] - ds.y[0])

# r_max: gap allowed between an object's edges in consecutive frames (metres)
r_max = 3 * max(dx, dy)  # starting guess, I need to tune it

# run the tracker
crest_mask_labelled, events = btrack(
    dmask.mask,
    dx,
    dy,
    structure,
    out_file="crest_mask_labelled.nc",
    event_out_file="breaking_events.json",
    r_max=r_max,
)

# we get infos for each object by 'inverting' the 'events' list
objects_infos = track_table(events, nt)

# ===========================================
print("Step 3: Compute the characteristics of a breaking")
print(" - length\n - speed\n - direction")

print(objects_infos)
print(crest_mask_labelled)

raise Exception

# add theses infos to 'object_infos'
for target in objects_infos.keys():
    istart = objects_infos[target]["birth"]
    iend = objects_infos[target]["death"]
    lengths = np.zeros((iend - istart, 2))
    speeds = np.zeros((iend - istart, 2))
    directions = np.zeros((iend - istart, 2))
    barys = np.zeros((iend - istart, 2))

    # print("id object", target)
    # print("object infos", objects_infos[target])
    # print("istart", istart)

    # this loop needs to be // later
    for it in range(istart, iend):
        snapshot = crest_mask_labelled.isel(time=it)
        obj_it = it - istart

        # step 3a: length
        lengths[obj_it, :] = length1(target, snapshot.values, dx, dy)

        # step 3b: speed and direction
        speeds[obj_it] = speed_direction(
            target, snapshot.values, ds["u.x"].isel(time=it), ds["u.y"].isel(time=it)
        )

        # step 3c: position
        barys[obj_it] = bbarycenter(target, snapshot)

    Lx, Ly = np.max(lengths, axis=0)
    avg_speed = np.mean(speeds, axis=0)
    print("object n° %d" % target)
    print(" max Lx,Ly = %0.1f, %0.1f" % (Lx, Ly))
    print(" avg_speed = %0.1f" % avg_speed[0])
    print(" avg_direction = %0.1f" % avg_speed[1])

raise Exception

# TEST =========================================
if True:
    print("Compute a length of a breaking")

    # Step 0: working setup
    target = 591  # the object
    itarget = 586  # the snapshot time index
    snapshot = crest_mask_labelled.isel(time=itarget)
    print("working at it=%d on object %d" % (itarget, target))

    sizes = length1(target, snapshot.values)

    print(f"Ny,Nx = {sizes}")
    print("Ly,Lx = %f m, %f m" % (sizes[0] * dy, sizes[1] * dx))

    # Test on a time serie for 1 breaker
    if True:
        istart = objects_infos[target]["birth"]
        iend = objects_infos[target]["death"]
        print("start = %d, end = %d" % (istart, iend))
        l_sizes = np.zeros((2, iend - istart))
        for it in range(istart, iend):
            sizes = length1(target, crest_mask_labelled.isel(time=it).values)
            l_sizes[0, it - istart] = sizes[0]
            l_sizes[1, it - istart] = sizes[1]

        lifetime = np.arange(istart, iend)
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), constrained_layout=True, dpi=200)
        ax.plot(lifetime, l_sizes[0, :] * dx, c="k", label=r"$L_y$")
        ax.plot(lifetime, l_sizes[1, :] * dy, c="b", label=r"$L_x$")
        ax.set_xlabel("timestep")
        ax.set_ylabel("size (m)")
        ax.legend()

    plt.show()


# TODO:
# - save the trajectory of barycenter of the breaking
# - compute characteristics of breaking (length, speed) and their evolution in the life of the object
#
# performance:
# - for each breaker, work on a small domain around it ?
#
# length of a snapshot:
# - count the number of pixel that have the id of the breaking ?
# - fit an ellipse and use its great axis ?
# - fit a polynomial ?

# speed:
# - mean of all pixel of the breaker ?

# 1 breaker = dict
# 'id'=int
# 'trajectory'= 3D numpy array
# 'size' = 1D numpy array (time)
# 'speed' = 1D numpy array

# THE length: max of the distribution of the breaker during its time evolution ?
# THE speed: same as for length

# to plot:
# - evolution in time of the size distribution (distribution of length on every snapshots, no time coherence needed)
# - same with speed

# - evolution of normalized length averaged over all breaker events

# TODO:
# - save a file with test of each function !


plt.show()
