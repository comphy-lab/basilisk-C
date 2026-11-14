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

import json
import os.path
import sys
from joblib import Parallel, delayed
import traceback
import warnings

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)

from breaking import (
    get_bins,
    btrack,
    simple_mapping_parallel,
    simple_mapping,
    breaker_characteristics,
)

# path = "/home/jacqhugo/basilisk/wiki/sandbox/hugoj/breaking_strat/N512_40_RE40000_TiniZl_remap1.0_cuda_smooth/"
path = "/home/jacqhugo/basilisk/wiki/sandbox/hugoj/test_windinput/test_real_field/PM_seastate_hold.cuda_dev/"
myfile = "out.nc"
seltime = slice(-4, -1, 1)
L0 = 200
kp = 5 * 2 * np.pi / L0
g = 9.81
omegap = np.sqrt(g * kp)
cp = omegap / kp
Tp = 2.0 * np.pi / np.sqrt(9.81 * kp)
lvl_surf = -1
threshold = -3 * kp
verbose = True

# =========================================
print("Step 0: get the data")
# chunking along 'time' but no 'x' or 'y' as we compute XY diags and we want to // on time
ds = xr.open_dataset(path + myfile, chunks={"time": 1, "x": -1, "y": -1})

# ds = ds.isel(time=slice(0, 5))  # work on subset for now
dssfx = ds.isel(level=-1)

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


NT_objects = np.amax(crest_mask_labelled).values
print(f"\tTotal number of objects = {NT_objects}")

# ===========================================
print("Step 3: Compute the characteristics of the breakings")
print(" - length\n - speed\n - direction")


file_name_breaking = "breaking_infos.npy"
method_length = 0


objects_infos = breaker_characteristics(
    events,
    crest_mask_labelled,
    dssfx,
    method_length=method_length,
    output=file_name_breaking,
    verbose=verbose,
)


# ===========================================
if True:
    print("Step 4: Compute the lambda(c) spectrum")

    NORMED = True

    def build_array(objects, length_method="mean", speed_method="mean", min_lifetime=0):
        """Return a dict of object with key=id, then dict L, S, D.

        length_method, speed_method: a key of REDUCERS, or any callable taking
            an (n, 2) array and returning a sequence of two values.
        min_lifetime: only objects with lifetime > min_lifetime are kept.
        """
        REDUCERS = {
            "first": lambda a: a[0],
            "last": lambda a: a[-1],
            "mean": lambda a: a.mean(axis=0),
            "median": lambda a: np.median(a, axis=0),
            "max": lambda a: a.max(axis=0),
            "min": lambda a: a.min(axis=0),
            "None": lambda a: a,
        }

        def get(m):
            return REDUCERS[m] if isinstance(m, str) else m

        reduce_length = get(length_method)
        reduce_speed = get(speed_method)

        rows = {}
        for id, info in objects.items():
            rows[id] = {}
            cond = (
                (info["lifetime"] > min_lifetime)
                # & (info["death"] is not None)
                & (("length") in info.keys())
            )
            if cond:
                rows[id]["L"] = reduce_length(np.asarray(info["length"]))
                rows[id]["S"] = reduce_speed(info["speeds_currents"])
                rows[id]["lifetime"] = info["lifetime"]

        return rows

    def nice_slope(ax, slope, X, Y0, Xt, Yt, fs, c, var):
        """
        Plot a slope on ax.
        X (tuple): move along x axis
        Y0: move vertically.
        (Xt,Yt): the position of the label, in data coordinates.
        'fs': fontsize
        'c': color of line and text
        """
        x0, x1 = X
        Y1 = Y0 * (X[1] / X[0]) ** slope
        ax.loglog([X[0], X[1]], [Y0, Y1], "--", lw=1, c=c)
        ax.text(Xt, Yt, rf"${var}^{slope}$", fontsize=fs, va="center", c=c)

    def slope_c(ax):
        return nice_slope(ax, -6, (0.36, 1), 0.27, 0.61, 0.041, 10, "k", "c")

    bdata = build_array(objects_infos, speed_method="None")

    bins = get_bins(kp)
    bins_center = bins[1:] - (bins[2] - bins[1]) / 2

    # ===============================
    # Lambda(c) the regular way
    #   here we do again the detection. If hist was saved we would reuse it !
    hists = np.zeros((nt, len(bins_center)))
    for it in range(nt):
        snapshot = dssfx.isel(time=it)

        hist, _ = simple_mapping(
            snapshot.eta,
            snapshot["u.x"],
            snapshot["u.y"],
            delta,
            bins=bins,
            method=0,
            threshold=threshold,
        )
        hists[it, :] = hist * delta / L0**2 / np.diff(bins)
    hist_m = np.mean(hists, axis=0)  # avg over all timesteps

    print(
        "Method 1: mean total crest length (m, sum(lambda(c)dc)*A)= %f"
        % (np.sum(hist_m * np.diff(bins)) * L0**2)
    )

    # ============================
    # Lambda(c) using object infos
    Lavg = np.zeros(len(bins_center))
    reduce_method_speed = None  # None, "mean"
    observed_S = L0**2
    observed_T = nt
    outside = 0
    speed_index = 0  # 0: norm, 1: direction
    length_index = 0  # 0: Y

    total_L = 0.0
    old_b = 0

    for io in range(1, len(bdata.keys()) + 1):
        if "L" in bdata[io].keys():
            length = bdata[io]["L"][length_index]
            total_L += length
            a_speed = bdata[io]["S"]
            lifetime = bdata[io]["lifetime"]
            if lifetime > 1:
                old_b += 1
            # loop over lifetime of object
            for it in range(lifetime):
                speed = a_speed[it][speed_index]
                ncells = len(np.atleast_1d(speed))
                if reduce_method_speed == "mean":
                    speed = np.atleast_1d(speed.mean())

                # loop over breaking speeds
                for ispeed in range(ncells):
                    # loop over the bins to add the length of this breaking to corresponding velocity
                    for ib in range(0, len(bins_center)):
                        if speed[ispeed] >= bins[ib] and speed[ispeed] < bins[ib + 1]:
                            Lavg[ib] += length / ncells / lifetime
                            break
                        if speed[ispeed] < bins[0] or speed[ispeed] > bins[-1]:
                            outside += 1
            if (outside > 1) & verbose:
                print(
                    f"WARNING: object {io} has {outside} pixels with velocity out of bins range [{bins[0]},{bins[-1]}"
                )
                outside = 0
    print("\tnumber of old (lifetime>1) breakers = %d" % old_b)
    # Normalize
    Lavg /= observed_S * observed_T * np.diff(bins)

    # for ik in range(len(bins_center)):
    #     dc = bins[ik + 1] - bins[ik]
    #     Lavg[ik] /= observed_S * observed_T * dc
    #     # Lavg[ik] /= observed_T * dc
    total_L /= nt

    print(
        "Method 2: mean total crest length (m)= %f"
        % (np.sum(Lavg * np.diff(bins)) * L0**2)
    )
    print("Method 2: mean total crest length (m, objects) = %f" % total_L)

    print("N_pixel_breakers*delta/nt=%.1f" % (dmask.mask.sum().values * delta / nt))
    mean_ratio = 0
    nobj = len(objects_infos.keys())
    counter = 0

    for io in range(1, nobj + 1):  # len(bdata.keys()) + 1)
        if "lifetime" in bdata[io].keys():
            No = objects_infos[io]["Npxl"][0]
            Lo = bdata[io]["L"][length_index]
            r = No * delta / Lo
            if False:
                print(
                    "object %d, No*delta=%.1f, No=%d, Lo=%.1f (ratio=%.1f)"
                    % (io, No * delta, No, Lo, r)
                )
            mean_ratio += r
            counter += 1
    mean_ratio /= counter
    print("mean ratio = %.3f" % mean_ratio)

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(5, 3), constrained_layout=True, dpi=200)
    if NORMED:
        ax.plot(bins_center / cp, hist_m * cp**3 / g, c="k", label="regular")
        ax.plot(bins_center / cp, Lavg * cp**3 / g, c="b", label="object")
        ax.set_xlabel(r"$c_b/c_p$")
        ax.set_ylabel(r"$\Lambda(c) c_p^3 g^{-1}$")
    else:
        ax.plot(bins_center / cp, hist_m, c="k", label="regular")
        ax.plot(bins_center / cp, Lavg, c="b", label="object")
        ax.set_xlabel(r"$c_b/c_p$")
        ax.set_ylabel(r"$\Lambda(c)$")
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim((0.1, 1.1))
    slope_c(ax)  # nice slope c**-6


# ===========================================
print("Step 5: Compute number of breakers with time")


L_nb = []
for it in range(nt):
    n = len(np.unique(crest_mask_labelled.isel(time=it)))
    L_nb.append(n)


fig, ax = plt.subplots(1, 1, figsize=(5, 3), constrained_layout=True, dpi=200)
ax.plot(ds.time / Tp, L_nb)
ax.set_ylabel("number of breakers")
ax.set_xlabel("t/Tp")
fig.savefig("number_of_breaker_with_time.png")


# TODO:
# [x] save the trajectory of barycenter of the breaking
# [x] compute characteristics of breaking (length, speed) and their evolution in the life of the object
#
# performance:
# - for each breaker, work on a small domain around it ?
# - compute characteristic of breakers in //
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
