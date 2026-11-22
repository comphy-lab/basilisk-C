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
    breaking_celerity_distribution,
    reducer,
    length_b,
)
from plot_tools import nice_slope

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


delta = (ds.x[1] - ds.x[0]).values
nt = len(ds.time)
bins = get_bins(kp)
bins_center = bins[1:] - (bins[2] - bins[1]) / 2

# =========================================
print("\nStep 1: detection of breaking fronts")

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
print("\nStep 2: track in time, birth, merge, death (time coherence)")

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
print("\nStep 3: Compute the characteristics of the breakings")
print(" - length\n - speed\n - direction")


file_name_breaking = "breaking_infos.npy"

# How to compute the length
# 0: fit a rectangle around the breaker
# 1: count the number of pixel
# 2: TBD, fit a segment
method_length = 0

objects_infos = breaker_characteristics(
    events,
    crest_mask_labelled,
    dssfx,
    method_length=method_length,
    output=file_name_breaking,
    verbose=verbose,
)

# TODO:
# - allow for a rotation of rectangle fit (method 0)
# - allow for parametric segment fit
#       this gives length, orientiation
# - allow to compute direction based on surface curvature
# - compute aspect ratio of breaker (is it defined for all length methods ?)
# - compute a non-dimensional ratio: breaker scale over curvature
# - allow to compute breaker velocity from position (not water current)

# ===========================================
print("\nStep 4: Compute the lambda(c) spectrum")


def slope_c(ax):
    return nice_slope(ax, -6, (0.36, 1), 0.12, 0.61, 0.041, 10, "k", "c")


reduce_method = {
    "speed": None,
    "length": None,
}
observed_S = L0**2
observed_T = nt

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
    " * Original method: mean total crest length (m, sum(lambda(c)dc)*A)= %f"
    % (np.sum(hist_m * np.diff(bins)) * L0**2)
)

# ============================
# Lambda(c) using object infos

Lavg, _, _ = breaking_celerity_distribution(
    objects_infos,
    bins,
    observed_S,
    observed_T,
    reduce_method=reduce_method,
    verbose=1,
    CHECK=True,
    crest_mask=dmask.mask,
)

# Plot
NORMED = True

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
fig.savefig("lambdaC.png")

# ===========================================
print("\nStep 5: Compute number of breakers with time")

length_index = 0
n_methods = 2

L_nb = np.zeros(nt)
Lb = np.zeros((nt, n_methods))
Lbstd = np.zeros((nt, n_methods))

for it in range(nt):
    snapshot = crest_mask_labelled.isel(time=it)
    unique, counts = np.unique(snapshot, return_counts=True)
    nobj = len(unique)

    Lbm = np.zeros(n_methods)
    # MEAN
    # loop over the object present in the current frame
    for k in range(nobj):
        target = unique[k]
        if target != 0:  # do not count background as object
            # loop over length methods
            for km in range(n_methods):
                Lbm[km] += length_b(target, snapshot, dx, dy, method=km)[length_index]
    Lbm /= nobj

    Lbs = np.zeros(n_methods)
    # VARIANCE
    for k in range(nobj):
        target = unique[k]
        if target != 0:  # do not count background as object
            # loop over length methods
            for km in range(n_methods):
                Lbs[km] += (
                    length_b(target, snapshot, dx, dy, method=km)[length_index]
                    - Lbm[km]
                ) ** 2
    Lbs /= nobj

    L_nb[it] = nobj
    Lb[it] = Lbm
    Lbstd[it] = Lbs


fig, ax = plt.subplots(1, 1, figsize=(5, 3), constrained_layout=True, dpi=200)
ax.plot(ds.time / Tp, L_nb, c="k")
ax.set_ylabel("number of breakers")
ax.set_xlabel("t/Tp")
fig.savefig("number_of_breaker_with_time.png")

fig, ax = plt.subplots(1, 1, figsize=(5, 3), constrained_layout=True, dpi=200)
ax.plot(ds.time / Tp, Lb[:, 0], c="k", label="method 0")
ax.plot(ds.time / Tp, Lb[:, 1], c="b", label="method 1")
ax.set_ylabel("<Lb> (m)")
ax.legend()
ax.set_xlabel("t/Tp")
fig.savefig("mean_length_breaker_with_time.png")


# ===========================================
print("\nStep 6: compare fluid velocity vs breaker velocity")

# -> for this I need a simulation with breaker that have lifetime > 1

# ===========================================
print("\nStep 7: Follow a breaker and plot vorticity")

# i) a tool to plot a breaker during its lifetime
# ii) get a unique grid for all breaker in their frame of reference
# iii) composite breaker

# ===========================================
print("\nStep 8: length, speed, direction during breaker lifetime")

# speed: both fluid based velocity AND dX/dt of each breaker object

# ===========================================
print("\nStep 9: distribution of lengths, velocity (=LambdaC) and direction")


# TODO:
# - a test for each function I use
# - what is the output frequency needed to follow breakers ?

plt.show()
