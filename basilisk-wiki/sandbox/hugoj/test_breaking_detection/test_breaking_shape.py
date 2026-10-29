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
    label_periodic,
    btrack,
    simple_mapping_parallel,
    track_table,
    length1,
    bbarycenter,
    speed_direction,
)


path = "./"
myfile = "surface_t589-591.nc"
L0 = 200
kp = 5 * 2 * np.pi / L0
lvl_surf = -1
threshold = -2 * kp

# =========================================
print("\n=> get the data")
# chunking along 'time' but no 'x' or 'y' as we compute XY diags and we want to // on time
ds = xr.open_dataset(path + myfile, chunks={"time": 1, "x": -1, "y": -1})

delta = (ds.x[1] - ds.x[0]).values
nt = len(ds.time)
bins = get_bins(kp)
bins_center = bins[1:] - (bins[2] - bins[1]) / 2
dx = delta
dy = delta
# r_max: gap allowed between an object's edges in consecutive frames (metres)
r_max = 3 * max(dx, dy)  # starting guess, I need to tune it

# =========================================
print("\n=> detection of breaking fronts")

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
print("\n=> label breakings (no time coherence)")
structure = ndi.generate_binary_structure(2, 2)

fig, ax = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True, dpi=200)
for it in range(nt):
    # load only the slice you need (lazy with xarray/dask)
    u = ds["u.x"].isel(time=it, level=lvl_surf)
    s = ax[it].pcolormesh(u.x, u.y, u, cmap="Greys_r", vmin=-4, vmax=4, rasterized=True)

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
fig.savefig("label_no_time_coherence.png")

# =========================================
print("\n=> track in time, birth, merge, death (time coherence)")
# Structure of the search around points
# A start is: ndi.generate_binary_structure(2, 2)
#
# which is equivalent to a kernel like
#
#   [[1,1,1],
#    [1,1,1],
#    [1,1,1]]
structure = ndi.generate_binary_structure(2, 2)


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

# TEST =========================================
print("\n=> Compute a length of a breaking")

# Step 0: working setup
target = 7  # the object
itarget = 1  # the snapshot time index
snapshot = crest_mask_labelled.isel(time=itarget)
print("working at it=%d on object %d" % (itarget, target))

sizes = length1(target, snapshot.values, dx, dy)

print(f"Nx,Ny = {sizes}")
print("Lx,Ly = %f m, %f m" % (sizes[0], sizes[1]))

# Test on a time serie for 1 breaker
if True:
    istart = 0
    iend = 3
    print("start = %d, end = %d" % (istart, iend))
    l_sizes = np.zeros((2, iend - istart))
    for it in range(istart, iend):
        sizes = length1(target, crest_mask_labelled.isel(time=it).values, dx, dy)
        l_sizes[0, it - istart] = sizes[0]
        l_sizes[1, it - istart] = sizes[1]

    lifetime = np.arange(istart, iend)
    fig, ax = plt.subplots(1, 1, figsize=(5, 3), constrained_layout=True, dpi=200)
    ax.plot(lifetime, l_sizes[0, :] * dx, c="k", label=r"$L_x$")
    ax.plot(lifetime, l_sizes[1, :] * dy, c="b", label=r"$L_y$")
    ax.set_xlabel("timestep")
    ax.set_ylabel("size (m)")
    ax.legend()
    fig.savefig("Length_evolution.png")

# TEST =========================================
print("\n => Compute barycenter of a breaking")

# 1 get the (X,Y) values for each point of a breaking on a snapshot
target = 7
itarget = 1  # int((istart + iend) / 2)  # working at mid age
print("working at it=%d" % itarget)
snapshot = crest_mask_labelled.isel(time=itarget)

# compute barycenter
bposition = bbarycenter(target, snapshot)
print(f"position of the barycenter for it=%d: {bposition}" % itarget)

# TEST =======================================
print("\n=> Compute speed and direction of the breaking")

target = 7
itarget = 1  # int((istart + iend) / 2)  # working at mid age
print("working at it=%d" % itarget)
snapshot = crest_mask_labelled.isel(time=itarget)

speed, direction = speed_direction(
    target, snapshot, ds["u.x"].isel(time=itarget), ds["u.y"].isel(time=itarget)
)

print(
    "speed = %.1f m/s, direction = %.1f° (>0 counterclockwise from x axis)"
    % (speed, direction)
)

# TEST =======================================
print("\n => Breaking-centered vue")
# Step 0: working setup
target = 7  # the object
itarget = 1  # the snapshot time index
snapshot = crest_mask_labelled.isel(time=itarget)
print("working at it=%d on object %d" % (itarget, target))

# Step 1: compute the size of the breaking
sizes = length1(target, snapshot.values, dx, dy)

# --- we now use this breaking length to get a adimensionned vue of the breaking
Lb = sizes.max()  # max between Lx and Ly
print("This breaking is of size Lb = %f m" % Lb)

# Step 2: compute barycenter
bposition = bbarycenter(target, snapshot)

# Step 3: extract values around +-0.5Lb in X,Y
offset = 0.5 * Lb
ds_local = ds.sel(
    x=slice(bposition[0] - offset, bposition[0] + offset),
    y=slice(bposition[1] - offset, bposition[1] + offset),
).isel(time=itarget, level=-1)
mask_local = snapshot.sel(
    x=slice(bposition[0] - offset, bposition[0] + offset),
    y=slice(bposition[1] - offset, bposition[1] + offset),
)

nX = (ds_local.x - bposition[0]) / Lb
nY = (ds_local.y - bposition[1]) / Lb
iy, ix = np.where(mask_local == target)
fig, ax = plt.subplots(1, 1, figsize=(4, 3), constrained_layout=True, dpi=200)
s = ax.pcolormesh(nX, nY, ds_local["w"], cmap="seismic", vmin=-1.5, vmax=1.5)
ax.scatter(nX[ix], nY[iy], marker="x", color="k", s=10, linewidth=0.5)
plt.colorbar(s, ax=ax, label="w (m/s)")
ax.set_xlabel("x/Lb")
ax.set_ylabel("y/Lb")
ax.set_aspect(1)
ax.set_title("breaking n°%d, it=%d" % (target, itarget))
fig.savefig("breaking_centered_vue.png")
plt.show()
