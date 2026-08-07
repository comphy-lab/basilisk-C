import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import gc

# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)

from breaking import get_bins, simple_mapping

path = "N512_nl30_0.000002_Tinizl/"
myfile = "out.nc"

path = "~/work/BasiLagrangian/models/N1024/"
myfile = "out954.nc"

seltime = -1
L0 = 200
step = 1
lvl_surf = -1
ds = xr.open_dataset(path + myfile, chunks={"t": -1, "x": -1, "y": -1})
# ds.isel(x=slice(0, 512, step), y=slice(0, 512, step))
ds = ds.isel(time=seltime)
#### Some metadata ####
delta = (ds.x[1] - ds.x[0]).values
kp = 10 * np.pi / L0
threshold = -2 * kp
bins = get_bins(kp)
bins_center = bins[1:] - (bins[2] - bins[1]) / 2
#### Compute stats ####
print("Computing stats for %s..." % (path + myfile))

hist, mask = simple_mapping(
    ds.eta,
    ds["u.x"].isel(level=lvl_surf),
    ds["u.y"].isel(level=lvl_surf),
    delta,
    bins=bins,
    method=0,
    threshold=threshold,
    EXTRA_FILTER=True,
    return_mask=True,
)
fig, ax = plt.subplots(1, 1, figsize=(10, 10), constrained_layout=True, dpi=100)
s = ax.pcolormesh(
    ds.x, ds.y, ds["u.x"].isel(level=lvl_surf), cmap="Greys_r", vmin=-1, vmax=2
)
plt.colorbar(s, ax=ax)
ax.pcolormesh(
    ds.x,
    ds.y,
    np.ma.masked_where(mask, np.ones(mask.shape)),
    cmap="cool",
    vmin=0,
    vmax=1,
)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_title("u.x and detected ridges")
plt.show()
hist = hist.assign_coords(c=bins_center)

#### Save to a separate file ####
filename = path + "breaking_hist.nc"
compression_settings = {
    "zlib": True,
    "complevel": 5,  # Compression level from 1 (fastest, least compression) to 9 (slowest, most compression)
}
hist.name = "hist"
hist.to_netcdf(filename, encoding={"hist": compression_settings})
print("Breaking stats saved!")
del hist  # Delete ds for memory
gc.collect()
