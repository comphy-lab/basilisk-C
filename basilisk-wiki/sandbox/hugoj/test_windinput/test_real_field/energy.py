"""
Just testing the output of the case of Rui with small CFL
CFL_H = 0.05
CFL = 0.05

(he indicated CFL_H=0.01, but it is very long so I tried with 0.05 first)
It seems that the simulation is still unstable, energy grows ...
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


res = "256"
case_list = ["decay", "hold", "grow"]

dfiles = {
    "decay": {"256": "N256_15_PM_decay.cuda/out.nc"},
    "hold": {"256": "PM_seastate_hold.cuda/out.nc"},
    "grow": {"256": "PM_seastate.cuda/out.nc"},
}

# decay_files = {
#     256: "N256_15_PM_decay.cuda/out.nc",
#     # "512": "N512_15_PM_decay.cuda/out.nc",
# }
#
# grow_files = {
#     256: "PM_seastate.cuda/out.nc",
# }
#
# # hold_files = {"256": "N256_15_PM_hold.cuda/out.nc"}
# hold_files = {256: "PM_seastate_hold.cuda/out.nc"}
#
# if case == "decay":
#     dfiles = decay_files
# elif case == "hold":
#     dfiles = hold_files
# elif case == "grow":
#     dfiles = grow_files
# else:
#     raise Exception(f"case {case} is not recognized")

fig, ax = plt.subplots(figsize=(6, 5))
for case in case_list:
    file = dfiles[case][res]
    ds = xr.open_dataset(file, chunks={"x": 128, "y": 128})

    N = len(ds.x)
    L0 = 200
    dx = ds.x[1] - ds.x[0]
    dy = dx
    kp = 5 * 2 * np.pi / L0
    omegap = np.sqrt(9.81 * kp)
    Tp = 2 * np.pi / omegap
    # print("computing Ek")
    skip = 2  # if = 1 just crash
    Ek = (
        0.5
        * (
            (ds["u.x"][::skip] ** 2 + ds["u.y"][::skip] ** 2 + ds["w"][::skip] ** 2)
            * dx
            * dy
            * ds.h[::skip]
        ).sum(dim=["level", "x", "y"])
        / L0
    ).values

    # print("computing Ep")
    Ep = (0.5 * 9.81 * (ds["eta"][::skip] ** 2 * dx).sum(dim=["x", "y"]) / L0).values
    E = Ek + Ep
    E0 = E[0]

    time = ds.time[::skip]

    # ax.plot(time / Tp, 2 * Ek / E0, c="b")  # , label="2Ek"
    # ax.plot(time / Tp, 2 * Ep / E0, c="g")  # , label="2Ep"
    ax.plot(time / Tp, E / E0, label=f"Et (case={case})")

ax.set_ylabel("E/E0")
ax.set_xlabel("t/Tp")
ax.legend()
fig.savefig(f"energy_{res}.pdf")

# fig, ax = plt.subplots(figsize=(8, 6))
# s = ax.pcolormesh(ds.x, ds.y, ds.eta.isel(time=0), label="eta t=0", cmap="Greys_r")
# plt.colorbar(s, ax=ax, label="eta (m)")
# ax.set_title("eta at t=0 ")
# ax.set_aspect(1)
#
# fig, ax = plt.subplots(figsize=(8, 6))
# s = ax.pcolormesh(ds.x, ds.y, ds.eta.isel(time=-1), label="eta t=-1", cmap="Greys_r")
# plt.colorbar(s, ax=ax, label="eta (m)")
# ax.set_title("eta at t/Tp=%d " % (ds.time[-1].values / Tp))
# ax.set_aspect(1)
#
# fig, ax = plt.subplots(figsize=(8, 6))
# s = ax.pcolormesh(
#     ds.x, ds.y, ds["u.x"].isel(time=0, level=-1), cmap="Greys_r", vmin=-2, vmax=2
# )
# plt.colorbar(s, ax=ax, label="u (m/s)")
# ax.set_title("surface u.x at t=0 ")
# ax.set_aspect(1)
#
# fig, ax = plt.subplots(figsize=(8, 6))
# s = ax.pcolormesh(
#     ds.x, ds.y, ds["u.x"].isel(time=-1, level=-1), cmap="Greys_r", vmin=-2, vmax=2
# )
# plt.colorbar(s, ax=ax, label="u (m/s)")
# ax.set_title("surface u.x at t/Tp=%d " % (ds.time[-1] / Tp))
# ax.set_aspect(1)

plt.show()
