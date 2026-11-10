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


SCATTER = False
PLOT_KE = True
PLOT_GPE = True
CASE_PLOT_KEGPE = "hold"
skip = 1
Hskp = 0.078223  # 0.247363
L0 = 200
kp = 2 * 5 * np.pi / L0
g = 9.81
P = 0.02
rho = 1000

E0th = g * (Hskp / kp / 4) ** 2 * L0**2

data = {
    "decay": np.loadtxt("PM_decay.cuda/out", skiprows=1),
    "hold": np.loadtxt("PM_seastate_hold.cuda/out", skiprows=1),
    # "constant": np.loadtxt("PM_seastate.cuda/out", skiprows=1),
}
colors = {"hold": "r", "decay": "k", "constant": "orange"}
colorke = "b"
colorgpe = "g"
lss = {"hold": "--", "decay": "-", "constant": "dashdot"}


time = {}
ke = {}
gpe = {}
E = {}
E0 = {}

for case in data.keys():
    time[case] = data[case][:, 0]
    ke[case] = data[case][:, 1]
    gpe[case] = data[case][:, 2]
    E[case] = ke[case] + gpe[case]
    E0[case] = E[case][0]

print("========================")
print("PM spectrum P = %f" % P)
print("kpHs initial = %f" % Hskp)
print("L0 = %f" % L0)
print("E0th/rho = %f" % (E0th / rho))
print("E0/rho (actual) = %f" % (E["hold"][0] / rho))
print("========================")


fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True, dpi=200)

ax.hlines(1.0, 0, 5, colors="gray", alpha=0.7)


if SCATTER:
    for case in data.keys():
        ax.scatter(
            time[case][::skip],
            E[case][::skip] / E0,
            marker="x",
            c=colors[case],
            ls=lss[case],
            label=f"E ({case})",
        )
else:
    for case in data.keys():
        if PLOT_KE and case == CASE_PLOT_KEGPE:
            ax.plot(
                time[case][::skip],
                2 * ke[case][::skip] / E0[case],
                color=colorke,
                ls=lss[case],
            )
        if PLOT_GPE and case == CASE_PLOT_KEGPE:
            ax.plot(
                time[case][::skip],
                2 * gpe[case][::skip] / E0[case],
                color=colorgpe,
                ls=lss[case],
            )

        ax.plot(
            time[case][::skip],
            E[case][::skip] / E0[case],
            color=colors[case],
            ls=lss[case],
            label=f"E ({case})",
        )
ax.legend(loc="lower left")
ax.set_ylabel("E/E0")
ax.set_xlabel("t/T0")

fig.savefig("energy.png")

# res = "256"
# case_list = ["decay", "hold", "grow"]
# # case_list = ["hold", "decay"]
# dfiles = {
#     # "decay": {"256": "N256_15_PM_decay.cuda/out.nc"},
#     "decay": {"256": "PM_decay.cuda/out.nc"},
#     "hold": {"256": "PM_seastate_hold/out.nc"},
#     "grow": {"256": "PM_seastate.cuda/out.nc"},
# }
#
# fig, ax = plt.subplots(figsize=(6, 5))
# for case in case_list:
#     file = dfiles[case][res]
#     ds = xr.open_dataset(file, chunks={"x": 64, "y": 64})
#     N = len(ds.x)
#     L0 = 200
#     dx = ds.x[1] - ds.x[0]
#     dy = dx
#     kp = 5 * 2 * np.pi / L0
#     omegap = np.sqrt(9.81 * kp)
#     Tp = 2 * np.pi / omegap
#     # print("computing Ek")
#     skip = 2  # if = 1 just crash
#     Ek = (
#         0.5
#         * (
#             (ds["u.x"][::skip] ** 2 + ds["u.y"][::skip] ** 2 + ds["w"][::skip] ** 2)
#             * dx
#             * dy
#             * ds.h[::skip]
#         ).sum(dim=["level", "x", "y"])
#         / L0**2
#     ).values
#
#     # print("computing Ep")
#     Ep = (0.5 * 9.81 * (ds["eta"][::skip] ** 2 * dx).sum(dim=["x", "y"]) / L0**2).values
#     E = Ek + Ep
#     E0 = E[0]
#
#     time = ds.time[::skip]
#
#     # ax.plot(time / Tp, 2 * Ek / E0, c="b")  # , label="2Ek"
#     # ax.plot(time / Tp, 2 * Ep / E0, c="g")  # , label="2Ep"
#     ax.plot(time / Tp, E / E0, label=f"Et (case={case})")
#
# ax.set_ylabel("E/E0")
# ax.set_xlabel("t/Tp")
# ax.legend()
# fig.savefig(f"energy_{res}.pdf")

plt.show()
