import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl

# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)

from data_reader import read_bas_data
from tools import build_z

folder = "ml_breaking_strat.gpu/"
file = "out.nc"
g = 9.81
betaT = 2e-4
N2 = 2e-6
T0 = 20.0
rho0 = betaT * T0
L0 = 200.0
kp = 10 * np.pi / L0
geoB = 1 / 3
D = 50
omegap = np.sqrt(g * kp)
Tp = 2 * np.pi / omegap
ds, grid = read_bas_data(
    folder + file, chunks={"time": 25, "x": 128, "y": 128, "level": -1}
)
# ds = ds.sel(time=ds.time[1:])  # it=0 is initial condition
nl, N = ds.sizes["level"], ds.sizes["x"]


def Ek(u, v, w, h):
    Ny = h.sizes["y"]
    Nx = h.sizes["x"]
    return (
        0.5
        * ((u**2 + v**2 + w**2) * h).sum(dim="level").sum(dim=["x", "y"])
        / (Nx * Ny)
    )


def Ep_wave(eta: xr.DataArray):
    # need times rho for energy
    Ny = eta.sizes["y"]
    Nx = eta.sizes["x"]
    return g * (eta**2).sum(dim=["x", "y"]) / (Nx * Ny)


def Ep_strat(drho: xr.DataArray, z: xr.DataArray, h: xr.DataArray):
    # need times rho for energy
    Ny = h.sizes["y"]
    Nx = h.sizes["x"]
    int_v = (drho * g * z * h).sum(dim="level")
    int_h = int_v.sum(dim=["x", "y"])
    return int_h / (Nx * Ny)


def EOT(betaT, T0, T: xr.DataArray) -> xr.DataArray:
    return betaT * (T0 - T)


def geo_law(rmin, nl):
    r = 1.0 + 2.0 * (1.0 / rmin - 1.0) / (nl - 1.0)
    hmin = (r - 1.0) / (pow(r, nl) - 1.0)
    beta = np.zeros(nl)
    for k in range(nl):
        beta[k] = hmin * r ** (nl - 1 - k)
    return beta


def layer_thickness(rmin, nl, D):
    return D * geo_law(rmin, nl)


h_ini1D = layer_thickness(geoB, nl, D)
h_ini = xr.DataArray(
    np.broadcast_to(h_ini1D[:, None, None], (nl, N, N)), dims=["level", "y", "x"]
)
zb = xr.DataArray(-D * np.ones((N, N)), dims=["y", "x"])
z_ini, z_lini = build_z(h_ini, zb)


Tini = N2 / g * z_ini + T0
drho_ini = EOT(betaT, T0, Tini)
drho = EOT(betaT, T0, ds.T)


Ep_s_t = Ep_strat(drho, ds.z, ds.h)
Ep_s_0 = Ep_strat(drho_ini, z_ini, h_ini)
Ep_w = Ep_wave(ds.eta)
Ep_s = Ep_s_t - Ep_s_0
Ep_total = Ep_w + Ep_s
Ec = Ek(ds["u.x"], ds["u.y"], ds["u.z"], ds.h)
Et = Ec + Ep_total

# profiles
# Tprof = ds.T.mean(dim=["x", "y"]).compute()
# Zm = ds.z.mean(dim=["x", "y"]).compute()
# cmap = mpl.colormaps["plasma"]
# colors = cmap(np.linspace(0, 1, len(ds.time)))
#
# fig, ax = plt.subplots(1, 1, figsize=(3, 7), constrained_layout=True, dpi=100)
# for it in range(len(ds.time)):
#     ax.plot(Tprof[it], Zm[it], c=colors[it])
# ax.set_xlabel("T(°C)")
# ax.set_ylabel("Z(m)")

fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
# ax.plot(ds.time / Tp, Ep_wave / Ep_wave[0], label="Ep_wave")
# ax.plot(ds.time / Tp, Ep_stratif / Ep_stratif[0], label="Ep_stratif")
ax.plot(ds.time / Tp, Ep_s, label="Ep_s")
ax.plot(ds.time / Tp, Ep_w, label="Ep_w")
ax.plot(ds.time / Tp, Ep_total, label="ep_total")
ax.plot(ds.time / Tp, Ec, label="Ec")
ax.set_ylabel("E")
ax.set_xlabel("t/Tp")
plt.legend()
fig.savefig("Energy.svg")


fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
ax.plot(ds.time / Tp, Ep_s / Ep_s[0], label="Ep_s")
ax.plot(ds.time / Tp, Ep_w / Ep_w[0], label="Ep_w")
ax.plot(ds.time / Tp, Ep_total / Ep_total[0], label="ep_total")
ax.plot(ds.time / Tp, Ec / Ec[0], label="Ec")
ax.set_ylabel("E/E0")
ax.set_xlabel("t/Tp")
plt.legend()
fig.savefig("Energy_normed_indiv.png")

fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
ax.plot(ds.time / Tp, Ep_s / Et[0], label="Ep_s")
ax.plot(ds.time / Tp, Ep_w / Et[0], label="Ep_w")
ax.plot(ds.time / Tp, Ep_total / Et[0], label="Ep_s+Ep_w")
ax.plot(ds.time / Tp, Ec / Et[0], label="Ec")
ax.plot(ds.time / Tp, Et / Et[0], label="Et")
ax.set_ylabel("E/E0")
ax.set_xlabel("t/Tp")
ax.set_title("bilan")
plt.legend()
fig.savefig("Energy_bilan.svg")


# plt.show()


# loading data
# rawT = np.loadtxt(folder + "T_profile.dat")
# rawU = np.loadtxt(folder + "u_profile.dat")
#
#
# def get_Nl_Nt(rawfile):
#     dt = 0.0
#     Nl = 0
#     while dt == 0:
#         Nl = Nl + 1
#         dt = rawfile[Nl][0]
#     Nt = int(rawfile[-1][0] / dt) + 1
#     return Nl, Nt
#
#
# def read_raw(raw, Ny, Nx):
#     data = np.zeros((Ny, Nx))
#     X = np.zeros(Nx)
#     Y = np.zeros(Ny)
#     for i in range(Nx):
#         X[i] = raw[i * Ny][0]
#         for k in range(Ny):
#             index = i * Ny + k
#             Y[k] = raw[k][1]
#             data[k, i] = raw[index][-1]
#     return X, Y, data
#
#
# # > T_profile
# Nl, Nt = get_Nl_Nt(rawT)
# time, layers, dataT = read_raw(rawT, Nl, Nt)
# T0 = dataT[-1, 0]
# # plot
# fig, ax = plt.subplots(1, 1, figsize=(7, 3), constrained_layout=True, dpi=100)
# # s = ax.pcolormesh(time, layers, dataT/T0, shading='nearest', cmap='magma')
# # plt.colorbar(s,ax=ax, label='$T/Ts_{0}$')
# s = ax.pcolormesh(time, layers, dataT, shading="nearest", cmap="magma")
# plt.colorbar(s, ax=ax, label="$T (°C)$")
# ax.set_xlabel("time (s)")
# ax.set_ylabel("layer")
# # fig.savefig('T_profile.png')
# plt.show()


#
# # > wT_profile
# Nl, Nt = get_Nl_Nt(rawU)
# time, layers, dataU = read_raw(rawU, Nl, Nt)
# #flx0 = 500/(1025*4.2e3)
# # plot
# fig, ax = plt.subplots(1,1,figsize = (7,3),constrained_layout=True,dpi=100)
# s=ax.pcolormesh(time,layers,dataU, shading="nearest", cmap='magma')
# plt.colorbar(s,ax=ax, label=r"U ($m.s^{-1}$)")
# ax.set_xlabel('time (s)')
# ax.set_ylabel('layer')
# fig.savefig('u_profile.png')

plt.show()
