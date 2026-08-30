import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import colorcet as cc

# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)
from data_reader import read_bas_data, build_grid
from diags import h_divergence, vorticity, grad_velocities, grad_dir
from L2_maker import make_L2_layer, make_L2_eulerian


# General parameters
L0 = 200
H0 = 50
N = 1024
kp = 10 * np.pi / L0
lambdaP = 2 * np.pi / kp
Ndeux = 2e-6
g = 9.81
betaT = 2e-4
Ts = 20.0
analysis = "Lagrange"  # Euler or Lagrange
hchunks = 256

# inpath = "N512_nl30_0.000002_Tinizl/"
# infile = "out.nc"
# inpath = "~/work/BasiLagrangian/models/N1024/"
# infile = "out954.nc"
inpath = "N1024_nl60_cuda/"
infile = "out551.nc"
outfile = "L2"
lsfx = -1
it = -1

# opening file
ds, grid = read_bas_data(
    inpath + infile, chunks={"level": -1, "x": hchunks, "y": hchunks}
)

ds = ds.isel(time=it)  # for this specific file I select a snapshot only

nl = len(ds.level)
dx = (ds.x[1] - ds.x[0]).values
dy = dx
X = np.linspace(-L0 / 2, L0 / 2 + dx, N) / L0
Y = np.linspace(-L0 / 2, L0 / 2 + dy, N) / L0
LEVELS = np.arange(0, nl)

wet_z = ds.z.isel(level=-1).min().values
print("Depth always under the surface: %f (m)" % wet_z)


# Semi-Lagrangian
if analysis == "Lagrange":
    dsL2 = make_L2_layer(
        ds, grid, inpath, outfile + "_L_" + infile, save=True, verb=True
    )
    tke_sfx = dsL2.tke.isel(level=lsfx).mean(dim=["x", "y"]).values
    Zm = ds.z.mean(dim=["x", "y"]).values
    zname = "level"
    zsfx = len(ds.level) - 1
    gridL2 = grid
    div = h_divergence(ds, gridL2, zvar=zname)
# Eulerian
elif analysis == "Euler":
    Nznew = 50
    znew = np.linspace(-8, 0, Nznew + 1)  # cartesian Z to interpolate on
    dsL2 = make_L2_eulerian(
        ds,
        grid,
        znew,
        inpath,
        outfile + "_E_" + infile,
        fill_value=np.nan,
        save=True,
        verb=True,
    )
    Zm = dsL2.znew.values
    zname = "znew"
    zsfx = 0.0
    gridL2 = build_grid(dsL2, zvar="znew")
    div = h_divergence(dsL2, gridL2, zvar=zname)
else:
    raise Exception(f"analysis value = {analysis} is not recognized")


# surface tke
tke_sfx0 = dsL2.tke.sel({zname: zsfx}).mean(dim=["x", "y"]).compute()

print(f"tke_sfx_avg = {tke_sfx0.values}")
print(f"lambda_p = {lambdaP}")

# ----2D Plots ------------------------
if False:
    fig, ax = plt.subplots(1, 1, figsize=(9, 7), constrained_layout=True, dpi=100)
    s = ax.pcolormesh(
        X, Y, div.sel({zname: zsfx}), cmap=cc.cm.coolwarm, vmin=-1.5, vmax=1.5
    )
    plt.colorbar(s, ax=ax)
    ax.set_aspect(1)
    ax.set_xlabel("X/L0")
    ax.set_ylabel("Y/L0")
    ax.set_title("divergence (s-1)")
    fig.savefig(f"div_sfx_{analysis}.png")

    fig, ax = plt.subplots(1, 1, figsize=(9, 7), constrained_layout=True, dpi=100)
    s = ax.pcolormesh(
        X,
        Y,
        dsL2.tke.sel({zname: zsfx}) / tke_sfx0,
        cmap=cc.cm.bmy,
        vmin=0,
        vmax=3,
    )
    plt.colorbar(s, ax=ax)
    ax.set_aspect(1)
    ax.set_xlabel("X/L0")
    ax.set_ylabel("Y/L0")
    ax.set_title("tke/tke_0")
    fig.savefig(f"tke_sfx_{analysis}.png")

# ---- Profiles Plots ------------------------
currents = ["u", "v", "w"]
ucolors = ["b", "g", "r"]
if False:
    znorm = lambdaP
    # TKE and flux
    fig, ax = plt.subplots(1, 2, figsize=(9, 7), constrained_layout=True, dpi=100)
    ax[0].plot(
        dsL2.tke.mean(dim=["x", "y"]) / tke_sfx0,
        Zm / znorm,
        c="k",
        label=r"$\overline{tke}$",
        marker="x",
    )
    for k, var in enumerate(currents):
        ax[0].plot(
            dsL2[var + var].mean(dim=["x", "y"]) / tke_sfx0,
            Zm / znorm,
            c=ucolors[k],
            label=rf"$\overline{{{var}'{var}'}}$",
            marker="x",
        )

    ax[0].set_xlabel(r"$\overline{tke}$/tke_0")
    ax[0].set_ylabel(r"z/$\lambda_p$")

    for k, var in enumerate(currents):
        ax[1].plot(
            dsL2[var + "w"].mean(dim=["x", "y"]) / tke_sfx0,
            Zm / znorm,
            c=ucolors[k],
            label=rf"$\overline{{{var}'w'}}$",
            marker="x",
        )

    ax[1].set_xlabel(r"$\overline{u_i' w'}$/tke_0")
    for axe in ax:
        axe.hlines(wet_z / znorm, 0, 2.0, color="gray")
        axe.set_ylim([-0.5, 0])
        axe.set_xlim([0, 1.1])
        axe.legend()

    fig.savefig(f"tke_prof_{analysis}.pdf")

    # HEAT flux
    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    ax.vlines(0, -H0 / znorm, 0, color="gray")
    # trick to get the TT flux in the legend
    ax.plot(
        dsL2["uT"].mean(dim=["x", "y"]),
        Zm / znorm,
        c="k",
        label=r"$\overline{T'T'}$",
        marker="x",
    )
    for k, var in enumerate(currents):
        ax.plot(
            dsL2[var + "T"].mean(dim=["x", "y"]),
            Zm / znorm,
            c=ucolors[k],
            label=rf"$\overline{{T'{var}'}}$",
            marker="x",
        )

    ax2 = ax.twiny()
    ax2.plot(
        dsL2["TT"].mean(dim=["x", "y"]),
        Zm / znorm,
        c="k",
        marker="x",
    )
    ax2.set_xlabel(r"$\overline{T'T'}$")

    ax.set_xlabel(r"$\overline{T'u_i'}$")
    ax.set_ylabel(r"z/$\lambda_P$")
    ax.set_ylim([-0.5, 0])
    ax.legend()
    xmin, xmax, ymin, ymax = ax.axis()
    ax.set_xlim([xmin, xmax])
    ax.hlines(wet_z / znorm, xmin, xmax, color="gray")
    fig.savefig(f"heat_flux_prof_{analysis}.pdf")

    # MEAN CURRENTS
    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    ax.vlines(0, -H0 / znorm, 0, color="gray")
    ax.plot(dsL2.u_m, Zm / znorm, c="b", label="u", marker="x")
    ax.plot(dsL2.v_m, Zm / znorm, c="g", label="v", marker="x")
    ax.plot(dsL2.w_m, Zm / znorm, c="r", label="w", marker="x")
    ax.set_xlabel("mean currents")
    ax.set_ylabel(r"z/$\lambda_P$")
    ax.set_ylim([-0.5, 0])
    ax.set_xlim([-0.01, 0.1])
    ax.hlines(wet_z / znorm, 19, 21, color="gray")
    ax.legend()
    fig.savefig(f"mean_currents_prof_{analysis}.pdf")

    # MEAN TEMPERATURE
    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    ax.plot(dsL2.T_m, Zm / znorm, c="k", label="T", marker="x")
    ax.plot(Ndeux / (g * betaT) * Zm + Ts, Zm / znorm, label="Tini", color="gray")
    ax.set_xlabel("mean temperature")
    ax.set_ylabel(r"z/$\lambda_P$")
    # ax.set_xlim([19.99, 20.0])
    # ax.set_ylim([-0.5, 0])
    ax.hlines(wet_z / znorm, 19.99, 20.0, color="gray")
    ax.legend()
    fig.savefig(f"mean_temp_prof_{analysis}.pdf")

# --- 3D identification ---
if False:
    ds, _ = grad_velocities(ds, grid, zvar=zname)
    ds, _ = vorticity(ds, grid)

    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    s = ax.pcolormesh(Y, LEVELS, ds.omegayp.isel(x=0), cmap="coolwarm")
    plt.colorbar(s, ax=ax, label=r"$\omega_y$")

    def compute_tensors(ds):
        T = ds["dudx"].values + ds["dvdy"].values + ds["dwdz"].values
        shape = (3, 3)
        J = np.zeros(shape)
        J[0, 0] = ds["dudx"]
        J[0, 1] = ds["dudy"]
        J[0, 2] = ds["dudz"]
        J[1, 0] = ds["dvdx"]
        J[1, 1] = ds["dvdy"]
        J[1, 2] = ds["dvdz"]
        J[2, 0] = ds["dwdx"]
        J[2, 1] = ds["dwdy"]
        J[2, 2] = ds["dwdz"]
        Jt = J.transpose()
        S = 0.5 * (J + Jt)
        O = 0.5 * (J - Jt)
        return T, J, Jt, S, O

    def compute_Q(T, S, O):
        return 0.5 * (T**2 + O**2 - S**2)

    def is_lambda2(S, O):
        eigens = np.linalg.eig(S**2 + O**2)
        if eigens.eigenvalues[2] < 0:
            return 1
        else:
            return 0

    # This is way too long, I need to make it faster
    is_vortex = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            for k in range(nl):
                ds2 = ds.isel(x=i, y=j, level=k)
                T, J, Jt, S, O = compute_tensors(ds2)
                is_vortex[i, j] = is_lambda2(S, O)
                Q = compute_Q(T, S, O)

# TODO : plot


huu = ds.h * ds["u.x"] * ds["u.x"]
huv = ds.h * ds["u.x"] * ds["u.y"]

dhuudx = grad_dir(huu, ds, grid, dir="X", zvar="z")
dhuvdy = grad_dir(huv, ds, grid, dir="Y", zvar="z")

term1_avg_of_grad = dhuudx.mean(dim=["x", "y"])
term2_avg_of_grad = dhuvdy.mean(dim=["x", "y"])
sum = term1_avg_of_grad + term2_avg_of_grad

fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
ax.plot(term1_avg_of_grad, Zm, c="b", label=r"proper $\overline{\partial_x huu}$")
ax.plot(term2_avg_of_grad, Zm, c="g", label=r"proper $\overline{\partial_y huy}$")
ax.plot(
    huu.diff("x").mean(dim=["x", "y"]),
    Zm,
    c="b",
    ls="--",
    label=r"dirty $\overline{\partial_x huu}$",
)
ax.plot(
    huv.diff("y").mean(dim=["x", "y"]),
    Zm,
    c="g",
    ls="--",
    label=r"dirty $\overline{\partial_y huv}$",
)
ax.plot(sum, Zm, c="k", label="sum")
ax.legend()
ax.set_ylabel("z")
ax.set_xlabel("advection terms")

plt.show()
