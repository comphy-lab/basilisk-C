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
from diags import h_divergence
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
analysis = "Euler"  # Euler or Lagrange
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

dx = (ds.x[1] - ds.x[0]).values
dy = dx
X = np.linspace(-L0 / 2, L0 / 2 + dx, N) / L0
Y = np.linspace(-L0 / 2, L0 / 2 + dy, N) / L0


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

else:
    raise Exception(f"analysis value = {analysis} is not recognized")


# surface tke
tke_sfx0 = dsL2.tke.sel({zname: zsfx}).mean(dim=["x", "y"]).compute()

# surface divergence
div = h_divergence(dsL2, gridL2, zvar="znew")

print(f"tke_sfx_avg = {tke_sfx0.values}")
print(f"lambda_p = {lambdaP}")

# ----2D Plots ------------------------
if True:
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
if True:
    znorm = lambdaP

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
        axe.set_xlim([0, 1.75])
        axe.legend()

    fig.savefig(f"tke_prof_{analysis}.pdf")

    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    ax.vlines(0, -H0 / znorm, 0, color="gray")
    for k, var in enumerate(currents):
        ax.plot(
            dsL2["T" + var].mean(dim=["x", "y"]),
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
        label=r"$\overline{T'T'}$",
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

    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    ax.vlines(0, -H0 / znorm, 0, color="gray")
    ax.plot(dsL2.u_m, Zm / znorm, c="b", label="u", marker="x")
    ax.plot(dsL2.v_m, Zm / znorm, c="g", label="v", marker="x")
    ax.plot(dsL2.w_m, Zm / znorm, c="r", label="w", marker="x")
    ax.set_xlabel("mean currents")
    ax.set_ylabel(r"z/$\lambda_P$")
    ax.set_ylim([-0.5, 0])
    ax.set_xlim([-0.01, 0.5])
    ax.hlines(wet_z / znorm, -1, 1, color="gray")
    ax.legend()
    fig.savefig(f"mean_currents_prof_{analysis}.pdf")

    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    ax.plot(dsL2.T_m, Zm / znorm, c="k", label="T", marker="x")
    ax.plot(Ndeux / (g * betaT) * Zm + Ts, Zm / znorm, label="Tini", color="gray")
    ax.set_xlabel("mean temperature")
    ax.set_ylabel(r"z/$\lambda_P$")
    ax.set_xlim([19.99, 20.0])
    ax.set_ylim([-0.5, 0])
    ax.hlines(wet_z / znorm, 19.99, 20.0, color="gray")
    ax.legend()
    fig.savefig(f"mean_temp_prof_{analysis}.pdf")


plt.show()
