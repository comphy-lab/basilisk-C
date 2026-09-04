"""
Figures for CSI report 2026
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import colorcet as cc
import matplotlib.ticker as ticker

# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)
from data_reader import read_bas_data, build_grid
from diags import compute_us_Kenyon1969
from L2_maker import make_L2_layer, make_L2_eulerian
from visu_3Dsnap import render_snapshot
from fftlib import get_spec_1D, get_wavenumber

# General parameters
L0 = 200
H0 = 100
N = 512
g = 9.81
kp = 10 * np.pi / L0
lambdaP = 2 * np.pi / kp
omegap = np.sqrt(g * kp)
Tp = 2 * np.pi / omegap
Ndeux = 2e-6
betaT = 2e-4
Ts = 20.0
analysis = "Lagrange"  # Euler or Lagrange
hchunks = 256

inpath = "N512_nl30_0.000002_Tinizl/"
# inpath = "N512_30_P0.02_RE40000_TiniZ/"
inpath = "N512_30_P0.02_RE40000_TiniZl/"
infile = "out.nc"
outfile = "L2_" + infile
lsfx = -1
att = 30 * Tp  # s, for snapshots

# opening file
ds, grid = read_bas_data(
    inpath + infile, chunks={"time": -1, "level": -1, "x": hchunks, "y": hchunks}
)
dsL2 = make_L2_layer(ds, grid, inpath, outfile, verb=True)

nl = len(ds.level)
nt = len(ds.time)
dx = (ds.x[1] - ds.x[0]).values
dy = dx
X = np.linspace(-L0 / 2, L0 / 2 + dx, N) / L0
Y = np.linspace(-L0 / 2, L0 / 2 + dy, N) / L0
LEVELS = np.arange(0, nl)

# wet_z = ds.z.isel(level=-1).min().values
# print("Depth always under the surface: %f (m)" % wet_z)

Zm = dsL2.z.mean(dim=["x", "y"])

colors = cc.fire[:: len(cc.fire) // nt]


# --- mean profiles evolution plots ---
if True:
    skip = 10

    fig, ax = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True, dpi=100)
    print("-> Tm profiles")
    for it in range(0, nt, skip):
        ax[0].plot(
            (dsL2.T_m.isel(time=it) - Ts) * 1000, Zm.isel(time=it), c="b", alpha=0.1
        )
    ax[0].plot(
        (dsL2.T_m.sel(time=att, method="nearest") - Ts) * 1000,
        Zm.sel(time=att, method="nearest"),
        c="k",
    )
    ax[0].set_ylim((-H0 / 5, 0))
    ax[0].set_xlim((-25, 1))
    ax[0].set_xlabel(r"$\overline{T}-T_s$ (mK)")
    ax[0].set_ylabel(r"\overline{z}")

    print("-> Um profiles")
    this_eta = ds.eta.sel(time=att, method="nearest").values
    kr, phi_k = get_spec_1D(this_eta, this_eta, L0 / N, averaging="radial")
    us_Kenyon = compute_us_Kenyon1969(
        phi_k / (2 * np.pi), kr * 2 * np.pi, Zm.sel(time=att, method="nearest").values
    )

    ax[1].vlines(0, -H0, 0, color="gray", ls="-")
    for it in range(0, nt, skip):
        ax[1].plot(dsL2.u_m.isel(time=it), Zm.isel(time=it), c="b", alpha=0.1)
    ax[1].plot(
        dsL2.u_m.sel(time=att, method="nearest"),
        Zm.sel(time=att, method="nearest"),
        c="k",
        label="t=%d Tp" % (att / Tp),
    )
    ax[1].plot(
        us_Kenyon, Zm.sel(time=att, method="nearest"), c="k", ls="--", label=r"$U_s$"
    )
    ax[1].vlines(0, -H0, 0, color="gray", ls="-")
    ax[1].set_xlim((-0.01, 0.2))
    ax[1].set_ylim((-H0 / 5, 0))
    ax[1].set_xlabel(r"$\overline{u}$ (m/s)")
    ax[1].set_ylabel(r"\overline{z}")

    print("-> Tke profiles")
    ax[2].vlines(0, -H0, 0, color="gray", ls="-")
    for it in range(0, nt, skip):
        ax[2].plot(
            dsL2["tke"].mean(dim=["x", "y"]).isel(time=it),
            Zm.isel(time=it),
            c="b",
            alpha=0.1,
        )
    ax[2].plot(
        dsL2["tke"].mean(dim=["x", "y"]).sel(time=att, method="nearest"),
        Zm.sel(time=att, method="nearest"),
        c="k",
    )
    ax[2].set_ylim((-H0 / 5, 0))
    ax[2].set_xlabel(r"$\overline{tke}$ (m2/s2)")
    ax[2].set_ylabel(r"\overline{z}")

    fig.savefig("CSI_mean_profile_evolution.pdf")

# TODO:
# - plot Stokes drift on U current plot


# --- spectrum evolution ---
if True:
    start = att - 20
    end = att + 20

    # find items in [start, end]
    subds = ds.sel(time=slice(start, end))

    _, _, _, k_sample = get_wavenumber(N, L0 / N)

    phi_k = np.zeros((len(subds.time), len(k_sample)))
    kr = np.zeros((len(subds.time), len(k_sample)))

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
    for it in range(0, len(subds.time)):
        this_eta = subds.eta.isel(time=it).values
        tmp = get_spec_1D(this_eta, this_eta, L0 / N, averaging="radial")
        phi_k[it] = tmp[1]
        kr[it] = tmp[0]

    kr *= 2 * np.pi

    phi_k = phi_k.mean(axis=0)
    kr = kr.mean(axis=0)
    ax.loglog(kr * L0, (phi_k / (2 * np.pi)) * kp**3, c=colors[it])
    ax.set_xlabel("kL0")
    ax.grid()
    ax.set_ylabel(r"$\phi(k) k_p^{3}$")
    fig.savefig("CSI_spectrum_evolution_t%d-%d.pdf" % (start, end))

if False:
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
    ax.plot((dsL2.T_m.isel(time=att) - Ts) * 1000, Zm.isel(time=att), c="k")
    # ax.set_xlim((0.999, 1))
    ax.set_ylim((-20, 0))
    ax.set_xlabel(r"$\overline{T}-T_s$ (mK)")
    ax.set_ylabel("z")
    fig.savefig(f"CSI_mean_T_at{att}.pdf")

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
    ax.plot(dsL2.u_m.isel(time=att), Zm.isel(time=att), c="k")
    ax.set_xlim((-0.1, 0.1))
    ax.set_ylim((-20, 0))
    ax.set_xlabel(r"$\overline{u}$ (m/s)")
    ax.set_ylabel("z")
    fig.savefig(f"CSI_mean_U_at{att}.pdf")

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
    ax.plot(dsL2.wT.isel(time=att).mean(dim=["x", "y"]), Zm.isel(time=att), c="k")
    # ax.set_xlim((-0.1, 0.1))
    ax.set_ylim((-20, 0))
    ax.set_xlabel(r"$\overline{w'T'}$ (K.m/s)")
    ax.set_ylabel("z")
    fig.savefig(f"CSI_mean_wT_at{att}.pdf")


# --- A snapshot of surface ---
def fmt(x, pos):
    a, b = "{:.5e}".format(x).split("e")
    b = int(b)
    return r"${} \times 10^{{{}}}$".format(a, b)


if False:
    fig, ax = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, dpi=100)
    s = ax.pcolormesh(
        ds.x,
        ds.y,
        # ds.T.isel(time=att, level=lsfx) / Ts,
        (ds.T.isel(time=att, level=lsfx) - Ts) * 1000,
        # vmin=19.996,
        # vmax=20.0,
        cmap=cc.m_fire,
    )
    ax.set_aspect(1)
    ax.set_xlabel("X(m)")
    ax.set_ylabel("Y(m)")
    # plt.colorbar(s, ax=ax, format=ticker.FuncFormatter(fmt), label="T-T0 (mK)")
    plt.colorbar(s, ax=ax, label="T-T0 (mK)")
    fig.savefig(f"CSI_sfx_T_at{att}.pdf")

# -- Diag de MLD ---
# TODO: !

# ===========
# 3D PLOTS
# ===========

# --- nice picture for 1rst page ---
if False:
    Tclim = (19.95, 20)
    render_snapshot(
        inpath + infile,
        "CSI_1rst_page",
        ttime=att,
        method="nearest",
        L0=L0,
        H0=H0,
        var_top="u.x",
        clim_top=(-2, 2.5),
        cmap_top="Greys_r",
        var_side="T",
        clim_side=Tclim,
        cmap_side="plasma",
        background="white",
        window_size=(1024, 1024),
        off_screen=False,
        verbose=True,
        xclip=None,
        yclip=None,
        zclip=None,
    )


# --- zoom on a breaking ---
xmin, xmax, ymin, ymax = -50, 50, -50, 50
if False:
    render_snapshot(
        inpath + infile,
        "zoom",
        ttime=att,
        L0=L0,
        H0=H0,
        var_top="u.x",
        clim_top=(-2.5, 2.5),
        cmap_top="Greys_r",
        var_side="T",
        clim_side=Tclim,
        cmap_side="plasma",
        window_size=(1024, 1024),
        off_screen=False,
        verbose=True,
        xclip=(xmin, xmax),
        yclip=(ymin, ymax),
    )

# --- slice ---
xmin, xmax, ymin, ymax = -50, 50, 0, 0
if False:
    render_snapshot(
        inpath + infile,
        "slice",
        ttime=att,
        L0=L0,
        H0=H0,
        var_top="u.x",
        clim_top=(-2.5, 2.5),
        cmap_top="Greys_r",
        var_side="T",
        clim_side=Tclim,
        cmap_side="plasma",
        window_size=(1024, 1024),
        off_screen=False,
        verbose=True,
        xclip=(xmin, xmax),
        yclip=(ymin, ymax),
    )

plt.show()
