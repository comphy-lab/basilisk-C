"""
Figures for CSI report 2026
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import colorcet as cc
import matplotlib.ticker as ticker
import time


# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)
from data_reader import read_bas_data, build_grid
from diags import (
    compute_us_Kenyon1969,
    grad_dir,
    compute_spectrum_trange,
    compute_Ep,
    compute_Ek,
    EOS,
)
from L2_maker import make_L2_layer, make_L2_eulerian
from visu_3Dsnap import render_snapshot
from fftlib import get_spec_1D, get_wavenumber

from dask.distributed import Client, LocalCluster


def main():

    if True:
        cluster = LocalCluster(
            n_workers=4,  # fewer workers, more memory each, tune to your machine
            threads_per_worker=2,
            memory_limit="4GB",  # per worker — dask will spill to disk before crashing
        )
        client = Client(cluster)
        print(client.dashboard_link)

    # General parameters
    L0 = 200
    H0 = 20
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
    hchunks = 64

    # inpath = "N512_nl30_0.000002_Tinizl/"
    # inpath = "N512_30_P0.02_RE40000_TiniZ/"
    # inpath = "N512_30_P0.02_RE40000_TiniZl/"
    inpath = "N512_40_RE40000_TiniZl_remap1.0/"
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

    Zm = dsL2.z.mean(dim=["x", "y"]).compute()

    colors = cc.fire[:: len(cc.fire) // nt]

    # ========================
    # --- Energy evolution ---
    # ========================
    if True:
        skip = 1

        drho = EOS(ds.T, Ts, betaT)

        # Potential energy
        # step 1 : at rest with linear stratification
        dz0 = H0 / nl
        z0 = np.arange(-H0 + dz0 / 2, 0, dz0)
        Tini = Ndeux / (g * betaT) * z0 + Ts
        drho0 = EOS(Tini, Ts, betaT)
        # step 2 : at every t
        Ep = compute_Ep(
            (1 + drho), g, ds, L0, Ep0=np.sum((1 + drho0) * g * z0 * dz0), skip=skip
        )
        # Note: we could split 'wave' and 'thermodynamic' potential energy
        # Ep_s0 = np.sum(drho0 * g * z0 * dz0)
        # Ep_s = compute_Ep(drho, ds.z, ds.h, Ep0=Ep_s0, skip=skip)
        # Ep_w = (0.5 * g * ds.eta**2 * dx**2).sum(dim=["x", "y"]) / L0**2

        # Kinetic energy
        Ek = compute_Ek(ds, L0, skip=skip)

        # Total energy
        Et = Ek + Ep

        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
        ax.plot(ds.time / Tp, Ep / Ep.isel(time=0), c="g", label="Ep")
        ax.plot(ds.time / Tp, Ek / Ek.isel(time=0), c="b", label="Ek")
        ax.plot(ds.time / Tp, Et / Et.isel(time=0), c="k", label="Et=Ek+Ep")
        ax.set_xlabel("t/Tp")
        ax.set_ylabel("E/E0")
        ax.legend()
        fig.savefig("CSI_energy_evolution.pdf")

    # =====================================
    # --- mean profiles evolution plots ---
    # =====================================

    if False:
        skip = 10

        fig, ax = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True, dpi=300)
        print("-> Tm profiles")
        for it in range(0, nt, skip):
            if it == 0:
                alpha = 1
                c = "gray"
            else:
                alpha = (it + 1) / (2 * nt)
                c = "b"
            ax[0].plot(
                (dsL2.T_m.isel(time=it).compute() - Ts) * 1000,
                Zm.isel(time=it),
                c=c,
                alpha=alpha,
                marker="x",
            )
        ax[0].plot(
            (dsL2.T_m.sel(time=att, method="nearest") - Ts) * 1000,
            Zm.sel(time=att, method="nearest"),
            c="k",
            marker="x",
        )
        ax[0].set_ylim((-H0, 0))
        ax[0].set_xlim((-25, 1))
        ax[0].set_xlabel(r"$\overline{T}-T_s$ (mK)")
        ax[0].set_ylabel(r"$\overline{z}$")

        print("-> Um profiles")
        this_eta = ds.eta.sel(time=att, method="nearest").values
        kr, phi_k = get_spec_1D(this_eta, this_eta, L0 / N, averaging="radial")
        us_Kenyon = compute_us_Kenyon1969(
            phi_k / (2 * np.pi),
            kr * 2 * np.pi,
            Zm.sel(time=att, method="nearest").values,
        )
        ax[1].vlines(0, -H0, 0, color="gray", ls="-", alpha=0.5)
        for it in range(0, nt, skip):
            if it == 0:
                alpha = 1
                c = "gray"
            else:
                alpha = (it + 1) / (2 * nt)
                c = "b"

            ax[1].plot(
                dsL2.u_m.isel(time=it),
                Zm.isel(time=it),
                c=c,
                alpha=alpha,
                marker="x",
            )
        ax[1].plot(
            dsL2.u_m.sel(time=att, method="nearest"),
            Zm.sel(time=att, method="nearest"),
            c="k",
            label="t=%d Tp" % (att / Tp),
            marker="x",
        )
        ax[1].plot(
            us_Kenyon,
            Zm.sel(time=att, method="nearest"),
            c="k",
            ls="--",
            label=r"$U_s$",
        )

        ax[1].set_xlim((-0.01, 0.2))
        ax[1].set_ylim((-H0, 0))
        ax[1].set_xlabel(r"$\overline{u}$ (m/s)")
        ax[1].set_ylabel(r"$\overline{z}$")

        print("-> Tke profiles")
        ax[2].vlines(0, -H0, 0, color="gray", ls="-", alpha=0.5)
        for it in range(0, nt, skip):
            if it == 0:
                alpha = 1
                c = "gray"
            else:
                alpha = (it + 1) / (2 * nt)
                c = "b"
            ax[2].plot(
                dsL2["tke"].mean(dim=["x", "y"]).isel(time=it),
                Zm.isel(time=it),
                c=c,
                alpha=alpha,
                marker="x",
            )
        ax[2].plot(
            dsL2["tke"].mean(dim=["x", "y"]).sel(time=att, method="nearest"),
            Zm.sel(time=att, method="nearest"),
            c="k",
            marker="x",
        )
        ax[2].set_ylim((-H0, 0))
        ax[2].set_xlabel(r"$\overline{tke}$ (m2/s2)")
        ax[2].set_ylabel(r"$\overline{z}$")

        fig.savefig("CSI_mean_profile_evolution.pdf")

    # TODO:
    # -

    # ==========================
    # --- spectrum evolution ---
    # ==========================

    if False:
        print("Plotting spectrum !")

        listt = [0, 2 * Tp, 5 * Tp, 10 * Tp]
        start = att - 5 * Tp
        end = att + 5 * Tp

        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=300)

        for it in range(len(listt)):
            this_eta = ds.eta.sel(time=listt[it], method="nearest")
            kr, phi_k = get_spec_1D(this_eta, this_eta, L0 / N, averaging="radial")

            ax.loglog(
                (kr * 2 * np.pi) * L0,
                (phi_k / (2 * np.pi)) * kp**3,
                c="b",
                alpha=(it + 1) / (2 * len(listt)),
                label="t/Tp=%d" % (listt[it] / Tp),
            )

        kr, phi_k = compute_spectrum_trange(ds, start, end)
        ax.loglog(
            kr * L0,
            (phi_k / (2 * np.pi)) * kp**3,
            c="k",
            label="t/Tp=[%d,%d]" % (start / Tp, end / Tp),
        )

        # nice slope k**-2.5
        x0, x1 = 1e1, 1e3  # adjust to sit in an empty-ish part of the plot
        slope = -2.5
        y0 = 8e-2  # adjust vertically so the line sits near your data
        y1 = y0 * (x1 / x0) ** slope
        ax.loglog([x0, x1], [y0, y1], "--", lw=1, c="g")
        # ax.text(x1 * 1.1, y1, r"$k^{-2.5}$", fontsize=10, va="center")
        ax.text(300, 5e-5, r"$k^{-2.5}$", fontsize=10, va="center", c="g")

        # nice slope k**-3
        x0, x1 = 1e1, 1e3  # adjust to sit in an empty-ish part of the plot
        slope = -3
        y0 = 3e-1  # adjust vertically so the line sits near your data
        y1 = y0 * (x1 / x0) ** slope
        ax.loglog([x0, x1], [y0, y1], "--", lw=1, c="coral")
        ax.text(400, 2e-5, r"$k^{-3}$", fontsize=10, va="center", c="coral")

        ax.set_xlabel("kL0")
        ax.grid()
        ax.legend()
        ax.set_xlim([10, 1000])
        ax.set_ylim([1e-7, 1e-2])
        ax.set_ylabel(r"$\phi(k) k_p^{3}$")
        fig.savefig("CSI_spectrum_t%d-%d.pdf" % (start, end))

    if False:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=300)
        ax.plot((dsL2.T_m.isel(time=att) - Ts) * 1000, Zm.isel(time=att), c="k")
        # ax.set_xlim((0.999, 1))
        ax.set_ylim((-20, 0))
        ax.set_xlabel(r"$\overline{T}-T_s$ (mK)")
        ax.set_ylabel("z")
        fig.savefig(f"CSI_mean_T_at{att}.pdf")

        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=300)
        ax.plot(dsL2.u_m.isel(time=att), Zm.isel(time=att), c="k")
        ax.set_xlim((-0.1, 0.1))
        ax.set_ylim((-20, 0))
        ax.set_xlabel(r"$\overline{u}$ (m/s)")
        ax.set_ylabel("z")
        fig.savefig(f"CSI_mean_U_at{att}.pdf")

        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=300)
        ax.plot(dsL2.wT.isel(time=att).mean(dim=["x", "y"]), Zm.isel(time=att), c="k")
        # ax.set_xlim((-0.1, 0.1))
        ax.set_ylim((-20, 0))
        ax.set_xlabel(r"$\overline{w'T'}$ (K.m/s)")
        ax.set_ylabel("z")
        fig.savefig(f"CSI_mean_wT_at{att}.pdf")

    # =============================
    # --- A snapshot of surface ---
    # =============================
    def fmt(x, pos):
        a, b = "{:.5e}".format(x).split("e")
        b = int(b)
        return r"${} \times 10^{{{}}}$".format(a, b)

    if False:
        fig, ax = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, dpi=300)
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

    # ==================
    # -- Diag de MLD ---
    # ==================
    # TODO: !
    if False:
        skip = 10

        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=300)
        print("-> MLD diag: dTm/dz")
        for it in range(0, nt, skip):
            dTmdz = grad_dir(
                dsL2.T_m.isel(time=it), dsL2.isel(time=it), grid, dir="Z", zvar="level"
            )
            ax.plot(
                dTmdz,
                Zm.isel(time=it),
                c="b",
                alpha=(it + 1) / (2 * nt),
                marker="x",
            )
        ax.set_ylim((-H0 / 5, 0))
        # ax[0].set_xlim((-25, 1))
        ax.set_xlabel(r"$\partial_z \overline{T}$")
        ax.set_ylabel(r"$\overline{z}$")

    # ===========
    # 3D PLOTS
    # ===========

    # --- nice picture for 1rst page ---
    if False:
        Tclim = (19.95, 20)
        render_snapshot(
            inpath + infile,
            "CSI_1rst_page",
            ttime=int(att),
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
            outpng=True,
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


if __name__ == "__main__":
    main()
