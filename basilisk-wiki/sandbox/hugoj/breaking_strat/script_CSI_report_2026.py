"""
Figures for CSI report 2026


you might need to change the temp dir for dask:
(if you get out of memory error because workers are spilling too much data on disk)
export TMPDIR=/home/jacqhugo/basilisk/wiki/sandbox/hugoj/breaking_strat/tmp
"""

import numpy as np
import matplotlib.pyplot as plt
import colorcet as cc
import time
import xarray as xr
from scipy.optimize import curve_fit

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
            n_workers=8,  # fewer workers, more memory each, tune to your machine
            threads_per_worker=1,
            memory_limit="2GB",  # per worker — dask will spill to disk before crashing
        )
        client = Client(cluster)
        print("dashboard :", client.dashboard_link)

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
    # inpath = "N512_40_RE40000_TiniZl_remap1.0/"
    inpath = "N512_40_RE40000_TiniZl_remap1.0_cuda/"
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
    if False:
        print("Energy evolution")
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

        ET0 = Et.isel(time=0)
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
        ax.plot(ds.time / Tp, 2 * Ep / ET0, c="g", label="2Ep")
        ax.plot(ds.time / Tp, 2 * Ek / ET0, c="b", label="2Ek")
        ax.plot(ds.time / Tp, Et / ET0, c="k", label="Et=Ek+Ep")
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

        kr, phi_k = compute_spectrum_trange(ds, N, L0, start, end)
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
        fig, ax = plt.subplots(1, 1, figsize=(6, 5), constrained_layout=True, dpi=300)
        s = ax.pcolormesh(
            ds.x,
            ds.y,
            # ds.T.isel(time=att, level=lsfx) / Ts,
            (ds.T.sel(time=att, method="nearest").isel(level=lsfx) - Ts) * 1000,
            vmin=-11,
            vmax=-4,
            cmap=cc.m_bmy,
        )
        ax.set_aspect(1)
        ax.set_xlabel("X(m)")
        ax.set_ylabel("Y(m)")
        # plt.colorbar(s, ax=ax, format=ticker.FuncFormatter(fmt), label="T-T0 (mK)")
        plt.colorbar(s, ax=ax, label="T-Ts (mK)")
        fig.savefig("CSI_sfx_T_at%d_sfx.pdf" % (att))

    if False:
        atlvl = 30
        zm = (
            ds.z.sel(time=att, method="nearest")
            .isel(level=atlvl)
            .mean(dim=["x", "y"])
            .values
        )
        fig, ax = plt.subplots(1, 1, figsize=(6, 5), constrained_layout=True, dpi=300)
        s = ax.pcolormesh(
            ds.x,
            ds.y,
            # ds.T.isel(time=att, level=lsfx) / Ts,
            (ds.T.sel(time=att, method="nearest").isel(level=atlvl) - Ts) * 1000,
            vmin=-11,
            vmax=-4,
            cmap=cc.m_bmy,
        )
        ax.set_aspect(1)
        ax.set_xlabel("X(m)")
        ax.set_ylabel("Y(m)")
        # plt.colorbar(s, ax=ax, format=ticker.FuncFormatter(fmt), label="T-T0 (mK)")
        plt.colorbar(s, ax=ax, label="T-Ts (mK)")
        fig.savefig("CSI_sfx_T_at%d_d%d.pdf" % (att, np.abs(zm)))

    # ==================
    # -- Diag de MLD ---
    # ==================
    if True:
        skip = 10

        dTdz_ini = Ndeux / (g * betaT)

        if False:
            fig, ax = plt.subplots(
                1, 1, figsize=(5, 5), constrained_layout=True, dpi=300
            )
            print("-> MLD diag: dTm/dz")
            for it in range(0, nt, skip):
                dTmdz = grad_dir(
                    dsL2.T_m.isel(time=it),
                    dsL2.isel(time=it),
                    grid,
                    dir="Z",
                    zvar="level",
                )
                ax.plot(
                    dTmdz,
                    Zm.isel(time=it),
                    c="b",
                    alpha=(it + 1) / (2 * nt),
                    marker="x",
                )
            ax.set_ylim((-H0, 0))
            ax.set_xlim((-dTdz_ini * 1.5, dTdz_ini * 1.5))
            ax.set_xlabel(r"$\partial_z \overline{T}$")
            ax.set_ylabel(r"$\overline{z}$")
            fig.savefig("CSI_gradT.pdf")

        """
        MLD = first level where change of dTdz is less than X% of dTdz_ini
        """
        if True:
            thrs = 0.25
            MLD = np.zeros(len(ds.time))

            dTmdz = grad_dir(
                dsL2.T_m,
                dsL2,
                grid,
                dir="Z",
                zvar="level",
            )

            dTmdz_diff = dTmdz.diff(dim="level")  # level dim now has nl-1 points...
            # dTmdz_diff[k] = dTmdz[k+1] - dTmdz[k]
            rel_change = np.abs(dTmdz_diff) / (
                np.abs(dTmdz.isel(level=slice(None, -1))) + 1e-8
            )

            # 2. Condition: gradient has stabilized (<10% change)
            cond = rel_change < 0.10  # dims: (time, level', ...) level'=nl-1

            # 3. Reverse along level so index 0 = surface-most point, then search
            cond_rev = cond.isel(level=slice(None, None, -1))

            first_idx_rev = cond_rev.argmax(
                dim="level"
            )  # first True from the surface down
            found = cond_rev.any(dim="level")

            first_idx_rev = first_idx_rev.compute()
            found = found.compute()

            # 4. Map the reversed index back to the original (bottom-to-top) level index
            n_levels = cond.sizes["level"]
            first_idx = (
                n_levels - 1
            ) - first_idx_rev  # index into the ORIGINAL cond/level' axis

            # 5. Get the depth at that level (Zm must share the same level' coordinate as cond,
            #    i.e. drop the last original level to match dTmdz_diff, OR the first — see note below)
            Zm_levels = Zm.isel(
                level=slice(None, -1)
            )  # align with dTmdz_diff's level axis
            MLD = xr.where(found, Zm_levels.isel(level=first_idx), np.nan)

            # 6. Fallback when no level ever stabilizes (e.g. fully mixed or noisy profile)
            MLD = MLD.ffill(dim="time").fillna(0)
            MLD = MLD.values

            MLD[0] = 0
            print(MLD)

            def h_conv(B0, Ndeux, offset, t):
                return -np.sqrt(2 * (B0) / Ndeux * t) + offset

            params, _ = curve_fit(h_conv, ds.time.values[1:], MLD[1:])
            y_smooth = h_conv(ds.time, *params)

            fig, ax = plt.subplots(
                1, 1, figsize=(5, 5), constrained_layout=True, dpi=300
            )
            ax.plot(ds.time / Tp, y_smooth, c="k", label=r"$\propto \sqrt{t}+c$")
            print(-y_smooth)
            ax.scatter(ds.time / Tp, MLD, c="b", marker="+", label="data")
            ax.set_xlabel("t/Tp")
            ax.set_ylabel("MLD (m)")
            ax.legend()
            fig.savefig("CSI_MLD.pdf")

            raise Exception
            # for it in range(1, nt):
            #     dTmdz = grad_dir(
            #         dsL2.T_m.isel(time=it),
            #         dsL2.isel(time=it),
            #         grid,
            #         dir="Z",
            #         zvar="level",
            #     )
            #     delta = dTmdz[-1] - dTmdz[-2]
            #     for kz in range(len(ds.level)):
            #         if delta > thrs * dTdz_ini:
            #             MLD[it] = Zm.isel(time=it, level=nl - kz)
            #         else:
            #             MLD[it] = MLD[it - 1]
            print(MLD)
            fig, ax = plt.subplots(
                1, 1, figsize=(5, 5), constrained_layout=True, dpi=300
            )
            ax.plot(ds.time / Tp, MLD)
            ax.set_xlabel("t/Tp")
            ax.set_ylabel("MLD (m)")
            fig.savefig("CSI_MLD.pdf")

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
