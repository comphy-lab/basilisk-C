import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os

g = 1.0
L0 = 1.0
ak = 0.35
k = 2 * np.pi / L0
lam = 2 * np.pi / k
Re = 20000
c = np.sqrt(g * k) / k
NT0 = 10

nu = c * lam / Re
T0 = 2 * np.pi / np.sqrt(g * k)

# ===========================
# compute energy evolution
# ===========================
if True:
    data = {
        "noforcing": np.loadtxt("stokes/out", skiprows=1),
        # "current_exactforced": np.loadtxt("../exact_forcing/out", skiprows=1),
        # "current_forced": np.loadtxt("../linear_wave_wind_input/out", skiprows=1),
    }
    time = {}
    ke = {}
    gpe = {}
    E = {}
    E0 = 1
    for case in data.keys():
        time[case] = data[case][:, 0]
        ke[case] = data[case][:, 1]
        gpe[case] = data[case][:, 2]
        E[case] = ke[case] + gpe[case]

    print("nu=%g, ak=%g, k=%g" % (nu, ak, k))
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.semilogy(
        time["noforcing"],
        E["noforcing"] / E0,
        color="pink",
        label="E (no forcing)",
    )
    # ax.semilogy(time["current_noforcing"], Eth / E0, color="r", label=r"$E(t)=E_0 e^{-2 \nu k^2 t}$")
    ax.set_xlabel("t/T0")
    ax.set_ylabel("E")
    # ax.set_xlim([0, 10])
    # ax.set_ylim([1.2e-6, 1.28e-6])
    ax.legend()
    ax.grid(axis="y", which="both")
    plt.savefig("energy.pdf", dpi=300)
    # plt.show()

# ==================================
# make a movie of surface elevation
# ==================================
if True:
    file = "stokes/out.nc"
    ds = xr.open_dataset(file)
    ds = ds.isel(y=0)
    fps = 30
    speed = 1
    t0 = float(ds.time.values[0])
    t1 = float(ds.time.values[-1])
    if t1 <= t0:
        raise ValueError(f"t_end ({t1}) must be greater than t_start ({t0})")

    duration_video = (t1 - t0) / speed
    nframes = int(np.round(duration_video * fps))
    if nframes < 1:
        raise ValueError("computed 0 frames — check fps / speed_factor / time range")
    target_times = t0 + np.arange(nframes) * (speed / fps)

    os.system("rm tmp/*")
    for it in range(len(target_times)):
        # print("t=%f" % target_times[it])
        pct = 100 * (it) / nframes
        print(f"  frame {it}/{nframes} ({pct:.0f}%)", end="\r", flush=True)
        fig, ax = plt.subplots(figsize=(6, 3), constrained_layout=True, dpi=100)
        ax.plot(ds.x, ds.eta.sel(time=target_times[it], method="nearest"))
        ax.set_xlabel("x (m)")
        ax.set_ylabel("z (m)")
        ax.set_ylim([-0.1, 0.1])
        ax.set_xlim([-L0 / 2, L0 / 2])
        ax.set_aspect(1)
        ax.set_title("t/T0 = %f" % (target_times[it] / T0))
        plt.savefig("tmp/eta_profile_t%03d.png" % it)
        plt.close(fig)

    os.system("rm movie.mp4")
    os.system(
        "ffmpeg -framerate 30 -i ./tmp/eta_profile_t%03d.png -c:v libx264 -pix_fmt yuv420p -r 30 movie.mp4"
    )
