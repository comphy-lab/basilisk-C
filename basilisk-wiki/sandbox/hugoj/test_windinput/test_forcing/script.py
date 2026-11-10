import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ====================
# PARAMETERS
CASE = "1D"
EXP = "linear"
PLOT_HPC = True
skip = 1  # plot every "skip"
PLOT_SCATTER = False

g = 1.0
L0 = 1.0
k = 2 * np.pi / L0
lam = 2 * np.pi / k
Re = 20000
NT0 = 10
rho = 1000
if EXP == "stokes":
    ak = 0.35
else:
    ak = 0.01
# ====================

c = np.sqrt(g * k) / k
T0 = 2 * np.pi / np.sqrt(g * k)
nu = c * lam / Re

p0th = 2 * rho * nu * k * c


def E_linwave_nostress(E0, nu, k, t):
    print("\nWave decay for linear wave (theory)")
    print(f"nu={nu}, k={k / np.pi}pi\n")
    return E0 * np.exp(-4 * nu * k**2 * t)


def E_linwave(E0, nu, k, t):
    return E0 * np.exp(-2 * nu * k**2 * t)


def load_cube(filename):
    blocks, run = [], []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):  # end of an experiment
                if run:
                    blocks.append(np.asarray(run, dtype=float))
                    run = []
            elif line.strip():
                run.append([float(x) for x in line.split()])
    if run:
        blocks.append(np.asarray(run, dtype=float))
    if len(blocks) > 1:
        return np.stack(blocks)
    else:
        return np.asarray(blocks)


if CASE == "1D":
    if EXP == "linear":
        data = {
            # "current_noforcing": np.loadtxt("no_forcing/out", skiprows=1),
            "current_noforcing": load_cube("no_forcing/out"),
            "current_exactforced": load_cube("exact_forcing/out"),
            "current_forced": load_cube("linear_wave_wind_input/out"),
        }
    else:  # EXP == "stokes"
        data = {
            "current_noforcing": np.loadtxt("stokes/out", skiprows=1),
            "current_forced": np.loadtxt("hold_stokes/out", skiprows=1),
        }
elif CASE == "2D":
    if EXP == "linear":
        data = {
            "current_noforcing": np.loadtxt("no_forcing2D.cuda/out", skiprows=1),
            "current_noforcing.hpc": np.loadtxt("no_forcing2D.hpc.cpu/out", skiprows=1),
            "current_exactforced": np.loadtxt("exact_forcing2D.cuda/out", skiprows=1),
            "current_exactforced.hpc": np.loadtxt("exact2D.hpc.cpu/out", skiprows=1),
            "current_forced": np.loadtxt(
                "linear_wave_wind_input2D.cuda/out", skiprows=1
            ),
            "current_forced.hpc": np.loadtxt("holding2D.hpc.cpu/out", skiprows=1),
        }
    else:  # EXP == "stokes"
        data = {
            "current_noforcing": np.loadtxt("stokes2D/out", skiprows=1),
            "current_forced": np.loadtxt("hold_stokes2D/out", skiprows=1),
        }

else:
    data = {}
    print(f"case {CASE} is not recognized")


time = {}
ke = {}
gpe = {}
E = {}
p0s = {}
for case in data.keys():
    if case == "current_forced":
        time[case], ke[case], gpe[case], E[case] = {}, {}, {}, {}
        for j in range(data[case].shape[0]):
            ampl = str(data[case][j, 0, -1])
            print(ampl)
            p0s[ampl] = data[case][j, :, -1]
            time[case][ampl] = data[case][j, :, 0]
            ke[case][ampl] = data[case][j, :, 1]
            gpe[case][ampl] = data[case][j, :, 2]
            E[case][ampl] = ke[case][ampl] + gpe[case][ampl]

    else:
        time[case] = data[case][0, :, 0]
        ke[case] = data[case][0, :, 1]
        gpe[case] = data[case][0, :, 2]
        E[case] = ke[case] + gpe[case]


print(f"CASE = {CASE}")
print(f"EXP = {EXP}")
print("T0 = %f s" % T0)
print("nu=%g, ak=%g, k=%g \n" % (nu, ak, k))

if EXP == "linear":
    E0th = 0.5 * rho * g * (ak / k) ** 2
    print(E0th)
    Eth = E_linwave(E0th, nu, k, time["current_noforcing"] * T0)
    E0 = E0th

    print("E0 theory = %g" % E0th)
    print(
        "E0 real = %g (need more resolution to match E0th)" % E["current_noforcing"][0]
    )
    print(
        "Initial energy missmatch = E0th - E0 = %g" % (E0th - E["current_noforcing"][0])
    )
    print("ratio is %g \n" % (E0th / E["current_noforcing"][0]))
    r = (
        (E["current_noforcing"][0] - E["current_noforcing"][-1])
        / E["current_noforcing"][0]
        * 100
    )
    print("decay: E decay rate = %g (percent, after %d T0)" % (r, NT0))
    print(
        "decay theory: E decay rate = %g (percent, after %d T0)"
        % ((1 - np.exp(-2 * nu * k**2 * (NT0 * T0))) * 100, NT0)
    )
    r = (
        (E["current_exactforced"][0] - E["current_exactforced"][-1])
        / E["current_exactforced"][0]
        * 100
    )
    print("exact: E decay = %g (percent, after %d T0)" % (r, NT0))
else:
    E0 = E["current_noforcing"][0]
    Eth = -1  # not defined yet


fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True, dpi=200)
ax.hlines(1.0, 0, 100, colors="gray", alpha=0.7)
ax2 = ax.twinx()
ax2.set_ylabel("p0 / p0th")
p0_colors = mpl.colormaps["Blues"](np.linspace(0.1, 1, len(p0s.keys())))
base_ampl = list(p0s.keys())[0]

if PLOT_SCATTER:
    if EXP == "linear":
        ax.plot(
            time["current_noforcing"],
            Eth / E0,
            color="r",
            label=r"$E(t)=E_0 e^{-2 \nu k^2 t}$",
        )
        ax.scatter(
            time["current_exactforced"][::skip],
            E["current_exactforced"][::skip] / E0,
            c="orange",
            label="exact forcing",
            marker="x",
        )
    ax.plot(
        time["current_forced"][base_ampl],
        E["current_forced"][base_ampl] / E0,
        c="k",
        marker="x",
        label=f"forced, p0={base_ampl}",
    )
    ax.scatter(
        time["current_noforcing"][::skip],
        E["current_noforcing"][::skip] / E0,
        marker="x",
        c="pink",
        label="noforcing",
    )
    if PLOT_HPC and CASE == "2D":
        ax.scatter(
            time["current_noforcing.hpc"][::skip],
            E["current_noforcing.hpc"][::skip] / E0,
            marker="o",
            c="pink",
            label="noforcing hpc",
        )
        ax.scatter(
            time["current_exactforced.hpc"][::skip],
            E["current_exactforced.hpc"][::skip] / E0,
            marker="o",
            c="orange",
            label="exact forcing hpc",
        )
        ax.scatter(
            time["current_forced.hpc"][::skip],
            E["current_forced.hpc"][::skip] / E0,
            marker="o",
            c="purple",
            label="current forced hpc",
        )

    for ck, p0value in enumerate(p0s.keys()):
        ax2.scatter(
            time["current_forced"][p0value][::skip],
            p0s[p0value][::skip] / p0th,
            c=p0_colors[ck],
            marker="x",
            label="p0_i/pth=%f.03f,
        )
else:
    if EXP == "linear":
        ax.plot(
            time["current_noforcing"],
            Eth / E0,
            color="r",
            label=r"$E(t)=E_0 e^{-2 \nu k^2 t}$",
        )
        ax.plot(
            time["current_exactforced"],
            E["current_exactforced"] / E0,
            color="orange",
            label="E (current exact forcing)",
            ls="-",
        )
    ax.plot(
        time["current_noforcing"],
        E["current_noforcing"] / E0,
        color="pink",
        label="E (current noforcing)",
        ls="-",
    )
    ax.plot(
        time["current_forced"][base_ampl],
        E["current_forced"][base_ampl] / E0,
        color="k",
        label="E (current forced)",
        ls="-",
    )
    if PLOT_HPC and CASE == "2D":
        ax.plot(
            time["current_noforcing.hpc"],
            E["current_noforcing.hpc"] / E0,
            color="pink",
            ls="dashdot",
            label="E (current noforcing hpc)",
        )
        ax.plot(
            time["current_exactforced.hpc"][::skip],
            E["current_exactforced.hpc"][::skip] / E0,
            color="orange",
            ls="dashdot",
            label="E (current exact forcing hpc)",
        )
        ax.plot(
            time["current_forced.hpc"][::skip],
            E["current_forced.hpc"][::skip] / E0,
            color="purple",
            ls="dashdot",
            label="E (current forced hpc)",
        )
    for ck, p0value in enumerate(p0s.keys()):
        ax2.plot(
            time["current_forced"][p0value][::skip],
            p0s[p0value][::skip] / p0th,
            color=p0_colors[ck],
            label="p0_i/pth=%0.2f" % (float(p0value) / p0th),
            ls="--",
        )

ax.set_xlabel("t/T0")
ax.set_ylabel("E/E0")
ax.set_xlim((0, NT0))
ax.set_ylim((0.95, 1.05))
ax2.set_ylim((0.0, 2.0))
# ax.set_ylim([1.16e-6, 1.28e-6])
ax.legend(loc="upper left")
ax2.legend(loc="upper right")
ax.grid(axis="y", which="both")
plt.savefig(f"energy_{CASE}_{EXP}.png", dpi=300)
plt.show()
