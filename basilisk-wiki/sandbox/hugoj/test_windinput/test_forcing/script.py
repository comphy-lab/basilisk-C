import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ====================
# PARAMETERS
CASE = "1D"
EXP = "linear"
PLOT_HPC = True
PLOT_KE_GPE = False

skip = 1  # plot every "skip"
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
            "current_noforcing": load_cube("no_forcing/out"),
            "current_exactforced": load_cube("exact_forcing/out"),
            "current_forced": load_cube("test_p0s/out"),
        }
    else:  # EXP == "stokes"
        data = {
            "current_noforcing": np.loadtxt("stokes/out", skiprows=1),
            "current_forced": np.loadtxt("hold_stokes/out", skiprows=1),
        }
elif CASE == "2D":
    if EXP == "linear":
        data = {
            # "current_noforcing": load_cube("no_forcing2D.cuda/out"),
            # "current_noforcing.hpc": load_cube("no_forcing2D.hpc.cpu/out"),
            # "current_exactforced": load_cube("exact_forcing2D.cuda/out"),
            # "current_exactforced.hpc": load_cube("exact2D.hpc.cpu/out"),
            "current_forced": load_cube("linear_wave_wind_input2D.cuda/out"),
            # "current_forced.hpc": load_cube("holding2D.hpc.cpu/out"),
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

    Eth = E_linwave(E0th, nu, k, time["current_noforcing"] * T0)
    E0 = E0th
    if False:
        print(E0th)
        print("E0 theory = %g" % E0th)
        print(
            "E0 real = %g (need more resolution to match E0th)"
            % E["current_noforcing"][0]
        )
        print(
            "Initial energy missmatch = E0th - E0 = %g"
            % (E0th - E["current_noforcing"][0])
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


if EXP == "linear":
    # subregion of the original image
    x1, x2, y1, y2 = 0.0, 10.0, 0.996, 1.001
    # position and size in data coordinates
    x0, y0, width, height = (
        7.5,
        0.98,
        2.0,
        0.015,
    )

fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True, dpi=300)

p0_colors = mpl.colormaps["brg"](np.linspace(0.1, 1, len(p0s.keys())))
base_ampl = list(p0s.keys())[0]

ax2 = ax.inset_axes(
    (x0, y0, width, height),
    xlim=[x1, x2],
    ylim=[y1, y2],
    transform=ax.transData,
    facecolor="white",
)
for axe in [ax, ax2]:
    axe.hlines(1.0, 0, 100, colors="gray", alpha=0.7)
    if EXP == "linear":
        axe.plot(
            time["current_noforcing"],
            Eth / E0,
            color="r",
            label=r"$E(t)=E_0 e^{-2 \nu k^2 t}$",
        )
        axe.plot(
            time["current_exactforced"],
            E["current_exactforced"] / E0,
            color="orange",
            label="forced, exact",
            ls="-",
        )
    axe.plot(
        time["current_noforcing"],
        E["current_noforcing"] / E0,
        color="pink",
        label="decay",
        ls="-",
    )
    # axe.plot(
    #     time["current_forced"][base_ampl],
    #     E["current_forced"][base_ampl] / E0,
    #     color="k",
    #     label="forced, p0/p0th=%0.1f" % (float(base_ampl) / p0th),
    #     ls="-",
    # )
    if PLOT_KE_GPE:
        axe.plot(
            time["current_forced"][base_ampl],
            2 * ke["current_forced"][base_ampl] / E0,
            color="b",
            label="2*Ek (forced)",
            ls="-",
        )
        axe.plot(
            time["current_forced"][base_ampl],
            2 * gpe["current_forced"][base_ampl] / E0,
            color="g",
            label="2*Ep (forced)",
            ls="-",
        )

    if PLOT_HPC and CASE == "2D":
        axe.plot(
            time["current_noforcing.hpc"],
            E["current_noforcing.hpc"] / E0,
            color="pink",
            ls="dashdot",
            label="E (current noforcing hpc)",
        )
        axe.plot(
            time["current_exactforced.hpc"][::skip],
            E["current_exactforced.hpc"][::skip] / E0,
            color="orange",
            ls="dashdot",
            label="E (current exact forcing hpc)",
        )
        axe.plot(
            time["current_forced.hpc"][::skip],
            E["current_forced.hpc"][::skip] / E0,
            color="purple",
            ls="dashdot",
            label="E (current forced hpc)",
        )
    for ck, p0value in enumerate(p0s.keys()):
        axe.plot(
            time["current_forced"][p0value][::skip],
            E["current_forced"][p0value][::skip] / E0,
            color=p0_colors[ck],
            label="forced, p0/p0th=%0.1f" % (abs(float(p0value) / p0th)),
            ls="--",
        )


ax.set_xlabel("t/T0")
ax.set_ylabel("E/E0")
ax.set_xlim((0, 10))  # NT0
ax.set_ylim((0.96, 1.005))
# ax.set_ylim([1.16e-6, 1.28e-6])
ax.legend(loc="lower left")
ax.grid(axis="both", which="both")
plt.savefig(f"energy_{CASE}_{EXP}.png", dpi=300)
plt.show()
