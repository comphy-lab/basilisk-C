import numpy as np
import matplotlib.pyplot as plt

CASE = "1D"

g = 1.0
L0 = 1.0
ak = 0.01
k = 2 * np.pi / L0
lam = 2 * np.pi / k
Re = 20000
c = np.sqrt(g * k) / k
NT0 = 10

nu = c * lam / Re
T0 = 2 * np.pi / np.sqrt(g * k)

print("T0 = %f" % T0)


def E_linwave_nostress(E0, nu, k, t):
    print("\nWave decay for linear wave (theory)")
    print(f"nu={nu}, k={k / np.pi}pi\n")
    return E0 * np.exp(-4 * nu * k**2 * t)


def E_linwave(E0, nu, k, t):
    return E0 * np.exp(-2 * nu * k**2 * t)


if CASE == "1D":
    data = {
        "current_noforcing": np.loadtxt("no_forcing/out", skiprows=1),
        "current_exactforced": np.loadtxt("exact_forcing/out", skiprows=1),
        "current_forced": np.loadtxt("linear_wave_wind_input/out", skiprows=1),
    }
elif CASE == "2D":
    data = {
        "current_noforcing": np.loadtxt("no_forcing2D.cuda/out", skiprows=1),
        "current_exactforced": np.loadtxt("exact_forcing2D.cuda/out", skiprows=1),
        "current_forced": np.loadtxt("linear_wave_wind_input2D.cuda/out", skiprows=1),
    }
else:
    print(f"case {CASE} is not recognized")


time = {}
ke = {}
gpe = {}
E = {}

for case in data.keys():
    time[case] = data[case][:, 0]
    ke[case] = data[case][:, 1]
    gpe[case] = data[case][:, 2]
    E[case] = ke[case] + gpe[case]

E0th = 0.5 * g * (ak / k) ** 2
print("\ntheoretical E0 = %g m3/s2" % E0th)
print("initial energy for current sim E0 = %g m3/s2" % E["current_noforcing"][0])
print("ratio is %g \n" % (E0th / E["current_noforcing"][0]))
Eth = E_linwave(E0th, nu, k, time["current_noforcing"] * T0)
E0 = E0th

print("Target for E=cst for this set of parameters:")
print("nu=%g, ak=%g, k=%g" % (nu, ak, k))
print("E0 theory = %g" % E0th)
print("E0 real = %g (need more resolution to match E0th)" % E["current_noforcing"][0])
print("Initial energy missmatch = E0th - E0 = %g" % (E0th - E["current_noforcing"][0]))
print(
    "noforcing: miss match Eth-E=%g (at %d T0) "
    % (Eth[-1] - E["current_noforcing"][-1], NT0)
)
print(
    "exact forcing: miss match Eth-E=%g (at %d T0) "
    % (Eth[-1] - E["current_exactforced"][-1], NT0)
)

fig, ax = plt.subplots(figsize=(8, 6))

# ax.plot(
#     time["current_noforcing"],
#     2 * gpe["current_noforcing"] / E0,
#     color="g",
#     label="2Ep (current)",
#     alpha=0.5,
# )
# ax.plot(
#     time["current_noforcing"],
#     2 * ke["current_noforcing"] / E0,
#     color="b",
#     label="2Ek (current_noforcing)",
#     alpha=0.5,
# )
ax.hlines(E0th / E0, 0, 100, colors="gray", alpha=0.7)
ax.semilogy(
    time["current_noforcing"], Eth / E0, color="r", label=r"$E(t)=E_0 e^{-2 \nu k^2 t}$"
)
ax.semilogy(
    time["current_noforcing"],
    E["current_noforcing"] / E0,
    color="pink",
    label="E (current noforcing)",
    ls="--",
)
ax.semilogy(
    time["current_forced"],
    E["current_forced"] / E0,
    color="purple",
    label="E (current forcing)",
    ls="--",
)
ax.semilogy(
    time["current_exactforced"],
    E["current_exactforced"] / E0,
    color="orange",
    label="E (current exact forcing)",
    ls="--",
)

ax.set_xlabel("t/T0")
ax.set_ylabel("E/E0")
ax.set_xlim([0, NT0])
# ax.set_ylim([1.16e-6, 1.28e-6])
ax.legend()
ax.grid(axis="y", which="both")
plt.savefig(f"energy_{CASE}.png", dpi=300)
plt.show()
