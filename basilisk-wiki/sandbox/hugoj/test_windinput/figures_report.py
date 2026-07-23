"""
This file plots the results from the data produced by 'run_test_instrumented'

run with
python3 figures_report.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
from pathlib import Path


results = "./results_instrumented/"


def E_linwave(nu, ak, k, t):
    # print("\nWave decay for linear wave (theory)")
    # print(f"nu={nu},ak={ak},k={k / np.pi}pi\n")
    return np.exp(-4 * nu * k**2 * t)


def read_data(path, skiprows=1):
    data = np.loadtxt(path, skiprows=skiprows)
    if Path(path).stat().st_size != 0:
        time = data[:, 0]
        ke = data[:, 1]
        gpe = data[:, 2]
        E = ke + gpe
        safe = 1
    else:
        time = ke = gpe = E = np.zeros(100)
        safe = 0
    return time, E, ke, gpe, safe


def read_data_limiters(path, skiprows=2):
    # time clip_count flux_count lim_zero lim_calls remap_flat remap_clip remap_calls w_clip w_calls

    if Path(path).stat().st_size != 0:
        # log = np.loadtxt(path, skiprows=skiprows)
        # matches a line that is only numbers (int/float/sci notation) separated by whitespace
        data_line_re = re.compile(
            r"^\s*-?\d+\.?\d*(?:[eE][-+]?\d+)?(\s+-?\d+\.?\d*(?:[eE][-+]?\d+)?)*\s*$"
        )
        rows = []
        with open(path, "r") as f:
            for line in f:
                if data_line_re.match(line):
                    rows.append([float(x) for x in line.split()])
        log = np.array(rows)

        ltime = log[:, 0]
        clip_count = log[:, 1]
        flux_count = log[:, 2]
        lim_zero = log[:, 3]
        lim_calls = log[:, 4]
        remap_flat = log[:, 5]
        remap_clip = log[:, 6]
        remap_calls = log[:, 7]
        w_clip = log[:, 8]
        w_calls = log[:, 9]

        CFL_clip = np.where(flux_count == 0, 0, clip_count / flux_count)
        gradient_limiter = np.where(lim_calls == 0, 0, lim_zero / lim_calls)
        remap_clip = np.where(remap_clip == 0, 0, remap_clip / remap_calls)
        remap_flat = np.where(remap_flat == 0, 0, remap_flat / remap_calls)
        w_clip = np.where(w_calls == 0, 0, w_clip / w_calls)
        safe = 1
    else:
        ltime = CFL_clip = gradient_limiter = remap_clip = remap_flat = w_clip = (
            np.zeros(100)
        )
        safe = 0
    return (
        ltime,
        CFL_clip,
        gradient_limiter,
        remap_flat,
        remap_clip,
        w_clip,
        safe,
    )


def parse_params(filename):
    with open(filename, "r") as f:
        line = f.readline()
    data = {}
    for match in re.finditer(r"(\w+)\s*=\s*(-?\d+\.?\d*(?:[eE][-+]?\d+)?)", line):
        key, val = match.groups()
        data[key] = int(val) if re.fullmatch(r"-?\d+", val) else float(val)
    data["T0"] = 2 * np.pi / np.sqrt(9.81 * data["k_"])
    return data


def myplot(file, ax, what, color):
    time, E, _, _, _ = read_data(file)
    ax.plot(time, E / E[0], c=color, label=f"{what}")


def func_plot(folders, thefile, outpath, what):
    timeth = np.linspace(0, 100, 100)

    nplots = len(folders)
    cmap = mpl.colormaps["plasma"]
    colors = cmap(np.linspace(0, 1, nplots))
    fig, ax = plt.subplots(figsize=(5, 5))

    count = 0
    for i, res in enumerate(folders.keys()):
        params = parse_params(folders[res] + "log")
        Eth = E_linwave(
            1 / params["RE"],
            params["ak"],
            params["k_"],
            timeth * params["T0"],
        )
        if what == "Re" or what == "ak":
            ax.plot(timeth, Eth, c=colors[i], label="Lamb (1932)", ls="--")
        else:
            if count == 0:
                ax.plot(timeth, Eth, c=colors[i], label="Lamb (1932)", ls="--")
                count = 1
        myplot(folders[res] + thefile, ax, f"{what}={res}", colors[i])
    ax.set_xlabel("t/T0")
    ax.set_ylabel("Total E/E0")
    ax.set_title(what)
    ax.legend()
    ax.grid()
    fig.savefig(outpath)


def list_folders(path):
    p = Path(path)
    return [f.name for f in p.iterdir() if f.is_dir()]


def plot_one_case(path, file):
    print(path)
    params = parse_params(path + "log")
    time, E, ke, gpe, safe = read_data(path + "out")
    ltime, CFL_clip, gradient_limiter, remap_flat, remap_clip, w_clip, safe = (
        read_data_limiters(path + "log")
    )

    if safe:
        E0 = E[0]
        Eth = E_linwave(
            1 / params["RE"],
            params["ak"],
            params["k_"],
            time * params["T0"],
        )
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot(time, 2 * ke / E0, color="b", label="2*ke")
        ax.plot(time, 2 * gpe / E0, color="g", label="2*gpe")
        ax.plot(time, E / E0, color="k", label="E")
        ax.plot(time, Eth, color="r", label=r"$E(t)=E_0 e^{-4 \nu k^2 t}$")
        ax.set_xlabel("t/T0")
        ax.set_ylabel("E/E0")
        ax.set_title(path)
        ax.set_ylim([0, 1.1])
        ax.grid()
        ax.legend(loc="lower left")
        plt.tight_layout()
        fig.savefig(path + "energy.png", dpi=150)

        ## LIMITER PLOT
        fig, ax = plt.subplots(figsize=(6, 6))
        # -> energy evolution
        ax.plot(time, E / E0, color="k", label="E", ls="--")
        ax.plot(time, Eth, color="k", label=r"$E(t)=E_0 e^{-4 \nu k^2 t}$")
        ax.set_ylabel("E/E0")
        ax.set_ylim([0.0, 1.0])
        ax.legend(loc="upper right")
        # -> limiters call evolution
        ax2 = ax.twinx()
        ax2.plot(ltime / params["T0"], 100 * CFL_clip, c="g", label="CFL clip")
        ax2.plot(
            ltime / params["T0"],
            100 * gradient_limiter,
            c="b",
            label="gradient_limiter",
        )
        ax2.plot(ltime / params["T0"], 100 * remap_clip, c="cyan", label="remap_clip")
        ax2.plot(ltime / params["T0"], 100 * remap_flat, c="yellow", label="remap_flat")
        ax2.plot(ltime / params["T0"], 100 * w_clip, c="r", label="w_clip")
        ax2.set_xlabel("t/T0")
        ax2.set_ylabel("limiter applied / limiter call (%)")
        ax2.set_title("Evolution of Energy and limiters calls (per 10 steps)")
        ax2.legend(loc="lower left")
        ax2.set_ylim([0, 10])
        ax.grid("both")
        plt.savefig(path + "limiters.png")
        plt.close("all")
    else:
        print("Code has crashed, not figure produced. Have a look at")
        print(f"the file {path}{file}")


# ===========================
# Plotting energy evolution
# ===========================
print("> plotting energy for each case")
basepath = "results_instrumented/"
folders = list_folders(basepath)
paths = [basepath + i + "/" for i in folders]
for folder in paths:
    plot_one_case(folder, "out")
# plt.close("all")

# ======================
# N
# ======================
thefile = "out"
folders = {
    "256": results + "N_256/",
    "512": results + "N_512/",
    "1024": results + "base/",
    "2048": results + "N_2048/",
    "4096": results + "N_4096/",
}

outpath = results + "change_of_N.png"
func_plot(folders, thefile, outpath, what="N")
# ========================
# DT
# ========================
thefile = "out"
folders = {
    "0.05": results + "DT_0.05/",
    "0.07": results + "DT_0.07/",
    "0.08": results + "DT_0.08/",
    "0.09": results + "DT_0.09/",
    "0.1": results + "DT_0.1/",
    "auto": results + "base/",
}

outpath = results + "change_of_DT.png"
func_plot(folders, thefile, outpath, what="DT")

# ========================
# CFL_H
# ========================
thefile = "out"
folders = {
    "0.3": results + "CFL_H_0.3/",
    "0.5": results + "CFL_H_0.5/",
    "0.7": results + "CFL_H_0.7/",
    "0.8": results + "CFL_H_0.8/",
    "0.9": results + "CFL_H_0.9/",
    "1.0": results + "base/",
}

outpath = results + "change_of_CFL_H.png"
func_plot(folders, thefile, outpath, what="CFL_H")


# ========================
# CFL
# ========================
thefile = "out"
folders = {
    "0.4": results + "CFL_0.4/",
    "0.5": results + "base/",
    "0.6": results + "CFL_0.6/",
    "0.7": results + "CFL_0.7/",
    "0.8": results + "CFL_0.8/",
    "0.9": results + "CFL_0.9/",
    "1.0": results + "CFL_1.0/",
}

outpath = results + "change_of_CFL.png"
func_plot(folders, thefile, outpath, what="CFL")

# ========================
# theta (minmod2)
# ========================
thefile = "out"
folders = {
    "1.0": results + "theta_1.0/",
    "1.3": results + "base/",
    "1.5": results + "theta_1.5/",
    "2.0": results + "theta_2.0/",
}

outpath = results + "change_of_theta.png"
func_plot(folders, thefile, outpath, what="theta")

# ========================
# theta_H
# ========================
thefile = "out"
folders = {
    "0.5": results + "base/",
    "0.52": results + "theta_h_0.52/",
    "0.54": results + "theta_h_0.54/",
    "0.56": results + "theta_h_0.56/",
}

outpath = results + "change_of_theta_h.png"
func_plot(folders, thefile, outpath, what="theta_h")

# ========================
# RE
# ========================
thefile = "out"
folders = {
    "400.": results + "base/",
    "1000.": results + "RE_1000./",
    "5000.": results + "RE_5000./",
    "10000.": results + "RE_10000./",
    "40000.": results + "RE_40000./",
}

outpath = results + "change_of_Re.png"
func_plot(
    folders,
    thefile,
    outpath,
    what="Re",
)


plt.show()
