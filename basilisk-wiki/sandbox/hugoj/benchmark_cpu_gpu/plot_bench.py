import os
import re
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Configuration ─────────────────────────────────────────────────────────────

PLOT_SPEEDUP = False

# Backends to plot (one subplot each); edit order freely
BACKENDS = ["cpu", "gpu", "cuda", "hip"] # , "hip"

MACHINES = ["sandbox","local","home","JZ","bigfoot"] # "bigfoot_NV100","JZ_NV100"
MACHINES = ["sandbox","local","JZ","bigfoot","ADAmi250x"]
GPU_NAMES = {"sandbox":"RTX4090",
            "local":"RTX1000",
             "home":"RTX4080",
             "bigfoot":"V100",
             "JZ":"A100",
             "ADAmi250x":"mi250x"}
CPU_NAMES = {"local":"i7x8",
             "home":"R7x8",
             "sandbox":"i7x8"}

# Colour palette (one colour per machine, cycles if more machines than colours)
PALETTE = plt.cm.tab10.colors

BASELINE_MACHINE = "sandbox"
BASELINE_BACKEND = "gpu" # == openGL

FILES = {'local':'results_bench_local_turbulence',
         'home':'results_bench_home_turbulence',
         'sandbox':'results_bench_sandbox_turbulence',
         'bigfoot':'result_bench_bigfoot_turbulence_cuda',
         'JZ':'result_bench_JZ_turbulence_cuda',
         'ADAmi250x':'results_bench_ADAmi250x_hip'}

CASE = 'turbulence'

# ── Parse files ───────────────────────────────────────────────────────────────

def check_line_speed(line: str) -> float | None:
    if line.startswith("#"):
        m = re.search(r"([\d.eE+\-]+)\s+points\.step/s", line)
        if m:
            return float(m.group(1))
        else:
            return None

def parse_file(filepath: str) -> dict:
    results = {}
    filename = filepath.split('/')[-1].split('_')
    with open(filepath, 'r') as f:
        lines = f.readlines()

    current_backend = None
    current_resolution = None

    # If the backend comes from the filename, set it once before the loop
    if filename[2] not in ['home', 'local', 'sandbox']:
        current_backend = filename[-1]
        if current_backend not in BACKENDS:
            raise Exception(f'Backend {current_backend} not recognized')
        results[current_backend] = {}

    for line in lines:
        line = line.strip()

        # Match backend from file content
        if filename[2] in ['home', 'local', 'sandbox']:
            if line in BACKENDS:
                current_backend = line
                if current_backend not in results:
                    results[current_backend] = {}  # ← only initialise once

        # Match resolution (standalone integer)
        if re.fullmatch(r'\d+', line):
            current_resolution = int(line)

        # Match the summary comment line with speed
        elif line.startswith('#') and 'points.step/s' in line:
            match = re.search(r'([\d.e+]+)\s+points\.step/s', line)
            if match and current_backend and current_resolution is not None:
                results[current_backend][current_resolution] = float(match.group(1))

    return results

def backend_to_axe(backend,axes):
    if backend=='cpu':
        return axes[0]
    elif backend=='gpu':
        return axes[1]
    elif backend=='cuda':
        return axes[2]
    elif backend=='hip':
        return axes[3]
    else:
        raise Exception(f'Backend {backend} not recognized')

def get_all_res(data):
    all_res = []
    for machine in data.keys():
        for backend in data[machine].keys():
            res=list(data[machine][backend].keys())
            if len(res)>len(all_res):
                all_res = res
    return all_res


data = {}
for machine in MACHINES:
    data[machine] = parse_file(FILES[machine])

if PLOT_SPEEDUP:
    nrow=2
    width=9
    height=6
else:
    nrow=1
    width=10
    height=3.5
fig, axes = plt.subplots(nrow, len(BACKENDS), figsize=(width, height), constrained_layout=True)
machine_colors = {m: PALETTE[i % len(PALETTE)] for i, m in enumerate(MACHINES)}
baseline = data[BASELINE_MACHINE][BASELINE_BACKEND]
baseline_speeds = list(baseline.values())

get_all_res(data)

all_res = get_all_res(data) #[float(value) for value in baseline.keys()]
n_m   = len(data.keys())
x     = np.arange(len(all_res))
width = 0.8 / n_m

maxspeed = 0.

if PLOT_SPEEDUP:
    base_ax_speed = axes[0,:]
    base_ax_speedup = axes[1,:]
else:
    base_ax_speed = axes

# log grid
for k in range(len(base_ax_speed)):
    base_ax_speed[k].set_axisbelow(True)
    base_ax_speed[k].grid(axis='y',which='both')

for i, machine in enumerate(data.keys()):
    for j, backend in enumerate(data[machine].keys()):
        print(machine, backend)
        ax_speed = backend_to_axe(backend,base_ax_speed) #axes[0, j]
        ax_speed.set_yscale('log', base=10)
        if PLOT_SPEEDUP:
            ax_speedup = backend_to_axe(backend,base_ax_speedup) #axes[1, j]
        # ── speed ──
        speeds = list(data[machine][backend].values())
        offset  = (i - n_m / 2 + 0.5) * width
        x = np.arange(len(data[machine][backend].keys()))
        if backend == "cpu":
            label=CPU_NAMES[machine]
        else:
            label=GPU_NAMES[machine]
        # The plot
        bars    = ax_speed.bar(
            x + offset, speeds, width,
            label=label,
            color=machine_colors[machine], linewidth=0.5,
        )
        ax_speed.set_title(backend)
        for speed in speeds:
            if speed>maxspeed:
                maxspeed = speed
        ax_speed.set_xticks(x)
        ax_speed.set_xticklabels([str(int(r)) for r in all_res[:len(x)]], rotation=45)
        if PLOT_SPEEDUP:
            # ── speedup ──
            speedups = np.zeros(len(speeds))
            for k in range(len(speeds)):
                if k>=len(baseline_speeds):
                    speedups[k]=0 # no data for baseline at this resolution
                else:
                    speedups[k] = speeds[k]/baseline_speeds[k]
            offset  = (i - n_m / 2 + 0.5) * width
            x = np.arange(len(data[machine][backend].keys()))
            if not (machine == BASELINE_MACHINE and backend == BASELINE_BACKEND): 
                bars    = ax_speedup.bar(
                    x + offset, speedups, width,
                    label=machine,
                    color=machine_colors[machine], edgecolor="white", linewidth=0.5,
                )
            if len(x)==len(all_res):
                ax_speedup.set_xticks(x)
                ax_speedup.set_xticklabels([str(int(r)) for r in all_res], rotation=45)

for backend in range(len(BACKENDS)):
    base_ax_speed[backend].set_ylim([1e5, 1.1*maxspeed])
    if PLOT_SPEEDUP:
        base_ax_speedup[backend].set_ylim([0,1.2])
base_ax_speed[0].set_ylabel('speed (points.step/s)')
base_ax_speed[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
          fancybox=True, shadow=True, ncol=2)
base_ax_speed[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
          fancybox=True, shadow=True, ncol=2)
base_ax_speed[3].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
          fancybox=True, shadow=True, ncol=2)
base_ax_speed[2].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
          fancybox=True, shadow=True, ncol=2)
if PLOT_SPEEDUP:
    base_ax_speedup[1].set_ylabel('speedup (vs RTX4090 OpenGL)')

fig.savefig(f"bench_{CASE}.png", dpi=150) #, bbox_inches="tight")
plt.show()


