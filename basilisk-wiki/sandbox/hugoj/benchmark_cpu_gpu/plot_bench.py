import os
import re
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Configuration ─────────────────────────────────────────────────────────────

# Backends to plot (one subplot each); edit order freely
BACKENDS = ["cpu", "gpu", "cuda", "hip"]

MACHINES = ["sandbox","local","home","JZ","bigfoot"] # "bigfoot_NV100","JZ_NV100"
GPU_NAMES = {"sandbox":"RTX4090",
            "local":"RTX1000",
             "home":"RTX4080",
             "bigfoot":"V100",
             "JZ":"A100"}
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
         'JZ':'result_bench_JZ_turbulence_cuda'}

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
        if current_backend not in ('cpu', 'gpu', 'cuda', 'hip'):
            raise Exception(f'Backend {current_backend} not recognized')
        results[current_backend] = {}

    for line in lines:
        line = line.strip()

        # Match backend from file content
        if filename[2] in ['home', 'local', 'sandbox']:
            if line in ('cpu', 'gpu', 'cuda', 'hip'):
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

data = {}
for machine in MACHINES:
    data[machine] = parse_file(FILES[machine])

print(data["sandbox"])

fig, axes = plt.subplots(2, len(BACKENDS), figsize=(9, 6), constrained_layout=True)
machine_colors = {m: PALETTE[i % len(PALETTE)] for i, m in enumerate(MACHINES)}
baseline = data[BASELINE_MACHINE][BASELINE_BACKEND]
baseline_speeds = list(baseline.values())

all_res = [float(value) for value in baseline.keys()]
n_m   = len(data.keys())
x     = np.arange(len(all_res))
width = 0.8 / n_m

maxspeed = 0.
for i, machine in enumerate(data.keys()):
    for j, backend in enumerate(data[machine].keys()):
        ax_speed = backend_to_axe(backend,axes[0,:]) #axes[0, j]
        ax_speedup = backend_to_axe(backend,axes[1,:]) #axes[1, j]
        # ── speed ──
        speeds = list(data[machine][backend].values())
        offset  = (i - n_m / 2 + 0.5) * width
        x = np.arange(len(data[machine][backend].keys()))
        if backend == "cpu":
            label=CPU_NAMES[machine]
        else:
            label=GPU_NAMES[machine]
        bars    = ax_speed.bar(
            x + offset, speeds, width,
            label=label,
            color=machine_colors[machine], edgecolor="white", linewidth=0.5,
        )
        ax_speed.set_title(backend)
        for speed in speeds:
            if speed>maxspeed:
                maxspeed = speed
        if len(x)==len(all_res):
            ax_speed.set_xticks(x)
            ax_speed.set_xticklabels([str(int(r)) for r in all_res], rotation=45)
        # ── speedup ──
        speedups = np.zeros(len(speeds))
        for k in range(len(speeds)):
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
    axes[0,backend].set_ylim([0, 1.1*maxspeed])
    axes[1,backend].set_ylim([0,1.1])
axes[0,0].set_ylabel('speed (points.step/s)')
axes[0,0].legend(fontsize=8, framealpha=0.8)
axes[0,1].legend(fontsize=8, framealpha=0.8)
axes[0,2].legend(fontsize=8, framealpha=0.8)
axes[0,3].legend(fontsize=8, framealpha=0.8)
axes[1,0].set_ylabel('speedup (vs RTX4090 OpenGL)')

fig.savefig(f"bench_{CASE}.png", dpi=150) #, bbox_inches="tight")
plt.show()


