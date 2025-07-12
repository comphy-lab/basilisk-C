#!/usr/bin/env python3
"""
Batch process multiple particle_*.dat datasets and produce smoothed ⟨vy⟩ and ⟨ρ⟩ plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from tqdm import tqdm

# ────────────── Settings ──────────────
data_dirs = [
    "plumes_y_periodic_buoyancy_NP_01_B_0.0",
    "plumes_y_periodic_buoyancy_NP_32_B_0.0",
    "plumes_y_periodic_buoyancy_NP_01_B_2.5",
    "plumes_y_periodic_buoyancy_NP_32_B_2.5"
]

base_data_path = "/home/fpicella/wiki/sandbox/fpicella/microswimmer_tHree_forces"
output_dir = os.path.join(base_data_path, "tHree_PLOTS")
os.makedirs(output_dir, exist_ok=True)

# ────────────── Plot appearance settings ──────────────
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

# Domain extents
x_min, x_max = -64, 64
y_min, y_max = -64, 64
x_bins = 256
y_bins = 256
dx = (x_max - x_min) / x_bins
dy = (y_max - y_min) / y_bins

# KDE kernel
kernel_radius_physical = 1.0
sigma_physical = kernel_radius_physical / 0.3
ker_rad_x = 10 * int(np.ceil(kernel_radius_physical / dx))
ker_rad_y = 10 * int(np.ceil(kernel_radius_physical / dy))
ker_size_x = 2 * ker_rad_x + 1
ker_size_y = 2 * ker_rad_y + 1
xk = np.linspace(-ker_rad_x * dx, ker_rad_x * dx, ker_size_x)
yk = np.linspace(-ker_rad_y * dy, ker_rad_y * dy, ker_size_y)
Xk, Yk = np.meshgrid(xk, yk)
gaussian_kernel = np.exp(-(Xk**2 + Yk**2) / (2 * sigma_physical**2)).astype(np.float32)

sliding_window = 500  # For smoothing

def sliding_average(data, window):
    if len(data) < window:
        return np.array(data)
    return np.convolve(data, np.ones(window)/window, mode='valid')

# ────────────── Main processing loop ──────────────
for folder in data_dirs:
    full_path = os.path.join(base_data_path, folder)
    file_pattern = os.path.join(full_path, "particle_*.dat")
    particle_files = sorted(glob.glob(file_pattern))
    if not particle_files:
        print(f"[Warning] No files found in: {folder}")
        continue

    print(f"Processing: {folder}")
    all_particle_data = [np.loadtxt(f) for f in particle_files]
    num_timesteps = all_particle_data[0].shape[0]

    avg_vy_per_timestep = []
    avg_kde_per_timestep = []
    timesteps = []

    for t_idx in tqdm(range(num_timesteps), desc=f"{folder}"):
        presence_map = np.zeros((y_bins, x_bins), dtype=np.float32)
        x_list, y_list, vy_list = [], [], []

        for pdata in all_particle_data:
            x, y  = pdata[t_idx, 1], pdata[t_idx, 2]
            vy    = pdata[t_idx, 5]

            if not (x_min <= x <= x_max and y_min <= y <= y_max):
                continue

            x_list.append(x)
            y_list.append(y)
            vy_list.append(vy)

            xi = int((x - x_min) / dx)
            yi = int((y - y_min) / dy)

            x0, y0 = xi - ker_rad_x, yi - ker_rad_y
            x1, y1 = x0 + ker_size_x, y0 + ker_size_y

            x_start = max(x0, 0)
            y_start = max(y0, 0)
            x_end   = min(x1, x_bins)
            y_end   = min(y1, y_bins)

            kx_start = x_start - x0
            ky_start = y_start - y0
            kx_end   = kx_start + (x_end - x_start)
            ky_end   = ky_start + (y_end - y_start)

            presence_map[y_start:y_end, x_start:x_end] += \
                gaussian_kernel[ky_start:ky_end, kx_start:kx_end]

        kde_values = []
        for x, y in zip(x_list, y_list):
            xi = int((x - x_min) / dx)
            yi = int((y - y_min) / dy)
            if 0 <= xi < x_bins and 0 <= yi < y_bins:
                kde_values.append(presence_map[yi, xi])

        kde_values = np.array(kde_values)
        vy_list = np.array(vy_list)

        avg_vy = np.mean(vy_list) if vy_list.size > 0 else 0.0
        avg_kde = np.mean(kde_values) if kde_values.size > 0 else 0.0

        avg_vy_per_timestep.append(avg_vy)
        avg_kde_per_timestep.append(avg_kde)
        timesteps.append(t_idx)

    avg_vy_smoothed  = sliding_average(avg_vy_per_timestep, sliding_window)
    avg_kde_smoothed = sliding_average(avg_kde_per_timestep, sliding_window)
    smoothed_timesteps = timesteps[sliding_window - 1:]

    # ────────────── Plotting ──────────────
    fig, ax = plt.subplots(2, 1, figsize=(6, 2.5), sharex=True)

    ax[0].plot(timesteps, avg_vy_per_timestep, color='crimson', alpha=0.25, linewidth=2)
    ax[0].plot(smoothed_timesteps, avg_vy_smoothed, color='crimson', linewidth=4)
    ax[0].axhline(0, color='black', linestyle='--', linewidth=1.5)
    ax[0].set_ylabel(r"$\left<U_y\right>, \overline{\left< U_y \right>}$")
    ax[0].grid(True)

    ax[1].plot(timesteps, avg_kde_per_timestep, color='teal', alpha=0.25, linewidth=2)
    ax[1].plot(smoothed_timesteps, avg_kde_smoothed, color='teal', linewidth=4)
    ax[1].set_ylabel(r"$\left<\rho\right>, \overline{\left<\rho\right>} $")
    ax[1].set_xlabel("Time")
    ax[1].grid(True)

    plt.tight_layout(pad=0.2)
    plt.tight_layout(pad=0.2)

    output_filename = f"summary_{folder}"
    output_path_svg = os.path.join(output_dir, output_filename + ".svg")
    output_path_pdf = os.path.join(output_dir, output_filename + ".pdf")
    
    plt.savefig(output_path_svg, bbox_inches='tight')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    plt.close()
    
    print(f"[✓] Saved: {output_path_svg}")
    print(f"[✓] Saved: {output_path_pdf}")
