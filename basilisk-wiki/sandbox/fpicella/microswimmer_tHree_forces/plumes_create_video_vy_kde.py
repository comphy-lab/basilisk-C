#!/usr/bin/env python3
"""
Render a frame-by-frame Gaussian-splat density map of particle positions,
with two subfigures per frame:
 - Left: vertical velocity (vy)
 - Right: local density via KDE (sampled from presence map)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
import glob
import os
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap

# ────────────────────── Custom Colormap Truncation ─────────────────────
def truncate_colormap(cmap, minval=0.25, maxval=0.75, n=256):
    new_cmap = LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval},{maxval})",
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

base_cmap = plt.get_cmap('seismic')
trunc_cmap = truncate_colormap(base_cmap, 0.25, 0.75)

# ───────────────────────────── Configuration ─────────────────────────────
folder        = "plumes_production"
file_pattern  = os.path.join(folder, "particle_*.dat")
frame_dir     = os.path.join(folder, "frames2D_vy_kde")
video_path    = os.path.join(folder, "particle_vy_kde.mp4")
os.makedirs(frame_dir, exist_ok=True)

# Domain extents
x_min, x_max = -64, 64
y_min, y_max = -64, 64

# Grid resolution
x_bins = 2 * 128
y_bins = 2 * 128
dx = (x_max - x_min) / x_bins
dy = (y_max - y_min) / y_bins

# Gaussian kernel in physical units
kernel_radius_physical = 1.0
sigma_physical = kernel_radius_physical / 0.3

# Expanded kernel for smoother blobs
ker_rad_x = 10 * int(np.ceil(kernel_radius_physical / dx))
ker_rad_y = 10 * int(np.ceil(kernel_radius_physical / dy))
ker_size_x = 2 * ker_rad_x + 1
ker_size_y = 2 * ker_rad_y + 1

# Gaussian kernel array
xk = np.linspace(-ker_rad_x * dx, ker_rad_x * dx, ker_size_x)
yk = np.linspace(-ker_rad_y * dy, ker_rad_y * dy, ker_size_y)
Xk, Yk = np.meshgrid(xk, yk)
gaussian_kernel = np.exp(-(Xk**2 + Yk**2) / (2 * sigma_physical**2)).astype(np.float32)

# ─────────────────────── Load particle data ────────────────────────
particle_files = sorted(glob.glob(file_pattern))
if not particle_files:
    raise FileNotFoundError("No particle_*.dat files found.")

all_particle_data = [np.loadtxt(f) for f in particle_files]
time = all_particle_data[0][:, 0]
num_steps = len(time)

# ────────────────────────── Frame generation ──────────────────────────
print("Rendering frames with vy and local KDE-based density overlay...")
for t_idx in tqdm(range(num_steps), desc="Rendering frames"):
    presence_map = np.zeros((y_bins, x_bins), dtype=np.float32)
    x_list, y_list, vy_list = [], [], []

    for p in all_particle_data:
        x, y  = p[t_idx, 1], p[t_idx, 2]
        vy    = p[t_idx, 5]

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

    # ────────────── Sample KDE value at each particle position ──────────────
    kde_values = []
    for x, y in zip(x_list, y_list):
        xi = int((x - x_min) / dx)
        yi = int((y - y_min) / dy)
        if 0 <= xi < x_bins and 0 <= yi < y_bins:
            kde_values.append(presence_map[yi, xi])
        else:
            kde_values.append(0.0)

    kde_values = np.array(kde_values)

    # ─────────────── Render side-by-side subplots ───────────────
    fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(16, 8), sharex=True, sharey=True)
    extent = [x_min, x_max, y_min, y_max]

    # Left: vy
    ax1.imshow(presence_map, origin='lower', extent=extent, cmap='gray_r', aspect='equal')
    sc1 = ax1.scatter(x_list, y_list, c=vy_list, cmap=trunc_cmap,
                      s=200, vmin=-0.5, vmax=0.5, alpha=1.0, edgecolors='none')
    ax1.set_title("Vertical Velocity")
    ax1.set_xticks([]); ax1.set_yticks([])
    fig.colorbar(sc1, ax=ax1, shrink=0.75)

    # Right: KDE value (local density)
    ax2.imshow(presence_map, origin='lower', extent=extent, cmap='gray_r', aspect='equal')
    sc2 = ax2.scatter(x_list, y_list, c=kde_values, cmap='cividis',
                      s=200, vmin=0.0, vmax=2.5,
                      alpha=1.0, edgecolors='none')
    ax2.set_title("Local Density (KDE)")
    ax2.set_xticks([]); ax2.set_yticks([])
    fig.colorbar(sc2, ax=ax2, shrink=0.75)

    fig.tight_layout()

    frame_path = os.path.join(frame_dir, f"frame_{t_idx:04d}.png")
    fig.savefig(frame_path, dpi=100)
    plt.close(fig)

print("All frames saved.")

# ───────────────────────────── Video Compilation ─────────────────────────────
print("Creating video …")
ffmpeg_cmd = (
    f"ffmpeg -y -framerate 30 -i {frame_dir}/frame_%04d.png "
    f"-c:v libx264 -pix_fmt yuv420p {video_path}"
)
os.system(ffmpeg_cmd)
print(f"Video saved: {video_path}")

