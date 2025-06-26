#!/usr/bin/env python3
"""
Render a frame-by-frame Gaussian-splat density map of particle positions
and stitch the frames into an MP4 video.

Kernel size is defined in *physical* units, not pixels.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')          # non-GUI backend for faster batch rendering
import matplotlib.pyplot as plt
import glob
import os
from tqdm import tqdm          # progress bar

# ─────────────────────────── Configuration ────────────────────────────
folder        = "plumes_production"
file_pattern  = os.path.join(folder, "particle_*.dat")
frame_dir     = os.path.join(folder, "frames2D")
video_path    = os.path.join(folder, "particle_presence_2D.mp4")
os.makedirs(frame_dir, exist_ok=True)

# Domain extents (physical units)
x_min, x_max  = -64, 64
y_min, y_max  = -64, 64

# Output grid resolution
x_bins        = 2*128
y_bins        = 2*128
dx            = (x_max - x_min) / x_bins
dy            = (y_max - y_min) / y_bins

# Gaussian kernel (physical units)
kernel_radius_physical = 1.0          # “radius” of kernel in physical units
sigma_physical         = kernel_radius_physical / 2.0  # std dev in same units
sigma_physical         = kernel_radius_physical / .3  # std dev in same units

# ────────────────────── Precompute Gaussian kernel ────────────────────
# Convert physical radius → pixel radius for each axis
#ker_rad_x = int(np.ceil(kernel_radius_physical / dx))
#ker_rad_y = int(np.ceil(kernel_radius_physical / dy))
# Made bigger, so to make less sharp blobs.
ker_rad_x = 10*int(np.ceil(kernel_radius_physical / dx))
ker_rad_y = 10*int(np.ceil(kernel_radius_physical / dy))
ker_size_x = 2 * ker_rad_x + 1
ker_size_y = 2 * ker_rad_y + 1

# Build anisotropic kernel in physical coordinates
xk = np.linspace(-ker_rad_x * dx, ker_rad_x * dx, ker_size_x)
yk = np.linspace(-ker_rad_y * dy, ker_rad_y * dy, ker_size_y)
Xk, Yk = np.meshgrid(xk, yk)
gaussian_kernel = np.exp(-(Xk**2 + Yk**2) / (2 * sigma_physical**2)).astype(np.float32)

# ───────────────────────── Load particle data ─────────────────────────
particle_files = sorted(glob.glob(file_pattern))
if not particle_files:
    raise FileNotFoundError("No particle_*.dat files found.")

all_particle_data = [np.loadtxt(f) for f in particle_files]
time = all_particle_data[0][:, 0]
num_steps = len(time)

# ────────────────────────── Frame generation ──────────────────────────
print("Rendering frames (physical-size Gaussian placement)…")
for t_idx in tqdm(range(num_steps), desc="Rendering frames"):
    presence_map = np.zeros((y_bins, x_bins), dtype=np.float32)

    for p in all_particle_data:
        x, y = p[t_idx, 1], p[t_idx, 2]

        # Skip particles outside the domain
        if not (x_min <= x <= x_max and y_min <= y <= y_max):
            continue

        # Grid indices of the particle centre
        xi = int((x - x_min) / dx)
        yi = int((y - y_min) / dy)

        # Kernel placement bounds in grid coordinates
        x0, y0 = xi - ker_rad_x, yi - ker_rad_y
        x1, y1 = x0 + ker_size_x, y0 + ker_size_y

        # Clip to grid
        x_start, y_start = max(x0, 0), max(y0, 0)
        x_end,   y_end   = min(x1, x_bins), min(y1, y_bins)

        # Corresponding kernel slice
        kx_start = x_start - x0
        ky_start = y_start - y0
        kx_end   = kx_start + (x_end - x_start)
        ky_end   = ky_start + (y_end - y_start)

        # Safe accumulation
        presence_map[y_start:y_end, x_start:x_end] += \
            gaussian_kernel[ky_start:ky_end, kx_start:kx_end]

    # ──────────────── Render and save current frame ───────────────────
    fig, ax = plt.subplots(figsize=(8, 8))       # square figure keeps aspect
    extent = [x_min, x_max, y_min, y_max]
    ax.imshow(presence_map, origin='lower', extent=extent,
              cmap='gray_r', aspect='equal')
    #ax.set_title(f"Time = {time[t_idx]:.2f}")
    #ax.set_xlabel("X")
    #ax.set_ylabel("Y")
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()

    frame_path = os.path.join(frame_dir, f"frame_{t_idx:04d}.png")
    fig.savefig(frame_path, dpi=100)
    plt.close(fig)

print("All frames saved.")

# ───────────────────────────── Video ──────────────────────────────────
print("Creating video …")
ffmpeg_cmd = (
    f"ffmpeg -y -framerate 30 -i {frame_dir}/frame_%04d.png "
    f"-c:v libx264 -pix_fmt yuv420p {video_path}"
)
os.system(ffmpeg_cmd)
print(f"Video saved: {video_path}")

