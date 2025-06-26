import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.ndimage import gaussian_filter

# Parameters
folder = "plumes_x_periodic_wider"
folder = "plumes_x_periodic_wider_NP_512"
folder = "plumes_x_periodic_SQUARE_NP_256"
file_pattern = os.path.join(folder, "particle_*.dat")
output_file = os.path.join(folder, "particle_presence_map.svg")

# Y cropping limits
y_min_crop = -32
y_max_crop = 24

# Gaussian resolution and kernel size
time_bins = 250
y_bins = 250
sigma_time = 10.0  # in bins
sigma_y = 10.0     # in bins

# Load files
particle_files = sorted(glob.glob(file_pattern))
if not particle_files:
    raise FileNotFoundError(f"No particle_*.dat files found in folder: {folder}")

# Load reference time
ref_data = np.loadtxt(particle_files[0])
time = ref_data[:, 0]
num_steps = len(time)

# Create 2D histogram
presence_map = np.zeros((y_bins, time_bins), dtype=np.float32)

time_min, time_max = time[0], time[-1]
time_edges = np.linspace(time_min, time_max, time_bins + 1)
y_edges = np.linspace(y_min_crop, y_max_crop, y_bins + 1)
time_centers = 0.5 * (time_edges[:-1] + time_edges[1:])
y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

# Fill 2D histogram
for f in particle_files:
    data = np.loadtxt(f)
    y_vals = data[:, 2]
    for t_idx in range(num_steps):
        y = y_vals[t_idx]
        if y_min_crop <= y <= y_max_crop:
            t_val = time[t_idx]
            t_bin = np.searchsorted(time_edges, t_val) - 1
            y_bin = np.searchsorted(y_edges, y) - 1
            if 0 <= t_bin < time_bins and 0 <= y_bin < y_bins:
                presence_map[y_bin, t_bin] += 1.0

# Apply Gaussian smoothing
presence_map_smoothed = gaussian_filter(presence_map, sigma=(sigma_y, sigma_time))

# Plot with contourf
T, Y = np.meshgrid(time_centers, y_centers)
plt.figure(figsize=(8, 5))
contour = plt.contourf(T, Y, presence_map_smoothed, levels=20, cmap='gray_r')
plt.colorbar(contour, label="Presence Density")
plt.xlabel("time")
plt.ylabel("y")
#plt.title("Particle Presence Map (Gaussian Contours)")
plt.tight_layout()
plt.savefig(output_file, format='svg')
print(f"Saved: {output_file}")

