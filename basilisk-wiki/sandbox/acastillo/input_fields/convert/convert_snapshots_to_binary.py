#!/usr/bin/env python

import numpy as np
import h5py
import sys
import argparse



def output_matrix(f, n, x, y, file):
  with open(file, 'wb') as fileID:
    # Write 'n' as single precision (float32)
    np.array([n], dtype=np.float32).tofile(fileID)

    # Write 'y' values as single precision
    for j in range(n):
      np.array([y[j]], dtype=np.float32).tofile(fileID)

    # Write 'x' values and the corresponding 'f' matrix values as single precision
    for i in range(n):
      np.array([x[i]], dtype=np.float32).tofile(fileID)
      for j in range(n):
        np.array([f[j, i]], dtype=np.float32).tofile(fileID)

def output_matrix_double(f, n, x, y, file):
  with open(file, 'wb') as fileID:
    # Write 'n' as double precision (float64)
    np.array([n], dtype=np.float64).tofile(fileID)

    # Write 'y' values as double precision
    for j in range(n):
      np.array([y[j]], dtype=np.float64).tofile(fileID)

    # Write 'x' values and the corresponding 'f' matrix values as double precision
    for i in range(n):
      np.array([x[i]], dtype=np.float64).tofile(fileID)
      for j in range(n):
        np.array([f[j, i]], dtype=np.float64).tofile(fileID)

def main(file_path, level, width):
  print(f"The file path is: {file_path}")

  # Open the HDF5 file
  with h5py.File(file_path, 'r') as hdf:
    # Access groups and datasets
    cells_group = hdf['/Cells']
    geometry_group = hdf['/Geometry']
    topology_dataset = hdf['/Topology']

    # Access datasets within groups
    cells_f_dataset = cells_group['f']
    cells_p_dataset = cells_group['p']
    cells_ux_dataset = cells_group['u.x']
    geometry_points_dataset = geometry_group['Points']

    # Read the data from datasets
    cells_f_data = cells_f_dataset[:]
    cells_p_data = cells_p_dataset[:]
    cells_ux_data = cells_ux_dataset[:]
    geometry_points_data = geometry_points_dataset[:]
    topology_data = topology_dataset[:]

    # Print the shape of each dataset
    print("Cells/f dataset shape:", cells_f_data.shape)
    print("Cells/p dataset shape:", cells_p_data.shape)
    print("Cells/u.x dataset shape:", cells_ux_data.shape)
    print("Geometry/Points dataset shape:", geometry_points_data.shape)
    print("Topology dataset shape:", topology_data.shape)

  L0 = width
  n = 2**level

  from scipy.interpolate import griddata
  cell_centers = np.mean(geometry_points_data[topology_data,:], axis=1)
  cells_x, cells_y, cells_z = np.split(cell_centers,3,axis=1)

  delta = np.max(np.max(geometry_points_data[topology_data,:], axis=1)-np.min(geometry_points_data[topology_data,:], axis=1), axis=1)
  cells_l_data = np.rint(np.log2(L0/delta))

  # Define a regular grid
  xc = ((np.arange(n)+0.5)/n - 0.5) * L0
  yc = ((np.arange(n)+0.5)/n - 0.5) * L0
  grid_x, grid_y = np.meshgrid(xc,yc)

  # Perform the interpolation
  grid_l = np.squeeze(griddata(cell_centers[:,0:2], cells_l_data, (grid_x, grid_y), method='nearest', fill_value=0))
  grid_f = np.squeeze(griddata(cell_centers[:,0:2], cells_f_data, (grid_x, grid_y), method='nearest', fill_value=0))
  grid_p = np.squeeze(griddata(cell_centers[:,0:2], cells_p_data, (grid_x, grid_y), method='nearest', fill_value=0))
  grid_u = np.squeeze(griddata(cell_centers[:,0:2], cells_ux_data[:,0], (grid_x, grid_y), method='nearest', fill_value=0))
  grid_v = np.squeeze(griddata(cell_centers[:,0:2], cells_ux_data[:,1], (grid_x, grid_y), method='nearest', fill_value=0))

  grid_l = np.rint(grid_l)

  output_matrix_double(grid_l, n, xc, yc, file_path[:-3]+'_l.bin')
  output_matrix_double(grid_f, n, xc, yc, file_path[:-3]+'_f.bin')
  output_matrix_double(grid_p, n, xc, yc, file_path[:-3]+'_p.bin')
  output_matrix_double(grid_u, n, xc, yc, file_path[:-3]+'_u.bin')
  output_matrix_double(grid_v, n, xc, yc, file_path[:-3]+'_v.bin')
  
  import matplotlib.pyplot as plt
  # Plot the array
  plt.imshow(grid_f, cmap='viridis')
  # Save the plot to a PNG file
  output_file = 'convert.png'
  plt.savefig(output_file)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Process a file given its path.")
  parser.add_argument("file_path", type=str, help="Path to the file to be processed")
  parser.add_argument("level", type=int, help="Path to the file to be processed")
  parser.add_argument("width", type=float, help="Path to the file to be processed")

  args = parser.parse_args()

  main(args.file_path, args.level, args.width)