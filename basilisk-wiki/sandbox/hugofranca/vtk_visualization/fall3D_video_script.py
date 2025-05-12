# This script generates the video in the example fall_vtk.c
# To run this script you need the following python packages:
# 1. numpy  2. matplotlib  3. scipy   4. pyvista
# You will also need ffmpeg for the final video creation
# 
# 
# The script generates a series of images in a temporary folder (__frames3D__/)
# And then calls ffmpeg to make a video out of the images

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import vtk
import pyvista as pv
import os


# You might need to uncomment the line below if you are using a headless display system
# pv.start_xvfb()

# Creating the temporary folder to store frames
frames_folder = "__frames3D__/"
if( not os.path.isdir(frames_folder) ):
  os.mkdir(frames_folder)
else:
  os.system("rm %s/*" % (frames_folder))


mesh_level = 7

simulation_step = -1
last_step = -1
time = 0.0
vtk_delta_time = 0.025
frame = 0
max_time = 5.0
while( time<=max_time ):

  groups = [(0, np.s_[:]), (1, 0), (1, 1)]
  plotter = pv.Plotter(shape=(2, 2), groups=groups, row_weights=[0.3, 0.7], window_size=(1600, 800), off_screen=True)

  f, ax = plt.subplots(tight_layout=True)

  # Plotting the graph Diameter vs Time
  base_folder = "../outputs/"
  simulation_folder = "fall3D_mesh%d_Newt/" % (mesh_level)
  file_name = base_folder + simulation_folder + "/log_file.txt"
  file_data = np.loadtxt(file_name)
  array_time = file_data[:, 0]
  array_diameter = file_data[:, 2]
  ax.plot(array_time, array_diameter, label=("Newtonian"))

  # Finding the next printed step for this simulation
  step = simulation_step + 1
  while( step<100000 and not os.path.exists(base_folder + simulation_folder + "/Interface-N%d.vtk" % (step)) ):
    step += 1
  simulation_step = step

  if( step<100000 ):
    last_step = step
  step = last_step


  # Loading the vtk files and slicing the mesh
  grid = pv.UnstructuredGrid(base_folder + simulation_folder + "/Mesh-N%d.vtk" % (step))
  slice_xy = grid.slice(normal=[0.0, 0.0, 1.0], origin=[0.0, 2.5, 0.0])
  slice_xz = grid.slice(normal=[0.0, 1.0, 0.0], origin=[0.0, 0.01, 0.0])
  interface = pv.PolyData(base_folder + simulation_folder + "/Interface-N%d.vtk" % (step))

  # Creating the gray plane that represents the wall
  wall = pv.Plane(center=[0.0, 0.0, 0.0], direction=[0.0, 1.0, 1.0], i_size=5.0, j_size=5.0)

  # Plotting the left image
  plotter.subplot(1, 0)
  plotter.set_background(color=[1.0, 1.0, 1.0])
  plotter.add_mesh(slice_xy, style='Wireframe', show_edges=True)
  plotter.add_mesh(interface, style='Surface', color="red", opacity=1.0)
  plotter.add_mesh(wall)
  plotter.camera_position = "xy"
  plotter.camera.zoom(1.3)

  # Plotting the right image
  plotter.subplot(1, 1)
  plotter.set_background(color=[1.0, 1.0, 1.0])
  plotter.add_mesh(slice_xz, style='Wireframe', show_edges=True)
  plotter.add_mesh(interface, style='Surface', color="red", opacity=1.0)
  plotter.camera_position = \
    [(0.0, 8.41958516411084, 12.490826626699343), \
    (0.0, 1.2080030487850308, 0.0), \
    (0.0, 1.0, 0.0)]
  plotter.camera.zoom(1.45)

  # Plotting the chart at the top
  plotter.subplot(0, 0)
  ax.grid(True)
  ax.legend()
  ax.plot([float(time), float(time)], [0, 10.0], "--")
  ax.set_ylim([0.95, 1.6])
  ax.set_xlim([0.0, 5.0])
  ax.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
  ax.set_title("Newtonian 3D droplet impact")
  ax.set_xlabel("Time")
  ax.set_ylabel("Diameter")
  h_chart = pv.ChartMPL(f)#, size=(1.0, 1.0), loc=(0.02, 0.06))
  h_chart.background_color = (1.0, 1.0, 1.0, 1.0)
  plotter.add_chart(h_chart)

  # Saving this frame
  plotter.screenshot("%s/frame_%06d" % (frames_folder, frame))
  plotter.close()
  plt.close()

  # Updating the simulation time
  print("Time: %g out of %g\n" % (time, max_time) )
  time += vtk_delta_time
  frame += 1

video_name = "video_fall_3D.avi"
framerate = 20
os.system('ffmpeg -framerate %s -i %s/frame_%%06d.png -b 5000k %s' % (str(framerate), frames_folder, video_name))