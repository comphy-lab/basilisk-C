# This script generates the video in the example fall_vtk.c
# To run this script you need the following python packages:
# 1. numpy  2. matplotlib  3. scipy   4. pyvista
# You will also need ffmpeg for the final video creation
# 
# 
# The script generates a series of images in a temporary folder (__frames__/)
# And then calls ffmpeg to make a video out of the images

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import vtk
import pyvista as pv
import os

# You might need to uncomment line below if you are using a headless display system
pv.start_xvfb()

# Creating the temporary folder to store frames
frames_folder = "__frames__/"
if( not os.path.isdir(frames_folder) ):
  os.mkdir(frames_folder)
else:
  os.system("rm %s/*" % (frames_folder))

vtk_number_frames = 101
list_Wi = [0, 0.5, 1, 5, 10, 100]

list_simulation_step = [-1, -1, -1, -1, -1, -1]
list_last_step = [-1, -1, -1, -1, -1, -1]

max_time = 5.0
time = 0.0
vtk_delta_time = 0.025
frame = 0
while( time<=max_time ):

  groups = [(0, np.s_[:]), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
  plotter = pv.Plotter(shape=(3, 3), groups=groups, row_weights=[0.5, 0.5, 0.5], window_size=(1600, 800), off_screen=True)

  f, ax = plt.subplots(tight_layout=True)

  for i_plot in range( len(list_Wi) ):
    Wi = float(list_Wi[i_plot])
    beta = 1 if (Wi==0) else 0.1
    mesh_level = 7

    # Plotting the graph Diameter vs Time
    base_folder = "../outputs/"
    simulation_folder = "fall_mesh%d_Wi%g/" % (mesh_level, Wi)
    file_name = base_folder + simulation_folder + "/log_file.txt"
    file_data = np.loadtxt(file_name)
    array_time = file_data[:, 0]
    array_diameter = file_data[:, 2]
    ax.plot(array_time, array_diameter, label=("Wi = %g" % (Wi)) )

    # Finding the next printed step for this simulation
    step = list_simulation_step[i_plot] + 1
    while( step<100000 and not os.path.exists(base_folder + simulation_folder + "/Interface-N%d.vtk" % (step)) ):
      step += 1
    list_simulation_step[i_plot] = step

    if( step<100000 ):
      list_last_step[i_plot] = step
    step = list_last_step[i_plot]

    # Loading the mesh vtk files and applying some filters
    grid = pv.UnstructuredGrid(base_folder + simulation_folder + "/Mesh-N%d.vtk" % (step))
    grid = grid.rotate_z(90, [0, 0, 0], inplace=False)
    grid = grid.clip_box(bounds=[-1.5, 0.0, 0.0, 2.5, 0.0, 0.0], invert=False, crinkle=True)
    grid += grid.reflect((1, 0, 0), point=(0, 0, 0))
    grid.active_scalars_name = "Axial-Velocity"

    # Loading the interface vtk files and applying filters
    interface = pv.PolyData(base_folder + simulation_folder + "/Interface-N%d.vtk" % (step))
    interface = interface.rotate_z(90, [0, 0, 0], inplace=False)
    interface += interface.reflect((1, 0, 0), point=(0, 0, 0))

    # Some visual settings for the colorbar
    minmax = grid.get_data_range()
    colorbar_args = dict(
      title="Axial Vel.",
        title_font_size=20,
        label_font_size=14,
        shadow=True,
        n_labels=5,
        italic=True,
        fmt="%.3f",
        font_family="arial",
        color="black",
        vertical=True,
        position_x=0.05,
        position_y=0.05,
        height=0.8,
    )

    # Adding the mesh and interface to the visualization window
    row_plot = 1 + int( i_plot/3 )
    column_plot = i_plot%3
    plotter.subplot(row_plot, column_plot)
    plotter.set_background(color=[1.0, 1.0, 1.0])
    plotter.camera_position = 'xy'
    plotter.reset_camera(bounds=[-0.5, 0.5, 0.1, 2.4, 0.0, 0.0])

    plotter.add_mesh(grid, style='Surface', show_edges=True, cmap="jet", scalar_bar_args=colorbar_args)
    plotter.add_mesh(interface, style='Wireframe', color="white", line_width=3)
    plotter.add_text("Wi = %g" % (Wi), position = "lower_right", font_size=10, color="black")


  plotter.subplot(0, 0)
  ax.grid(True)
  ax.legend()
  ax.plot([float(time), float(time)], [0, 10.0], "--")
  ax.set_ylim([1.0, 2.5])
  ax.set_xlim([0.0, 5.0])
  ax.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
  ax.set_title("Viscoelastic droplet impact")
  ax.set_xlabel("Time")
  ax.set_ylabel("Diameter")
  h_chart = pv.ChartMPL(f)
  h_chart.background_color = (1.0, 1.0, 1.0, 1.0)
  plotter.add_chart(h_chart)


  plotter.screenshot("%s/frame_%06d" % (frames_folder, frame))
  plotter.close()
  plt.close()

  # Updating the simulation time
  print("Time: %g out of %g\n" % (time, max_time) )
  time += vtk_delta_time
  frame += 1


video_name = "video_fall.avi"
framerate = 20
os.system('ffmpeg -framerate %s -i %s/frame_%%06d.png -b 5000k %s' % (str(framerate), frames_folder, video_name))