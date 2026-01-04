/**
# Verification of Curvature Correction

This test verifies the implementation of the curvature correction for a vortex filament.
We use a simple ring geometry where analytical solutions for the geometry are available.

We check:
1. The Frenet-Serret frame and local coordinates against analytical values.
2. The divergence-free condition of the computed velocity and vorticity fields ** which is not met here **. 


~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName('./slice_y.pvtu')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
c_array = vtk_to_numpy(data.GetCellData().GetArray('rho')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
z = centers[:, 2]
c = c_array[:]

# Plot with matplotlib
plt.figure(figsize=(8, 6))
im = plt.tripcolor(z, x, c, cmap='magma', rasterized=True)
cs = plt.tricontour(z, x, c, levels=8, colors='white', linewidths=0.5)
plt.xlabel('z')
plt.ylabel('x')
plt.axis('image')
plt.xlim(-0.65, 0.65)
plt.ylim(-0.0, 1.30)
plt.colorbar(im, label=r'$\rho(x,y,z)$')
plt.clabel(cs, inline=True, fontsize=10, fmt=r'$\phi$=%1.2f')
plt.title('Local radial coordinate')
plt.savefig('slice_y_field_rho.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName('./slice_y.pvtu')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
c_array = vtk_to_numpy(data.GetCellData().GetArray('phi')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
z = centers[:, 2]
c = c_array[:]*180./np.pi;

# Plot with matplotlib
plt.figure(figsize=(8, 6))
im = plt.tripcolor(z, x, c, cmap='twilight', rasterized=True)
cs = plt.tricontour(z, x, c, levels=8, colors='white', linewidths=0.5)
plt.xlabel('z')
plt.ylabel('x')
plt.axis('image')
plt.xlim(-0.65, 0.65)
plt.ylim(-0.0, 1.30)
plt.colorbar(im, label=r'$\phi(x,y,z)$')
plt.clabel(cs, inline=True, fontsize=10, fmt=r'$\phi$=%1.0f°')
plt.title('Local angular coordinate')
plt.savefig('slice_y_field_phi.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName('./slice_y.pvtu')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
u_array = vtk_to_numpy(data.GetCellData().GetArray('omega0.x')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
z = centers[:, 2]
u = u_array[:, :]
# Plot with matplotlib
fig, ax = plt.subplots(1,3,figsize=(8*3, 6))

im = ax[0].tripcolor(z, x, u[:,0], cmap='magma', rasterized=True)
cs = ax[0].tricontour(z, x, u[:,0], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[0], label=r'$\omega_\rho^{(0)}(x,y,z)$')
ax[0].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[0].set_title('Local radial vorticity (axi-symmetric)')

im = ax[1].tripcolor(z, x, u[:,1], cmap='magma', rasterized=True)
cs = ax[1].tricontour(z, x, u[:,1], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[1], label=r'$\omega_\phi^{(0)}(x,y,z)$')
ax[1].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[1].set_title('Local azimuthal vorticity (axi-symmetric)')

im = ax[2].tripcolor(z, x, u[:,2], cmap='magma', rasterized=True)
cs = ax[2].tricontour(z, x, u[:,2], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[2], label=r'$\omega_\xi^{(0)}(x,y,z)$')
ax[2].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[2].set_title('Local axial vorticity (axi-symmetric)')

for i in range(3):
    ax[i].set_xlabel('z')
    ax[i].set_ylabel('x')
    ax[i].axis('image')
    ax[i].set_xlim(-0.65, 0.65)
    ax[i].set_ylim(-0.0, 1.30)
    
plt.savefig('slice_y_field_omega0.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName('./slice_y.pvtu')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
u_array = vtk_to_numpy(data.GetCellData().GetArray('omega1.x')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
z = centers[:, 2]
u = u_array[:, :]
# Plot with matplotlib
fig, ax = plt.subplots(1,3,figsize=(8*3, 6))

im = ax[0].tripcolor(z, x, u[:,0], cmap='magma', rasterized=True)
cs = ax[0].tricontour(z, x, u[:,0], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[0], label=r'$\omega_\rho^{(1)}(x,y,z)$')
ax[0].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[0].set_title('Local radial vorticity (curvature correction)')

im = ax[1].tripcolor(z, x, u[:,1], cmap='magma', rasterized=True)
cs = ax[1].tricontour(z, x, u[:,1], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[1], label=r'$\omega_\phi^{(1)}(x,y,z)$')
ax[1].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[1].set_title('Local azimuthal vorticity (curvature correction)')

im = ax[2].tripcolor(z, x, u[:,2], cmap='magma', rasterized=True)
cs = ax[2].tricontour(z, x, u[:,2], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[2], label=r'$\omega_\xi^{(1)}(x,y,z)$')
ax[2].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[2].set_title('Local axial vorticity (curvature correction)')

for i in range(3):
    ax[i].set_xlabel('z')
    ax[i].set_ylabel('x')
    ax[i].axis('image')
    ax[i].set_xlim(-0.65, 0.65)
    ax[i].set_ylim(-0.0, 1.30)
    
plt.savefig('slice_y_field_omega1.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName('./slice_y.pvtu')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
u_array = vtk_to_numpy(data.GetCellData().GetArray('omega.x')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
z = centers[:, 2]
u = u_array[:, :]
# Plot with matplotlib
fig, ax = plt.subplots(1,3,figsize=(8*3, 6))

im = ax[0].tripcolor(z, x, u[:,0], cmap='magma', rasterized=True)
cs = ax[0].tricontour(z, x, u[:,0], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[0], label=r'$\omega_x(x,y,z)$')
ax[0].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[0].set_title('Local x-vorticity')

im = ax[1].tripcolor(z, x, u[:,1], cmap='magma', rasterized=True)
cs = ax[1].tricontour(z, x, u[:,1], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[1], label=r'$\omega_y(x,y,z)$')
ax[1].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[1].set_title('Local y-vorticity')

im = ax[2].tripcolor(z, x, u[:,2], cmap='magma', rasterized=True)
cs = ax[2].tricontour(z, x, u[:,2], levels=5, colors='white', linewidths=0.5)
fig.colorbar(im, ax=ax[2], label=r'$\omega_z(x,y,z)$')
ax[2].clabel(cs, inline=True, fontsize=10, fmt='%1.2f')
ax[2].set_title('Local z-vorticity')

for i in range(3):
    ax[i].set_xlabel('z')
    ax[i].set_ylabel('x')
    ax[i].axis('image')
    ax[i].set_xlim(-0.65, 0.65)
    ax[i].set_ylim(-0.0, 1.30)
    
plt.savefig('slice_y_field_omega.png', dpi=150, bbox_inches='tight')
plt.close()
~~~


~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName('./slice_y.pvtu')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
c_array = vtk_to_numpy(data.GetCellData().GetArray('div_u')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
z = centers[:, 2]
c = c_array[:];

# Plot with matplotlib
plt.figure(figsize=(8, 6))
im = plt.tripcolor(z, x, c, cmap='twilight', rasterized=True)

plt.xlabel('z')
plt.ylabel('x')
plt.axis('image')
plt.xlim(-0.65, 0.65) 
plt.ylim(-0.0, 1.30)
plt.colorbar(im, label=r'$\text{div}(\mathbf{u})$')
plt.title('Divergence of velocity field')
plt.savefig('slice_y_div_u.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName('./slice_y.pvtu')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
c_array = vtk_to_numpy(data.GetCellData().GetArray('div_omega')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
z = centers[:, 2]
c = c_array[:];

# Plot with matplotlib
plt.figure(figsize=(8, 6))
im = plt.tripcolor(z, x, c, cmap='twilight', rasterized=True)

plt.xlabel('z')
plt.ylabel('x')
plt.axis('image')
plt.xlim(-0.65, 0.65) 
plt.ylim(-0.0, 1.30)
plt.colorbar(im, label=r'$\text{div}(\mathbf{\omega})$')
plt.title('Divergence of vorticity field')
plt.savefig('slice_y_div_omega.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

*/

#include "grid/octree.h"
#include "view.h"
#include "acastillo/input_fields/filaments.h"
#include "acastillo/output_fields/vtu/output_vtu.h"

int minlevel = 4;
int maxlevel = 7;

int main() {
  L0 = 1.50;
  X0 = -L0 / 16;
  Y0 = -L0 / 16;
  Z0 = -L0 / 2;
  N = 1 << minlevel;
  init_grid(N);

  scalar v_mag[], v_ang[];
  vector u[], omega[];
  
  // Initialize fields
  foreach() {
    v_mag[] = 0;
    v_ang[] = 0;
    foreach_dimension() {
      u.x[] = 0;
      omega.x[] = 0;
    }
  }

  // --- 1. Setup Ring Filament ---
  int turns = 4; // Single ring
  int nseg_per_turn = 256; // High resolution for geometry
  int nseg = nseg_per_turn * turns + 1;
  double R = 1.0;
  double U_c = 0.1;
  double dphi = 2 * pi / ((double)nseg_per_turn);
  double phi_arr[nseg];
  double a1[nseg];
  coord C1[nseg];

  for (int i = 0; i < nseg; i++) {
    phi_arr[i] = dphi * (double)i - 3 * pi; 
    C1[i].x = R * cos(phi_arr[i]);
    C1[i].y = R * sin(phi_arr[i]);
    C1[i].z = 0.0;
    a1[i] = 0.10;
  }

  coord xshift = {0, 0, 0}, dxshift = {0, 0, 0};
  struct vortex_filament filament1;
  allocate_vortex_filament_members(&filament1, nseg);
  initialize_filaments(filament1, nseg, dphi, phi_arr, a1, C1, xshift, dxshift);

  // --- 2. Compute Fields ---
  // Refine grid
  scalar dmin[];
  for (int i = (maxlevel - minlevel - 1); i >= 0; i--) {
    foreach () {
      struct vortex_filament params1 = filament1;
      params1.pcar = (coord){x, y, z};
      dmin[] = (get_min_distance(spatial_period=1, max_distance=4*L0, vortex=&params1) < (i + 1) * a1[0]) * noise();
    }
    adapt_wavelet((scalar *){dmin}, (double[]){1e-12}, maxlevel - i, minlevel);
  }

  scalar phi[], rho[];
  vector e_r[], e_th[];
  vector local_coords[];
  vector omega0[], omega1[];

  foreach () {
    struct vortex_filament params1 = filament1;
    params1.pcar = (coord){x, y, z};
    struct local_filament vortex1 = get_local_coordinates(spatial_period=1, max_distance=4*L0, vortex=&params1);

    if (vortex1.near == 1) {
      integration_results corr;
      compute_A_with_derivatives(vortex1.rho, vortex1.a, U_c, &corr);

      batchelor_vortex sol;
      sol = calculate_vortex_flow(&vortex1, U_c, &corr);

      rho[] = vortex1.rho;
      phi[] = vortex1.theta - vortex1.theta0;
      foreach_dimension() {
        e_r.x[]  =  cos(phi[]) * vortex1.Nvec.x + sin(phi[]) * vortex1.Bvec.x;
        e_th.x[] = -sin(phi[]) * vortex1.Nvec.x + cos(phi[]) * vortex1.Bvec.x;
      }

      local_coords.x[] = vortex1.Mcar.x;
      local_coords.y[] = vortex1.Mcar.y;
      local_coords.z[] = vortex1.Mcar.z;

      omega0.x[] = sol.w0_r;
      omega0.y[] = sol.w0_th;
      omega0.z[] = sol.w0_xi;

      omega1.x[] = sol.w1_r;
      omega1.y[] = sol.w1_th;
      omega1.z[] = sol.w1_xi;

      foreach_dimension() {
        u.x[] += (sol.u0_r * e_r.x[]) + (sol.u0_th * e_th.x[]) + (sol.u0_xi * vortex1.Tvec.x);
        u.x[] += (sol.u1_r * e_r.x[]) + (sol.u1_th * e_th.x[]) + (sol.u1_xi * vortex1.Tvec.x);

        omega.x[] += (sol.w0_r * e_r.x[]) + (sol.w0_th * e_th.x[]) + (sol.w0_xi * vortex1.Tvec.x);
        omega.x[] += (sol.w1_r * e_r.x[]) + (sol.w1_th * e_th.x[]) + (sol.w1_xi * vortex1.Tvec.x);
      }
    }
  }
  boundary({u, omega});

  scalar div_u[], div_omega[];
  foreach() {
    div_u[] = 0.;
    div_omega[] = 0.;
    foreach_dimension() {
      div_u[] += center_gradient(u.x);
      div_omega[] += center_gradient(omega.x);
    }
  }

  stats s_du = statsf(div_u);
  stats s_dw = statsf(div_omega);
  fprintf(stderr, "div_u: max %g, std %g\n", s_du.max, s_du.stddev);
  fprintf(stderr, "div_omega: max %g, std %g\n", s_dw.max, s_dw.stddev);

  output_slice_vtu({rho, phi, div_u, div_omega}, {u, omega, e_r, e_th, local_coords, omega0, omega1}, "slice_x", (coord){1,0,0}, 0);
  output_slice_vtu({rho, phi, div_u, div_omega}, {u, omega, e_r, e_th, local_coords, omega0, omega1}, "slice_y", (coord){0,1,0}, 0);
  output_slice_vtu({rho, phi, div_u, div_omega}, {u, omega, e_r, e_th, local_coords, omega0, omega1}, "slice_z", (coord){0,0,1}, 0);

  free_vortex_filament_members(&filament1);
  return 0;
}
