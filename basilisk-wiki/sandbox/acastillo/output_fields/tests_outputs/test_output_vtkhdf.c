/** 
# Testing `output_vtkhdf()` and `output_vtkhdf_slice()`

Tests VTKHDF output functions using standard analytical test functions:

- Level set function (signed distance to circle)
- Gaussian scalar field with polynomial modulation  
- Taylor-Green vortex (divergence-free velocity field)

~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the VTKHDF format
reader = vtk.vtkHDFReader()
reader.SetFileName('./domain.hdf')
reader.Update()

data = reader.GetOutput()

# Extract the field data (Cell Data)
u_array = vtk_to_numpy(data.GetCellData().GetArray('u.x')) 

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]
u = u_array[:, 0]

# Plot with matplotlib
plt.figure(figsize=(8, 6))
plt.tricontourf(x, y, u, levels=14, cmap='magma')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar(label='u.x')
plt.title('Cell Data Function u.x')
plt.axis('image')
plt.savefig('field_u.x_cell_data.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

*/

#include "acastillo/output_fields/vtkhdf/output_vtkhdf.h"
#define MAXLEVEL 8
#define r2 (sq(x) + sq(y))

int main(){
  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << (MAXLEVEL-1);
  init_grid(N);

 #if TREE
  double outer_radius = 0.25;
  double inner_radius = 0.1 ;
  refine(((r2 < sq(outer_radius)) && (r2 > sq(inner_radius))) && level < MAXLEVEL);
#endif

  // Create test fields with common analytical functions
  scalar f[], p[];
  vector u[];
  foreach(){
    // Level set function: signed distance to a circle
    f[] = sqrt(r2) - 0.2;
    
    // Scalar field: Gaussian with polynomial modulation
    p[] = exp(-5*r2) * (1 + x*y);
    
    // Vector field: Taylor-Green vortex (divergence-free)
    u.x[] = sin(2*pi*x) * cos(2*pi*y);
    u.y[] = -cos(2*pi*x) * sin(2*pi*y);
  }

  // And write a vtkhdf file
  output_vtkhdf({f,p}, {u}, "domain.hdf");

#if dimension > 2
  output_vtkhdf_slice({f,p}, {u}, "slice_x.hdf", (coord){1,0,0}, 0);
  output_vtkhdf_slice({f,p}, {u}, "slice_y.hdf", (coord){0,1,0}, 0);
  output_vtkhdf_slice({f,p}, {u}, "slice_z.hdf", (coord){0,0,1}, 0);
#endif

}
