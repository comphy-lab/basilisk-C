/** 
# Testing `output_facets_vtu()` 

In this example use the sample results to test the output routines. 

~~~pythonplot Sample Visualization
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np

# Use the correct reader for the .pvtu (parallel) format
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName('./Interface.vtu')
reader.Update()

data = reader.GetOutput()
points_vtk = data.GetPoints()
points = vtk_to_numpy(points_vtk.GetData())
kappa = vtk_to_numpy(data.GetCellData().GetArray('kappa'))

# Compute Cell Centers using a VTK filter
cell_centers_filter = vtk.vtkCellCenters()
cell_centers_filter.SetInputData(data)
cell_centers_filter.Update()

# Extract the coordinates of the cell centers
centers_vtk = cell_centers_filter.GetOutput().GetPoints()
centers = vtk_to_numpy(centers_vtk.GetData())

x = centers[:, 0]
y = centers[:, 1]

# Plot with matplotlib
plt.figure(figsize=(8, 6))
plt.scatter(x, y, c=kappa, cmap='magma', s=5)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar(label='kappa')
plt.title('Cell Data Function kappa')
plt.axis('image')
plt.savefig('field_kappa_cell_data.png', dpi=150, bbox_inches='tight')
plt.close()
~~~

*/

#define MAXLEVEL 8
vector h[];

#include "navier-stokes/centered.h"
#include "vof.h"

#include "acastillo/output_fields/vtu/output_vtu.h"
#include "view.h"
#include "curvature.h"

scalar f[], * interfaces = {f};

#define r2 (sq(x) + sq(y) + sq(z))

int main(){
  L0 = 2.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << (MAXLEVEL-1);
  init_grid(N);

 #if TREE
  double outer_radius = 0.25;
  double inner_radius = 0.1 ;
  refine(((r2 < sq(outer_radius)) && (r2 > sq(inner_radius))) && level < MAXLEVEL);
#endif

  // Create a test field with common analytical functions
  // Cardioid: (x^2 + y^2 - 2*a*x)^2 = 4*a^2*(x^2 + y^2), simplified for level set
  double a = 0.15;  // Size parameter for cardioid
  fraction (f, sq(sq(x) + sq(y) - 2*a*x) - 4*sq(a)*(sq(x) + sq(y)));

  scalar kappa[];
  curvature (f, kappa);
  output_facets_vtu(f, kappa, "Interface");
}
