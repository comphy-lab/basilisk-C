/** 
# Testing `output_vtkhdf_box()`

Tests VTKHDF output functions using standard analytical test functions:

- Level set function (signed distance to circle)
- Gaussian scalar field with polynomial modulation  
- Taylor-Green vortex (divergence-free velocity field)

~~~pythonplot Sample Visualization
import pyvista as pv

# --- 1. Read Data and Extract Scalar ---
# PyVista automatically handles .pvtu and loads data into a mesh object.
mesh = pv.read('./domain.hdf')
scalar_name = 'u.x'

# --- 2. Prepare Scalar for Plotting (Extract X-component) ---
if mesh.cell_data[scalar_name].ndim > 1:
    # If 'u.x' is a vector, extract the X-component and name the new array
    mesh[f'{scalar_name}_component'] = mesh.cell_data[scalar_name][:, 0]
    active_scalar = f'{scalar_name}_component'
else:
    active_scalar = scalar_name

# --- 3. Visualize and Save 3D Mesh ---
# Plotter handles 3D rendering and visualization automatically.
plotter = pv.Plotter(off_screen=True, window_size=(1000, 800))
plotter.add_mesh(
    mesh,
    scalars=active_scalar,
    cmap='magma',
    scalar_bar_args={'title': 'Cell Data u.x'}
)
plotter.view_isometric()
plotter.screenshot('field_u.x_pyvista_3d.png', scale=2)
plotter.close()
~~~

*/

#include "grid/octree.h"
#include "acastillo/output_fields/vtkhdf/output_vtkhdf_box.h"
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
    u.z[] = 0.;
  }

  // And write a vtkhdf file, but only inside a region defined by box
  coord box[2] = {{-0.25, -0.25, -0.025}, {0.25, 0.25, 0.025}};
  output_vtkhdf_box({f,p}, {u}, "domain.hdf", box);

}
