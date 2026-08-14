/** 
# Testing `output_facets_vtu()` 

In this example use the sample results to test the output routines. 

~~~pythonplot
import sys, os
if os.path.exists('test_output_vtu.py'):
  sys.path.append('.')
else:
  sys.path.append('../')
import test_output_vtu
data = test_output_vtu.print_vtu_info('Interface.vtu')
~~~

*/

#define MAXLEVEL 8
#define ASPECTRATIO 8
vector h[];

#include "navier-stokes/centered.h"
#include "vof.h"

#include "acastillo/output_fields/vtu/output_vtu.h"
#include "view.h"
#include "curvature.h"

scalar f[], * interfaces = {f};

#define r2 (sq(x) + sq(y) + sq(z))

int main(){

  #if !TREE
    #if dimension == 2
      dimensions (nx = 1, ny = ASPECTRATIO);
    #else
      dimensions (nx = 1, ny = ASPECTRATIO, nz = 1);
    #endif
  #endif

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

  /**
  Verify the interface that was just written: every facet vertex must lie on
  the cardioid used to initialise the volume fraction, to within one cell. The
  report goes to stderr, i.e. to the `log` diffed against
  `test_output_vtu_facets.ref`. */
  if (pid() == 0)
    system ("python3 ../test_output_vtu.py --cardioid Interface.vtu 1>&2");
}
