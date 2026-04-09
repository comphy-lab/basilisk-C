/** 
# Testing `output_facets_vtkhdf()` 

In this example use the sample results to test the output routines. 
*/

#define MAXLEVEL 9
#define ASPECTRATIO 8
vector h[];

#include "navier-stokes/centered.h"
#include "vof.h"

#include "acastillo/output_fields/vtkhdf/output_vtkhdf.h"
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

  double a = 0.15;
  fraction (f, sq(sq(x) + sq(y) - 2*a*x) - 4*sq(a)*(sq(x) + sq(y)));

  scalar kappa[];
  curvature (f, kappa);
  output_facets_vtkhdf(f, kappa, "Interface_vtkhdf.hdf");
}
