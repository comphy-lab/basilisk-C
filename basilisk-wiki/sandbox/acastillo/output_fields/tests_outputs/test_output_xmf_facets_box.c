/**
# Testing `output_facets_xmf_box()`

In this example use the sample results to test the output routines for the
VOF interface restricted to a box region.

*/

#define MAXLEVEL 8
#define ASPECTRATIO 8
vector h[];

#include "navier-stokes/centered.h"
#include "vof.h"

#include "acastillo/output_fields/xdmf/output_xdmf.h"
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

  // Restrict the export to the half of the domain with x >= 0: the cardioid
  // is not symmetric about x=0, so this actually drops facets, unlike a box
  // that happened to contain the whole shape.
  #if dimension == 2
    coord box[2] = {{0., -L0/2.}, {L0/2., L0/2.}};
  #else
    coord box[2] = {{0., -L0/2., -L0/2.}, {L0/2., L0/2., L0/2.}};
  #endif
  output_facets_xmf_box(f, {kappa}, "Interface_xmf_box", box);

  /**
  Verify the interface that was just written: every facet vertex must lie on
  the cardioid used to initialise the volume fraction, to within one cell,
  every attribute must be finite, and -- the property specific to this
  writer -- no emitted facet's cell may lie outside the requested box. The
  report goes to stderr, i.e. to the `log` diffed against
  `test_output_xmf_facets_box.ref`. */
  if (pid() == 0)
    system ("python3 ../test_output_xmf.py --cardioid "
            "--box=0,-1,1,1 Interface_xmf_box.xmf 1>&2");
}
