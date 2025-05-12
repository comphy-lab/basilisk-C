/** 
# Testing `output_facets_vtu()` 

In this example use the sample results to test the output routines. 

*/

#define MAXLEVEL 7
vector h[];

#include "navier-stokes/centered.h"
#include "vof.h"

#include "./output_vtu.h"
#include "view.h"
#include "curvature.h"

scalar f[], * interfaces = {f};

#define r2 (sq(x) + sq(y) + sq(z))

// Function to generate a triangle wave
double triangle_wave(double t, double amplitude, double period) {  
  return (4 * amplitude / period) * fabs(t - (period / 2) * floor(2 * t / period + 0.5)) - amplitude/2;
}

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

  fraction (f, triangle_wave(x,0.33,0.5) > y);

  draw_vof("f");
  squares("f", linear = false);
  save ("init_f.png");

  scalar kappa[];
  curvature (f, kappa);

  output_facets_vtu(f, kappa, "Interface");
}
