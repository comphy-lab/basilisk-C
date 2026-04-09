/** 
# Testing `output_catalyst()` 

In this example use the sample results to test the output routines and verify
that the zero-copy arrays are mapped correctly to Catalyst.

To run this, you must have ParaView Catalyst installed and linked.
*/

#ifndef MAXLEVEL
#define MAXLEVEL 5
#endif
#define ASPECTRATIO 8
vector h[];

#include "navier-stokes/centered.h"
#include "vof.h"

// Note: you may need to define HAVE_HDF5 for the internal helpers
#define HAVE_HDF5 1
#include "acastillo/output_fields/catalyst/output_catalyst.h"
#include "curvature.h"

scalar f[], * interfaces = {f};

#define r2 (sq(x) + sq(y) + sq(z))

int main(int argc, char* argv[]){
  
  if (pid() == 0) {
      fprintf(stderr, "Binary Booting. Total Arguments Recevied: %d\n", argc);
      for(int i=0; i<argc; i++) {
          fprintf(stderr, "  Arg[%d]: %s\n", i, argv[i]);
      }
      fflush(stderr);
  }

  // Parse any catalyst scripts provided in CLI
  output_catalyst_init(argc, argv);

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

  // Initial conditions 
  double a = 0.15;  
  fraction (f, sq(sq(x) + sq(y) - 2*a*x) - 4*sq(a)*(sq(x) + sq(y)));

  scalar p[];
  vector u[];
  foreach() {
    p[] = sq(x);
    u.x[] = x;
    u.y[] = y;
    #if dimension == 3
    u.z[] = z;
    #endif
  }

  // Execute Catalyst across multiple cycles to trigger the Extractors!
  for (int cycle = 0; cycle <= 2; cycle++) {
      output_catalyst({f, p}, {u}, cycle, cycle * 0.1);
  }

  // Clean-up
  output_catalyst_finalize();
}
