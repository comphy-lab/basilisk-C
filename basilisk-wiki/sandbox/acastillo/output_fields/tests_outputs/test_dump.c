#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

#define MAXLEVEL 8
#define ASPECTRATIO 2
#define r2 (sq(x) + sq(y))

int main(){

  #if !TREE
    #if dimension == 2
      dimensions (nx = 1, ny = ASPECTRATIO);
    #else
      dimensions (nx = 1, ny = ASPECTRATIO, nz = 1);
    #endif
  #endif

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
  scalar f[];
  foreach(){    
    f[] = x;
    p[] = y;    
    u.x[] = sin(2*pi*x) * cos(2*pi*y);
    u.y[] = -cos(2*pi*x) * sin(2*pi*y);
  }

  dump();
}
