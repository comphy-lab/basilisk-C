/** 
# Testing `output_xmf_box()`

Tests XDMF output functions constrained by a box region.
*/

#include "acastillo/output_fields/xdmf/output_xdmf.h"
#define MAXLEVEL 8
#define ASPECTRATIO 8
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

  // Create test fields 
  scalar f[], p[];
  vector u[];
  foreach(){
    f[] = point.level;
    p[] = pid();
    u.x[] = x;
    u.y[] = y;
    #if dimension == 3
      u.z[] = z;
    #endif
  }

  // Write an xmf file, but only inside a region defined by box
  #if dimension == 2
    coord box[2] = {{-0.25, -0.25}, {0.25, 0.25}};
  #else
    coord box[2] = {{-0.25, -0.25, -0.25}, {0.25, 0.25, 0.25}};
  #endif
  output_xmf_box({f,p}, {u}, "domain_xmf", box);

}
