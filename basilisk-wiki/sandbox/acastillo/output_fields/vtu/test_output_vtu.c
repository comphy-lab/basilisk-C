/** 
# Testing `output_vtu()` and `output_vtu_slice()`

In this example use the sample results to test the output routines. 

*/


#include "output_vtu.h"
#define MAXLEVEL 8
#define r2 (sq(x) + sq(y))

// Function to generate a triangle wave
double triangle_wave(double t, double amplitude, double period) {  
  return (4 * amplitude / period) * fabs(t - (period / 2) * floor(2 * t / period + 0.5)) - amplitude/2;
}

int main(){
  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << (MAXLEVEL-4);
  init_grid(N);

 #if TREE
  double outer_radius = 0.25;
  double inner_radius = 0.1 ;
  refine(((r2 < sq(outer_radius)) && (r2 > sq(inner_radius))) && level < MAXLEVEL);
#endif

  // We create some fields
  scalar f[], p[];
  vector u[];
  foreach(){
    f[] = (triangle_wave(x,0.33,0.25) > y);
    p[] = (y + 0.1*noise())*f[];    
    u.x[] = (  sq(y) + sin(2*pi*x)*cos(2*pi*y) + 0.1*noise())*f[] + (sin(3*pi*x)*cos(3*pi*y) + 0.1*noise())*(1-f[]);
    u.y[] = (cube(y) + sin(4*pi*x)*cos(4*pi*y) + 0.1*noise())*f[] + (sin(5*pi*x)*cos(5*pi*y) + 0.1*noise())*(1-f[]);
  }

  // And write a vtu file
  output_vtu({f,p}, {u}, "domain");

#if dimension > 2
  output_slice_vtu({f,p}, {u}, "slice_x", (coord){1,0,0}, 0);
  output_slice_vtu({f,p}, {u}, "slice_y", (coord){0,1,0}, 0);
  output_slice_vtu({f,p}, {u}, "slice_z", (coord){0,0,1}, 0);
#endif

}
