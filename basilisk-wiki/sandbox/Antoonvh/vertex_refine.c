/**
# Vertex-value prolongation

Using injection if possible...

~~~gnuplot
set xr [-0.1:1.1]
set yr [-0.1:1.1]
set key off
set size square
plot 'out' w l lw 2, 'log' u ($1+0.02):($2+0.015):3 w labels
~~~
*/
#include "utils.h"
int Maxlevel = 2;
#define FUNC (x + y + z)

void refine_vertex (Point point, scalar s) {
  for (int i = 0; i < 2; i++) 
    for (int j = 0; j < 1 + (dimension > 1); j++) 
      for (int k = 0; k < 1 + (dimension > 2); k++) { 
	fine (s, 2*i, 2*j, 2*k) = s[i, j, k];
	fine (s, 2*i + 1, 2*j    , 2*k)     = nodata; //1
#if (dimension > 1)
	fine (s, 2*i    , 2*j + 1, 2*k)     = nodata; //1
	fine (s, 2*i + 1, 2*j + 1, 2*k)     = nodata; //2
#if (dimension > 2)
	fine (s, 2*i    , 2*j    , 2*k + 1) = nodata; //1
	fine (s, 2*i + 1, 2*j    , 2*k + 1) = nodata; //2
	fine (s, 2*i    , 2*j + 1, 2*k + 1) = nodata; //2
	fine (s, 2*i +1 , 2*j + 1, 2*k + 1) = nodata; //3
#endif 
#endif
      }
}

int main() {
  init_grid (1 << Maxlevel);
  vertex scalar phi[];
  phi.prolongation = phi.refine = refine_vertex;
  foreach_vertex() 
    phi[] = FUNC;
  refine (level == Maxlevel);
  output_cells();
  foreach_vertex()
    if (phi[] != nodata)
      fprintf (stderr, "%g %g %g\n", x, y , phi[]);
}
