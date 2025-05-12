/**
# `foreach_segment()` misses intersection

~~~gnuplot Example
reset
set xr [-1.3:-0.6]
set yr [-0.2:0.5]
plot 'out' w l t 'cells' ,				\
'Total_segment' w l lw 4 lc 'red' ,			\
'log' w l t 'segments' lw 2 lc 'green
~~~
*/

#include "utils.h"

int main() {
  L0 = 2.5;
  X0 = Y0 = -L0/2.0123;
  init_grid (16);
  output_cells();
  coord S[2] = {{sin(2*pi*8./11), cos(2*pi*8./11)},
		{sin(2*pi*9./11), cos(2*pi*9./11)}};
  FILE * fp = fopen ("Total_segment", "w");
  fprintf (fp, "%g %g\n%g %g\n", S[0].x, S[0].y, S[1].x, S[1].y);
  fclose (fp);
  foreach_segment (S, r) {
    fprintf (stderr, "%g %g\n%g %g\n\n",
	     r[0].x, r[0].y, r[1].x, r[1].y);
  }
}
