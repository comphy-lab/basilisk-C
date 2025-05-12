/**
# Domain decompoisition on trees

See the movie for the grid structure and the `npe()=5` decomposition.

![](dcomp/mov.mp4)

The number of non-local neighbors varies:

~~~gnuplot
set size square
set xlabel 'Refinement central coordinate'
set ylabel 'Non local neighbors' 
plot 'out' t 'Horizontal' , 'log' t 'vertical'
~~~

 */
#include "view.h"
scalar pids[];

int minlevel = 3, maxlevel = 7;
double W = 0.1;

int main() {
  for (double yh = 0.8; yh >= 2*W; yh -= 0.01) {
    init_grid (1 << minlevel);
    unrefine (level > minlevel);
    refine (fabs(y - yh) < W && level < maxlevel);
    foreach() 
      pids[] = pid();
    boundary ({pids});
    view (tx = -0.5, ty = -0.5);
    cells();
    squares ("pids", min = 0, max = npe() - 1);
    save ("mov.mp4");
    int nc = 0;
    foreach(reduction(+:nc)) { 
      if (!is_local(neighbor(1,0,0)))
	nc++;
      if (!is_local(neighbor(0,1,0)))
	  nc++;
    }
    if (pid() == 0)
      printf ("%g %d\n", yh, nc);
  }

  for (double xh = 0.8; xh >= 2*W; xh -= 0.01) {
    init_grid (1 << minlevel);
    unrefine (level > minlevel);
    refine (fabs(x - xh) < W && level < maxlevel);
    foreach() 
      pids[] = pid();
    boundary({pids});
    cells();
    squares ("pids", min = 0, max = npe() - 1);
    save ("mov.mp4");
    int nc = 0;
    foreach(reduction(+:nc)) { 
      if (!is_local(neighbor(1,0,0)))
	nc++;
      if (!is_local(neighbor(0,1,0)))
	nc++;
    }
    if (pid() == 0)
      fprintf (stderr,"%g %d\n", xh, nc);
  }
}
