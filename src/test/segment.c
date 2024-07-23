/**
# Segment traversal

This can be used to compute quantities (such as fluxes) or display
profiles along cross-sections defined by a segment. 

~~~gnuplot Endpoints of segments (symbols) and associated cells (arrows)
unset key
set size ratio -1
plot 'out' w  l, '< grep a log' u 2:3 pt 7 ps 1, 'seg' u 1:2 w l lw 4, \
  '< grep a log' u 4:5:6 w labels, '' u 4:5:($2-$4):($3-$5) w vectors
~~~
*/

#include "utils.h"

int main()
{
  init_grid (8);
  periodic (right);
  periodic (left);
  vector u[];
  foreach() {
    u.x[] = sin(2.[-1]*pi*x);
    u.y[] = sin(6.[-1]*pi*x);
  }
  output_cells();

  coord S[2][2] = {{{0.1, 0.1}, {0.76, 0.76}},
		   {{0.8, 0.1}, {0.2, 0.7}}};
  for (int i = 0; i < 2; i++) {
    fprintf (stderr, "%g %g\n%g %g\n\n",
	     S[i][0].x, S[i][0].y,
	     S[i][1].x, S[i][1].y);
    foreach_segment (S[i],r)
      for (int i = 0; i < 2; i++)
	fprintf (stderr, "a %g %g %g %g %d\n", r[i].x, r[i].y, x, y, _n);
  }
}
