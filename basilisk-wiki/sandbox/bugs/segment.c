/**
# Segment traversal bad name

The global scalar field *p* used often in solvers collides with the
coord *p* of *foreach_segment*. Proposal: use in
*foreach_segment* *r* instead of *p*.
*/

#include "utils.h"

scalar p[];

int main()
{
  init_grid (8);
  periodic (right);
  periodic (left);
  vector u[];
  foreach() {
    u.x[] = sin(2.*pi*x);
    u.y[] = sin(6.*pi*x);
  }
  output_cells();

  coord S[2][2] = {{{0.1, 0.1}, {0.76, 0.76}},
		   {{0.8, 0.1}, {0.2, 0.7}}};
  for (int i = 0; i < 2; i++) {
    fprintf (stderr, "%g %g\n%g %g\n\n",
	     S[i][0].x, S[i][0].y,
	     S[i][1].x, S[i][1].y);
    foreach_segment (S[i], r)
      for (int i = 0; i < 2; i++)
	fprintf (stderr, "a %g %g %g %g %d\n", r[i].x, r[i].y, x, y, _n);
  }
}
