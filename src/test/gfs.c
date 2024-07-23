#include "utils.h"

int main()
{
  int depth = 6;
  origin (-0.5, -0.5, -0.5);
  init_grid (8);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y + z*z)));
  
  scalar a[];
  vector u[];
  double k = 2.*pi;
  foreach()
    u.x[] = u.y[] = a[] = sin(k*x)*cos(k*y);

  FILE * fp = fopen ("gfsi.gfs", "w");
  output_gfs (fp);
  fclose (fp);
}
