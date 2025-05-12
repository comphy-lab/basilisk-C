/**
# The patch 'Automatic boundary conditions' seems to be incompatible
  with some of the procedure below...*/

#include "embed.h"
#include "run.h"

int lvl;
int main()
{
  size(5.);
  origin(-2.5, -2.5);
  for (lvl = 5; lvl <= 5; lvl++) {
    N = 1 << (lvl);
    init_grid (N);
    run();
  }
}  

static double bilinear_corrected_embed_gradient (Point point,
						 coord * o, coord * n)
{
  assert (cs[] > 0. && cs[] < 1.); //only in mixed cells
  double area = 0.;
  area = embed_geometry (point, o, n);
  return area;
}


event defaults (t = 0) {
vertex scalar phiv[];
  foreach_vertex()
    phiv[] = sq(x) + sq(y) - sq(1.);
  fractions (phiv, cs, fs);
  fractions_cleanup (cs, fs);

  // this is not a good way to do things, but this should not crash qcc
  coord o, n, np;
  foreach() 
    if(cs[] > 0. && cs[] < 1.) {
      // coord o, n, np; // this works
      bilinear_corrected_embed_gradient (point, &o, &n);
      foreach_dimension()
        np.x = sign(-n.x);
      printf(" %g %g  %g  %g \n", o.x, o.y, np.x, np.y);
    }
}