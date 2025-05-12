/**
# LS_reinit() test case

Extreme test case */

#define BICUBIC 1
#define BGHOSTS 2
#include "popinet/distance_point_ellipse.h"

#include "embed.h"
#include "alimare/alex_functions.h"
#include "alimare/LS_reinit.h"
#include "view.h"

double perturb (double x, double y, double eps, coord center)
{
  return eps + sq(x - center.x) + sq(y - center.y);
}

void draw_isolines(scalar s, double smin, double smax, int niso, int w)
{
  scalar vdist[];
  cell2node(s,vdist);
  boundary ({vdist});
  for (double sval = smin ; sval <= smax; sval += (smax-smin)/niso)
    isoline ("vdist", sval, lw = w);
}

scalar dist[];
scalar * level_set = {dist};

int main() 
{
  origin (-5., -5.);
  L0 = 10;

  int MAXLEVEL = 9;  
  init_grid (1 << MAXLEVEL);
  
  coord  center_perturb = {-3.5,-2.};
  foreach(){
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    dist[] = (1.0 + 0.15*cos(6.*theta) - r)*perturb(x,y, 0.1, center_perturb)/3.;
  }
  boundary({dist});

  view (fov = 15.);
  squares ("dist", map = cool_warm, min = -1, max = 1);
  draw_isolines (dist, -1., 1., 10, 1);
  save ("dist_init.png");

  LS_reinit (dist, it_max = 200);
  squares ("dist", map = cool_warm, min = -1, max = 1);
  draw_isolines (dist, -1., 1., 10, 1);
  save ("dist_first_reinit.png");

  squares ("dist", map = cool_warm, min = -1, max = 1);
  draw_isolines (dist, -1, 1, 10, 1);
  save ("dist_final.png");
}

/**
We show here the initial and final level-set for the same isovalues.

![Initial level-set](reinit_circle/dist_init.png) 

![first reinit level-set](reinit_circle/dist_first_reinit.png) 

![Final level-set](reinit_circle/dist_final.png)

## References

~~~bib
@article{russo_remark_2000,
  title = {A remark on computing distance functions},
  volume = {163},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Russo, Giovanni and Smereka, Peter},
  year = {2000},
  pages = {51--67}
}
~~~
*/
