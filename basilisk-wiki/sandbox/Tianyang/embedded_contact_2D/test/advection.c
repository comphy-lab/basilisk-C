#include "grid/multigrid.h"
#include "embed.h"
#include "run.h"
#include "timestep.h"

vector u[];
face vector uf[];

event defaults (i = 0)
{
  CFL = 0.8;
  
#if TREE
  uf.x.refine = refine_face_solenoidal;
#if EMBED
  uf.x.refine = refine_face;
  foreach_dimension()
    uf.x.prolongation = refine_embed_face_x;
  for (scalar s in {u}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
#endif // EMBED
#endif // TREE

  foreach()
    foreach_dimension()
      dimensional (u.x[] == Delta/t);
}

double dtmax;

double t0 = 1. [0, 1];

event init (i = 0)
{
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  
  dtmax = DT;
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

event vof (i++,last);

#include "contact_ebm.h"
#include "vof_ebm.h"

double radius, distance;
double theta0 = 30.;

scalar f[], * interfaces = {f};

vector h[];

scalar contact_angle[];

int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << 7);

  f.height = h;

  f.angle = contact_angle;
  
  radius = 0.1409;
  distance = (sqrt(2*(1-cos(pi/6.)))*radius);

  run();
}

event init (i = 0)
{
  foreach() {
    u.x[] =  2.*pi*y/(1. [0,1]);
    u.y[] = -2.*pi*x/(1. [0,1]);

    contact_angle[] = pi*theta0/180.;
  }

  solid_ebm (cs, fs, (sq(x) + sq(y) - sq(radius)));
  fraction_ebm (f, - (sq(x) + sq(y-distance) - sq(radius)));
  intersect_vof_solid (f, cs, fs);

  update_cfg (f, cs, fs, f.cfg);
  update_cet (f, cs, fs, f.cfg, f.cet, f.angle);
}

event logfile (t = 1) {}

event end (t = end)
{
  char name[90];
  sprintf (name, "profile_%d_%d.dat", N, (int) theta0);
  FILE *fp = fopen(name, "w");
  output_facets_ebm (f, cs, f.cet, fp);
  fclose (fp);

  output_facets_ebm (f, cs, f.cet, stderr);
}

/**
## Results
~~~gnuplot Shapes of the interface.
reset
set size ratio -1
set xrange [-0.4:0.4]
set yrange [-0.4:0.4]
set object 1 circle at 0,0 size 0.1409 lw 3 lc rgb "#808080"
plot 'profile_128_30.dat' w l lw 3 lc 1 t "theta = 30"
~~~

## See also

*/