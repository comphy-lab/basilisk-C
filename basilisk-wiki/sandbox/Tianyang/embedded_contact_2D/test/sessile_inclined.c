#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "contact_ebm.h"
#include "vof_ebm.h"
#include "tension_ebm.h"

scalar f[], * interfaces = {f};

vector h[];

scalar contact_angle[];

double theta0 = 30;
int level_init = 7;

FILE *fp_stat, *fp_rad;

int main()
{
  size (2 [1]);
  init_grid (1 << level_init);
  
  const face vector muc[] = {0.04, 0.04};
  mu = muc;

  f.height = h;

  f.angle = contact_angle;
  
  f.sigma = 0.1;

  fclose (fopen ("ustats", "w"));
  fclose (fopen ("radius", "w"));

  double thetas[2] = {15., 165.};
  for (int i = 0; i <= 1; i++) {
    theta0 = thetas[i];
    run();
  }
}

event init (t = 0)
{
  foreach()
    contact_angle[] = pi*theta0/180.;

  solid_ebm (cs, fs, (y + x - 1.4));
  fraction_ebm (f, - (sq(x - 0.7) + sq(y - 0.7) - sq(0.28)));
  intersect_vof_solid (f, cs, fs);

  update_cfg (f, cs, fs, f.cfg);
  update_cet (f, cs, fs, f.cfg, f.cet, f.angle);

  fp_stat = fopen ("ustats", "a");
  fp_rad = fopen ("radius", "a");
}

scalar fpre[];

event logfile (i++, t <= 100)
{
  double df = change (f, fpre);
  double uxmax = normf(u.x).max;
  double uymax = normf(u.y).max;
  fprintf (fp_stat, "%d %g %g %g %g\n", i, t, uxmax, uymax, df);

  if (i > 1 && df < 1e-6)
    return 1;
}

event end (t = end)
{
  scalar kappa[];
  curvature_ebm (f, cs, kappa);

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf_no_cm(f).sum;

  fprintf (fp_rad, "%d %g %g %.3g\n", N, theta0, R/sqrt(V/pi), s.stddev);

  char name[90];
  sprintf (name, "profile_%d_%d.dat", N, (int) theta0);
  FILE *fp = fopen(name, "w");

  output_facets_ebm (f, cs, f.cet, fp);
  output_facets_ebm (f, cs, f.cet, stderr);

  fclose(fp_stat);
  fclose(fp_rad);
}


/**
## Results

~~~gnuplot Shapes of the interface.
reset
set size ratio -1
set xrange [0.:2.]
set yrange [0.:2.]
set object 1 circle at -1.5065665,-1.5065665 size 3.23064 lw 2 dt 2 lc rgb "black" front
set object 2 circle at 0.8354843,0.8354843 size 0.19836 lw 2 dt 2 lc rgb "black" front
plot 'profile_128_15.dat' w l lw 2 lc rgb "red" t "theta = 15", \
  'profile_128_165.dat' w l lw 2 lc rgb "blue" t "theta = 165", \
  1.4-x w l lw 3 lc rgb "black" t "", \
  3 w l lw 2 dt 2 lc rgb "black" t "Theoretical"
~~~

## See also

*/