#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "contact_ebm.h"
#include "vof_ebm.h"
#include "tension_ebm.h"

double alpha_cylinder = 0.8305581443, alpha_drop = 0.7402381825;
double r_ratio = 1.094659125 [0];

double r_cylinder = 0.1409 [1];
double ypos = -0.3359375;
double r_drop, dis_cyc;

scalar f[], * interfaces = {f};

vector h[];

scalar contact_angle[];

double theta0 = 30;
int level_init = 7;

FILE *fp_stat, *fp_rad;

int main()
{
  origin (0., ypos);
  init_grid (1 << level_init);
  
  const face vector muc[] = {0.04, 0.04};
  mu = muc;

  f.height = h;

  f.angle = contact_angle;
  
  f.sigma = 0.1;

  r_drop = r_cylinder*r_ratio;
  dis_cyc = r_cylinder*cos(alpha_cylinder) + r_drop*cos(alpha_drop);

  fclose (fopen ("ustats", "w"));
  fclose (fopen ("radius", "w"));

  theta0 = 30.;
  run();

}

event init (t = 0)
{
  foreach()
    contact_angle[] = pi*theta0/180.;

  solid_ebm (cs, fs, (sq(x) + sq(y) - sq(r_cylinder)));
  fraction_ebm (f, - (sq(x) + sq(y - dis_cyc) - sq(r_drop)));
  intersect_vof_solid (f, cs, fs);

  update_cfg (f, cs, fs, f.cfg);
  update_cet (f, cs, fs, f.cfg, f.cet, f.angle);

  fp_stat = fopen ("ustats", "a");
  fp_rad = fopen ("radius", "a");
}

scalar fpre[];

event logfile (i++; t <= 10)
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
set xrange [-0.3:0.3]
set yrange [-0.2:0.4]
set object 1 circle at 0,0 size 0.1409 lw 3 lc rgb "#808080"
set object 2 circle at 0,0.0982831 size 0.190553 lw 2 dt 2 lc rgb "black" front
plot 'profile_128_30.dat' w l lw 2 lc rgb "red" t "theta = 30", \
  'profile_128_30.dat' using (-$1):($2) w l lc rgb "red" lw 2 t "", \
  0.5 w l lw 2 dt 2 lc rgb "black" t "Theoretical"
~~~

## See also

*/