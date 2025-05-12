/**
# Droplet spreading on a cylinder of the same size

![Animation of the droplet (blue) and cylinder (grey).](circle-droplet/movie.mp4)

~~~gnuplot Time evolution of the normalised total volume
set xlabel 'Time'
set ylabel 'V/V_0'
plot 'log' u 1:2 w l t ''
~~~
*/

#include "Chongsen/src/EBM_VOF/myembed.h"
#include "navier-stokes/centered.h"
#include "Chongsen/src/EBM_VOF/embed_contact.h"
#include "Chongsen/src/EBM_VOF/embed_two-phase.h"
#include "Chongsen/src/EBM_VOF/embed_tension.h"
#include "navier-stokes/perfs.h"
#include "profiling.h"

vector tmp_h[], o_interface[], ncc[], hnew1[];

double thetac = 30.;

//This doesn't take much time.
#define tend 0.001

double csTL = max(1e-2, VFTL);

#define MAXLEVEL  6
#define r0        0.2/200
#define l0        5*r0
#define ed        l0/pow(2,MAXLEVEL)

#define hsolid    0.2/200
#define hcircle   0.5826223/200
#define rr0       0.201191/200

#define rho01     1.
#define rho02     1.
#define mu01      7.5e-4
#define mu02      7.5e-4
#define sigma0    0.1

u.t[embed]  = dirichlet(0.);
u.n[embed]  = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

u.n[top]    = dirichlet(0.);
u.t[top]    = dirichlet(0.);

u.n[right]  = dirichlet(0.);
u.t[right]  = dirichlet(0.);


u.n[left]   = dirichlet(0.);
u.t[left]   = neumann(0.);

f[left] = neumann(0.);
cs[left] = neumann(0.);
tmp_c[left] = neumann(0.);

cs[bottom] = 1.;
f[bottom] = 0.;
tmp_c[bottom] = 0.;

int main()
{
  size(l0);
  origin(0,0);
  rho1 = rho01;
  rho2 = rho02;
  mu1 = mu01;
  mu2 = mu02;
  f.sigma = sigma0;
  tmp_c.height = tmp_h;
  tmp_c.hnew1 = hnew1;
  tmp_c.oxyi = o_interface;
  tmp_c.nc = ncc;
  N = 1 << MAXLEVEL;
  run();
}

event init (t = 0)
{
  solid(cs,fs, (sq(y - hsolid) + sq(x) - sq(r0)));
  cleansmallcell (cs, fs, csTL);

  fraction(f, - (sq(x) + sq(y - hcircle) - sq(rr0)));
    
  foreach(){
    f[]*=cs[];
    contact_angle[] = thetac;
  }
}

event volume (i = 0; t<=tend; i+=10)
{  
  double v = 0;
  foreach(reduction(+:v))
    if (dv() > 0)
      v += dv()/cm[]*f[];
  static double v0;
  if (i == 0) v0 = v;
  fprintf (stderr, "%.16f %.16f\n", t, v/v0);
}

#include "view.h"
event movie (t += tend/300.)
{
  view (quat = {0.000, 0.000, 0.000, 1.000},
	fov = 30, near = 0.01, far = 1000,
	tx = 0, ty = -0.358, tz = -1.8,
	width = 600, height = 600);
  draw_vof (c = "cs", s = "fs", filled = -1, fc = {0.5,0.5,0.5});
  draw_vof (c = "tmp_c", fc = {0.447,0.717,0.972}, filled = 1);
  mirror (n = {1,0}) {
    draw_vof (c = "cs", s = "fs", filled = -1, fc = {0.5,0.5,0.5});
    draw_vof (c = "tmp_c", fc = {0.447,0.717,0.972}, filled = 1);
  }  
  save ("movie.mp4");
}

#if TREE
event adapt (i++) {
  scalar sf1[];
  foreach() {
    sf1[] = (8. * tmp_c[] +
	     4. * (tmp_c[-1] + tmp_c[1] +
		   tmp_c[0, 1] + tmp_c[0, -1] +
		   tmp_c[0, 0, 1] + tmp_c[0, 0, -1]) +
	     2. * (tmp_c[-1, 1] + tmp_c[-1, 0, 1] + tmp_c[-1, 0, -1] + tmp_c[-1, -1] +
		   tmp_c[0, 1, 1] + tmp_c[0, 1, -1] + tmp_c[0, -1, 1] + tmp_c[0, -1, -1] +
		   tmp_c[1, 1] + tmp_c[1, 0, 1] + tmp_c[1, -1] + tmp_c[1, 0, -1]) +
	     tmp_c[1, -1, 1] + tmp_c[-1, 1, 1] + tmp_c[-1, 1, -1] + tmp_c[1, 1, 1] +
	     tmp_c[1, 1, -1] + tmp_c[-1, -1, -1] + tmp_c[1, -1, -1] + tmp_c[-1, -1, 1]) / 64.;
    sf1[] += cs[];
  }
  adapt_wavelet ({sf1}, (double[]){1e-5}, minlevel = max(3, MAXLEVEL - 7), maxlevel = MAXLEVEL);
}
#endif
