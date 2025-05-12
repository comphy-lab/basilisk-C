/**
# Waves refraction on a shallow elliptical reef platform */

#include "grid/multigrid.h"
#include "green-naghdi.h"

#define HS 1.5
#define TS 12
#define DEPTH 1.5
#define XYR 0.70
#define XW 708
#define RROTATE 0

int main()
/** The domain is 1536m long (reef flat is 1000m x 700m) 
This is a low resolution run (N=256 gives 6m cells) 
for higher resolution change N to 512 (3m) or 1024 (1.5m) */
{
  L0 = 1536;
  X0 = -L0/2.;
  Y0 = -L0/2.;
  G = 9.81;
  N = 256;
  breaking = 1;
  run();
}

/**
monochromatic wave field is generated at the left boundary */

u.x[left]  = - radiation ((1*HS)*sin(2.*pi*t/TS));
u.x[right] = + radiation (0);

/**
Define average and maximum stats variables */

scalar maxa[];
scalar SH[];
scalar SHc[];
scalar AH[];
scalar Vel[];
scalar SVel[];
scalar AVel[];
scalar vmax[];
scalar vmin[];
scalar Sux[];
scalar Suy[];
scalar Aux[];
scalar Auy[];

event init (i = 0)
{
  
  /**
  bathymetry is an elongate reef platform with ~35 degree slope and flat 1.5m deep reef flat 
  */

  double h0 = 100;
  double cosa = cos (RROTATE*pi/180.), sina = sin (RROTATE*pi/180.);
  foreach() {
    double xr = x*cosa - y*sina, yr = x*sina + y*cosa;
    double z0 = 0; /** xr >= -5.82 ? (5.82 + xr)/50. : 0.; */
    double zs = sq(xr/XW) + sq(yr/(XW*XYR)) <= 1. ?
      1 + 50*(sqrt(1 - sq(xr/XW) - sq(yr/(XW*XYR)))) : 0.;
    zb[] = min(-h0+(max((z0 + (zs*5.5) - h0), 0)),0) - DEPTH;

    h[] = max (0., -zb[]);
    eta[] = h[] > dry ? h[] + zb[] : 0;
    maxa[] = 0.;
  }
}

/**
quadratic bottom firction with extra friction added near right 
boundary to prevent reflection */

event friction (i++) {
  foreach() 
       {
      double a = x < 650 ? h[] < dry ? HUGE : 1. + 0.01*dt*norm(u)/h[] :
		           h[] < dry ? HUGE : 1. + 2*(x - 650.)*dt*norm(u)/h[];
      foreach_dimension()
	u.x[] /= a;
	Vel[] = norm(u);
    }
}

event logfile ( t+= 1) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);
}

/**
MOVIES */
event movies (t += 0.5) {
  static FILE * fp = popen ("ppm2mpeg > Vel.mpg", "w");
  scalar m[], etam[];
  foreach() {
    etam[] = Vel[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, fp, mask = m, min = 0, max = 2, n = 512, linear = true);
}

/**
max and time-averaged stats are calculated after the wave field is 
developed on the reef platform */

event maximum (t = 600; i++) {
  foreach() {
 /** Set amax as highest wave amplitude */
    if (fabs(eta[]) > maxa[])
      maxa[] = fabs(eta[]);
 /** Set SH as the sum water level */
    if (h[] > dry)
      SH[] = SH[] + (h[] + zb[]);
 /** Set SHc as the number of sum water level */
    if (h[] > dry)
      SHc[] = SHc[] + 1;
 /** Set AH as the average water level */
    if (h[] > dry)
      AH[] = SH[] / SHc[];
 /** Set SVel as the sum velocity */
    if (h[] > dry)
      SVel[] = SVel[] + (Vel[]);
 /** Set AH as the average velocity */
    if (h[] > dry)
      AVel[] = SVel[] / SHc[];
/** Set vmax as highest velocity */
    if (h[] > dry && Vel[] > vmax[])
      vmax[] = Vel[];
/** Set vmax as highest velocity */
    if (h[] > dry && Vel[] < vmin[])
      vmin[] = Vel[];
/** Set Sux as the sum u.x */
    if (h[] > dry)
      Sux[] = Sux[] + u.x[];
/** Set Suy as the sum u.y */
    if (h[] > dry)
      Suy[] = Suy[] + u.y[];
/** Set Aux as the average u.x */
    if (h[] > dry)
      Aux[] = Sux[] / SHc[];
/** Set Auy as the average u.y */
    if (h[] > dry)
      Auy[] = Suy[] / SHc[];
 }
}

/**
grid data is output every 300 seconds and at the end */

event out (t += 300) {
  char name[20]; sprintf (name, "outzxy-%g.out", t);
  FILE * fp = fopen (name, "w");
  output_field ({eta,zb,maxa,AH,AVel,Vel,vmax,vmin,u,Aux,Auy}, fp, linear = true);
}

event end (t = 1209) {
  FILE * fp = fopen ("end", "w");
  output_field ({eta,zb,maxa,AH,AVel,Vel,vmax,vmin,u,Aux,Auy}, fp, linear = true);

}

/**
  ~~~gnuplot Bathymetry
  set term @PNG enhanced size 640,640 font ",8"
  set pm3d map
  set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
                        0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392,	  \
                        0.625 1 0.9333 0, 0.75 1 0.4392 0,		  \
                        0.875 0.9333 0 0, 1 0.498 0 0 )
  set size ratio -1
  set xlabel 'x (m)'
  set ylabel 'y (m)'
  splot [-768:768][-768:768]'end' u 1:2:4 w pm3d t ''
  ~~~

~~~gnuplot Instantaneous wave field
a0 = 1
splot [-768:768][-768:768]'end' u 1:2:($3/a0) w pm3d t ''
~~~

~~~gnuplot Maximum wave amplitude
a0 = 1
splot [-768:768][-768:768]'end' u 1:2:($5/a0) w pm3d t ''
~~~

~~~gnuplot mean water level
a0 = 1
splot [-768:768][-768:768]'end' u 1:2:($6/a0) w pm3d t ''
~~~

~~~gnuplot mean velocity
a0 = 1
splot [-768:768][-768:768]'end' u 1:2:($7/a0) w pm3d t ''
~~~

~~~gnuplot velocity
a0 = 1
splot [-768:768][-768:768]'end' u 1:2:($8/a0) w pm3d t ''
~~~

*/

