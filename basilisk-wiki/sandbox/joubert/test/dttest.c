/**
# Time derivative computation

A bubble rises due to buoyancy. We want to chek different way of computing
time derivative for example vertical velocity here.
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "view.h"

int main (int argc, char * argv[])
{
  origin (-L0/2., 0);
  rho2 = 1./100.;
  mu1 = 1e-6;
  mu2 = mu1/rho1*rho2;
  f.sigma = 0.5;
  TOLERANCE = 1e-4;
  run();
}

event acceleration (i++)
{
  face vector av = a;
  foreach_face(y)
    av.y[] -= 9.8;
}

event init (t = 0)
{
  fraction (f, sq(x) + sq(y-0.15) - sq(0.125));
}

event logfile (i+=10) {
  static double xb0=0.15, xb20= 0.15, t0= 0.;
  double xb0n= 0.15, xb20n= 0.15, t0n= 0.;
  xb0n= xb20n= 0.15, t0n= 0.;
  double xb = 0., vb = 0., sb = 0., xb2 = 0., vb2 = 0., sb2 = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb) reduction(+:xb2) reduction(+:vb2) reduction(+:sb2)) {
//v1
    double dv = (1. - f[])*dv();
    xb += y*dv;
    vb += u.y[]*dv;
    sb += dv;
//v2
    if (f[] <= 1e-10) {
      xb2 += y*dv();
      vb2 += u.y[]*dv();
      sb2 += dv();
    }
  }
  if (i > 0) {
  double xdotb = (xb/sb - xb0/sb)/(t - t0),
	 xdotb2 = (xb2/sb2 - xb20/sb2)/(t - t0);
  double vdotb = fabs(xdotb),
	 vdotb2 = fabs(xdotb2);
  double xdotbn = (xb/sb - xb0n/sb)/(t - t0n),
	 xdotb2n = (xb2/sb2 - xb20n/sb2)/(t - t0n);
  double vdotbn = fabs(xdotbn),
	 vdotb2n = fabs(xdotb2n);
  fprintf (ferr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n", 
	t, sb, xb/sb, vb/sb, xdotb, sb2, xb2/sb2, vb2/sb2, xdotb2, xdotbn,
	xdotb2n, vdotbn, vdotb2n, vdotb, vdotb2, dt, perf.t, perf.speed);
    xb0=xb, xb20=xb2, xb0n=xb, xb20n=xb2, t0=t, t0n=t;
  }
}

event end (t = 0.5)
{
view (fov = 24, quat = {0,0,0,1}, tx = -0.00240032, ty = -0.472865, bg = {1,1,1}, width = 600, height = 600);
  box();
  draw_vof ("f");
  squares ("f", spread = -1);
  save ("f.png");
}

#if 1 // for debugging
event dumping (t +=0.1; t <= 0.5)
{
  char name[80];
  sprintf (name, "snap-%g", t);
  dump(file = name);
}
#endif

/**
## Results

![Bubble at end of simulation.](dttest/f.png)

~~~gnuplot Rise velocity as a function of time.
reset
set grid
set xlabel 'Time'
set key bottom right
plot 'log' u 1:4 w l t 'vy v1',\
'log' u 1:8 w l t 'vy v2'
~~~

~~~gnuplot Vertical position as a function of time .
reset
set grid
set xlabel 'Time'
set key bottom right
plot 'log' u 1:3 w l t 'y v1',\
'log' u 1:7 w l t 'y v2'
~~~

~~~gnuplot Rise velocity as a function of time .
reset
set grid
set xlabel 'Time'
set key bottom right
set logscale y
plot 'log' u 1:4 w l t ' vy',\
'log' u 1:14 w l t '|dy/dt| v1 static',\
'log' u 1:15 w l t '|dy/dt| v2 static'
~~~

~~~gnuplot Rise velocity as a function of time .
reset
set grid
set xlabel 'Time'
set key bottom right
set logscale y
plot 'log' u 1:4 w l t ' vy',\
'log' u 1:12 w l t '|dy/dt| v1 normal',\
'log' u 1:13 w l t '|dy/dt| v2 normal'
~~~
*/
