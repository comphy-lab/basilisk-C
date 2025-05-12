/**
# Layered flow entry 

This is adapted from [overflow.c](/src/test/overflow.c) to simulate a flow of variable density layers into water, using the [Boussinesq buoyancy](/src/layered/dr.h) module of the [non-hydrostatic](/src/layered/nh.h) [layered scheme](/src/layered/hydro.h).


<video controls>
<source src="https://video.wixstatic.com/video/e0f4ed_25ed64db037049a79769475c0f7b14e5/360p/mp4/file.mp4" type="video/mp4">
</video>

*/



#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"

/**
Setting $\Delta \rho$ directly.
*/

#define drho(T) T
#include "layered/dr.h"
#include "layered/remap.h"
//#include "layered/perfs.h"

double rho_0 = 0.7;
double rho_1 = 0.5;
double rho_2 = 0.;
double rho_3 = 0.4;

double nu_H = 1000.;

int main()
{
  size (16e3);
  N = 1024;
  breaking = 0.12;
  nl = 40;
  nu = 1e-4;
  G = 9.81;
  DT = 0.05;
  cell_lim = mono_limit;
  system ("rm -r gnuplot");
  system ("mkdir gnuplot");
  system ("mkdir gnuplot/all");
  system ("mkdir gnuplot/closeup");
  run();
}

/**
Constant slope down to 150m depth, flow initialised as block 20m tall up slope split into three densities.
*/

event init (i = 0)
{
  foreach() {
    zb[] = x < 3125 ? 50 - (x*0.064) : -150; // Constant slope down to minimum
    foreach_layer() {
      if (x < 625) // Flow
	h[] = point.l < nl/2 ? 20/(nl/2) : 0.0002;
      else if (x < 787.5)
        h[] = 0.0002; // Gap
      else
	h[] = - zb[]/nl; //Water
    }

    foreach_layer() {
      if (x < 625) {
        if(point.l < nl/6)
          T[] = rho_0; // Bottom third of flow
        else if(point.l < nl/3)
          T[] = rho_1; // Middle third of flow
        else if(point.l < nl/2)
          T[] = rho_2; // Top third of flow
      }
      else
        T[] = rho_3; // Water
    }
  }
}



scalar lambda_q[];
event viscous_term1 (i++)
{
  double Cf = 1e-2;
  lambda_b = lambda_q;
  foreach() {
    double au = norm(u);
    lambda_q[] = au < 1e-6 ? HUGE : nu*h[]/(Cf*au);
  }
}

event viscous_term2 (i++)
{
  scalar d2u[];
  foreach_layer() {
    foreach()
      d2u[] = (u.x[1] + u.x[-1] - 2.*u.x[])/sq(Delta);
    foreach()
      u.x[] += dt*nu_H*d2u[];
  }
  boundary ((scalar *){u});
}



/**
### Outputs */

event logfile (i += 10)
{
  fprintf (stderr, "%g %g %d\n", t, dt, mgp.i);
}



void setup (FILE * fp)
{
  fprintf (fp,
	   "set pm3d map\n"
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [0:0.7]\n"
	   "set xlabel 'x (km)'\n"
	   "set ylabel 'depth (m)'\n"
	   "set xrange [0:12.5]\n"
	   "set yrange [-150:75]\n"
	   );
}

#define minute 60.
void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f minutes'\n"
	   "sp '-' u ($1/1e3):2:4\n",
	   t/minute);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event gnuplot (t += 1)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    setup (fp);
/* Uncomment for animation during run.
  if (getenv ("DISPLAY")) {
    fprintf (fp, "set term x11\n");
    plot (fp);
  }
*/
  int timeint = (int)t; // This is dodgy
  if (timeint%6 == 0) {
  	fprintf (fp,
		   "set term pngcairo font \",10\" size 800,500\n"
		   "set xrange [0:12.5]\n"
		   "set yrange [-150:75]\n"
		   "set output 'gnuplot/all/plot-%06d.png'\n", i);
  	plot (fp);
  }
  if (t < 180) {
  	fprintf (fp,
		   "set term pngcairo font \",10\" size 1000,300\n"
		   "set xrange [0.7:1.5]\n"
		   "set yrange [-50:50]\n"
		   "set output 'gnuplot/closeup/plot-%06d.png'\n", i);
	  plot (fp);
  }
}


event end (t = 30*minute)
{
  fprintf (stderr, "\n\nDone\n");
}