/**
# Dense flow entry 

Similar to [stratified-flow.c](stratified-flow.c) but using [isopycnal.h](/src/layered/isopycnal.h) to model a single denser flow's entry into water.



<div class="figure">
<video controls>
<source src="https://video.wixstatic.com/video/e0f4ed_21d387d26864400fa89f0d3da419f1d7/480p/mp4/file.mp4" type="video/mp4">
</video>
<p class="caption">
Whole simulation.
</p>
</div>




*/



#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/isopycnal.h"

double rho_0 = 0.1;
double rho_1 = 0;
double * drho = (double []){ 0.1, 0 };

double nu_H = 1000.; // 1e-2;

int main()
{
  size (16e3);
  N = 256;
  breaking = 0.12;

  nl = 2;
  nu = 1e-4;

  G = 9.81;
  DT = 0.05;

  system ("rm -f plot-*.png");
  run();
}

event init (i = 0)
{
  foreach() {
    zb[] = x < 3.125e3 ? 50 - (x*0.064) : -150;
    foreach_layer() {
      if (x < 0.625e3)
	h[] = point.l <  nl/2 ? 20/(nl/2) : 0.002;
      else if (x < 0.7875e3)
        h[] = 0.002;
      else
	h[] = point.l >= nl/2 ? - zb[]/(nl/2) : 0.002;
    }
  }
}


scalar lambda_q[];
event viscous_term1 (i++)
{
  // Quadratic bottom friction,
  // see also: gerris-snapshot/doc/figures/diffusion.tm
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
#if ISOPYCNAL
	   "set pm3d map corners2color c2\n"
#else
	   "set pm3d map\n"
#endif
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [0:0.05]\n"
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
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], rho_0);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], rho_1);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event gnuplot (t += 6)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    setup (fp);
  //if (getenv ("DISPLAY")) {
  //  fprintf (fp, "set term x11\n");
  //  plot (fp);
  //}
  fprintf (fp,
	   "set term pngcairo font \",10\" size 800,500\n"
	   "set xrange [0:12.5]\n"
	   "set output 'plot-%05d.png'\n", i);
  plot (fp);
}

event end (t = 43000)
{
  fprintf (stderr, "done");
}
