/**
# Layered flow entry 

Similar to [stratified-flow.c](stratified-flow.c) but adjusted to be more similar to a pyroclastic density current. Thin dense layer at base of flow, overlayed by lighter.

This is adapted from [overflow.c](/src/test/overflow.c) to simulate a flow of variable density layers into water, using the [Boussinesq buoyancy](/src/layered/dr.h) module of the [non-hydrostatic](/src/layered/nh.h) [layered scheme](/src/layered/hydro.h).




<div class="figure">
<video controls>
<source src="https://video.wixstatic.com/video/e0f4ed_f21751dc0af142399cbc0e2238c94d9c/360p/mp4/file.mp4" type="video/mp4">
</video>
<p class="caption">
Close-up of entry, displaying relative density.
</p>
</div>




<div class="figure">
<video controls>
<source src="https://video.wixstatic.com/video/e0f4ed_74fd72659982463fa6070be94bc30cf2/360p/mp4/file.mp4" type="video/mp4">
</video>
<p class="caption">
Close-up of entry, displaying u.x.
</p>
</div>




<div class="figure">
<video controls>
<source src="https://video.wixstatic.com/video/e0f4ed_813e764cbb95432a8f174db1d766bc29/480p/mp4/file.mp4" type="video/mp4">
</video>
<p class="caption">
Whole simulation.
</p>
</div>




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

double rho_0 = 0.4;
double rho_1 = -0.2;
double rho_2 = -0.6;
double rho_3 = 0.;

double nu_H = 1000.; // 1e-2;

int main()
{
  size (16e3);
  N = 1024;//2048;
  breaking = 0.12;

  nl = 40;
  nu = 0.001;//0.1;//1e-3;

  G = 9.81;
  DT = 0.05;//0.025;

  cell_lim = mono_limit;
  system ("rm -r gnuplot");
  system ("mkdir gnuplot");
  system ("mkdir gnuplot/all");
  system ("mkdir gnuplot/closeup");
  system ("mkdir gnuplot/closeup/density");
  system ("mkdir gnuplot/closeup/ux");
  //system ("rm -f plot-*.png");
  run();
}


/**
Constant slope down to 150m depth, flow initialised as block 20m tall up slope split into three densities.
*/

event init (i = 0)
{
  foreach() {
    zb[] = x < 2000 ? 50 - (x*0.05) : -50; // Constant slope down to minimum
    foreach_layer() {
      if (x < 625) // Flow
	//h[] = 20/nl;
        h[] = point.l < nl/2 ? 20/(nl/2) : 0.0002;
      else if (x < 1000)
        h[] = 0.0002; // Gap
      else
	h[] = - zb[]/nl; //Water
    }

    foreach_layer() {
      if (x < 425) {
        if(point.l < nl/12) {
          T[] = rho_0; // Bottom sixth of flow
          //u.x[] = 0.;
        }
        else if(point.l < nl/3) {
          T[] = rho_1; // Middle third of flow
          //u.x[] = 5.;
        }
        else if(point.l < nl/2) {
          T[] = rho_2; // Top third of flow
          //u.x[] = 10.;
        }
      }
      else if (x < 525) { // between front and main flow
        if(point.l < nl/4)
          T[] = rho_1;
        else
          T[] = rho_2;
      }
      else if (x < 625) // front
        T[] = rho_2;
      else		// Main water body
        T[] = rho_3; // Water
    }
  }
}


scalar lambda_q[];
event viscous_term1 (i++)
{
  double Cf = 0.05;//1e-2;
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




/**
 Animations including density, u.x...
*/


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
	   "set xlabel 'x (km)'\n"
	   "set ylabel 'depth (m)'\n"
	   //"set xrange [0:12.5]\n"
	   //"set yrange [-50:75]\n"
	   );
}

#define minute 60.
void plot_ux (FILE * fp)
{
  fprintf (fp,
	   "set cbrange [-0.1:30]\n" // u.x
	   "set title 't = %.2f minutes'\n"
	   "sp '-' u ($1/1e3):2:3\n",
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
void plot_density (FILE * fp)
{
  fprintf (fp,
	   "set cbrange [-0.6:0.4]\n" // Density
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
	   "set term pngcairo font \",10\" size 1000,400\n"
	   "set xrange [0:12.5]\n"
	   "set yrange [-50:75]\n"
	   "set output 'gnuplot/all/plot-%06d.png'\n", i);
    plot_density (fp);
  }
  if (t < 180) {
    fprintf (fp,
	   "set term pngcairo font \",10\" size 1000,300\n"
	   "set xrange [0.7:1.5]\n"
	   "set yrange [-50:50]\n"
	   "set output 'gnuplot/closeup/density/plot-%06d.png'\n", i);
    plot_density (fp);
  }
  if (t < 180) {
    fprintf (fp,
	   "set term pngcairo font \",10\" size 1000,300\n"
	   "set xrange [0.7:1.5]\n"
	   "set yrange [-50:50]\n"
	   "set output 'gnuplot/closeup/ux/plot-%06d.png'\n", i);
    plot_ux (fp);
  }
}


event end (t = 30*minute)
{
  fprintf (stderr, "\n\nDone\n");
}