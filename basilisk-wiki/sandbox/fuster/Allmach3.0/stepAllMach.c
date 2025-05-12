/**
# A Mach-3 wind tunnel with a step

This test included in the review paper of [Woodward & Colella, 1984](#woodward1984)
shows the interaction between a shock wave and a step.

In this case we use the Bell-Colella-Glaz advection scheme using the minmod2 slope limiter

This script obtain the isocontours below that can be compared with [Woodward & Colella, 1984](#woodward1984) 
*/

#include "two-phase-compressible.h"

# define LEVEL 10

double rho0 = 1.4;
double p0 = 1.;
double tend = 2.;

/** Boundary conditions */

q.n[left]  	= dirichlet(3*rho0);
q.t[left]  	= dirichlet(0.);

q.n[right] 	= neumann(0.);
fE1[right]	= neumann(0.);
frho1[right] 	= neumann(0.);
  
q.n[bottom] 	= dirichlet(0.);
q.t[bottom] 	= neumann(0.);
	
q.n[top] 	= dirichlet(0.);
q.t[top] 	= neumann(0.);
	
frho1[left]	= dirichlet(rho0);

/** Main program */
int main() {

  gamma1 = 1.4; 

  size (3. [0]);
  DT = HUGE [0];
  init_grid (1 << LEVEL);
  run(); 
}

event init (t = 0)
{
  mask (y > 1. ? top : none); 
  mask ( (x > 0.6 && y < 0.2) ? bottom : none );

  f.gradient = minmod2;
  theta = 1.3;

  foreach() {
    f[]         = 1.;
    p[]		= p0;
    frho1[]	= rho0;
    frho2[] 	= 0.;
    q.x[] 	= 3*frho1[];
    q.y[] 	= 0.;
    fE1[]	= p[]/(gamma1 - 1.) + 0.5*pow(q.x[],2)/frho1[];
    fE2[] 	= 0.;
  }
}

/** We test the grid adaptation */
event adapt (i++) {
  adapt_wavelet((scalar *){p},(double[]){0.1},maxlevel = LEVEL);
}
 
/**
## Output */

event logfile (i++) {
  if (i == 0)
    fprintf (ferr, "t dt grid->tn perf.t perf.speed\n");
  fprintf (ferr, "%g %g %ld %g %g\n", t, dt, grid->tn, perf.t, perf.speed);
}

/**
## Movies */

event movies (i += 5)
{
  static FILE * fp1 = popen ("ppm2mp4 rho.mp4", "w");
  output_ppm (frho1, fp1, box = {{0.,0.},{3.,1.}},
	      linear = true, spread = 2, n = 512);

  scalar l[];
  foreach()
    l[] = level;

  static FILE * fp2 = popen ("ppm2mp4 level.mp4", "w");
  output_ppm (l, fp2, box = {{0.,0.},{3.,1.}},
	      linear = false, min = 0, max = LEVEL, n = 512); 
}

/**
![Evolution of the density field](stepAllMach/rho.mp4)

![Evolution of the level of refinement](stepAllMach/level.mp4)

We plot the evolution of iso-density contours. 

![Evolution of iso-density contours](stepAllMach/isovalue.gif)
*/

event isovalue (t += tend/20.; t <= tend)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (t == 0.)
    fprintf (fp,
	     "set term gif animate crop\n"
	     "set output 'tmp.gif'\n"
	     "set xrange [0:3]\n"
	     "set yrange [0:1]\n"
	     "unset surface\n"
	     "set view map\n"
	     "unset key\n"
	     "set size ratio -1\n"
	     "set contour base\n"
 	     // "set cntrlabel onecolor\n"
 	     "unset clabel\n"
	     "set cntrparam levels incremental 0.267,(6.75-0.267)/30,6.75\n"
	     );
  fprintf (fp, "splot '-' w l lt -1\n");
  output_field ({rho}, fp, n = 512, linear = true);
  fputs ("e\n", fp);
  fflush (fp);
}

event end (t = tend)
{
  FILE * fp = fopen ("rho.end", "w");
  output_field ({rho}, fp, n = 512, linear = true);
  fclose (fp);

  system ("gifsicle --optimize --delay 10 --loopcount=0 tmp.gif > isovalue.gif "
	  "&& rm -f tmp.gif");
}

/**
## References

~~~bib
@article{woodward1984,
  title={The numerical simulation of two-dimensional 
  fluid flow with strong shocks},
  author={Woodward, Paul and Colella, Phillip},
  journal={Journal of computational physics},
  volume={54},
  number={1},
  pages={115--173},
  year={1984},
  publisher={Elsevier},
  url={http://seesar.lbl.gov/anag/publications/colella/A3_1984.pdf}
}
~~~
*/
