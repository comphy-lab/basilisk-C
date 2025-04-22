/**
# Double Mach reflection of a Mach 10 shock from a wall

This test, included in the review paper of [Woodward & Colella,
1984](#woodward1984), shows the interaction between an oblique shock
wave and a wall.

In this case we use the Bell-Colella-Glaz advection scheme with the
minmod2 slope limiter.

The isocontours below can be compared with Figure 4 of [Woodward
& Colella, 1984](#woodward1984) using the same effective resolution
$\Delta x = 1/120$.

~~~gnuplot
set term svg enhanced size 640,280 font ",10"
set xrange [0:3]
set yrange [0:1]
unset surface
set view map
unset key
set size ratio -1
set contour base
# set cntrlabel onecolor (for gnuplot >= 4.7)
unset clabel
set cntrparam levels incremental 1.731,(20.92 - 1.731)/30,20.92
splot 'rho.end' w l lt -1
~~~

![Evolution of the density field](inclined-shock/rho.mp4)

![Evolution of the level of refinement](inclined-shock/level.mp4)
*/

#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"

/**
Problem parameters */

double rhoL, rhoR = 1.4;
double pL, pR = 1.;
double tend = 0.2;
double uL, vL;
double Ldomain = 4.26667 [1];
int LEVEL = 9;

/**
Boundary Conditions */

q.n[right] 	= neumann(0.);
fE1[right] 	= neumann(0.);
frho1[right] 	= neumann(0.);

q.n[left] 	= dirichlet(uL*rhoL);
q.t[left] 	= dirichlet(vL*rhoL);
p[left]         = dirichlet(pL);
frho1[left] 	= dirichlet(rhoL); 
f[left] 	= dirichlet(1.);
frho2[left] 	= dirichlet(0.);

q.t[top] 	= dirichlet(uL*rhoL);
q.n[top] 	= dirichlet(vL*rhoL);
p[top]          = dirichlet(pL); 
frho1[top] 	= dirichlet(rhoL);
f[top] 		= dirichlet(1.);
frho2[top] 	= dirichlet(0.);

q.n[bottom] 	= x < 1./6. ? neumann(0.) : dirichlet(0.);
q.t[bottom] 	= neumann(0.);
fE1[bottom]     = neumann(0.);
frho1[bottom] 	= neumann(0.);

/**
## Main program */

int main()
{
  gamma1 = 1.4; 
  TOLERANCE = 1.e-6 [*];
  size (Ldomain);
  init_grid (1 << LEVEL);
  run();
}

/**
## Initial conditions */

event init (i = 0)
{   
  double Ms = 10., gr = (gamma1 +1.)/(gamma1 -1.);

  /** Pre-shocked state: */
  
  double cR = sqrt(pR/rhoR*gamma1);   
   
  /** Post-shocked state: */
  
  pL = pR*(2.*gamma1*sq(Ms) - (gamma1 - 1.))/(gamma1 + 1.);
  double VL = Ms*cR*(1. - (((gamma1 - 1.)*sq(Ms) + 2.)/((gamma1 + 1.)*sq(Ms))));
  uL =  VL*sqrt(3.)/2.; // V*cos(30)
  vL = - VL/2.; // V*sin(30)

  scalar m[];
  fraction (m, -sqrt(3.)*(x - 1./6.) + y);
  
  rhoL = rhoR*(gr*pL/pR + 1.)/(gr + pL/pR);

  foreach() {
    f[]     = 1.;
    p[]     = m[]*pL + (1. - m[])*pR;
    frho1[] = m[]*rhoL + (1. - m[])*rhoR;
    frho2[] = 0.;
    q.x[]   = m[]*frho1[]*uL;
    q.y[]   = m[]*frho1[]*vL;
    fE1[]   = p[]/(gamma1 - 1.) + 0.5*(sq(q.x[]) + sq(q.y[]))/frho1[];
    fE2[]   = 0.;
  }
}

/**
## Grid adaptation */

event adapt (i++) {
  adapt_wavelet ((scalar *){p}, (double[]){0.01}, maxlevel = LEVEL);
}
 
/**
## Output */

event logfile (i++)
{
  if (i == 0) {
    fprintf (stdout, "grid->tn perf.t perf.speed\n");
    fprintf (stderr, "t dt pmax pmin rhomax rhomin\n");
  }
  fprintf (stdout, "%ld %g %g\n", grid->tn, perf.t, perf.speed);
  stats sp = statsf(p), sr = statsf(rho);
  fprintf (stderr, "%g %g %g %g %g %g\n", t, dt,
	   sp.min, sp.max, sr.min, sr.max);
}

/**
## Movies */

event movies (i += 5)
{
  output_ppm (frho1, file = "rho.mp4", box = {{0.,0.},{3.,1.}},
	      linear = true, spread = -1, n = 512);

  scalar l[];
  foreach()
    l[] = level;

  output_ppm (l, file = "level.mp4", box = {{0.,0.},{3.,1.}},
	      linear = false, min = 0, max = LEVEL, n = 512); 
}

event end (t = tend)
{
  FILE * fp = fopen ("rho.end", "w");
  output_field ({rho}, fp, n = 512, linear = true);
  fclose (fp);
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
