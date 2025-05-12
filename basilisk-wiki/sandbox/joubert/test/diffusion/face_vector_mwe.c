/**
# Test several way to compute a face vector. 

This checks the way that a face vector wich value depend on the $f$
is computed in simple configuration. The $Diff.x[]$ value should be zero in every cells that
contains $finside$ phase (ie. centered circle). Adaptative mesh and three-phase
case can be swithched. */

#define ADAPT 1
#define REFINE 0
#define UNREFINE 0

#include "fractions.h"

#define VOFTHR 1e-6

#define dx (1/128.)
#define x0 (50*dx)
#define D1 1.
#define D2 0.
#ifndef D
# define D(f)  (clamp(f,0.,1.)*(D1 - D2) + D2)
#endif

FILE * fp,* fp1,* fp2, * fp3;
int minlevel=6;
int maxlevel=7;

int main() {

  size (1);
  init_grid (1 << minlevel);

  /**
    Three fractions can be defined using vertex scalar in case we use adaptivity:

   * $finside$ corresponding to a circle center in the middle of the domain;
   * $fbis$ corresponding to a half circle in the bottom middle of the domain;
   * $f$ corresponding to the interception of $finside$ and $fbis$;
   */
  scalar f[];
  vertex scalar phi[];
  face vector s[];

  foreach_vertex() {
    phi[] = x-x0;
  }
  boundary ({phi});
  fractions (phi, f, s);

  fp=fopen("averagediff.dat", "w");
  fp1=fopen("noaveragediff.dat", "w");
  fp2=fopen("interfacellsf.dat", "w");
  fp3=fopen("nointerfacellsf.dat", "w");

  /**
    The cell grid and values of the fractions are saved for regular or
    adaptative grid. */

#if TREE
  f.refine = f.prolongation = fraction_refine;
#if REFINE
  refine ((x > x0-x0/10. && x < x0+x0/10.) && level < maxlevel);
#endif
#if UNREFINE
  unrefine ((x > x0-x0/10. && x < x0+x0/10.) && level < maxlevel);
#endif
#if ADAPT
  adapt_wavelet ({f}, (double[]){1e-2}, maxlevel);
#endif
#endif
  output_cells (stdout);
  output_facets (f, stderr, s);

  /**
    Here we compute  the Diff value for two phase case with average and
    conditions on $f$ the face vector value. */

  face vector Diff[];
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    Diff.x[]=D(ff)*fm.x[];
    fprintf (fp, "%g %g %2.1g\n", x, y, Diff.x[]);
  }
  foreach_face() {
    Diff.x[]=0.;
    if (f[] > 1.- VOFTHR && f[-1] >1.-VOFTHR)
      Diff.x[]=D1*fm.x[];
    else 
      Diff.x[]=D2*fm.x[];
    fprintf (fp1, "%g %g %2.1g\n", x, y, Diff.x[]);
  } 
  /**
    Here we want to detect the interface cells and their $f$ values.*/

  foreach() {
    if (f[] <=1.-VOFTHR && f[] > VOFTHR)
      fprintf (fp2, "%g %g %2.1g\n", x, y, f[]);
    else
      fprintf (fp3, "%g %g %2.1g\n", x, y, f[]);
  }
  fclose (fp);
  fclose (fp1);
  fclose (fp2);
  fclose (fp3);
}

/**
## Ouputs

~~~gnuplot Reconstruction of $\phi(x,y)$ with $f[]$ value outside interface
set term @SVG size 900,900
unset border
unset key
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'nointerfacellsf.dat' u 1:2:3 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with $Diff.x[]$ value computed with average
set term @SVG size 900,900
unset key
unset border
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'averagediff.dat' u 1:2:3 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with $Diff.x[]$ value imposed with conditions
set term @SVG size 900,900
unset key
unset border
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'noaveragediff.dat' u 1:2:3 with labels
~~~
*/
