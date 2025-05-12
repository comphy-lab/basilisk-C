/**
# Test tracers diffusion with several way to compute face vector. 

This checks the result of a way to compute face vector wich value depend on the $f$
on the diffusion of tracers in simple configuration. The $Tf[]$ value should stay at zero in every cells that
contains $finside$ (ie. centered circle) phase (no diffusion) and diffusion elsewhere.
Adaptative mesh and three-phase case can be swithched. */

#define ADAPT 1
#define LINEAR 1
#define LIMITER 1
#define VOL_CORRECTION 1
#define REFINE 0
#define UNREFINE 0

#include "fractions.h"
#include "run.h"
#include "diffusion.h"

#define VOFTHR 1e-6
#define MAXTIME (5*DT)

#define dx (1/128.)
#define x0 (50*dx)
#define D1 1.
#define D2 0.
#ifndef D
#define D(f)  (clamp(f,0.,1.)*(D1 - D2) + D2)
#endif

FILE * fp,* fp1,* fp2, * fp3, * fp4, * fp5, * fp6, * fp7, * fp8;
int minlevel=6;
int maxlevel=7;

double dt, initial_Tf;
scalar f[], Tfbegin[], Tf[];
vertex scalar phi[];
face vector s[], Diff[];

int main() {

  size (1);
  init_grid (1 << minlevel);
  DT = 0.005; 
  TOLERANCE = 1e-6;

  /**
    Three fractions can be defined using vertex scalar in case we use adaptivity:
   * $f$ corresponding to the outsied of circle centered in the middle of the domain;
   */
  fp=fopen("tracers_cons.dat", "w");
  fp1=fopen("interfacellsf.dat", "w");
  fp2=fopen("nointerfacellsf.dat", "w");
  fp3=fopen("tracers.dat", "w");
  fp4=fopen("averagediff.dat", "w");
  fp5=fopen("noaveragediff.dat", "w");
  fp6=fopen("facetsbegin.dat", "w");
  fp7=fopen("cellsbegin.dat", "w");
  fp8=fopen("values_begin.dat", "w");

  /**
    The cell grid and values of the fractions are saved for regular or
    adaptative grid. */

#if TREE
  f.refine = f.prolongation = fraction_refine;
#if LIMITER
  Tf.gradient = minmod2;
  Tfbegin.gradient = minmod2;
#endif
#if LINEAR
  Tf.refine = refine_linear;
  Tf.restriction = restriction_volume_average;
  Tfbegin.refine = refine_linear;
  Tfbegin.restriction = restriction_volume_average;
#endif
#endif
  run();
}

event init (t=0) {
  foreach_vertex() {
    phi[] = x-x0;
  }
  boundary ({phi});
  fractions (phi, f, s);
#if ADAPT
  adapt_wavelet ({f}, (double[]){1e-2}, maxlevel);
#endif
#if UNREFINE
  unrefine ((x > x0-x0/10. && x < x0+x0/10.) && level < maxlevel);
#endif
#if REFINE
  refine ((x > x0-x0/10. && x < x0+x0/10.) && level < maxlevel);
#endif
  foreach()
    Tf[]= f[];
  boundary ((scalar *){Tf});
#if ADAPT
  adapt_wavelet ({Tf}, (double[]){1e-2}, maxlevel);
#endif
  foreach() {
    Tfbegin[]=Tf[];
    fprintf (fp8, "%g %g %2.1g %2.1g %2.1g\n", x, y, Tfbegin[], Tf[], f[]);
  }
  boundary ((scalar *){Tfbegin});
  output_facets (f, fp6, s);
  output_cells (fp7);
  initial_Tf = statsf(Tf).sum;
}

event tracer_diffusion (i++) {
  /**
    Here we compute  the Diff value for two phase case with average and
    conditions on $f$ the face vector value. */

  face vector Diff[];
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    Diff.x[]=D(ff)*fm.x[];
    if (t == MAXTIME)
      fprintf (fp4, "%g %g %2.1g\n", x, y, Diff.x[]);
  }
  foreach_face() {
    Diff.x[]=0.;
    if (f[] >1.-VOFTHR && f[-1] >1.-VOFTHR)
      Diff.x[]=D1*fm.x[];
    else 
      Diff.x[]= D2*fm.x[];
    if (t == MAXTIME)
      fprintf (fp5, "%g %g %2.1g\n", x, y, Diff.x[]);
  }
  boundary ((scalar *){Diff});
  dt = dtnext (DT);
#if VOL_CORRECTION
  scalar volume_correction[];
#if TREE
  volume_correction.prolongation = volume_correction.refine = fraction_refine;
#endif
  foreach() {
    f[] = clamp(f[], 0., 1.);
    volume_correction[] = cm[]*max(f[], VOFTHR);
  }
  boundary ({f, volume_correction});
  foreach()
    Tf[] = (f[] > VOFTHR ? Tf[]/f[] : 0.);
  boundary({Tf});
  diffusion (Tf, dt, Diff, theta = volume_correction); 
  foreach()
    Tf[] *= f[];
  boundary({Tf});
#else
  diffusion (Tf, dt, Diff); 
#endif
}

event tracer_statistic (i++) {
  double current_Tf = statsf(Tf).sum;
  double variation_Tf = (initial_Tf-current_Tf)/initial_Tf;
  double Tff = 0., Tfg = 0., volg = 0., voll = 0.;
  foreach(reduction(+:Tff) reduction(+:Tfg) reduction(+:volg) reduction(+:voll)) {
    if (dv() > 0.) {
      if (f[] == 0.) {
	Tfg += Tf[]*dv();
	volg += dv();
      }
      else if (f[] >=1.-VOFTHR) {
	Tff += Tf[]*dv();
	voll += dv();
      }
    }
  }
  fprintf (fp, "%g %g %g %g %g %g %g %g %g\n", 
      t, initial_Tf, current_Tf, variation_Tf, volg, voll, Tff/voll, Tfg/volg, interface_area(f));
  fflush (fp);
}

event output (t = MAXTIME) {
  scalar Tfdiff[];
  /**
    Here we want to detect the interface cells and their $f$ values.*/
  foreach() {
    if (f[] <=1.-VOFTHR && f[] > VOFTHR)
      fprintf (fp1, "%g %g %2.1g\n", x, y, f[]);
    else
      fprintf (fp2, "%g %g %2.1g\n", x, y, f[]);
  }
  foreach() {
    Tfdiff[] = Tf[]-Tfbegin[];
    fprintf (fp3, "%g %g %2.1g %2.1g %2.1g\n", x, y, Tfbegin[], Tf[], Tfdiff[]);
  }
  boundary ((scalar *){Tfdiff});
  output_facets (f, stderr, s);
  output_cells (stdout);
  fclose (fp1);
  fclose (fp2);
  fclose (fp3);
  dump();
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,Tf}, (double[]){1e-2, 1e-2}, maxlevel);
}
#endif

/**
## Ouputs

~~~gnuplot Conservation of $Tf$ 
set xlabel 'Diffusion time step [-]'
set ylabel '(Tf0-Tf)/Tf0 [-]'
DT=0.005
plot 'tracers_cons.dat' u ($1/DT):4 w lp t ''
~~~

~~~gnuplot Conservation of $Tf$ 
set xlabel 'Diffusion time step [-]'
set ylabel 'Tf in f[] [-]'
DT=0.005
plot 'tracers_cons.dat' u ($1/DT):7 w lp t ''
~~~

~~~gnuplot Conservation of $Tf$ 
set xlabel 'Diffusion time step [-]'
set ylabel 'Tf in 1-f[] [-]'
DT=0.005
plot 'tracers_cons.dat' u ($1/DT):8 w lp t ''
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with $f[]$ value at start of simulation.
reset
set term @SVG size 900,900
unset border
unset key
unset tics
plot [0.25:0.5][0.25:0.5][:] 'cellsbegin.dat' w l, 'facetsbegin.dat' w l, 'values_begin.dat' u 1:2:5 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with $f[]$ value outside interface at the end of simulation.
reset
set term @SVG size 900,900
unset border
unset key
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'nointerfacellsf.dat' u 1:2:3 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with $Diff.x[]$ value computed with average at the end of simulation.
reset
set term @SVG size 900,900
unset key
unset border
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'averagediff.dat' u 1:2:3 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with $Diff.x[]$ value imposed with conditions at the end of simulation.
reset
set term @SVG size 900,900
unset key
unset border
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'noaveragediff.dat' u 1:2:3 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with values of $Tf$ at the begining of simulation
reset
set term @SVG size 900,900
unset key
unset border
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'tracers.dat' u 1:2:3 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with values of $Tf$ at the end of simulation
reset
set term @SVG size 900,900
unset key
unset border
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'tracers.dat' u 1:2:4 with labels
~~~

~~~gnuplot Reconstruction of $\phi(x,y)$ with values of difference between begining and end of simulation of $Tf$ 
reset
set term @SVG size 900,900
unset key
unset border
unset tics
plot [0.25:0.5][0.25:0.5][:] 'out' w l, 'log' w l, 'tracers.dat' u 1:2:5 with labels
~~~
*/
