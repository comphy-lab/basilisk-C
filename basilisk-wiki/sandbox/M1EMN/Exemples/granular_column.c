/**
# Collapse of granular columns
 
## Problem
Geophysics, moutains...
 

## Equations
We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology. This is the counterpart in Basilisk  of the [gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/column.html#htoc9) case.
 
It consists to solve for three aspect ratios the collapse of three columns. Each collapse is compared to the contact dynamics solution.

## Code
Includes and definitions, note `granular.h`, and here `multigrid`
*/
#include "grid/multigrid.h"
#include "granular.h"
// Domain extent
#define LDOMAIN 5.
// heap definition
double  H0,R0,D;
//film size
double  Lz;
// Maximum refinement
#define LEVEL 7
char s[80];
FILE * fpf;

/**
Boundary conditions for granular flow, pressure must be zero at the surface */

p[top] = dirichlet(-RHOF*LDOMAIN);
pf[top] = dirichlet(-RHOF*LDOMAIN);
u.n[top] = neumann(0);
u.t[bottom] = dirichlet(0);

/**
the three cases in the main */

int main() {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 1e-2;
  TOLERANCE = 1e-3;
  
// Initial conditions a=.5
  H0=1;
  R0=2.0001;

  
  const face vector g[] = {0.,-1.};
  a = g;
    
  fpf = fopen ("out05", "w");
  Lz = LDOMAIN;
  run();
  fclose (fpf);
  system("gnuplot comp05.gnu");

// Initial conditions a=1.42
  H0=1;
  R0=1./1.42; 
  D=1./65;
  fpf = fopen ("out142", "w");
  Lz = 3;
  run();
  fclose (fpf);
  system("gnuplot comp142.gnu");
  
// Initial conditions a=6.26
  H0=1;
  R0=1./6.26; 
  D=1./65;
  fpf = fopen ("out626", "w");
  Lz = 2.5;
  run();
  fclose (fpf); 
  system("gnuplot comp626.gnu");  
  
 
}
/**
initial heap, a rectangle
*/

event init (t = 0) {
 

  scalar phi[];
  foreach_vertex()
    phi[] = min(H0 - y, R0 - x);
  fractions (phi, f);
/*
initialisation of hydrostatic pressure for granular phase  
*/
  foreach()
    p[] = f[]*(H0-y);
}
/**
convergence outputs
*/
void mg_print (mgstats mg)
{
   if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0);
}
/**
convergence outputs
*/
event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}
/**
convergence outputs
*/
event interface (t = {0, 1. , 2., 3., 4.}) {
  output_facets (f, fpf); 
  char s[80];
  sprintf (s, "field-%g-%5.3g.txt", t, R0);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf}, fp, linear = true);
  fclose (fp);
}
/**
film output
*/

event movie (t += 0.05) {
  //static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, file = "level.mp4", min = 0, max = LEVEL,
	      n = 2048, box = {{0,-1},{Lz,Lz}});

  foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
//  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, file = "velo.mp4", min = 0, max = 1.5, linear = true,
	      n = 1024, box = {{0,-1},{5,1.5}});

 // static FILE * fp3 = popen ("ppm2mpeg > f.mpg", "w");
  foreach()
    l[] = f[]*p[];
  output_ppm (l, file = "f.mp4", min = 0, linear = true,
	      n = 2048, box = {{0,-1},{Lz,Lz}});
}
 event pictures (t=3) {
  output_ppm (f, file = "f.png", min=0, max = 2,  spread = 2, n = 256, linear = true,
  box = {{0,0},{4,2}});
}


/**

## Run

 to run without `make`


~~~bash
qcc -g -O2 -DTRASH=1 -Wall  -o granular_column granular_column.c -lm
./granular_column > out
~~~

to run with `make`
 
~~~bash
make granular_column.tst
make granular_column/plots
make granular_column.c.html;
~~~

 
# Results

Exemples of comparisons Discrete Contacts Dynamics vs Continuum Model for 3 aspect ratio
 
~~~gnuplot collapse Contacts Dynamics vs Continuum Model a =0.5
d=0.005
h0=0.149
set xlabel "x"
set ylabel "h(x,t)"
p[0:4][0:1.5]'../../granular_column/ShapeTime.A-01.dat' u (($1*d)/h0):(($2*d)/h0) t'DCM'w l,\
'../granular_column/out05' w l t 'h'
~~~

 
~~~gnuplot collapse Contacts Dynamics vs Continuum Model a =1.42
set xlabel "x"
set ylabel "h(x,t)"
d=0.005
h0=0.324
p[0:4][0:1.5]'../Img/ShapeTime.A-04.dat' u (($1*d)/h0):(($2*d)/h0) t'DCM'w l,\
'../granular_column/out142' w l t 'h'
~~~


~~~gnuplot collapse Contacts Dynamics vs Continuum Model a =6.26
set xlabel "x"
set ylabel "h(x,t)"
d=0.005
h0=0.67
p[0:4][0:1.5]'../Img/ShapeTime.A-03.dat' u (($1*d)/h0):(($2*d)/h0) t'DCM'w l,\
'../granular_column/out626' w l t 'h'
set term png
set output 'a.png'
 replot
~~~
 
  
[![](./granular_column/f.png)](./granular_column/velo.mp4)
  velocity  (click on image for animation)



 
# Related examples
 
* [The basic Bagnold flow](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
* [dry column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c)
 
* [silos](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c)
 
 
 

# Bibliography
* Lagrée, Staron, Popinet 
["The granular column collapse as a
continuum: validity of a two–dimensional
Navier–Stokes model with a μ(I)-rheology"](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf) J. Fluid Mech. 2011 


Montpellier 07/14
*/
