/**
# Flow in a Sandglass / Hourglass

 
## Problem

A silo is a structure for storing bulk materials like grains,  coal, cement, gravels,  woodchips, food products...
An hourglass or sandglass (marine term), or sand timer, or sand clock, or even egg timer, is a device used to measure the passage of time.
 
 In both cases, a granular material (sand, grains...) is going out from a reservoir   through a small aperture.
 
 The fact that flow rate is constant has been used for time measurement for a long time.
Hagen (the same Hagen from Hagen-Poiseuille) documented this and found the scaling law with aperture size,this was rediscoverd by Beverloo in the 60s.
 
 
~~~gnuplot silo
 set arrow 1 from 4,0 to 4,1.1 front
 set arrow 2 from 4,1 to 4,0 front
 set arrow 3 from 0.1,.2 to 4.9,0.2 front
 set arrow 4 from 4.9,.2 to 0.1,0.2 front
 set arrow 5 from 2,.1 to 3.,0.1 front
 set arrow 6 from 3,.1 to 2.,0.1 front
 set label 1 "L0" at 3,.25 front
 set label 2 "2 W" at 2.5,.05 front
 set label 3 "height \nof grains" at 4.1,.5 front
 set label 4 "grains" at 1.2,.9 front
 set label 5 "air" at 1.2,1.2 front
 set xlabel "x"
 set ylabel "y"
 p [0:5]0 not,1+.05*(x-2.5)**2 w filledcurves x1 linec 3 t'free surface'
 unset arrow 1; unset arrow 2; unset arrow 3
 unset arrow 4; unset arrow 5; unset arrow 6
 unset label 1; unset label 2; unset label 3
 unset label 4; unset label 5
 
 ~~~
 
## Equations
 
 We propose an implementation of the Jop Pouliquen Forterre $\mu(I)$ rheology for the flow in a hourglass.
 
 
We find that the flow through the orifice follows the Beverloo-Hagen discharge law.
 
## Code
Includes and definitions
*/
#include "grid/multigrid.h"
#include "granular.h"
// Domain extent
#define LDOMAIN 5.
// heap definition 2W is the size of the gate of the silo
double  H0,R0,D,W,tmax,Q,Wmin,DW;
//film size
double  Lz;
// Maximum refinement
#define LEVEL 6
///6 here to speed up for web site but 7 better
char s[80];
FILE * fpf,*fwq;

/**
Boundary conditions for granular flow, pressure must be zero at the surface.
The pressure is zero in the hole (of width $2W$),
 but the lithostatic gradient is given elswhere
on the bottom.
No slip boundary conditions on the sides.
*/
p[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[bottom] =  fabs(x-LDOMAIN/2)<= W ? neumann(0):  dirichlet(0);
u.n[bottom] =  fabs(x-LDOMAIN/2)<= W ? neumann(0):  dirichlet(0);
p[bottom]   =  fabs(x-LDOMAIN/2)<= W ? dirichlet(0): neumann(0); 
u.n[right] = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.t[left] = dirichlet(0);
f[left]= neumann(0);
/**
the three cases in the main */
int main() {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 0.01;
  TOLERANCE = 1e-3; 
// Initial conditions a=.5
  H0=4.5;
  R0=20.000;
  // size of the hole
  Wmin = 0.157;
  DW = LDOMAIN/N;
  
  const face vector g[] = {0.,-1.};
  a = g;
  
  fwq = fopen ("outWQ", "w");
  fclose(fwq);
  Lz = LDOMAIN;
  tmax = 10.;
 
#if 1
  for(W = Wmin;  W <= .6 ; W+= 2*DW ){
   Q = 0; 
   tmax = 4.; 
   fpf = fopen ("interface.txt", "w");
   run();
   fclose (fpf);
   fprintf (stdout,"\n");
   fwq = fopen ("outWQ", "a");
   fprintf(fwq," %lf %lf \n", W, Q);
   fclose (fwq);
   }
#endif
  Q = 0; 
  tmax = 20;
  W = .5;
  fpf = fopen ("interface.txt", "w");
  run();
  fclose (fpf);
}
/**
initial heap, a rectangle
*/
event init (t = 0) {
  scalar phi[];
  foreach_vertex()
    phi[] = min(H0 - y, R0 - x);
  fractions (phi, f);
/**
lithostatic pressure, with a zero pressure near the hole
to help 
*/
   foreach()
     p[] = (fabs(x-LDOMAIN/2)<= W && fabs(y)<= .1) ?  0 : max(H0 - y,0) ;    
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
convergence stats
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
save some interfaces
*/
//event interface (t = {0, 1. , 2., 3., 4.}) {
 event interface ( t = 0; t += 1 ; t <= tmax) {
  output_facets (f, fpf); 
  char s[80];
  sprintf (s, "field-%g.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf}, fp, linear = true);
  fclose (fp);
}
/**
Rate of flowing materials across the hole
 $$-\frac{dV}{dt} \text{ with } V = \int f dv $$
*/
event debit (t += 0.05 ) {
  static double Vold,V=1,Qinst=0;
  Vold=V; 
  V=0;
  foreach()
    V = V + f[]* Delta * Delta;
  Qinst = -(V-Vold)/.05;  
  if(Qinst > Q) Q = Qinst;
  if(t>=.1) fprintf (stdout,"%lf %lf %lf %lf \n",t,V/L0/H0,W,Q); 
  fflush (stdout);
}  
/**
film output
*/
#if 0
event movie (t += 0.05) {
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = 0, max = LEVEL, 
	      n = 2048, box = {{0,0},{Lz,Lz}});

  foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, fp2, min = 0, max = 2., linear = true, 
	      n = 2048, box = {{0,0},{Lz,Lz}});

  static FILE * fp3 = popen ("ppm2mpeg > pressure.mpg", "w");
  foreach()
    l[] = f[]*p[];
  output_ppm (l, fp3, min = 0, linear = true,
	      n = 2048, box = {{0,0},{Lz,Lz}});
}
#else


event pictures (t=3) {
  scalar l[];
  foreach()
     l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  output_ppm (l, file = "f.png", min=0, max = 2,  spread = 2, n = 256, linear = true,
  box = {{0,0},{Lz,Lz}});
}

event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = level;
    output_ppm (l, file = "level.mp4", min = 0, max = LEVEL,
                n = 2048, box = {{0,-1},{Lz,Lz}});
    
    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    output_ppm (l, file = "velo.mp4", min = 0, max = 2., linear = true,
                n = 2048, box = {{0,-1},{Lz,Lz}});
    foreach()
    l[] = f[]*p[];
    output_ppm (l, file = "pressure.mp4", min = 0,max = 2., linear = true,
                n = 2048, box = {{0,-1},{Lz,Lz}});
}/**

event pictures (t=3) {
    output_ppm (f, file = "f.png", min=0, max = 2,  spread = 2, n = 512, linear = true,
                box = {{0,0},{2,2}});
    output_ppm (p, file = "p.png", min=0, max = 1,  spread = 2, n = 512, linear = true,
                box = {{0,0},{2,2}});
}
  */
#endif
/**
If gfsview is installed on your system you can use this to visualise the simulation as it runs.
 If you manage to install `Bview`it is better !
*/
#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D -s  ", "w");
  output_gfs (fp, t = t);
}
#endif
/**
if "grid/multigrid.h" is not included, then this is the quadtree with possibility of adaptation
*/
#if QUADTREE
event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, (double[]){5e-3,0.001,0.001}, LEVEL, LEVEL-2,
     list = {p,u,pf,uf,g,f});
}
#endif
/**
## Run

to run without `make`

~~~bash
qcc -g -O2 -DTRASH=1 -Wall  -o granular_sandglass granular_sandglass.c -lm
./granular_sandglass > out
~~~

 
to run with `make`
 
~~~bash
make granular_sandglass.tst
make granular_sandglass/plots
make granular_sandglass.c.html;
~~~

## Results

Volume as function of time for different orifices


~~~gnuplot Beverloo 
 set xlabel "t"
 set ylabel "V(t)"
 set key right
 p[0:]'out' t 'V' w l 
~~~


Fitting the Beverloo law, $2W$ is the size of the hole of the silo
 
~~~gnuplot Beverloo 
 set xlabel "W"
 set ylabel "Debit"
 f(x) = a*(x +b)**1.5 
 fit f(x) 'outWQ' via a,b
 set key left
 p[0:]'outWQ' t'debit' w lp,4.254*x**(1.5),f(x) t'fit' ,'../REFCASES/outWQ_granular_sandglass.dat' w p t'for control'
~~~

 note that the `'../REFCASES/outWQ_granular_sandglass.dat'`  was computed with
 `#define RHOF 1e-4`
 `#define mug  1e-5`
  and with  `LEVEL 7`
 
 

Plot of the interface every at time
~~~gnuplot interface of grain during the flow t=0, 1, 2, 3 ... 
 reset  
 p[0:][0:]'interface.txt' not w l
~~~

Plot of pressure at time 4

~~~gnuplot pressure
reset
set pm3d; set palette rgbformulae 22,13,-31;unset surface;
set ticslevel 0;
unset border;
unset xtics;
unset ytics;
unset ztics;
unset colorbox;
#set xrange[-3:3];set yrange[-3:3];
set view 0,0
sp'field-4.txt' u 1:2:($3>.9 ? $4 :0) not
~~~


Plot of velocity at time 4

~~~gnuplot velocity
reset
set pm3d; set palette rgbformulae 22,13,-31;unset surface;
set ticslevel 0;
unset border;
unset xtics;
unset ytics;
unset ztics;
unset colorbox;
#set xrange[-3:3];set yrange[-3:3];
set view 0,0
sp'field-1.txt' u 1:2:($3>.9 ? sqrt($7*$7+$6*$6) :0) not
~~~ 
 
Film of pressure  
  
  [![](granular_sandglass/p.png)](granular_sandglass/pressure.mp4)
  velocity  (click on image for animation)

Film of velocity
  
  [![](granular_sandglass/f.png)](granular_sandglass/velo.mp4)
  velocity  (click on image for animation)

if quadtree

  [![](granular_sandglass/f.png)](granular_sandglass/level.mp4)
  level of mesh (click on image for animation)


~~~gnuplot velocity
reset
set border 4095 lt -1 lw 1.000
set format cb "%.01t*10^{%T}"
set samples 31, 31
set isosamples 31, 31
unset surface
set ticslevel 0 
set xlabel "X" 
set xrange [  :  ] noreverse nowriteback
set ylabel "Y"
set yrange [  :  ] noreverse nowriteback
set cblabel "the colour gradient"
set pm3d  
set palette positive nops_allcF maxcolors 0 gamma 1.5 gray
set view 0,0
sp'field-3.txt' u 1:2:($3>.9 ? sqrt($7*$7+$6*$6) :0) not
~~~
 
## Links
 
 * [granular include](http://basilisk.fr/sandbox/M1EMN/Exemples/granular.h)
 
 * [The basic Bagnold flow](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
 * [dry column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c)
 
 * [dry cohesive column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column_cohesif.c)
 
 * [silos](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c)
 
 * [silos](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c)
 
 * [Hele Shaw silos with friction](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass_muw.c)
 
## Bibliography
 
* L. Staron, P.-Y. Lagrée & S. Popinet (2012)
"The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra" 
Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390

* L. Staron, P.-Y. Lagrée, & S. Popinet (2014)
"Continuum simulation of the discharge of the granular silo, A validation test for the μ(I) visco-plastic flow law" 
Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6
           
* Y. Zhou P.-Y. Lagrée S. Popinet, P. Ruyer and P. Aussillous 
"Experiments on, and discrete and continuum simulations of, the discharge of granular media from silos with a lateral orifice"
Journal of Fluid Mechanics Volume 829 25 October 2017 , pp. 459-485
http://doi.org/10.1017/jfm.2017.543 

* Y. Zhou, P.-Y. Lagrée, S. Popinet, P. Ruyer, and P. Aussillous
"Gas-assisted discharge flow of granular media from silos "
Phys. Rev. Fluids 4, 124305 – Published 18 December 2019 DOI: 10.1103/PhysRevFluids.4.124305 

* Luke Fullard, Daniel J. Holland, Petrik Galvosas, Clive Davies, P.-Y. Lagrée, and Stéphane Popinet (2019)
"Quantifying silo flow using MRI velocimetry for testing granular flow models"
Phys. Rev. Fluids 4, 074302, DOI: 10.1103/PhysRevFluids.4.074302

*/
