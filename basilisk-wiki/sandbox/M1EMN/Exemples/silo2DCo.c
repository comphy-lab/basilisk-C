/**
# Flow in a Sandglass / Hourglass

 
## Problem

 This is the problem of flow in a silo as the code [granular_sandglass](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c).
  Now the grains are cohesive, with cohesion as in [dry cohesive column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column_cohesif.c)
 
This is the same cohesive problem as [http://basilisk.fr/sandbox/M1EMN/Exemples/siloaxiCo.c]()
but in 2D, not in axi.
 
 
## Equations
 
 We propose an implementation of the Jop Pouliquen Forterre $\mu(I)$ rheology for the flow in the cohesive  hourglass.
 $$\tau = \tau_0 + \mu(I) p$$
 
We find that the flow through the orifice follows the Beverloo-Hagen discharge law with a modification due to cohesion $\tau_0$ and friction ($\mu_s$).
 
 $W$ is the width of the aperture of the silo (by symmetry $W/2$ is used). This $W$ is the fundamental length scale (in the code $W=1$).
 
## Code
Includes and definitions
*/

#include "granular.h"
// Domain extent
#define LDOMAIN 10.
// heap definition
double  H0,Largeur,D,W,tmax,Q,Wmin,DW;
//film size
double  Lz;
// Maximum refinement
#define LEVEL 6
///6 here to speed up but 7 better
char s[80];
FILE * fpf,*fwq;
/**
Boundary conditions for granular flow, pressure must be zero at the surface.
The pressure is zero in the hole, but the lithostatic gradient is given elswhere 
on the bottom.
No slip boundary conditions on the sides.
*/
p[top] =  dirichlet(0.);
u.n[top] = neumann(0);
u.t[bottom] =  fabs(x )<= W/2 ? neumann(0):  dirichlet(0);
u.n[bottom] =  fabs(x )<= W/2 ? neumann(0):  dirichlet(0);
p[bottom]   =  fabs(x )<= W/2 ? dirichlet(0): neumann(0);
p[right]   = neumann(0);
p[left]   =   neumann(0);
u.n[right] = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.t[left] = neumann(0);
/**
the three cases in the main */
int main() {
  L0 = LDOMAIN;
  X0 =0;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 0.01;
  TOLERANCE = 1e-5;
/**
  $W$ is the fundamental length scale (in the code $W=1$).
*/
  H0=L0;
  Lz = LDOMAIN;
  W = 1. ;
  Largeur= 10.;
    
  fwq = fopen ("outWQ", "w");
  fclose(fwq);
  const face vector g[] = {0.,-1.};
  a = g;
  mus=0.3;

/**
     when local define  full 1`
*/
#define full 0
    
#ifdef full
 for(mus = 0.1;  mus <= .5 ; mus+= .1 ){
    for(Co = 0.;  Co <= .7 ; Co+= .1 ){
#else
    for(Co = 0.;  Co <= .6 ; Co+= .25 ){  
#endif

   Q = 0;
   tmax = 4;
   fpf = fopen ("interface.txt", "w");
   run();
   fclose (fpf);
   fprintf (stdout,"\n");
   fwq = fopen ("outWQ", "a");
   fprintf(fwq," %lf %lf %lf %lf \n", W, Q, Co, mus);
   fclose (fwq);
  }
    fwq = fopen ("outWQ", "a");
    fprintf(fwq,"\n");
    fclose (fwq);
#ifdef full
 }
#endif
}
/**
initial fill
*/
event init (t = 0) {
    mask (x >   Largeur/2.  ? right : none);
  scalar phi[];
  foreach_vertex()
    phi[] = (4.5 - y);
  fractions (phi, f);
/**
lithostatic pressure, with a zero pressure near the hole
to help 
*/
   foreach()
     p[] = (fabs(x)<= W/2 && fabs(y)<= .1) ?  0 : max(4.5 - y,0) ;
}
// use the NZ technology: impose always p=0 at the top!
// but stop it
/*
 event go (t = 1) {
    scalar phi[];
    foreach_vertex()
    phi[] = (L0-.1 - y);
    fractions (phi, f);
}*/
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
 event interface ( t = 0; t += 1 ; t <= tmax) {
  output_facets (f, fpf); 
  char s[80];
  sprintf (s, "field-%g.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,eta}, fp, linear = true);
  fclose (fp);
}
/**
Rate of flowing materials across the hole (mind  the 2 by symmetry)
 $$-\frac{dV}{dt} \text{ with } V = 2 \int f d v $$
 and $-2 \int_0^{W/2} u_y ds$
*/
event debit (t += 0.05 ) {
  static double Vold,V=1,Qinst=0;
  Vold=V; 
  V=0;
    // mind  the 2
  foreach()
    V = V + 2*f[]*Delta*Delta;
  Qinst = -(V-Vold)/.05;
  Q=0;
    // mind  the 2
  foreach_boundary (bottom)
    Q -= 2*(Delta)*u.y[] ;
  if(t>=.1) fprintf (stdout,"%lf %lf %lf %lf  %lf \n",t,V/L0/H0,W,Qinst,Q);
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
    output_ppm (l, file = "pressure.mp4", min = 0, linear = true,
                n = 2048, box = {{0,-1},{Lz,Lz}});
}
#endif
event pictures (t=3) {
    scalar l[];
    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    output_ppm (l, file = "f.png", min=0, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{Lz,Lz}});
    output_ppm (p, file = "p.png", min=0, max = 1,  spread = 2, n = 256, linear = true,
                box = {{0,0},{Lz,Lz}});
}
/**
If `gfsview` is installed on your system you can use this to visualise the simulation as it runs.
 If you manage to install `Bview`it is better !
*/
#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D vue.gfs ", "w");
  output_gfs (fp, t = t);
}
#endif

/**
## Run

to run without `make`

~~~bash
qcc -g -O2 -DTRASH=1 -Wall  -o silo2DCo silo2DCo.c -lm
./silo2DCo > out
~~~

 
to run with `make`
 
~~~bash
make silo2DCo.tst
make silo2DCo/plots
make silo2DCo.c.html;
~~~

## Results


 Increasing the cohesion $Co=\frac{\tau_0}{\rho g W}$, note that here $L_0=W$,
 $Co=\frac{\tau_0}{\rho g W}$ decreases the flow rate.
 
In red for $\mu_s=0.3$ LEVEL 6, in green for $\mu_s=0.1, 0.2,...,0.5$ LEVEL 7
~~~gnuplot
reset
set xlabel "Co"
set ylabel "Q"
p [0:.8]'outWQ'u ($3):2 w l  t'mus=.3 LEVEL=6',\
'../REFCASES/outWQmussilo2DCo.dat'u ($3):2 w l t'mus=.1--0.5  LEVEL=7'
~~~ 

 at the limit of the flow,
 $$\tau = \tau_0 + \mu_s p$$
 
 as $p$ scales with $\rho g W$ hence
 $$\tau  \sim  \rho g W (\frac{\tau_0}{\rho g W }   + \mu_s  )$$
 so we guess that $(\frac{\tau_0}{\rho g W }   + \mu_s  )=Co + \mu_s$ is a pertinent quantity, indeed
 the previous flow rate $Q$ curves collapse when plotted as a fundtion of $Co + \mu_s$.
 
 In red for $\mu_s=0.3$, in green for $\mu_s=0.1, 0.2,...,0.5$
 ~~~gnuplot
 reset
 set xlabel "Co+mus"
 set ylabel "Q"
 p [0:.8]'outWQ'u ($3+$4):2 w l t'mus=.3',\
 '../REFCASES/outWQmussilo2DCo.dat'u ($3+$4):2 w l not
 ~~~

 
 
 
Film of pressure  
  
  [![](silo2DCo/p.png)](silo2DCo/pressure.mp4)
  velocity  (click on image for animation)

Film of velocity
  
  [![](silo2DCo/f.png)](silo2DCo/velo.mp4)
  velocity  (click on image for animation)

 

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
sp'field-3.txt' u 1:2:(($3>.9)&&($3<100) ? sqrt($7*$7+$6*$6) :0) not
~~~
 
## Links
 
 * [dry column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c)
 
 * [dry cohesive column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column_cohesif.c)
 
 * [silos](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c)
 
 * axi cohesive problem as [http://basilisk.fr/sandbox/M1EMN/Exemples/siloaxiCo.c]()

 
## Bibliography
 
* L. Staron, P.-Y. Lagrée, & S. Popinet (2014)
"Continuum simulation of the discharge of the granular silo, A validation test for the μ(I) visco-plastic flow law" 
Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6

* L. Staron, P.-Y. Lagrée & S. Popinet (2012)
"The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra" 
Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390
 
* Y. Zhou,   <a href="http://www.lmm.jussieu.fr/%7Elagree/TEXTES/PDF/PhysRevFluids.4.124305.pdf">P.-Y. Lagr&eacute;e</a>,  S. Popinet, P. Ruyer, and P. Aussillous<br>"Gas-assisted discharge flow of granular media from silos "<br>
Phys. Rev. Fluids 4, 124305 – Published 18 December 2019
DOI: 10.1103/PhysRevFluids.4.124305
 
* Abramian
 
  
* Luke Fullard, Daniel J. Holland, Petrik Galvosas, Clive Davies,    
 <a href="http://www.lmm.jussieu.fr/%7Elagree/TEXTES/PDF/PhysRevFluids.4.074302.pdf">P.-Y. Lagr&eacute;e</a>, and St&eacute;phane Popinet  <br>
"Quantifying silo flow using MRI velocimetry for testing granular flow models"<br>
Phys. Rev. Fluids 4, 074302, DOI: 10.1103/PhysRevFluids.4.074302<br>
  

 * Y. Zhou  <a href="http://www.lmm.jussieu.fr/%7Elagree/TEXTES/PDF/ZhouJFM2017.pdf">P.-Y. Lagr&eacute;e</a>  S. Popinet, P. Ruyer and P. Aussillous  <br>
"Experiments on, and discrete and continuum simulations of, the discharge of granular media from silos with a lateral orifice"<br>
Journal of Fluid Mechanics
Volume 829 25 October 2017 , pp. 459-485<br>
http://doi.org/10.1017/jfm.2017.543 <br>

 * L. Staron,  <a href="http://www.lmm.jussieu.fr/%7Elagree/TEXTES/PDF/epje130141.pdf">P.-Y. Lagr&eacute;e</a>,   &amp; S. Popinet  (2014)<br>
"Continuum simulation of the discharge of the granular silo,
 A validation test for the &#956;(I) visco-plastic flow law"
<br>
 Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6<br>
 
 
 * L. Staron,   <a href="http://www.lmm.jussieu.fr/%7Elagree/TEXTES/PDF/PhysFluids_24_103301.pdf">P.-Y. Lagr&eacute;e</a>  &amp; S. Popinet (2012)<br>
"The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra"
<br>Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390 <br>
 
*/
