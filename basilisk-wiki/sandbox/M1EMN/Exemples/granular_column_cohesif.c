/**
# Collapse of cohesive granular columns 

## Problem

We look at a 2D column of cohesive grains collapsing on a flat surface.


![Snapshot of the heap. t=0](granular_column_cohesif/f0.png)
![Snapshot of the heap. t=1](granular_column_cohesif/f1.png)
![Snapshot of the heap. t=2](granular_column_cohesif/f2.png)
 
## Equations
 
 We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology plus cohesion. It is just an extra term:
$$ \tau = (Co + \mu(I) P ) \frac{D}{|D|}$$
 the equivalent viscosity
 $$\eta =\frac{Co + \mu(I) P}{\sqrt{2} D_2}$$
is included in Navier Stokes in [http://basilisk.fr/sandbox/M1EMN/Exemples/granular.h]()

## Code
 
Includes and definitions, note `granular.h`
*/
#include "granular.h"
#include "heights.h"
// Domain extent
#define LDOMAIN 5.
// heap definition
double  H0,R0,xfront;
//film size
double  Lz;
// Maximum refinement
#define LEVEL 8
// 8 is OK
#define LEVELmin 2
char s[80];
FILE * fpf;

/**
 Boundary conditions for granular flow, pressure must be zero at the surface */

p[top]      = dirichlet(-RHOF*LDOMAIN);
pf[top]     = dirichlet(-RHOF*LDOMAIN);
u.n[top]    = neumann(0);
u.t[bottom] = dirichlet(0);


/**
the three cases in the main */

int main() {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 25.e-3;
  TOLERANCE = 1e-3;

/** Initial conditions aspect ratio a=.5 or .25 by symmetry
 Grain size $D$, value of the $\mu_s$ $\Delta \mu$ and $I_0$ of $\mu(I)$ law,
 and gravity $\overrightarrow{g}=-g \overrightarrow{e}_y$ */

  H0=1;
  R0=2.0001;
  D=1./32;
  //mus=.32;
  mus=.4;
  dmu=.28;
  I0=0.4;

  const face vector gravity[] = {0.,-1.};
  a = gravity;
 
  FILE * f = fopen ("tauYxmax.txt", "w");

  fpf = fopen ("out000", "w");
  Lz = LDOMAIN;
  Co=0.00;
  run();
  fprintf (f, "%g %g  \n", Co, xfront );
  fclose (fpf); 
  system("cp velo.mpg velo0.00.mpg");

#define loc
    
#ifdef loc
  fpf = fopen ("out005", "w");
  Co=0.05;
  run();
  fprintf (f, "%g %g  \n", Co, xfront );
  fclose (fpf); 
  system("cp velo.mpg velo0.05.mpg");
    

  fpf = fopen ("out010", "w");
  Co=0.10;
  run();
  fprintf (f, "%g %g  \n", Co, xfront );
  fclose (fpf); 
  system("cp velo.mpg velo0.10.mpg"); 

  Co=0.15;
  fpf = fopen ("out015", "w");
  run();
  fprintf (f, "%g %g  \n", Co, xfront );
  system("cp velo.mpg velo0.15.mpg"); 
  Co=0.10; 
#endif

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
## outputs
 
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
surface output
*/
event interface (t = {0, 1. , 2., 3., 4., 5.}) {
  output_facets (f, fpf); 
  char s[80];
  sprintf (s, "field-%g-%g.txt", t, Co);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf}, fp, linear = true);
  fclose (fp);
}

/**
#if QUADTREE
  if no #include "grid/multigrid.h" then adapt */

event adapt(i++){

 scalar K1[],K2[];
 foreach()
   K1[]=(f[0,1]+f[0,-1]+f[1,0]+f[-1,0])/4*noise();
 boundary({K1});
 
 for(int k=1;k<8;k++)
 {
 foreach()
   K2[]=(K1[0,1]+K1[0,-1]+K1[1,0]+K1[-1,0])/4;
 boundary({K2});
 foreach()
   K1[]=K2[]*noise();
}
 adapt_wavelet({K1,f},(double[]){0.001,0.01}, maxlevel = LEVEL, minlevel = LEVELmin);

}


/**
film output
*/
#if 0
event movie (t += 0.05) {
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level*(1-f[]);
    boundary ({l});
  output_ppm (l, fp1, min = 0, max = LEVEL,
	      n = 2048, box = {{0,-1},{Lz,Lz}});

  foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, fp2, min = 0, max = 1.5, linear = true, 
              n = 1024, box = {{0,-1},{4,1.5}});

  static FILE * fp3 = popen ("ppm2mpeg > f.mpg", "w");
  foreach()
    l[] = f[]*p[];
    boundary ({l});
  output_ppm (l, fp3, min = 0, max =1, linear = true,
              n = 2048, box = {{0,0},{Lz,Lz}});
}
#else
event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = level*(1-f[]);
    boundary ({l});
    output_ppm (l, file = "level.mp4", min = 0, max = LEVEL,
                n = 2048, box = {{0,-1},{Lz,Lz}});
    
    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    output_ppm (l, file = "velo.mp4", min = 0, max = 1.5, linear = true,
                n = 1024, box = {{0,-1},{4,1.5}});
    
    foreach()
    l[] = f[]*p[];
    boundary ({l});
    output_ppm (l, file = "f.mp4", min = 0, max =1, linear = true,
                n = 2048, box = {{0,0},{Lz,Lz}});
}

#endif

vector h[];
event timeseries (t += 0.1 ) {
    heights (f, h);
    double maxy = - HUGE,maxx = - HUGE;;
    foreach()
    if ((h.y[] != nodata) && (h.x[] != nodata)) {
        double yi = y + height(h.y[])*Delta;
        double xi = x + height(h.x[])*Delta;
        if (yi > maxy)
            maxy = yi;
        if (xi > maxx)
            maxx = xi;
    }
    char s[80];
    sprintf (s, "hmax-%g",Co);
    static FILE * fp0 = fopen (s, "w");
    fprintf (fp0, "%g %g %g\n", t, maxx, maxy);
    fflush (fp0);
    xfront=maxx;
}

event pictures (t = {0, 1. , 2., 3., 4.} ) {
    if(Co<0.0000001){
    scalar l[];
    foreach()
      l[] = f[]*p[];
    boundary ({l});
    if(t==0.)
    output_ppm (f, file = "f0.png", min=-1, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{4,2}});
    if(t==1.)
    output_ppm (f, file = "f1.png", min=-1, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{4,2}});
    if(t==2.)
    output_ppm (f, file = "f2.png", min=-1, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{4,2}});
    if(t==3.)
    output_ppm (f, file = "f3.png", min=0, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{4,2}});}

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
make granular_column_cohesif.tst
make granular_column_cohesif/plots
make granular_column_cohesif.c.html;
~~~

## Results

 First we check that we have exactly the same results than in the pure
 [dry column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c) test
 (small differences are maybe due to mesh adaptation) and
 with superposition of contact dynamics method DCM:
 
~~~gnuplot
 d=0.005
 h0=0.149
 set xlabel "x"
 set ylabel "h(x,t)"
 
 p[0:4][0:1.5]'../../granular_column/ShapeTime.A-01.dat' u (($1*d)/h0):(($2*d)/h0) t'DCM'w l,\
 'out000' w l t 'this code','../REFCASES/out05_granular_column.dat' t'for control' w l
~~~
 

Comparison for several values of $Co$ of the change of chape at time 1,2,3 and 4

~~~gnuplot collapse  Continuum Model various Co
set xlabel "x"
set ylabel "h(x,t)"
p[0:4][0:1.5]'out000' w l t 'Co=0.0',\
             'out005' w l t 'Co=0.05',\
             'out010' w l t 'Co=0.10','out015' w l t 'Co=0.15'
~~~


 
Iso plot of velocity and pressure at $t=2$
 

~~~gnuplot velocity and pressure
set pm3d map
set palette rgbformulae 22,13,-31;
unset colorbox
set multiplot layout 2,1
set xlabel "x  iso u"
set ylabel "y"
 splot [:4][:1] 'field-3-0.15.txt' u 1:2:($3>.999? $5*$3 :NaN)   not
set xlabel "x  iso p"
 splot [:4][:1] 'field-3-0.15.txt' u 1:2:($3>.999? $4*$3 :NaN)   not
unset multiplot
reset

~~~

a zoom on the solid corner... at time 1,2,3 and 4

~~~gnuplot
set view 0,0
set contour;unset surface
splot [0:3][0:1.5] 'field-1-0.1.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)  w l not, \
                  'field-2-0.1.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)   w l not,\
                  'field-3-0.1.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)   w l not,\
                  'field-4-0.1.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)   w l not
~~~

Looking at the rigid corner
 
~~~gnuplot
set view 0,0
set contour;unset surface
splot [0:3][0:1.5] 'field-2-0.15.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)  w l not, \
 'field-2-0.15.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)   w l not,\
 'field-2-0.15.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)   w l not,\
 'field-2-0.15.txt' u 1:2:(sqrt($5*$5+$6*$6)*$3)   w l not
~~~
 
Comparaison of runout for different cohesions as function of time
 
~~~gnuplot
 reset
 set key left
 set xlabel "t"
 set ylabel "h_m(t),x_m(t)"
 p'hmax-0' u 1:2 w l t 'Co=0.00','' u 1:3 w l t 'Co=0.00',\
 'hmax-0.05' u 1:2 w l t 'Co=0.05','' u 1:3 w l t 'Co=0.05',\
 'hmax-0.1' u 1:2 w l t 'Co=0.10','' u 1:3 w l t 'Co=0.10',\
 'hmax-0.15' u 1:2 w l t 'Co=0.15','' u 1:3 w l t 'Co=0.15'
~~~
 
 
 Comparaison of runout for 1D Savage Hutter,   2D Multilayer (RNSP), and 2D NS
 as function of cohesion
 
~~~gnuplot
 reset
 set xlabel "tauY"
 set ylabel "(x max- x init)"
 p'tauYxmax.txt' u ($1):($2-2)  w lp  t'NS',\
  '../cohesive_muI_collapse_ML/tauYxmax.txt'u ($1):($2-2)  t '2D RNSP' w l,\
  '../cohesive_savagehutter/tauYxmax.txt' u ($1):($2-2) t'1D SH' w l
~~~

 
 
 
 
Some films,
Velocity


[![](./granular_column_cohesif/f0.png)](./granular_column_cohesif/velo.mp4)
  velocity  (click on image for animation)

 level
 
 [![](./granular_column_cohesif/f1.png)](./granular_column_cohesif/level.mp4)
 velocity  (click on image for animation)
 
collapse
 
 [![](./granular_column_cohesif/f2.png)](./granular_column_cohesif/f.mp4)
 velocity  (click on image for animation)
 

## Related examples, links
 

 * [The basic Bagnold flow](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic_cohesif.c]()

 * [dry column](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c)

 * [silos](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c)

 
 * le 1D [http://basilisk.fr/_edit/sandbox/M1EMN/Exemples/cohesive_savagehutter.c]()
 
 * le 2D RNSP en multilayer
 [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapse_ML.c]()


## Bibliography

* Lagrée, Staron, Popinet 
["The granular column collapse as a
continuum: validity of a two–dimensional
Navier–Stokes model with a μ(I)-rheology"](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf) J. Fluid Mech. 2011 

 
 * Anaïs Abramian,    Lydie Staron, Pierre-Yves Lagrée
 "[http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/8.0000049.pdf](The Slumping of a Cohesive Granular Column: Continuum and Discrete Modelling)"
 Journal of Rheology (Vol.64, Issue 5).
 https://doi.org/10.1122/8.0000049
 
 
 
 
CISM 2019
*/
