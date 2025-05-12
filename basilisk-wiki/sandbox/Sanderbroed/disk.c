/**
#Numerical setup for the fruit-frost based example

*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"  
#include "navier-stokes/perfs.h"
#include "fractions.h"
#include "view.h"
#include "../Antoonvh/slicer.h"
#include "run.h"

u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);

scalar b[], * tracers = {b};
face vector av[];

double Uo = 5, Ro = 0.75; //SI units
double to, nu, V, Nbv;
double Pi1 = 7, Re = 4000, Pi3 = 500000, Pi4 = 0.02; // =to*Uo/R;

double height  = 10, width = 0.25; //cannon mounting height
double dir = 0.4; // injector angle in x-y plane
scalar cannon[]; // Cannon volume fraction field

b[bottom] = dirichlet(0.);
b[top] = dirichlet(sq(Nbv)*L0);

int maxlevel = 8;

int main() {
  periodic (left);
  periodic (back);
  //cannon is a fraction field
  cannon.refine = cannon.prolongation = fraction_refine;
  foreach_dimension()
    u.x.refine = refine_linear;
  p.refine  = pf.refine  = refine_linear; //exact for linear stratifications
  L0 = 20.;
  X0 = Z0 = -L0/2;
  to = Pi1*Ro/Uo;
  nu = Uo*Ro/Re;
  Nbv = Pi4*Uo/Ro;
  a = av; //gravity
  const face vector muc[] = {nu, nu, nu};
  mu = muc;
  V  = Uo/(Pi3*Ro);
  run();
}

//Cannon:

event init (t = 0) {
  refine (sq(x) + sq(y - height) + sq(z) < sq(3*Ro) && level < maxlevel -2);
  refine (sq(x) + sq(y - height) + sq(z) < sq(2*Ro) && level < maxlevel - 1);
  refine (sq(x) + sq(y - height) + sq(z) < sq(Ro) && level < maxlevel);
  scalar sph[], planeup[], planedown[];
  fraction(sph, -sq(x) - sq(y - height) - sq(z) + sq(Ro));
  fraction(planeup, -cos(dir)*(y - height - cos(dir)*width/2.) - sin(dir)*(x - sin(dir)*width/2.));
  fraction(planedown, cos(dir)*(y - height + cos(dir)*width/2.) + sin(dir)*(x + sin(dir)*width/2.));
  foreach ()
    cannon[] = sph[]*planeup[]*planedown[];
  boundary ({cannon});
  foreach()
    u.z[] = y*V;
  boundary({u.z});
  foreach()
    b[] = sq(Nbv)*y;
  boundary({b});
  DT = 0.1;
}

event acceleration (i++) {
  coord grav = {0, 1, 0};
  boundary({b});
  foreach_face()
    av.x[] = grav.x*face_value(b,0);
}

event vortex_generator (i++; t <= to) {
  coord cdir = {-sin(dir), -cos(dir), 0};
    foreach() {
      foreach_dimension()
	u.x[] = u.x[]*(1. - cannon[]) + cannon[]*Uo*cdir.x;
    }
    boundary ((scalar*){u})
}


event done_generating (t = to) { //do not draw unresolved facets
  foreach()
    cannon[] = 0;
}

event adapt (i++) {
  double ue = Uo/20;
  adapt_wavelet ({b, u}, (double[]){sq(Nbv)*Ro*2, ue, ue, ue}, maxlevel);
}

#include "lambda2.h"
event mov (t += 0.05) {
  view (phi = 0.4, ty = -0.2);
  scalar l2[];
  lambda2(u,l2);
  isosurface("l2", -0.1);                   //Vortex
  draw_vof ("cannon", fc = {0.8, 0.6, 0.5});//Cannon
  squares ("x + z/2.", n = {0,1,0}, min = -3*L0 , max = 3*L0); // A surface
  translate(z = -L0/2){
  cells();
  squares ("b", min = 0, max = sq(Nbv)*height, linear = true);
  }
  save ("mov7.10.001.mp4");
}
/**
![](disk/mov7.10.001.mp4)

## Output

Before generating slices and compute area fractions, the buoyancy is
converted to a temperature field in $K$elvin.

$$T = T_{ref} + \frac{T_{ref}*b}{g},$$, 

with $T_{ref}$ the surface value. 

This gives easier plots. 

~~~gnuplot Tempeature at $t = 6s$ at $2m$ above the surface 
set xlabel 'z [px]'
set ylabel 'x [px]'
set cbrange [270:272]
set xrange [0:255]
set yrange [0:255]
set size ratio 1
set cblabel 'Temperature [K]'
plot 'slice016' matrix with image
~~~


~~~gnuplot Tempeature at $t = 10s$ at $2m$ above the surface 
set xlabel 'z [px]'
set ylabel 'x [px]'
set cbrange [270:272]
set xrange [0:255]
set yrange [0:255]
set size ratio 1
set cblabel 'Temperature [K]'
plot 'slice0110' matrix with image
~~~

~~~gnuplot Heated area evolution. 
  reset
  set size square
  set xlabel 't [s]'
  set ylabel 'A [m^2]
  set grid
  plot 'areas01' u 1:5 w l lw 3 t '>271.5 K',\
       '' u 1:6 w l lw 3 t '>271 K',\
       '' u 1:7 w l lw 3 t '>272.5 K'
~~~
 */

double T_ref = 270;
event slices (t += 1) {
  scalar T[];
  foreach()
    T[] = T_ref + T_ref*b[]/9.81;
  boundary ({T});
  double yp = 2;
  char fname[99];
  sprintf (fname, "slice01%g", t);
  sliceXZ (fname, T, yp, maxlevel);
}
#define NR_VALS (6) 

event areas (t += 0.1) {
  scalar T[];
  foreach()
    T[] = T_ref + T_ref*b[]/9.81;
  boundary ({T});
  double yp = 2;
  static FILE * fp = fopen ("areas01", "w");
  
  double values[NR_VALS] = {270, 270.5, 271, 272.5, 273, 274}; //values in array
  
  //header:
  if (t == 0) {
    fprintf (fp, "#t");
    for (int i = 0; i < NR_VALS; i++)
      fprintf (fp, "\t>%g", values[i]);
    fprintf (fp, "\n");
  }

  //values:
  fprintf (fp, "%g", t);
  for (int i = 0; i < NR_VALS; i++) {
    double A = 0;
    foreach(reduction (+:A)) {
      if (fabs(y - yp) < Delta/2.)                  //slice is in cell...
	if (interpolate (T, x, yp, z) > values[i])  // sufficient value
	  A += sq(Delta);                           // Area of slice in sufficient cell
    }
    fprintf (fp, "\t%g", A);
  }
  fprintf (fp, "\n");
}

event stop (t = 60); //one minute

