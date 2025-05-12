/**
# A vortex ring injector mounted above a surface

This page examplifies how to setup a vortex-ring cannon above a
surface. The numerical cannon is a [puck](/sandbox/Antoonvh/puck.c),
targeting fluid in a `dir`ection. We choose `y` as the height
coordinate.

![The ring's fluid mixes with the atmosphere at the surface](diskinjector/mov.mp4)
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "view.h"

u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);

double Uo = 5, Ro = 0.75; //SI units
double to, nu; 
double Pi1 = 7, Re = 4000; //= to*Uo/Ro; 

double height = 10, width = 0.25;  // Cannon mounting height
double dir = 0.4; // injector angle in x-y plane
scalar cannon[];  // Cannon volume fraction field

int maxlevel = 8;

int main() {
  periodic (left);
  periodic (back);
  //cannon is a fraction field
  cannon.refine = cannon.prolongation = fraction_refine;
  foreach_dimension()
    u.x.refine = refine_linear;
  L0 = 20.;
  X0 = Z0 = -L0/2;
  to = Pi1*Ro/Uo;
  nu = Uo*Ro/Re;
  const face vector muc[] = {nu, nu, nu};
  mu = muc;
  run();
}
/**
## Vortex-ring cannon

During initialization, we compute the cell-volume fraction of the
cannon. `cannon[] = 1` for cells inside the cannon, `cannon[] = 0`
outside, and a fraction (`0 > cannon[] < 1`) for cells on the
cannon boundary.
 */

event init (t = 0) {
  refine (sq(x) + sq(y - height) + sq(z) < sq(3*Ro) && level < maxlevel - 2);
  refine (sq(x) + sq(y - height) + sq(z) < sq(2*Ro) && level < maxlevel - 1);
  refine (sq(x) + sq(y - height) + sq(z) < sq(Ro) && level < maxlevel);
  scalar sph[], planeup[], planedown[];
  fraction(sph, -sq(x) - sq(y - height) - sq(z) + sq(Ro));
  fraction(planeup,  -cos(dir)*(y - height - cos(dir)*width/2.) - sin(dir)*(x - sin(dir)*width/2.));
  fraction(planedown, cos(dir)*(y - height + cos(dir)*width/2.) + sin(dir)*(x + sin(dir)*width/2.));
  foreach ()
    cannon[] = sph[]*planeup[]*planedown[];
  boundary ({cannon});
}
/**
Depending on the cannon fraction, the fluid is forced in the
chosen `dir`ection at speed `Uo` for a time `to`. 
 */
event vortex_generator (i++; t <= to) {
  coord cdir = {-sin(dir), -cos(dir), 0};
  foreach() {
    foreach_dimension() 
      u.x[] = u.x[]*(1. - cannon[]) + cannon[]*Uo*cdir.x; 
  }
  boundary ((scalar*){u})
}

event done_generating (t = to) { // Do not draw unresolved facets
  foreach()
    cannon[] = 0;
}

event adapt (i++) {
  double ue = Uo/20;
  adapt_wavelet ((scalar*){u}, (double[]){ue, ue, ue}, maxlevel);
}

#include "lambda2.h"
event mov (t += 0.05) {
  view (phi = 0.4, ty = -0.2);
  scalar l2[];
  lambda2(u, l2);
  isosurface("l2", -0.1);                                     // Vortex
  draw_vof ("cannon", fc = {0.8, 0.6, 0.5});                  // Cannon
  squares ("x + z/2.", n = {0,1,0}, min = -3*L0, max = 3*L0); // A surface 
  cells();
  save ("mov.mp4");
}

event stop (t = 60); //one minute

/**
## To do

* Add relevant fruit-frost physics  
* Diagnose heating  
 */
