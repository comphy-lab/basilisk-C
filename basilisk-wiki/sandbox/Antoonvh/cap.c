/**
# Why do plants transpire?

To pump nutrients through their capillaries:

![Without transpiration, a balance at rest emerges. With evaporation:
 the capillary pumps fluid in a quasi-steady balance.](cap/movie.mp4)
*/
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#define BVIEW 1
#include "particles.h"

/**
The contact angle is setup. 
 */
bid capillary;
vector h[];
double theta0 = 50.;
h.t[left] = contact_angle (theta0*pi/180.);
h.t[capillary] = contact_angle (x < 2 ? theta0*pi/180.: pi/2.);

face vector av[];

int main() {
  L0 = 11.;
  mu1 = mu2 = 0.1;
  a = av;
  rho1 = 2.;
  rho2 = 1.;
  f.sigma = 0.0015;
  f.height = h;
  N = 32;
  run();
}

event init (t = 0) {
  foreach()
    f[]  = (y < 4);
  /**
Contant angles do work on straight `mask`ed cells.
   */
  mask (y > 2 && y < 8 && x > 1 && x < 3 ? capillary : none);
  /**
Nutrients are suspended in the water.
   */
  foreach() 
    if (f[] >= 0.9)
      n_part++;
  loc = (coord*) malloc (sizeof(coord)*n_part);
  n_part = 0;
  foreach()
    if (f[] >= 0.9) {
      coord new = {x, y, z};
      loc[n_part++] = new;
    }
}

/**
Gravity is an important part of the column's force balance.
 */
event acceleration (i++) {
  foreach_face(y)
    av.y[] = -0.001;
}


/**
Evaporation of the fluid in the cappilary is hacked in via the event
below. The aim is to implement a fixed evaportation at the surface of
the column.
 */
bool evaporate = false;
event turn_on_evaporation (t = 2000)
  evaporate = true;

event evaporation (i++) {
  if (evaporate) {
    double evap_rate = 0.01;
    foreach() {
      if (f[] > 1e-5 && f[] < 1. - 1e-5 && x < 2.) {
	coord n = interface_normal (point, f), p;
	double alpha = plane_alpha (f[], n);
	double l = plane_area_center (n, alpha, &p);
	double flux = dt*l*evap_rate;
	f[] = max (f[] - flux, 0);
      }
    }
  }
}

event movie (t += 10; t <= 12500) {
  char str[99];
  scalar white[];
  foreach()
    white[] = 0;
  view (tx = -0.5, ty = -0.5, bg = {0.9, 0.9, 0.9});
  if (evaporate)
    sprintf (str, "With evaporation");
  else
    sprintf (str, "Without evaporation");
  draw_string (str, size = 30);
  draw_vof ("f");
  draw_vof ("f", filled = 1, fc = {0.2, 0.4, 0.9});
  scatter (loc);
  translate (z = -0.1)
    squares ("white", map = cool_warm, min = -1., max = 1);
  save ("movie.mp4");
}
 
event adapt (i++) 
  adapt_wavelet ({f}, (double[]){0.001}, 8, 6);

