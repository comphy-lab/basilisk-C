/**
# Marangoni flow vs Laplace Pressure

We test here the ability of the multilayer solver to reproduce correctly the 
Marangoni flow given a surface tension gradient. 
This is a qualitative test controlled by surface tension and under zero-gravity.

A periodic tank of finite dimensions containing a liquid film of mean height 
$h_0$ is considered. A periodic surface tension is imposed on the surface of 
the liquid (phisically, such a gradient could be obtained for a thin film with 
a temperature gradient). Marangoni flows tend to accumulate liquid on the low 
surface tension areas whereas Laplace pressure tends to flatten the surface.
 
![Surface Profil](marangoni_stress_3D/movie.mp4)
*/  

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid.h"
#include "view.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/implicit.h"
#include "layered/remap.h"

#define h_0 0.5 //adimensionned (epsilon)
#define T_END 20.
#define DELTA_T (T_END/200.)
const double k =1.;

/**
## Main function
The boundary conditions are set to periodic. The test is done with zero-gravity. */
int main()
{
  L0 = 2[1]*pi; 
  origin (-L0/2., -L0/2.);
  foreach_dimension()
    periodic (right);
  N = 64;
  nl = 2;
  G = 0;
  nu = 10.;
  CFL = 0.5;
  TOLERANCE = 1e-6; 
  linearised = true;       
  run();
}

/**
## Initialisation */
scalar sig[];
event init (i = 0)
{
  foreach() {
    zb[] = - h_0;
    double H = h_0/k;
    foreach_layer() {
      h[] = H/nl;
      sig[] = 1.[3,-2] * (1. - 0.8 * cos(sqrt(sq(x)+sq(y))*k));
    }
  }
  sigma = sig;
}

/**
# Movie 
This is how we export the headline movie. */
event movie (t += DELTA_T; t <= T_END)
{
  view (fov = 20, quat = {0.475152,0.161235,0.235565,0.832313},
	tx = 0, ty = -0.015, width = 600, height = 400);
  char s[80];
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 80);
  squares ("eta", linear = true, z = "eta", min = -0.05, max = 0.05);
  save ("movie.mp4");
}
