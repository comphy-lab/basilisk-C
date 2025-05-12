/**
# Marangoni flow vs Laplace Pressure

We test here the ability of the multilayer solver to reproduce correctly the Marangoni flow given a surface tension gradient. 
This is a qualitative test controlled by surface tension and under zero-gravity.

A periodic tank of finite dimensions containing a liquid film of mean height $h_0$ is considered. 
A periodic surface tension is imposed on the surface of the liquid (phisically, such a gradient could be obtained for a thin film with a temperature gradient).
Marangoni flows tend to accumulate liquid on the low surface tension areas whereas Laplace pressure tends to flatten the surface.
 
![Surface Profil](test_Marang3D/movie.mp4)
*/ 

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid.h"
#include "view.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

#define h_ 0.5
#define T_END 5.
#define DELTA_T (T_END/500.)

/**
## Main function
The boundary conditions are set to periodic. The test is done with zero-gravity. */
int main()
{
  origin (-L0/2., -L0/2.);
  foreach_dimension()
    periodic (right);
  N = 64;
  nl = 2;
  G = 0;
  nu = 10;
  run();
}

/**
## Initialisation 
A wave relaxe slowly.*/
event init (i = 0)
{
  foreach() {
    zb[] = - h_;
    double H = h_;
    foreach_layer() {
      h[] = H/nl;
      sigma[] = 1. - 0.5*cos(2*pi*sqrt(sq(x)+sq(y)));
    }
  }
}

/**
# Movie 
This is how we export the headline movie. */
event movie (t = 0.; t += DELTA_T; t <= T_END)
{
  view (fov = 20, quat = {0.475152,0.161235,0.235565,0.832313},
	tx = 0, ty = -0.015, width = 600, height = 400);
  char s[80];
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 80);
  squares ("eta", linear = true, z = "eta", min = -0.1, max = 0.05);
  save ("movie.mp4");
}