/**
# Relaxation due to Laplace pressure in 3D
Here, we test qualitatively the ability of the multilayer solver to capture Laplce pressure.
A periodic wave profil is initialized. Surface tension is supposed to flatten the surface. 

![Surface profil relaxation](relaxation3D/movie.mp4)
*/ 
 
/**
## Include*/
#include "grid/multigrid.h"
#include "view.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

/**
## Geometry and resolution*/
#define L 2.
#define h_ 0.5
#define ah h_/5
#define LEVEL 6
#define layers 3
#define T_END 2.
#define DELTA_T (T_END/100.)

/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_ 0.
#define Re 0.01

/**
## Main function
The test is done without any viscosity or gravity. */
int main()
{
  L0 = L;
  origin (-L0/2., -L0/2.);
  foreach_dimension()
    periodic (right);
  N = 1 << LEVEL;
  nl = layers;
  G = g_;
  nu = 1/Re;
  run();
}

/**
## Initialisation 
A wave relaxe slowly.*/
event init (i = 0)
{
  foreach() {
    double H =  h_ + ah*cos(k_*x)*cos(k_*y);
    zb[] = - h_;
    foreach_layer() {
      h[] = H/nl;
      sigma[] = 1.;
    }
  }
}

/**
## Outputs*/
event outputs(t = 0.; t += DELTA_T; t <= T_END) {
  if (i < 10) {
    static FILE * fpp = fopen ("phi.txt", "w");
    foreach() {
      fprintf (fpp, "phi_th : %.4g", kappa_c(point));
//      foreach_layer ()
 //       fprintf (fpp, "      %.4g       ", phi[]);
      fprintf (fpp, "\n");
    }
    fprintf (fpp, "\n");
  }
}


/**# Movie
To see the movement of the wave*/
event movie (t = 0.; t += DELTA_T; t <= T_END)
{
  view (fov = 20, quat = {0.475152,0.161235,0.235565,0.832313},
	tx = 0.02, ty = -0.015, width = 600, height = 400);
  char s[80];
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 80);
  squares ("eta", linear = true, z = "eta", min = -0.15, max = 0.15);
  save ("movie.mp4");
}