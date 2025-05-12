/**
# MWE for vertical solute diffusion.
Solute diffusion answers to a Fick diffusion equation of a specie with a concentration $c$.
$$
\frac{\partial c}{\partial t} =\nabla^2 c 
$$
In this minimum working example, this equation is solved vertically with a Dirichlet boundary condition at the bottom and a Neumann boundary condition at the top.*/

#include "grid/multigrid1D.h"
#include "../hydro_c.h"
#include "layered/nh.h"
#include "layered/remap.h"


double ak = 0.35;
double RE = 100.;
double Delta;

/**
## Geometry and resolution*/
#define L 10.
#define LEVEL 1

#define T_END 5.
#define DT_MAX 0.5
#define DELTA_T (T_END/10.)

/**
## Constants*/
#define k_  (2.*pi)
#define h_   0.5
#define g_   1.
#define T0  (k_/sqrt(g_*k_))

/**
## Physical parameters
In the bulk:
* Diffusion coefficient of the solute in the liquid
* Initial concentration (molar fraction) of the solute*/
#define Diff 1
#define solute0 0

scalar c_bb[];

int main()
{
  L0 = L;
  origin (-L0/2.);
  periodic (right);
  N = 1 << LEVEL;
  Delta = L0/N;
  nl = 100;
  G = g_;
  nu = 1./RE;
  D = Diff;
  bottom_boundary_neumann = false;
  c_b = c_bb;
  run();
}


event defaults (i = 0)
{
  foreach()
    c_bb[] = 1.;
  assert (cl == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar c = new scalar;
    cl = list_append (cl, c);
  }
  reset (cl, 0.);
  int l = 0;
  for (scalar c in cl) {
    tracers[l] = list_append (tracers[l], c);
    l++;
  }
}

event init (i = 0)
{
  foreach() {
    double H = 10.;
    double z = 0;
    vector u;
    scalar h, w, c;
    for (h,u,w,c in hl,ul,wl,cl) {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = 0;
      w[] = 0;
      z += h[]/2.;
	  	c[] = solute0;
    }
  }
}


/**
## Results*/
event profiles (t += DELTA_T; t = DELTA_T; t <= T_END) {
  static FILE * fpx = fopen("datax.txt", "w");
    foreach(){
      double y = 0;
      scalar h,c;
      for (c, h in cl, hl){
        y += h[]/2;
        fprintf (fpx, "%g %g %g %g %g \n", y, t, y/(2*sqrt(t)), c[], c[] - 1 + erf(y/(2*sqrt(t))));
        y += h[]/2;
      }   
  fprintf (fpx, "\n");
    }
}

/**
We compare the exact and computed solutions for a cross-section at
$x=0$.

~~~gnuplot Concentration profil over $y$
set terminal @PNG enhanced size 640,640 font ",8"
set output 'profil.png'
L0=10
set xlabel "y"
set ylabel "concentration"
set key left
plot \
  './datax.txt' u 1:4 t 'c' w l
~~~

~~~gnuplot Rescaled concentration profil over $y$
set output 'rescaled.png'
set xlabel "y/2sqrt(t)"
set ylabel "concentration"
plot \
  './datax.txt' u 3:4 t 'c' w l, \
  1-erf(x) t 'erf'
~~~

~~~gnuplot Comparison with erfc
set output 'plot3.png'
set xlabel "y/2sqrt(t)"
set ylabel "Error"
plot \
  './datax.txt' u 3:5 t 'c' w l
~~~
There is a strong error (10 %) on the first time step that reduce quickly in the next steps. This error do not depend on the number of layers.
*/
