/**
# Transcritical flow over a bump with multiple layers

We want to reproduce the transcritical test case of [Audusse et al,
2011](/src/references.bib#audusse2011), section 5.6.2. */

#include "grid/cartesian1D.h"
#include "saint-venant.h"

/**
We need a field to store the variable bottom friction. */

scalar lambda[];

/**
The constant inflow flow rate (m^2/s). */

const double Q = 1.;

/**
The outlet water level. */

const double Ho = 0.6;

int main() {
  X0 = 0.;
  L0 = 21. [1];
  G = 9.81;
  N = 256;
  
  /**
  The viscosity is set to $\nu = 0.01 m^2/s$ and the bottom friction
  is variable. */
  
  nu = 0.01;
  lambda_b = lambda;
  
  /**
  We vary the number of layers. */

  nl = 2;  run();
  nl = 5;  run();
  nl = 15; run();
}

/**
We impose the outlet water level. */

h[right]   = dirichlet (Ho);
eta[right] = dirichlet (Ho);

/**
## Initialisation

We initialise the topography, the initial water depth *h* and we
create a field *hc* to check convergence on *h*. */

scalar hc[];

event init (i = 0) {
  double a = 0.2, b = 5.75/2.;
  foreach() {
    zb[] = max(0., a*(1.[0] - sq((x - 10.)/b)));
    hc[] = h[] = Ho - zb[];
  }

  /**
  ## Boundary conditions on velocity
  
  We impose a constant flow rate `Q` at the inlet and a Neumann
  condition at the outlet.*/
  
  for (vector u in ul) {
    u.n[left] = dirichlet(h[left] ? Q/h[left] : 0.);
    u.n[right] = neumann(0.);
  }
}

/**
## Bottom friction

We use the Strickler relation:
$$
k(h,\mathbf{U}) = \frac{g}{S^2h^{1/3}}|\mathbf{U}| 
$$
with $S = 25 m^{1/3}/s$ the Strickler coefficient, $h$ the water depth
and $\mathbf{U}$ the depth-averaged velocity. Note that we have to use
a lower Strickler coefficient (i.e. larger friction) to get results
comparable to those of [Audusse et al,
2011](/src/references.bib#audusse2011). */

event friction (i++) {
  foreach() {
    double U = 0.;
    int l = 0;
    for (vector u in ul)
      U += u.x[]*layer[l++];
    double S = 25., k = G/(sq(S)*pow(h[],1./3.))*fabs(U);
    lambda[] = k > 0. ? nu/k : 0.;
  }
}

/**
We check for convergence. */

event logfile (t += 0.1; i <= 100000) {
  double dh = change (h, hc);
  if (i > 0 && dh < 1e-4)
    return 1;
}

/**
Uncomment this part if you want on-the-fly animation. */

#if 0
event output (i++) {
  static FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
           "set title 'nl=%d, t=%f'\n"
           "set xl 'x'\nset yl 'h'\n"
           "plot [0:21][] '-' u 1:2 w l t 'eta', '-' u 1:3 w l t 'zb'\n",
	   nl, t); 
  foreach()
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n");
  fflush (fp);
}
#endif

/**
## Outputs

At the end of the simulation we save the profiles. */

event output (t = end) {
  char name[80];
  sprintf (name, "end-%d", nl);
  FILE * fp = nl == 15 ? stderr : fopen (name, "w");
  foreach() {
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
    if (nl == 15) {
      double z = zb[];
      int l = 0;
      printf ("%g %g %g\n", x, z, u.x[]);
      for (vector u in ul) {
	z += layer[l++]*h[];
	printf ("%g %g %g\n", x, z, u.x[]);
      }
      printf ("\n");
    }
  }
}

/**
## Results

~~~gnuplot Free surface and topography. This can be compared to figure 9 of [Audusse et al, 2011](/src/references.bib#audusse2011).
set xr [0:21]
set yr [0:1]
set xlabel 'x'
set ylabel 'z'
plot 'end-2' u 1:3 w l t 'topography', \
     'end-2' u 1:2 w l t '2 layers', \
     'end-5' u 1:2 w l t '5 layers', \
     'log'   u 1:2 w l t '15 layers'
~~~

~~~gnuplot Horizontal velocity field (15 layers).
set term PNG enhanced font ",10"
set output 'vel.png'
set pm3d
set pm3d map interpolate 10,1
unset key
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,	\
0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,	\
0.875 0.9333 0 0, 1 0.498 0 0 )

splot 'out' u 1:2:3
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/shock.html#layered)
*/
