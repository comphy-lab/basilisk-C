/**
# Spurious velocity at the intersections between non-periodic faces in a rotational stokes flow

We compute here a rotational flow with $\omega = [0,0,1]$. The code
generates a spurious *u.z* velocity at the intersections between the
non-periodic faces when the grid is not initialiazed at the maximum
level of refinement. The amplitude of the spurious *u.z* velocity is
independent of the minimum and maximum refinement levels and is equal
to $|8|$. This spurious velocity triggers mesh adaptation,
unnecessarily increasing the cell load and therefore the computational
time. This result can be reproduced for other rotation directions (*x*
and *y*).

We will solve the Navier-Stoles equations on an adaptive grid. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"

int lmin = 3; // Min mesh refinement level
int lmax = 7; // Max mesh refinement level
double cmax = 1e-2; // Refinement criteria for the velocity field

int l, lm, lp;
double uzavg, uzmin, uzmax;

int main ()
{  
  /**
  The domain is $128\times 128 \times 128$. */

  L0 = 128.;
  size (L0);
  origin (-L0/2., -L0/2., -L0/2.);

  /**
  We set front periodic boundary conditions. */
  
  periodic (front);
  
  DT = 1e-1; // The characteristic time is *D^2/nu = 1*.
 
  /**
  We tune the Poisson solver. Note that the *TOLERANCE* of the
  multigrid solver is divided by $\Delta t^2$. */
    
  stokes = true;
  TOLERANCE = 1e-4*sq (DT);
  
  /**
  We initialize the grid at different refinement levels. */

  FILE * fp = fopen ("results.dat", "w");
  
  for (lm = lmin; lm < 5; lm++) {
    for (lp = 5; lp <= lmax; lp++) {
      for (l = lm; l <= lp; l++) {

	N = 1 << l;
	init_grid (N);
	run();
	
	fprintf (fp, "%d %d %g %g %g %g\n",
		 lm, lp, ((double) (l - lm)/(lp - lm)), uzavg, uzmin, uzmax);
	fflush (fp);
      }
    }
  }
  fclose (fp);
}

/**
## Boundary conditions

We apply rotational flow boundary conditions $\omega \times x$. */

u.n[left] = dirichlet (-y);
u.t[left] = dirichlet (x);
u.r[left] = dirichlet (0);
/* uf.n[left] = dirichlet (-y); */
/* uf.t[left] = dirichlet (x); */
/* uf.r[left] = dirichlet (0); */
p[left] = neumann (0);
/* pf[left] = neumann (0); */

u.n[right] = dirichlet (-y);
u.t[right] = dirichlet (x);
u.r[right] = dirichlet (0);
/* uf.n[right] = dirichlet (-y); */
/* uf.t[right] = dirichlet (x); */
/* uf.r[right] = dirichlet (0); */
p[right] = neumann (0);
/* pf[right] = neumann (0); */

u.n[bottom] = dirichlet (x);
u.t[bottom] = dirichlet (0);
u.r[bottom] = dirichlet (-y);
/* uf.n[bottom] = dirichlet (x); */
/* uf.t[bottom] = dirichlet (0); */
/* uf.r[bottom] = dirichlet (-y); */
p[bottom] = neumann (0);
/* pf[bottom] = neumann (0); */

u.n[top] = dirichlet (x);
u.t[top] = dirichlet (0);
u.r[top] = dirichlet (-y);
/* uf.n[top] = dirichlet (x); */
/* uf.t[top] = dirichlet (0); */
/* uf.r[top] = dirichlet (-y); */
p[top] = neumann (0);
/* pf[top] = neumann (0); */

/**
## Initial conditions
*/

event init (t = 0)
{  
  foreach() {
    u.x[] = -y;
    u.y[] = x;
    u.z[] = 0.;
  }
  boundary ((scalar *) {u});
}

/**
## Movie */

event movie (i++)
{
  if (lm == 3 && lp == 7 && l == 3) {
    clear ();
    view (fov = 30,
  	  tx = 0., ty = 0., bg = {1,1,1},
  	  width = 400, height = 400,
  	  camera = "iso");
    cells (n = {1,0,0}, alpha = -64);
    cells (n = {0,1,0}, alpha = -64);
    cells (n = {0,0,1}, alpha = -64);
    squares ("u.z", n = {1,0,0}, alpha = -64);
    squares ("u.z", n = {0,1,0}, alpha = -64);
    squares ("u.z", n = {0,0,1}, alpha = -64);
    save ("mesh-coarse.mp4", opt = " -r 1");
  }
  if (lm == 3 && lp == 7 && l == 7) {
    clear ();
    view (fov = 30,
  	  tx = 0., ty = 0., bg = {1,1,1},
  	  width = 400, height = 400,
  	  camera = "iso");
    cells (n = {1,0,0}, alpha = -64);
    cells (n = {0,1,0}, alpha = -64);
    cells (n = {0,0,1}, alpha = -64);
    squares ("u.z", n = {1,0,0}, alpha = -64);
    squares ("u.z", n = {0,1,0}, alpha = -64);
    squares ("u.z", n = {0,0,1}, alpha = -64);
    save ("mesh-fine.mp4", opt = " -r 1");
  }
}

/**
## Adaptive mesh refinement */

event adapt (i++)
{
  adapt_wavelet ({u.x,u.y,u.z}, (double[]){cmax,cmax,1e-1}, maxlevel = (lp), minlevel = (lm));
}

event stop (t = 0.7) {
  uzavg = normf(u.z).avg;
  uzmin = statsf(u.z).min;
  uzmax = statsf(u.z).max;
  return 1;
}

/**
## Results

We observe that when the mesh is initialized at the maximum level of
refinement, there is no spurious *u.z* velocity.

![Time evolution of the initial coarse mesh](rotationz/mesh-coarse.mp4)
![Time evolution of the initial fine mesh](rotationz/mesh-fine.mp4)

~~~gnuplot Average, min and max values of *u.z* at the final time step as a function of the initial grid refinement level *l*
reset
set xlabel '(l-lm)/(lp-lm)'
set ylabel 'u_z'
plot 'results.dat' u 3:4 w p ps 1.5 lc rgb 'blue' t 'avg', 'results.dat' u 3:5 w p ps 1.5 lc rgb 'red' t 'min', 'results.dat' u 3:6 w p ps 1.5 lc rgb 'green' t 'max'
~~~

~~~
*/
