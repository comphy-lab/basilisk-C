/**
# Poiseuille flow in a corner
I simulate a poiseuille flow in a corner geometry. I impose a pressure boundary condition on the right and the top. 
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

u.n[right] = neumann(0.);
u.n[top] = neumann(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);

p[right] = dirichlet((double) N);

p[top] = dirichlet(0.);
double EPS;
int main()
{

  size (1. [0]);
  
  origin (0., 0.);
  
  stokes = true;
  
  TOLERANCE = 1e-5;
  for (int lvl = 7; lvl <= 7; lvl++){
    N = pow(2, lvl);
    run();
  }
}

scalar un[];

#define WIDTH 0.5

event init (t = 0) {

  /**
  The viscosity is unity. */
  
  mu = fm;

  /**
  The geometry of the system is a corner flow */
  EPS = L0/N;

  /**
  We can now initialize the volume fractions in the domain. */
  vertex scalar phi[], phi_vert[], phi_hor[];
  foreach_vertex(){
    phi_vert[] = -(x - 0.2 * L0);
    phi_hor[] = -(y - 0.2 * L0);
    phi[] = union(phi_hor[], phi_vert[]);
  }

    fractions (phi, cs, fs);
    boundary ({cs, fs});
    fractions_cleanup (cs, fs);
    boundary ({cs, fs});
    restriction ({cs, fs});

  /**
  # Embedded boundary conditions
  
  No slip boundary condition */
  
  u.x[embed] = y > 0. ? dirichlet(0.) :  dirichlet(0.);
  u.y[embed] = y > 0. ? dirichlet(0.) :  dirichlet(0.);

}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-12)
    return 1; /* stop */
}

/**
We compute error norms and display the horizontal velocity, pressure and
error fields using bview. */


#define poiseuille(y, N) -N/2. * (y*y - pow(L0/8. + EPS/4., 2))
scalar mag_u[];
event profile (t = end)
{
  scalar e[];
  foreach() {
    if (cs[] > 0.) {
      e[] = u.x[] - poiseuille (y, N);
      mag_u[] = sqrt(sq(u.x[]) + sq(u.y[]));
    }
  }

  norm n = normf (e);
  foreach(){
    if (cs[] > 0.){
      fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	     N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
    }
  }

  char file_v_y[80];
  char file_v_x[80];

  sprintf(file_v_x, "vel_along_x.dat");
  sprintf(file_v_y, "vel_along_y.dat");
  FILE * fp_y = fopen(file_v_y, "w");
  FILE * fp_x = fopen(file_v_x, "w");
  foreach(){
    if (cs[] > 0.){
      if (point.i == 90) fprintf(fp_x, "%g %g \n", y, u.x[]);
      if (point.j == 90) fprintf(fp_y, "%g %g \n", x, u.y[]);
    }
  }

  dump();
  
  draw_vof ("cs", "fs", filled = -1, fc = {0, 0, 0});
  squares ("u.x", spread = -1);
  save ("ux.png");
  
  draw_vof ("cs", "fs", filled = -1, fc = {0, 0, 0});
  squares ("u.y", spread = -1);
  save ("uy.png");
  
  draw_vof ("cs", "fs", filled = -1, fc = {0, 0, 0});
  squares ("mag_u", spread = -1);
  save ("magu.png");
  
  char file[80];
  sprintf(file, "data_couette_%d", N);
  FILE * fp = fopen(file, "w");
  
  foreach() {
    if (cs[] > 0.){
      fprintf (fp, "%g %g %g %g %g\n",
          y, u.x[], u.y[], p[], e[]);
    }
  }
}

/**
## Results

![Horizontal velocity](pressure_corner/ux.png)

![Vertical velocity](pressure_corner/uy.png)

![velocity magnitude](pressure_corner/magu.png)


The velocity is the same in the horizontal and vertical part of the channel but I am not sure what the theoretical should be.

~~~gnuplot Velocity profile in the horizontal and vertical part of the channel
L0 = 1.
N = 128.
set xlabel 'y'
set ylabel 'u_x'
plate_loc = (L0/8. + (L0/N/4.))
poiseuille(y, N) = - N/2 * (y*y - plate_loc * plate_loc);
set grid
# set arrow from 0.25, graph 0 to 0.25, graph 1 nohead
# set arrow from 0.5, graph 0 to 0.5, graph 1 nohead
plot [-0.02:0.22][0:0.4]'vel_along_y.dat' title 'vertical channel', 'vel_along_x.dat' u 1:(-1*$2) title 'horizontal channel'
~~~
*/
