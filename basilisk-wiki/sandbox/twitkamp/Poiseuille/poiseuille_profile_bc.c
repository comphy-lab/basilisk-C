/**
# Poiseuille embed with velocity boundary condition

The Developed Poiseuille Velocity-profile dictates a consistent (and non-zero) pressure-gradient forcing, also at the boundaries!

in 2D:
$$u_x(y) = -\frac{\partial p}{\partial x} \frac{1}{2\mu}y(H-y)$$

![$u_x$ looks good](poiseuille_profile_bc/ux.png)

![$p$ looks good](poiseuille_profile_bc/p.png)

![Error field shows some inlet adjustment effect](poiseuille_profile_bc/e.png)

*/

double dpdx = -1;
double muv = 1;
double H = 0.3;

#define U_X(y) ( (y) < 0 ? 0 : ((y) > H ? 0 :(-dpdx/(2*muv)*(y)*(H - (y))))) 

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

double EPS = pi*1e-2; // This is an imporantant value
// u.n[left] = dirichlet(1.);
u.n[left] = dirichlet(U_X((y + EPS)));

// uf.n[left] = dirichlet(1.);
u.n[right] = neumann(0.);
/**
We should set a consistent pressure drop
 */

p[left] = dirichlet(-dpdx*L0);
// p[right] = dirichlet(0.); // This is fine
p[right] = neumann(dpdx);    // this too
int main()
{
  
  size (1. [0]);
  origin (-L0/2., -L0/2.);

  stokes = true;
  // quadratic interpolation for quadratic profile
  u.x.third = true;
  TOLERANCE = 1e-5;
  for (int lvl = 4; lvl <= 7; lvl++){ // lvl = 7 can be problematic due to small fluid-fraction cells
    N = pow(2, lvl);
    EPS = L0/(2*N);
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
  The geometry of the embedded boundary is defined as two plates located at $y = 0$ and $y = H$. */
  //EPS = L0/N;
  solid (cs, fs, intersection( H - (y + EPS), 0 + (y + EPS)));
  multigrid_restriction ({cs, fs});
  /**
   # Embedded boundary conditions

  */
  
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  /**
  We initialize the reference velocity field. */
  
  // foreach()
  //   un[] = u.y[];
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-8)
    return 1; /* stop */
}

/**
We compute error norms and display the horizontal velocity, pressure and
error fields using bview. */



event profile (t = end)
{
  scalar utheta[], e[];
  foreach() {
    if (cs[] > 0.) {
      e[] = u.x[] - U_X(y + EPS);
    }
    else
      e[] = p[] = utheta[] = nodata;
  }

  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  
  dump();
  
  draw_vof ("cs", "fs", filled = -1, fc = {0, 0, 0});
  squares ("u.x", spread = -1);
  save ("ux.png");

  draw_vof ("cs", "fs", filled = -1, fc = {0, 0, 0});
  squares ("p", spread = -1);
  save ("p.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("e", spread = -1);
  save ("e.png");
  
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
Not sure why Velocity is no showing but the Mass flux $Q =\int u dy$ is not constant in space. This is an issue,
![Horizontal velocity](ux.png)


*/
