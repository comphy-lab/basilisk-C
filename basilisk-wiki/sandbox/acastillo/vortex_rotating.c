/**
# Two-dimensional dynamics of co-rotating vortices

This test is similar to [vortex.c]() but uses a non-inertial rotating frame with
a rotation rate $\Omega = (\Gamma_1+ \Gamma_2)/(2 \pi b^2)$, where $\Gamma_1$
and $\Gamma_2$ are the circulation of the two vortices which are separated by a
distance $b$. Here, the idea is to obtain a base state to be used for a linear
stability analysis.
*/


#define GAMMA 1.0
#define OMEGA (GAMMA/pi)
double EndTime=4*sq(pi);
double ad = 0.125, bd = 0.5;

/**
## Model equations
In this frame the incompressible Navier-Stokes equations read as
$$
\frac{d}{dt}\vec{u} + (\vec{u}\cdot\nabla)\vec{u}
= -\frac{1}{\rho}\nabla p + \nu\Delta\vec{u} - 2\vec{\Omega}\times\vec{u}
- \vec{\Omega}\times(\vec{\Omega}\times\vec{r})
$$
where the last two terms correspond to the Coriolis and centrifugal acceleration
respectively.

Let us consider a modified pressure
$$
p^* = p - \frac{1}{2} (\vec{\Omega}\times\vec{r})^2
$$
such that only the Coriolis term is required
$$
\frac{d}{dt}\vec{u} + (\vec{u}\cdot\nabla)\vec{u}
= -\frac{1}{\rho}\nabla p^* + \nu\Delta\vec{u} - 2\vec{\Omega}\times\vec{u}
$$

As in [vortex.c]() we use centered Navier–Stokes solver (without viscosity)
with the Coriolis term included in the acceleration term.

*/

#include "navier-stokes/centered.h"
#include "view.h"


event acceleration (i++) {
  coord cor = {-2.0*OMEGA, 2.0*OMEGA};
  face vector av = a;
  foreach_face(x)
    av.x[] += cor.x*uf.y[] + sq(OMEGA)*x;

  foreach_face(y)
    av.y[] += cor.y*uf.x[] + sq(OMEGA)*y;
}


/**
## Velocity at the boundaries
In the rotating frame, we expect velocity at the boundaries to match the
background velocity
$$
\vec{u}_b = (\nabla\times\psi) = \Omega (y \hat{e}_x - x \hat{e}_y)
$$
where
$$
\psi = -\frac{1}{2} \Omega (x^2 + y^2)
$$
*/

u.n[right]   = dirichlet( OMEGA*y);
u.n[left]    = dirichlet( OMEGA*y);
u.n[top]     = dirichlet(-OMEGA*x);
u.n[bottom]  = dirichlet(-OMEGA*x);

u.t[right]   = dirichlet(-OMEGA*x);
u.t[left]    = dirichlet(-OMEGA*x);
u.t[top]     = dirichlet( OMEGA*y);
u.t[bottom]  = dirichlet( OMEGA*y);

uf.n[right]   = dirichlet( OMEGA*y);
uf.n[left]    = dirichlet( OMEGA*y);
uf.n[top]     = dirichlet(-OMEGA*x);
uf.n[bottom]  = dirichlet(-OMEGA*x);

uf.t[right]   = dirichlet(-OMEGA*x);
uf.t[left]    = dirichlet(-OMEGA*x);
uf.t[top]     = dirichlet( OMEGA*y);
uf.t[bottom]  = dirichlet( OMEGA*y);

scalar psi[];
psi[right]  = dirichlet(-0.5*OMEGA*(sq(x) + sq(y)));
psi[left]   = dirichlet(-0.5*OMEGA*(sq(x) + sq(y)));
psi[top]    = dirichlet(-0.5*OMEGA*(sq(x) + sq(y)));
psi[bottom] = dirichlet(-0.5*OMEGA*(sq(x) + sq(y)));


/**
## Domain size and spatial resolution
The domain is centered on $(0,0)$ and the maximum level of refinement
is N i.e. the initial grid has $2^N$ grid points per dimension.
*/

#define MAXLEVEL 10
int main() {
  size(16);
  origin (-8, -8);
  init_grid (1 << MAXLEVEL);
  DT = 0.5;
  run();
}

#if TREE
event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){5e-5,5e-5}, MAXLEVEL);
}
#endif

/**
## Initial conditions
We include two Gaussian vortices
$$
\omega_v = \frac{\Gamma}{\pi a^2} \exp{\left(-(r^2/a^2)\right)}
$$
in addition to the background vorticity.
*/
event init (t = 0)
{
  a = new face vector;

  scalar omega[];
  foreach(){
    omega[] =  (GAMMA/(pi * sq(ad))) * (exp(-(sq(x - bd) + sq(y))/sq(ad))); // Vortex 1
    omega[] += (GAMMA/(pi * sq(ad))) * (exp(-(sq(x + bd) + sq(y))/sq(ad))); // Vortex 2
    omega[] -= 2*OMEGA; // Background rotation
    psi[] = -0.5*OMEGA*(sq(x) + sq(y));
  }
  boundary ({psi, omega});
  poisson (psi, omega);

  coord f = {-1.,1.};
  foreach()
  foreach_dimension()
    u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar*){u});
}

/**
# Outputs
We make animations of the vorticity.
![We follow the rotation of the vortex pair, we should see the elliptical deformation but little else](vortex_rotating/vortex_rotating_zoomed.mp4)
*/

event movie (t += 0.1; t <= EndTime){
  scalar omega[];
  vorticity (u, omega);
  view (fov=25);
  squares ("omega", linear = true);
  isoline ("omega", n = 20, min=0.05/(pi*sq(ad)), max=1/(pi*sq(ad)));
  box();
  save ("vortex_rotating.mp4");

  view (fov=8);
  squares ("omega", linear = true);
  isoline ("omega", n = 20, min=0.05/(pi*sq(ad)), max=1/(pi*sq(ad)));
  box();
  save ("vortex_rotating_zoomed.mp4");
}


/**
We also use [tag.h]() to compute the circulation, vortex center and second
moment of vorticity to validate our results.

~~~gnuplot
reset
set xlabel 'x'
set ylabel 'y'
set grid
plot 'vortex_rotating/vortex_rotating.asc' u 5:6
~~~

*/

#include "tag.h"
event logs (t += 0.1; t<=EndTime) {
  scalar omega[];
  vorticity (u, omega);
  stats s = statsf (omega);
  fprintf (ferr, "%g %d %g %g %g %d\n", t, i, dt, s.sum, s.max, mgp.i);

  scalar m[];
  foreach()
  m[] = omega[] > (GAMMA/(pi * sq(ad)))*0.001;
  int n = tag (m);

  double v[n], c[n];
  double b_x[n], b_y[n], b_z[n]; // Using this instead to coord, I have and issue with pointers.
  for (int j = 0; j < n; j++)
    v[j] = b_x[j] = b_y[j] = b_z[j] = c[j] = 0.;

  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*omega[];
      coord p = {x,y,z};
      foreach_dimension()
        b_x[j] += dv()*omega[]*p.x;
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b_x, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b_y, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b_z, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  foreach_leaf()
  if (m[] > 0) {
    int j = m[] - 1;
    coord p = {x,y,z};
    foreach_dimension()
    c[j] += dv()*omega[]*sq(b_x[j] - p.x);
  }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, c, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  FILE * fp = fopen("vortex_rotating.asc", "a");
  for (int j = 0; j < n; j++){
    fprintf (fp, "%d %g %d %g %g %g %g\n", i, t, j, v[j], b_x[j]/v[j], b_y[j]/v[j], sqrt(c[j]/v[j]));
    fprintf (stderr, "%d %g %d %g %g %g %g\n", i, t, j, v[j], b_x[j]/v[j], b_y[j]/v[j], sqrt(c[j]/v[j]));
  }
  fclose (fp);
}


/**
We may also take a snapshot of the velocity field and the streamfunction at
different times. Here, we store the fields as a Paraview-compatible `.vtu` file
and as a `basilisk` `.dump` file in case we want to restart the simulation.
*/
#include "output_fields/output_vtu_foreach.h"
#include "streamfunction.h"
void backup_fields (vector u, int nf)
{
  FILE * fp ;
  char name[80], stamp[1024];
  sprintf(stamp, "%6.6i", nf);

  scalar omega[];
	vorticity (u, omega);
	boundary ({psi, omega});
	poisson (psi, omega);

  sprintf(name, "vortex_fixed_%3.3d.vtu", nf);
  fp = fopen(name, "w");
  output_vtu_bin_foreach ((scalar *) {psi}, (vector *) {u}, N, fp, false);
  fclose (fp);

  sprintf(name, "vortex_fixed_%6.6d.dump", nf);
  dump (name);
}

event image (t += 5.0; t<=EndTime) {
  static int nf=0;
  backup_fields(u,nf);
  nf++;

  scalar omega[];
	vorticity (u, omega);
	boundary ({psi, omega});
	poisson (psi, omega);

  stats s0 = statsf (omega);
  fprintf (stderr, "Final condition: %.9g %.9g %.9g \n", s0.min, s0.max, s0.sum);

  scalar * list = {omega, psi, p, u};
  for (scalar s in list){
    view(fov=25);
    squares (s.name, linear = false);
    isoline (s.name, n=11);
    box();

    char name[80];
    sprintf(name, "vortex_rotating_%s.png", s.name);
    save (name);
  }
}

event wrapup (t = end) {
	backup_fields(u,999);
}
