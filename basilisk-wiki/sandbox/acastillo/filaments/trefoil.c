/**
# Simulated Navier-Stokes trefoil

The evolution of a trefoil vortex knot is simulated as in Kerr (2015). The goal
of the physical space initialisation is to map an analytically defined trefoil
vortex onto an Eulerian (static) numerical mesh. The trefoil trajectory
discussed in here is defined by:
$$
x = \sin{(t)} + 2 \sin{(2t)}
$$
$$
y = \cos{(t)} - 2 \cos{(2t)}
$$
$$
z = -3 \sin{(3t)}
$$
where $t$ ranges between 0 and $2\pi$.

For this example, we use the (compressible) Navier-Stokes equations inside a triple periodic box.

![Space-curve of a trefoil vortex relative to the domain size](trefoil/trefoil0.png)

Which results in something like this
![The movie shows a $\lambda_2=0$ isosurface](trefoil/trefoil1.mp4)

![The movie shows a $\lambda_2=0$ isosurface on top of a slice of $\omega_z$](trefoil/trefoil3.mp4)

*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "PointTriangle.h"
#define MINLEVEL 4

int n_seg = 128;
double as = 0.1;

int main() {
  L0 = 16;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.01;
  N = 1<<MINLEVEL;
  periodic(left);
  periodic(top);
  periodic(front);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-5, 1e-5, 1e-5}, MAXLEVEL, MINLEVEL);
}

/**
## Initialisation
For the initialisation process, we proceed as follows.

1. The position of the vortex segments is given by a space-curve $c$
parametrized as function of t0. We also require the first and second derivatives
of $c$, as to define a local Serret-Frenet frame $(\vec{t},\vec{n},\vec{b})$.
$$
\vec{t} = \frac{d\vec{c}}{dt_0}/\|\frac{d\vec{c}}{dt_0}\|
$$
$$
\vec{n} = \frac{d^2\vec{c}}{dt_0^2}/\|\frac{d^2\vec{c}}{dt_0^2}\|
$$
$$
\vec{b} = \vec{t} \times \vec{n}
$$

2. Each position $\vec{p}$ is projected into the local Frenet-Serret frame
to obtain a set of local coordinates, such that:
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{t} = 0
$$
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{n} = x_n
$$
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{b} = x_b
$$
which is done through a minization process.

3. We use the local coordinates $(x_n, x_b)$ to define a radial coordinate
$\rho$ required to compute the vorticity of a Lamb-Oseen vortex as
$$
\vec{\omega} = \Gamma/(\pi a^2) \exp(-\rho^2/a^2) \cdot \vec{t}
$$
where $\Gamma$ is the circulation and $a$ the core size.

4. We consider the flow is expressed as a vector potential,
$$
\vec{u} = \nabla \times \vec{\psi}
$$
such that
$$
\nabla^2 \vec{\psi} = -\vec{\omega}
$$
which requires solving a Poisson problem for each vorticity component.

5. We use the vector potential to evaluate the velocity field.
*/
#include "view.h"
#include "filaments.h"
#include "filaments.c"
event init (t = 0) {
  refine  (sqrt(sq(x) + sq(y) + sq(z)) < 5 && level < MAXLEVEL);

  // Step 1
  double delta_t0 = 6*pi/((double)n_seg-1);
  double t0[n_seg];
  coord c[n_seg], dc[n_seg], d2c[n_seg];

  for (int i = 0; i < n_seg; i++){
    t0[i] = delta_t0 * (double)i - 2*pi;

    c[i].x = sin(t0[i]) + 2*sin(2*t0[i]);
    c[i].y = cos(t0[i]) - 2*cos(2*t0[i]);
    c[i].z = -sin(3*t0[i]);

    dc[i].x =  cos(t0[i]) + 4*cos(2*t0[i]);
    dc[i].y = -sin(t0[i]) + 4*sin(2*t0[i]);
    dc[i].z = -3*cos(3*t0[i]);

    d2c[i].x = -sin(t0[i]) - 8*sin(2*t0[i]);
    d2c[i].y = -cos(t0[i]) + 8*cos(2*t0[i]);
    d2c[i].z = 9*sin(3*t0[i]);
  }
  view (theta = pi/6, phi = pi/6, fov=40);
  box();
  draw_space_curve(n_seg, c);
  save ("trefoil0.png");

  coord tvec[n_seg], nvec[n_seg], bvec[n_seg];
  for (int i = 0; i < n_seg; i++){
    foreach_dimension(){
      tvec[i].x =  dc[i].x/sqrt(vecdot( dc[i],  dc[i]));
      nvec[i].x = d2c[i].x/sqrt(vecdot(d2c[i], d2c[i]));
    }
    bvec[i] = vecdotproduct(tvec[i], nvec[i]);
  }

  // Steps 2 and 3
  vector omega[];
  foreach(){
    coord pcar = {x,y,z};
    coord val_omega = get_vorticity_filament(pcar, n_seg, as, t0, c, tvec, nvec, bvec, 1);
    foreach_dimension()
      omega.x[] = val_omega.x;
  }

  // Step 4
  vector psi[];
  foreach()
    foreach_dimension()
      psi.x[] = 0.;
  boundary ((scalar*){psi,omega});
  poisson (psi.x, omega.x);
  poisson (psi.y, omega.y);
  poisson (psi.z, omega.z);

  // Step 5
  foreach(){
    u.x[] = ((psi.z[0,1,0] - psi.z[0,-1,0]) - (psi.y[0,0,1] - psi.y[0,0,-1]))/(2.*Delta);
    u.y[] = ((psi.x[0,0,1] - psi.x[0,0,-1]) - (psi.z[1,0,0] - psi.z[-1,0,0]))/(2.*Delta);
    u.z[] = ((psi.y[1,0,0] - psi.y[-1,0,0]) - (psi.x[0,1,0] - psi.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});

  FILE * fp = fopen("trefoil_x0.asc", "w");
  fputs ("[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e\t [14]maxvor \n", fp);
  fclose(fp);

  fp = fopen("trefoil_y0.asc", "w");
  fputs ("[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e\t [14]maxvor \n", fp);
  fclose(fp);

  fp = fopen("trefoil_z0.asc", "w");
  fputs ("[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e\t [14]maxvor \n", fp);
  fclose(fp);
}


/**
## Results
For this example, we track the evolution of the kinetic energy, the viscous
enery dissipation rate, as well as the total helicity
$$
H = \int \vec{u} \cdot \vec{\omega} dV
$$
*/
#include "lambda2.h"
#include "3d/ellipticity.h"
event logfile (t += 0.1) {
  scalar l2[];
  lambda2 (u, l2);

  vector omega[];
  vorticity3d(u, omega);

  scalar m[];
  foreach()
    m[] = l2[] < 0;

  FILE * fp = fopen("trefoil_x0.asc", "a");
  vorticity_moments_plane(omega.z, m, fp, (coord){1,0,0}, 0.);
  fclose(fp);

  fp = fopen("trefoil_y0.asc", "a");
  vorticity_moments_plane(omega.z, m, fp, (coord){0,1,0}, 0.);
  fclose(fp);

  fp = fopen("trefoil_z0.asc", "a");
  vorticity_moments_plane(omega.z, m, fp, (coord){0,0,1}, 0.);
  fclose(fp);
}

event movie (t += 0.2) {
  scalar l2[];
  lambda2 (u, l2);

  vector omega[];
  vorticity3d(u, omega);

  view (theta = pi/6, phi = pi/6, fov=35);
  isosurface ("l2", 0);
  squares ("omega.z", linear = false, alpha=-L0/2);
  box();
  save ("lambda2.mp4");
}


#include "../output_fields/output_vtu_foreach.h"
event snapshots (t += 20.0) {

  scalar l2[];
  lambda2 (u, l2);
  stats f = statsf (l2);

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  static int nf = 0;
  char name[80];
  sprintf(name, "trefoil_%3.3d", nf);
  output_vtu ((scalar *) {l2}, (vector *) {u, omega}, name);
  nf++;
}

event slices (t += 1.0) {
  scalar l2[];
  lambda2 (u, l2);

  vector omega[];
  vorticity3d(u, omega);

  static int ns = 0;
  char name[80];
  sprintf(name, "trefoil_x0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){1,0,0}, 0.0);

  sprintf(name, "trefoil_y0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,1,0}, 0.0);

  sprintf(name, "trefoil_z0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,0,1}, 0.0);

  view (theta = pi/6, phi = pi/6, fov=35);
  isosurface ("l2", 0);
  squares ("omega.z", linear = false, alpha=-L0/2);
  box();
  sprintf(name, "omegaz_%3.3d.png", ns);
  save (name);
  ns++;

}

event stop (i = 5)
  dump();
