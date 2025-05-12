/**
# Simulated Navier-Stokes generalized helical vortex pair

The evolution of a generalized helical vortex pair is simulated as in
[helical.c](helical.c), but the trajectories of the vortices are read from an
external file. Here, we consider consider first a simple configuration
composed of two co-axial helices of opposite circulation and same core size a.
The internal and external helices are defined by their radii $R_{int}$ and
$R_{ext}$ and pitches $h_{int}$ and $h_{ext}$, as illustrated in the figure

![Configuration of two undeformed co-axial helical vortices, taken from
E. Durán Venegas and S. Le Dizès (2019)](helical_pair/duran1.png)

The goal of the physical space initialisation is to map an analytically defined
trefoil vortex onto an Eulerian (static) numerical mesh. The helical trajectories
discussed in here have no analytical form, but are obtained from a previous
calculation as in E. Durán Venegas and S. Le Dizès (2019):
With each solution is associated a moving frame where the solution is steady.
For a given vortex parametrized by its radial position $r(z)$ and angular
position $\phi(z)$ as a function of the $z$ coordinate, the condition of
steadiness reads as
$$
\frac{dr}{dz} = \frac{V_r^{ind}}{V_z^{ind} - V_z^F}
$$
$$
\frac{d\phi}{dz} = \frac{\Omega^{ind}-\Omega^F}{V_z^{ind} - V_z^F}
$$
where $\vec{V}^{ind}$ indicates the induced velocity, while $\Omega^F$ and
$V_z^F$ are the frame velocties. Solutions are obtained using a
filamentary code developed by E. Durán Venegas during this thesis,
stored inside `save_R0.5_h1_alpha1.25_epsilon0.03c.mat`. These results are
used to generate `hext.asc` and `hint.asc`, see `preprocess.py`.

For this example, we use the (compressible) Navier-Stokes equations inside a
periodic box corresponding to $R_{int}/R_{ext}=0.5$, $h_{ext}/R_{ext}=1.0$,
$h_{int}/h_{ext}=1.25$, and $a/R_{ext}=0.03$.

![Space-curve of a helical vortex relative to the domain size](helical_pair/helical_pair0.png)

Which results in something like this
![The movie shows a $\lambda_2=0$ isosurface](helical_pair/helical_pair1.mp4)

![The movie shows a $\lambda_2=0$ isosurface on top of a slice of $u_z$](helical_pair/helical_pair2.mp4)

*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "PointTriangle.h"

int n_seg;
int minlevel = 5;
int maxlevel = 10;
double as;
double Hs1, Rs1, n_turns1;
double Hs2, Rs2, n_turns2;

int main() {
  L0 = 5;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.01;
  N = 1<<minlevel;

  periodic(back);
  run();
}

event adapt (i++)
  adapt_wavelet ((scalar*){u}, (double[]){0.001, 0.001, 0.001}, maxlevel, minlevel);

/**
## Initialisation
For the initialisation process, we proceed as follows.

1. The position of the vortex segments is given by a space-curve $c$
parametrized as function of t0.

2. We require the first and second derivatives
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

3. Each position $\vec{p}$ is projected into the local Frenet-Serret frame
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

4. We use the local coordinates $(x_n, x_b)$ to define a radial coordinate
$\rho$ required to compute the vorticity of a Lamb-Oseen vortex as
$$
\vec{\omega} = \Gamma/(\pi a^2) \exp(-\rho^2/a^2) \cdot \vec{t}
$$
where $\Gamma$ is the circulation and $a$ the core size.

5. We consider the flow is expressed as a vector potential,
$$
\vec{u} = \nabla \times \vec{\psi}
$$
such that
$$
\nabla^2 \vec{\psi} = -\vec{\omega}
$$
which requires solving a Poisson problem for each vorticity component.

6. We use the vector potential to evaluate the velocity field.
*/
#include "view.h"
#include "filaments.h"
#include "filaments.c"
#define RAD (sqrt(sq(x) + sq(y)))
event init (t = 0) {

  // Step 1a: Load external vortex
  FILE * fp = fopen("hext.asc", "r");
  fscanf(fp, "%d",  &n_seg);
  fscanf(fp, "%lf", &Rs1);
  fscanf(fp, "%lf", &as);
  fscanf(fp, "%lf", &Hs1);

  n_turns1 = L0/Hs1;
  refine  (RAD < Rs1 + 3*as   && RAD > Rs1 - 3*as   && level < maxlevel - 1);
  refine  (RAD < Rs1 + 1.5*as && RAD > Rs1 - 1.5*as && level < maxlevel);

  double t0[n_seg];
  coord c1[n_seg], dc1[n_seg], d2c1[n_seg];
  coord c2[n_seg], dc2[n_seg], d2c2[n_seg];

  for (int i = 0; i < n_seg; i++)
    fscanf(fp, "%lf %lf %lf %lf", &t0[i], &c1[i].x, &c1[i].y, &c1[i].z);
  fclose (fp);

  // Step 1b: Load internal vortex
  fp = fopen("hint.asc", "r");
  fscanf(fp, "%d",  &n_seg);
  fscanf(fp, "%lf", &Rs2);
  fscanf(fp, "%lf", &as);
  fscanf(fp, "%lf", &Hs2);

  n_turns2 = L0/Hs2;
  refine  (RAD < Rs2 + 3*as   && RAD > Rs2 - 3*as   && level < maxlevel - 1);
  refine  (RAD < Rs2 + 1.5*as && RAD > Rs2 - 1.5*as && level < maxlevel);

  for (int i = 0; i < n_seg; i++)
    fscanf(fp, "%lf %lf %lf %lf", &t0[i], &c2[i].x, &c2[i].y, &c2[i].z);
  fclose (fp);

  view (camera="iso", fov=40);
  box();
  draw_space_curve(n_seg, c1);
  draw_space_curve(n_seg, c2);
  save ("helical_pair0.png");



  // Step 2
  double delta_t0 = t0[1] - t0[0];
  coord x1shift = {0, 0, (3*n_turns1)*Hs1}, dx1shift = {0, 0, 0};
  coord x2shift = {0, 0, (3*n_turns2)*Hs2}, dx2shift = {0, 0, 0};
  fd_derivative(n_seg, delta_t0,  x1shift,  c1,  dc1);
  fd_derivative(n_seg, delta_t0,  x2shift,  c2,  dc2);
  fd_derivative(n_seg, delta_t0, dx1shift, dc1, d2c1);
  fd_derivative(n_seg, delta_t0, dx2shift, dc2, d2c2);

  coord tvec1[n_seg], nvec1[n_seg], bvec1[n_seg];
  coord tvec2[n_seg], nvec2[n_seg], bvec2[n_seg];
  for (int i = 0; i < n_seg; i++){
    foreach_dimension(){
      tvec1[i].x =  dc1[i].x/sqrt(vecdot( dc1[i],  dc1[i]));
      nvec1[i].x = d2c1[i].x/sqrt(vecdot(d2c1[i], d2c1[i]));

      tvec2[i].x =  dc2[i].x/sqrt(vecdot( dc2[i],  dc2[i]));
      nvec2[i].x = d2c2[i].x/sqrt(vecdot(d2c2[i], d2c2[i]));
    }
    bvec1[i] = vecdotproduct(tvec1[i], nvec1[i]);
    bvec2[i] = vecdotproduct(tvec2[i], nvec2[i]);
  }

  // Steps 3 and 4
  vector omega[];
  foreach(){
    coord pcar = {x,y,z};
    coord val_omega1 = get_vorticity_filament(pcar, n_seg, as, t0, c1, tvec1, nvec1, bvec1, 0);
    coord val_omega2 = get_vorticity_filament(pcar, n_seg, as, t0, c2, tvec2, nvec2, bvec2, 0);
    foreach_dimension()
      omega.x[] = val_omega1.x - val_omega2.x;
  }

  // Step 5
  vector psi[];
  foreach()
    foreach_dimension()
      psi.x[] = 0.;
  boundary ((scalar*){psi,omega});
  poisson (psi.x, omega.x);
  poisson (psi.y, omega.y);
  poisson (psi.z, omega.z);

  // Step 6
  foreach(){
    u.x[] = ((psi.z[0,1,0] - psi.z[0,-1,0]) - (psi.y[0,0,1] - psi.y[0,0,-1]))/(2.*Delta);
    u.y[] = ((psi.x[0,0,1] - psi.x[0,0,-1]) - (psi.z[1,0,0] - psi.z[-1,0,0]))/(2.*Delta);
    u.z[] = ((psi.y[1,0,0] - psi.y[-1,0,0]) - (psi.x[0,1,0] - psi.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});

  fp = fopen("helical_pair.asc", "w");
  fclose (fp);

  scalar omega_mag[];
  foreach()
    omega_mag[] = sqrt(sq(omega.x[]) + sq(omega.y[]) + sq(omega.z[]));

  stats f = statsf (omega_mag);

  view (camera="iso", fov=40);
  box();
  isosurface ("omega_mag", f.stddev);
  isosurface ("omega_mag", -f.stddev);
  draw_space_curve(n_seg, c1);
  draw_space_curve(n_seg, c2);
  save ("helical_pair1.png");
}


/**
## Results
For this example, we track the evolution of the kinetic energy, the viscous
enery dissipation rate, as well as the total helicity
$$
H = \int \vec{u} \cdot \vec{\omega} dV
$$
*/
event logs (i++){
  vector du[], dv[], dw[]; // Would be nice to use a tensor here.
	foreach(){
    du.x[] = (u.x[1,0] - u.x[-1,0])/2./Delta;
    dv.x[] = (u.y[1,0] - u.y[-1,0])/2./Delta;
    dw.x[] = (u.z[1,0] - u.z[-1,0])/2./Delta;

    du.y[] = (u.x[0,1] - u.x[0,-1])/2./Delta;
    dv.y[] = (u.y[0,1] - u.y[0,-1])/2./Delta;
    dw.y[] = (u.z[0,1] - u.z[0,-1])/2./Delta;

    du.z[] = (u.x[0,0,1] - u.x[0,0,-1])/2./Delta;
    dv.z[] = (u.y[0,0,1] - u.y[0,0,-1])/2./Delta;
    dw.z[] = (u.z[0,0,1] - u.z[0,0,-1])/2./Delta;
  }

  vector omega[];
  foreach(){
    omega.x[] = dw.y[] - dv.z[];
    omega.y[] = du.z[] - dw.x[];
    omega.z[] = dv.x[] - du.y[];
  }
  boundary ((scalar*){omega});

  double kin=0., eps=0., vor=0., hel=0.;
	foreach(reduction(+:kin), reduction(+:eps), reduction(+:vor), reduction(+:hel)){
		foreach_dimension(){
      kin += dv() * sq(u.x[]);
		  eps += dv() * (sq(du.x[]) + sq(dv.y[]) +  sq(dw.x[]));
      vor += dv() * sq(omega.x[]);
      hel += dv() * u.x[] * omega.x[];
    }
	}

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[1] - u.x[];
    div[] /= Delta;
  }
  stats s0 = statsf (div);

  FILE * fp = fopen("helical_pair.asc", "a");
  fprintf (fp, "%d %f %.9g %.9g %.9g %.9g %.9g %.9g \n", i, t, s0.sum/s0.volume, s0.max, kin, eps, vor, hel);
  fclose (fp);
}

#include "lambda2.h"
event movie (t += 0.005) {
  scalar l2[];
  lambda2 (u, l2);
  foreach()
    l2[] = l2[] > -0.01 ? nodata : l2[];
  stats f = statsf (l2);
  foreach()
    l2[] = l2[] == nodata ? 0 : l2[];
  boundary ({l2});

  view (theta = pi/6, phi = pi/6, fov=30);
  isosurface ("l2", -f.stddev);
  save ("helical_pair1.mp4");

  isosurface ("l2", -f.stddev);
  squares ("u.z", linear = true);
  save ("helical_pair2.mp4");

  vector du[], dv[], dw[];
	foreach(){
    du.x[] = (u.x[1,0] - u.x[-1,0])/2./Delta;
    dv.x[] = (u.y[1,0] - u.y[-1,0])/2./Delta;
    dw.x[] = (u.z[1,0] - u.z[-1,0])/2./Delta;

    du.y[] = (u.x[0,1] - u.x[0,-1])/2./Delta;
    dv.y[] = (u.y[0,1] - u.y[0,-1])/2./Delta;
    dw.y[] = (u.z[0,1] - u.z[0,-1])/2./Delta;

    du.z[] = (u.x[0,0,1] - u.x[0,0,-1])/2./Delta;
    dv.z[] = (u.y[0,0,1] - u.y[0,0,-1])/2./Delta;
    dw.z[] = (u.z[0,0,1] - u.z[0,0,-1])/2./Delta;
  }

  vector omega[];
  foreach(){
    omega.x[] = dw.y[] - dv.z[];
    omega.y[] = du.z[] - dw.x[];
    omega.z[] = dv.x[] - du.y[];
  }
  boundary ((scalar*){omega});

  isosurface ("l2", -f.stddev);
  squares ("omega.z", linear = true);
  save ("helical_pair3.mp4");
}

event stop (t = 5)
  dump();
