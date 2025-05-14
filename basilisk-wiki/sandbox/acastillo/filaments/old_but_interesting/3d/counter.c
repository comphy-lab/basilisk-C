/**
# The dynamics of a viscous vortex dipole

The evolution of a viscous vortex dipole is simulated as in [../2d/counter.c]()
but in 3D.

![Vorticity field](counter/omega.png)

## Initial conditions and simulation parameters
The initial condition used in the direct numerical simulation (DNS) is chosen to
be a superposition of two Lamb–Oseen vortices with circulation $\pm \Gamma$
located at $(\pm b/2, 0)$ in the ($x,y$)-plane. Velocity field is initialized
from the vorticity field which requires solving a Poisson problem for each
component. Additionally, we follow the vortex pair moving with
velocity $\Gamma/2\pi$.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define MINLEVEL 4
#define RAD (sqrt(sq(x) + sq(y)))

int main(){
  L0 = 16;
  X0 = Y0 = Z0 = -L0/2;
  init_grid (1 << MINLEVEL);
  periodic(top);
  periodic(left);
  periodic(front);
  double reynolds= 2500;
  const face vector muc[] = {1./reynolds,1./reynolds};
  mu = muc;
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-5, 1e-5, 1e-5}, MAXLEVEL, MINLEVEL);
}

event init (t = 0){
  refine  (RAD < 1 && level < MAXLEVEL);

  vector psi[], omega[];
  double ad = 0.1, bd=0.5, Gamma= 1.0;

  foreach_cell() {
    omega.x[] = 0.;
    omega.y[] = 0.;
    omega.z[] =  Gamma/(pi*sq(ad)) * exp(-(sq(x - bd)/sq(ad) + sq(y)/sq(ad)));
    omega.z[] -= Gamma/(pi*sq(ad)) * exp(-(sq(x + bd)/sq(ad) + sq(y)/sq(ad)));
  }

  foreach()
    foreach_dimension()
      psi.x[] = 0.;
  boundary ((scalar*){psi,omega});
  poisson (psi.x, omega.x);
  poisson (psi.y, omega.y);
  poisson (psi.z, omega.z);

  foreach(){
    u.x[] = -((psi.z[0,1,0] - psi.z[0,-1,0]) - (psi.y[0,0,1] - psi.y[0,0,-1]))/(2.*Delta);
    u.y[] = -((psi.x[0,0,1] - psi.x[0,0,-1]) - (psi.z[1,0,0] - psi.z[-1,0,0]))/(2.*Delta);
    u.z[] = -((psi.y[1,0,0] - psi.y[-1,0,0]) - (psi.x[0,1,0] - psi.x[0,-1,0]))/(2.*Delta);

    u.y[] -= Gamma/(2*pi);
  }
  boundary ((scalar *){u});

  FILE * fp = fopen("vortex.asc", "w");
  fputs ("[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e\t [14]maxvor \n", fp);
  fclose(fp);
}

/**
## Ellipticity of the vortex dipole
The goal is to evaluate several quantities that are pertinent to describe the
flow field such as circulation, sizes, and ellipticity. When a vorticity field
possesses an elliptical shape, these quantities are sufficient to define its
geometry. In this case, we measure the ellipticity from the intersection of
each vortex with the plane $z=0$.
*/

#include "ellipticity.h"
event logfile (t += 0.1) {

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  scalar m[];
  foreach()
    m[] = abs(omega.z[]) > s.stddev/100;

  FILE * fp = fopen("vortex.asc", "a");
  vorticity_moments_plane(omega.z, m, fp, (coord){0,0,1}, 0.);
  fclose(fp);
}

/**
~~~pythonplot Vertical position $\\mu_y$ (in the fixed reference frame), core size $a$ and ellipticity $e$ as function of $t$.
Vorticity threshold needs to be tweaked...
import numpy as np
import matplotlib.pyplot as plt

file1 = 'vortex.asc'
data = np.loadtxt(file1, skiprows=1)

ix1 = np.where(data[:,2] > 0)[0];
ix2 = np.where(data[:,2] < 0)[0];

vortex1 = data[ix1,:]
vortex2 = data[ix2,:]

params = {'backend': 'ps',
          'axes.labelsize': 12,
          'font.size': 12,
          'legend.fontsize': 12,
          'font.size': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family': 'serif',
          'figure.figsize': (4.5*3, 4.5*4/5),
          'contour.negative_linestyle':'dashed'}
plt.rcParams.update(params)

f, axes = plt.subplots(1,3, sharex=False, sharey=False)
ax= np.ravel(axes)

ax[0].plot(vortex1[:,0], vortex1[:,4] - vortex1[:,0]/(2*np.pi), vortex2[:,0], vortex2[:,4] - vortex2[:,0]/(2*np.pi), '--')
ax[0].set_xlabel('$t$')
ax[0].set_ylabel('$\\mu_y(t)$')
ax[0].set_xlim([0,20])
ax[0].legend(['Vortex 1', 'Vortex 2'])

ax[1].plot(vortex1[:,0], vortex1[:,8], vortex2[:,0], vortex2[:,8], '--')
ax[1].set_xlabel('$t$')
ax[1].set_ylabel('$a(t)$')
ax[1].set_xlim([0,20])
ax[1].legend(['Vortex 1', 'Vortex 2'])

ax[2].plot(vortex1[:,0], vortex1[:,11], vortex2[:,0], vortex2[:,11], '--')
ax[2].legend(['Vortex 1', 'Vortex 2'])
ax[2].set_xlabel('$t$')
ax[2].set_ylabel('$e(t)$')
ax[2].set_xlim([0,20])

plt.tight_layout()
plt.savefig('series.png')
~~~
*/

/**
# Additional outputs

The vorticity and $\lambda_2$ field look something like this

![Vorticity field](counter/omega.mp4) ![Vorticity field](counter/lambda2.mp4)

*/

#include "lambda2.h"
#include "view.h"
event movie (t += 0.2) {
  scalar l2[];
  lambda2 (u, l2);
  stats f = statsf (l2);

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  view (theta = pi/6, phi = pi/6, fov=35);
  isosurface ("l2", -f.stddev);
  squares ("l2", linear = false, alpha=-L0/2);
  box();
  save ("lambda2.mp4");

  isosurface ("omega.z",  s.stddev);
  isosurface ("omega.z", -s.stddev);
  squares ("omega.z", linear = false, alpha=-L0/2);
  box();
  save ("omega.mp4");
}

#include "../../output_fields/output_vtu_foreach.h"
event snapshots (t += 50.0) {

  scalar l2[];
  lambda2 (u, l2);
  stats f = statsf (l2);

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  static int nf = 0;
  char name[80];
  sprintf(name, "counter_%3.3d", nf);
  output_vtu ((scalar *) {l2}, (vector *) {u, omega}, name);
  nf++;
}

event slices (t += 10.0) {

  scalar l2[];
  lambda2 (u, l2);
  stats f = statsf (l2);

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  static int ns = 0;
  char name[80];
  sprintf(name, "counter_x0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){1,0,0}, 0.0);

  sprintf(name, "counter_y0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,1,0}, 0.0);

  sprintf(name, "counter_z0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,0,1}, 0.0);

  view (theta = pi/6, phi = pi/6, fov=35);
  isosurface ("l2", -f.stddev);
  squares ("omega.z", linear = false, alpha=-L0/2);
  box();
  sprintf(name, "omegaz_%3.3d.png", ns);
  save (name);
  ns++;
}


event output (i = 5)
  dump();
