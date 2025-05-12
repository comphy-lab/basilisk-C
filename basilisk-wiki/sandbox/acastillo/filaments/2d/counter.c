/**
# The dynamics of a viscous vortex dipole

The evolution of a viscous vortex dipole is simulated as in Delbende & Rossi,
Phys. Fluids 21, 073605 (2009).

![Vorticity field with streamfunction isolines](counter/omega.png)

## Initial conditions and simulation parameters
The initial condition used in the direct
numerical simulation (DNS) is chosen to be a superposition of two Lamb–Oseen
vortices with circulation $\pm \Gamma$ located at $(\pm b/2, 0)$ in the
($x,y$)-plane. Velocity field is initialized from the vorticity field as in
[vortex.c](). Additionally, we follow the vortex pair moving with velocity
$\Gamma/2\pi$.
*/

#include "navier-stokes/centered.h"
#define MAXLEVEL 11
#define MINLEVEL 4
#define RAD (sqrt(sq(x) + sq(y)))
int main(){
  L0 = 16;
  X0 = Y0 = -L0/2;
  init_grid (1 << MINLEVEL);
  periodic(top);
  periodic(left);
  double reynolds= 8000;
  const face vector muc[] = {1./reynolds,1./reynolds};
  mu = muc;
  run();
}

event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){1e-5,1e-5}, MAXLEVEL, MINLEVEL);
}

event init (t = 0){

  refine  (level < MAXLEVEL);

  scalar psi[], omega[];

  double Gamma=1.0, ad = 0.1, bd=0.5;
  foreach() {
    omega[]  =  Gamma/(pi*sq(ad)) * exp(-(sq(x - bd)/sq(ad) + sq(y)/sq(ad)));
    omega[] -=  Gamma/(pi*sq(ad)) * exp(-(sq(x + bd)/sq(ad) + sq(y)/sq(ad)));
    psi[] = 0.;
  }
  boundary ({psi,omega});

  poisson (psi, omega);
  coord f = {-1.,1.};
  foreach(){
    foreach_dimension(){
      u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
    }
    u.y[] += Gamma/(2*pi);
  }
  boundary ((scalar *){u});

  FILE * fp = fopen("vortex.asc", "w");
  fputs ("[1]t\t [2]j\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]M20\t [7]M02\t [8]M11\t [9]a\t [10]b\t [11]c\t [12]e\t [13]maxvor \n", fp);
  fclose(fp);
}

/**
## Ellipticity of the vortex dipole
The goal is to evaluate several quantities that are pertinent to describe the
flow field such as circulation, sizes, and ellipticity. When a vorticity field
possesses an elliptical shape, these quantities are sufficient to define its
geometry.
*/

#include "ellipticity.h"
event logfile (t += 0.05) {
  scalar omega[];
  vorticity (u, omega);
  stats s = statsf (omega);

  scalar m[];
  foreach()
    m[] = abs(omega[]) > s.stddev/100;

  FILE * fp = fopen("vortex.asc", "a");
  vorticity_moments(omega, m, fp);
  fclose(fp);
}

/**
~~~pythonplot Size $a$ and ellipticity $e$ as function of $t$
import numpy as np
import matplotlib.pyplot as plt

file1 = 'vortex.asc'
data = np.loadtxt(file1, skiprows=1)

ix1 = np.where(data[:,3] > 0)[0];
ix2 = np.where(data[:,3] < 0)[0];

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
          'figure.figsize': (4.5*2, 4.5*4/5*2),
          'contour.negative_linestyle':'dashed'}
plt.rcParams.update(params)

f, axes = plt.subplots(2,2, sharex=False, sharey=False)
ax= np.ravel(axes)

ax[0].plot(vortex1[:,0], vortex1[:,12])
ax[0].plot(vortex2[:,0], vortex2[:,12], linestyle='dashed')
ax[0].set_yscale('log')
ax[0].set_xlabel('$t$')
ax[0].set_ylabel('$\omega(t)$')
ax[0].set_xlim([0,175])
ax[0].set_ylim([1e0,1e2])

ax[1].plot(vortex1[:,0], vortex1[:,2])
ax[1].plot(vortex2[:,0], vortex2[:,2], linestyle='dashed')
ax[1].set_xlabel('$t$')
ax[1].set_ylabel('$\Gamma(t)$')
ax[1].set_xlim([0,175])


ax[2].plot(vortex1[:,0], vortex1[:,8])
ax[2].plot(vortex2[:,0], vortex2[:,8], linestyle='dashed')
ax[2].set_xlabel('$t$')
ax[2].set_ylabel('$a(t)$')
ax[2].set_xlim([0,175])

ax[3].plot(vortex1[:,0], vortex1[:,11])
ax[3].plot(vortex2[:,0], vortex2[:,11], linestyle='dashed')
ax[3].set_xlabel('$t$')
ax[3].set_ylabel('$e(t)$')
ax[3].set_xlim([0,10])
ax[3].set_ylim([0,0.5])

plt.tight_layout()
plt.savefig('series.png')
~~~
*/

/**
# Additional outputs

The vorticity field looks like this

![Vorticity field](counter/omega.mp4)

*/

#include "view.h"
event movie (t += 0.02) {
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});

  squares ("omega", linear = false);
  box();
  save ("omega.mp4");
}

#include "../../output_fields/output_vtu_foreach.h"
event snapshots (t += 10.0) {

  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});

  static int nf = 0;
  char name[80];
  sprintf(name, "counter_%3.3d", nf);
  output_vtu ((scalar *) {omega}, (vector *) {u}, name);

  squares ("omega", linear = false);
  box();
  sprintf(name, "omega_%3.3d.png", nf);
  save (name);
  nf++;
}

event output (t = 200)
  dump();
