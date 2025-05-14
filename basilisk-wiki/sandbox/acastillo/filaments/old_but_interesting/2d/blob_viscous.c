/**
# The dynamics of a blob of vorticity

## Initial conditions and simulation parameters

The initial conditions are read from a gnuplot-compatible binary file, obtained
from `vorticity_18deg.mat`. We preprocess this field using `blob.m`, which
essentially extends the domain, centers the vorticity patches and applies
a Hanning window to ensure that vorticity vanishes smoothly away from the
center. The resulting vorticity field looks something like this

![Initial vorticity field](blob_viscous/omega.mp4)

The goal is to indentify a secondary growth mechanism due to the initial
vortiity distribution. Initially, we neglect the viscous diffusion.
*/


#include "navier-stokes/centered.h"
#define MAXLEVEL 11
#define MINLEVEL 4

int main(){
  L0 = 20;
  X0 = Y0 = -L0/2;
  init_grid (256);
  periodic(top);
  periodic(left);
  const face vector muc[] = {1./100, 1./100};
  mu = muc;
  run();
}

event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){1e-4,1e-4}, MAXLEVEL, MINLEVEL);
}

#include "view.h"
#include "../../input_fields/auxiliar_input.h"
#define RAD (sqrt(sq(x) + sq(y)))
event init (t = 0){
  refine  (RAD < 10 && level < MAXLEVEL-3);
  refine  (RAD <  5 && level < MAXLEVEL-2);

  scalar psi[], omega[];
  FILE * fp = fopen("../blob/omega.bin", "r");
	if (!fp) printf("Binary file not found");
	input_matrix(omega,fp,N,X0,Y0,L0);
	fclose (fp);

  foreach() {
    psi[] = 0.;
  }
  boundary ({psi,omega});

  poisson (psi, omega);
  coord f = {-1.,1.};
  foreach()
    foreach_dimension()
      u.x[] = -f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});

  fp = fopen("vortex.asc", "w");
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
event logfile (t += 0.01) {
  scalar omega[];
  vorticity (u, omega);
  stats s = statsf (omega);

  scalar m[];
  foreach()
    m[] = abs(omega[]) > s.stddev/10;

  FILE * fp = fopen("vortex.asc", "a");
  vorticity_moments(omega, m, fp);
  fclose(fp);
}

/**

~~~pythonplot Circulation $\Gamma$, size $a$ and ellipticity $e$ as function of $t$
import numpy as np
import matplotlib.pyplot as plt

file1 = 'vortex.asc'
data = np.loadtxt(file1, skiprows=1)

ix1 = np.where(data[:,3] > 130)[0];
ix2 = np.where((data[:,3] <=130) & (data[:,3] > 100))[0];
ix3 = np.where(data[:,3] <= 100)[0];



vortex1 = data[ix1,:]
vortex2 = data[ix2,:]
vortex3 = data[ix3,:]


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

ax[0].plot(vortex1[:,0], vortex1[:,3])
ax[0].plot(vortex2[:,0], vortex2[:,3])
ax[0].scatter(vortex3[:,0], vortex3[:,3], color='C3', s=1)
ax[0].set_xlim([0,10])
ax[0].set_xlabel('$t$')
ax[0].set_ylabel('$\Gamma(t)$')
ax[0].legend(['Vortex1+', 'Vortex2+', 'Other'])

ax[1].plot(vortex1[:,0], vortex1[:,9])
ax[1].plot(vortex2[:,0], vortex2[:,9])
ax[1].set_xlim([0,10])
ax[1].set_ylim([0,1])
ax[1].set_xlabel('$t$')
ax[1].set_ylabel('$a(t)$')
ax[1].legend(['Vortex1+', 'Vortex2+'])

ax[2].plot(vortex1[:,0], np.sqrt(1 - (vortex1[:,11]**2/vortex1[:,10]**2)))
ax[2].plot(vortex2[:,0], np.sqrt(1 - (vortex2[:,11]**2/vortex2[:,10]**2)))
ax[2].set_xlim([0,10])
ax[2].set_ylim([0,1])
ax[2].set_xlabel('$t$')
ax[2].set_ylabel('$e(t)$')
ax[2].legend(['Vortex1+', 'Vortex2+'])

plt.tight_layout()
plt.savefig('series.png')
~~~

*/

/**
# Additional outputs
*/

event movie (t += 0.02) {
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  stats s = statsf (omega);

  scalar m[];
  foreach()
    m[] = abs(omega[]) > s.stddev/10;
  int n = tag (m);

  squares ("omega", linear = false);
  box();
  save ("omega.mp4");

  squares ("m", min=0, max=n);
  box();
  save ("tag.mp4");
}

#include "../../output_fields/output_vtu_foreach.h"
event snapshots (t += 0.5) {

  scalar omega[];
  vorticity (u, omega);
  poisson (omega);
  stats s = statsf (omega);

  scalar m[];
  foreach()
    m[] = abs(omega[]) > s.stddev/10;
  int n = tag (m);

  static int nf = 0;
  char name[80];
  sprintf(name, "blob_%3.3d", nf);
  output_vtu ((scalar *) {omega, m}, (vector *) {u}, name);

  squares ("omega", linear = false);
  box();
  sprintf(name, "omega_%3.3d.png", nf);
  save (name);

  squares ("m", min=0, max=n);
  box();
  sprintf(name, "tag_%3.3d.png", nf);
  save (name);
  nf++;
}

event output (t = 5.0)
  dump();
