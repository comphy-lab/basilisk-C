/**
# The dynamics of a viscous vortex dipole

The evolution of a viscous vortex dipole is simulated as in
[Le Dizes and Verga., 2002](#ledizes2002)).
.

![Vorticity field with streamfunction isolines](corotating/omega.png)

## Initial conditions and simulation parameters
The initial condition used in the direct
numerical simulation (DNS) is chosen to be a superposition of two Lamb–Oseen
vortices with circulation $\Gamma$ located at $(\pm b/2, 0)$ in the
($x,y$)-plane. Velocity field is initialized from the vorticity field as in
[vortex.c]().
*/

#include "navier-stokes/centered.h"
#define MAXLEVEL 11
#define MINLEVEL 4
#define RAD (sqrt(sq(x) + sq(y)))
double reynolds = _REYNOLDS;
int main(){
  L0 = 16;
  X0 = Y0 = -L0/2;
  init_grid (1 << MINLEVEL);
  periodic(top);
  periodic(left);
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
    omega[] +=  Gamma/(pi*sq(ad)) * exp(-(sq(x + bd)/sq(ad) + sq(y)/sq(ad)));
    psi[] = 0.;
  }
  boundary ({psi,omega});

  poisson (psi, omega);
  coord f = {-1.,1.};
  foreach()
    foreach_dimension()
      u.x[] = -f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});

  char name[80];
  sprintf(name, "vortex_re%g.asc", reynolds);
  FILE * fp = fopen(name, "w");
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

  char name[80];
  sprintf(name, "vortex_re%g.asc", reynolds);
  FILE * fp = fopen(name, "a");
  vorticity_moments(omega, m, fp);
  fclose(fp);
}

/**
~~~pythonplot Size $a$ and ellipticity $e$ as function of $t$. Vorticity threshold must be tweaked.
import numpy as np
import matplotlib.pyplot as plt

# Read the time series and separate each vortex using the minimal distance between
# consecutive time steps
file1 = 'vortex_re8000.asc'
data = np.loadtxt(file1, skiprows=1)
data = data[~np.isnan(data).any(axis=1)]
dt = 0.05
t = 0.0
from_index = np.zeros([np.shape(data)[0],], int)
to_index = np.zeros([np.shape(data)[0],], int)

for it, t in enumerate(np.arange(t, t + 40, dt)):
    index1 = np.where( np.abs(data[:,0] - t) < 1e-8 )[0]
    for i in index1:
        mindist = 1e30
        index2 = np.where( np.abs(data[:,0] - (t+dt)) < 1e-8 )[0]
        for j in index2:
            dist = np.sqrt((data[i,3]-data[j,3])**2 + (data[i,4]-data[j,4])**2 + (data[i,5]-data[j,5])**2)

            if (dist < mindist):
                mindist = dist
                from_index[i] = i
                to_index[i]   = j

index = 0
ix1 = 0
while 1 :
    current = to_index[index]
    index = np.where(from_index[:] == current)[0]
    if index.size == 0:
        break
    ix1 = np.append(ix1, index[0])

index = 1
ix2 = 1
while 1 :
    current = to_index[index]
    index = np.where(from_index[:] == current)[0]
    if index.size == 0:
        break
    ix2 = np.append(ix2, index[0])

vortex1a = data[ix1,:]
vortex2a = data[ix2,:]


# Read the time series and separate each vortex using the minimal distance between
# consecutive time steps
file1 = 'vortex_re2000.asc'
data = np.loadtxt(file1, skiprows=1)
data = data[~np.isnan(data).any(axis=1)]
dt = 0.05
t = 0.0
from_index = np.zeros([np.shape(data)[0],], int)
to_index = np.zeros([np.shape(data)[0],], int)

for it, t in enumerate(np.arange(t, t + 40, dt)):
    index1 = np.where( np.abs(data[:,0] - t) < 1e-8 )[0]
    for i in index1:
        mindist = 1e30
        index2 = np.where( np.abs(data[:,0] - (t+dt)) < 1e-8 )[0]
        for j in index2:
            dist = np.sqrt((data[i,3]-data[j,3])**2 + (data[i,4]-data[j,4])**2 + (data[i,5]-data[j,5])**2)

            if (dist < mindist):
                mindist = dist
                from_index[i] = i
                to_index[i]   = j

index = 0
ix1 = 0
while 1 :
    current = to_index[index]
    index = np.where(from_index[:] == current)[0]
    if index.size == 0:
        break
    ix1 = np.append(ix1, index[0])

index = 1
ix2 = 1
while 1 :
    current = to_index[index]
    index = np.where(from_index[:] == current)[0]
    if index.size == 0:
        break
    ix2 = np.append(ix2, index[0])

vortex1b = data[ix1,:]
vortex2b = data[ix2,:]


def correct_angles(y):

    i = np.where(y[1:]-y[0:-1] > np.pi)[0]
    for k in np.arange(0,np.size(i)):
        y[i[k]+1:] = y[i[k]+1:] - 2*np.pi

    i = np.where(y[1:]-y[0:-1] < -np.pi)[0]
    for k in np.arange(0,np.size(i)):
        y[i[k]+1:] = y[i[k]+1:] + 2*np.pi

    return y


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

ax[0].plot(vortex1a[:,0], vortex1a[:,8], label='Vortex 1, $Re=8000$')
ax[0].plot(vortex2a[:,0], vortex2a[:,8], label='Vortex 2, $Re=8000$', linestyle='dashed')
ax[0].plot(vortex1b[:,0], vortex1b[:,8], label='Vortex 1, $Re=2000$')
ax[0].plot(vortex2b[:,0], vortex2b[:,8], label='Vortex 2, $Re=2000$', linestyle='dashed')
ax[0].set_xlim([0,23])
ax[0].set_ylim([0,0.25])
ax[0].legend(fontsize=10)
ax[0].set_ylabel('$a(t)$')
ax[0].set_xlabel('$t$')

ba = np.sqrt((vortex1a[:,3]-vortex2a[:,3])**2 + (vortex1a[:,4]-vortex2a[:,4])**2)
bb = np.sqrt((vortex1b[:,3]-vortex2b[:,3])**2 + (vortex1b[:,4]-vortex2b[:,4])**2)

ax[1].plot(vortex1a[:,0], ba, label='$Re=8000$')
ax[1].plot(vortex1b[:,0], bb, label='$Re=2000$', color='C2')
ax[1].set_ylim([0,1.5])
ax[1].set_xlim([0,23])
ax[1].legend(fontsize=10)
ax[1].set_ylabel('$b(t)$')
ax[1].set_xlabel('$b$')

phi1a = correct_angles(np.arctan2(vortex1a[:,4], vortex1a[:,3]))
phi1b = correct_angles(np.arctan2(vortex1b[:,4], vortex1b[:,3]))
phi2a = correct_angles(np.arctan2(vortex2a[:,4], vortex2a[:,3]))
phi2b = correct_angles(np.arctan2(vortex2b[:,4], vortex2b[:,3]))

omega1a = (phi1a[1:]-phi1a[0:-1])/(t[1:]-t[0:-1])
omega2a = (phi2a[1:]-phi2a[0:-1])/(t[1:]-t[0:-1])
omega1b = (phi1b[1:]-phi1b[0:-1])/(t[1:]-t[0:-1])
omega2b = (phi2b[1:]-phi2b[0:-1])/(t[1:]-t[0:-1])

t = vortex1a[:,0]
ax[2].plot(t[1:], omega1a, label='Vortex 1, $Re=8000$')
ax[2].plot(t[1:], omega2a, label='Vortex 2, $Re=8000$', linestyle='dashed')
ax[2].plot(t[1:], omega1b, label='Vortex 1, $Re=2000$')
ax[2].plot(t[1:], omega2b, label='Vortex 2, $Re=2000$', linestyle='dashed')
ax[2].set_xlim([0,23])
ax[2].set_ylim([0,0.5])
ax[2].legend(fontsize=10)
ax[2].set_ylabel('$\dot{\\theta}(t)$')
ax[2].set_xlabel('$t$')

ax[3].plot(vortex1a[:,0], vortex1a[:,11], label='Vortex 1, $Re=8000$')
ax[3].plot(vortex2a[:,0], vortex2a[:,11], label='Vortex 2, $Re=8000$', linestyle='dashed')
ax[3].plot(vortex1b[:,0], vortex1b[:,11], label='Vortex 1, $Re=2000$')
ax[3].plot(vortex2b[:,0], vortex2b[:,11], label='Vortex 2, $Re=2000$', linestyle='dashed')
ax[3].set_xlim([0,10])
ax[3].set_ylim([0,0.6])
ax[3].legend(fontsize=10)
ax[3].set_ylabel('$e(t)$')
ax[3].set_xlabel('$t$')

plt.tight_layout()
plt.savefig('series.png')

~~~
*/

#include "view.h"
#include "../../output_fields/output_vtu_foreach.h"
event snapshots (t += sq(pi)/4.0) {

  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});

  static int nf = 0;
  char name[80];
  sprintf(name, "corotating_re%g_%3.3d", reynolds, nf);
  output_vtu ((scalar *) {omega}, (vector *) {u}, name);

  squares ("omega", linear = false);
  box();
  sprintf(name, "omega_re%g_%3.3d.png", reynolds, nf);
  save (name);
  nf++;
}

event output (t = 4*sq(pi))
  dump();

/**
# Additional outputs

The vorticity field looks like this

![Vorticity field](corotating/omega.mp4)

*/


event movie (t += 0.2) {
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});

  squares ("omega", linear = false);
  box();

  char name[80];
  sprintf(name, "omega_re%g.mp4", reynolds);
  save (name);
}

/**
## References

~~~bib
@article{ledizes2002,
  title={Viscous interactions of two co-rotating vortices before merging},
  author={Le Dizes, St{\'e}phane and Verga, Alberto},
  journal={Journal of Fluid Mechanics},
  volume={467},
  number={1},
  pages={389--410},
  year={2002},
  publisher={Cambridge; New York: Cambridge University Press, 1956-}
}
~~~
*/
