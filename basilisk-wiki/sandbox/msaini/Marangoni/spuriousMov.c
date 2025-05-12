#define JACOBI 1

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

#if CLSVOF
# include "two-phase-clsvof.h"
#include "integral.h"
scalar sigmaf[];
#elif HF2D
# include "two-phase-HF.h"
#include "integral.h"
scalar sigmaf[];
#else
#include "vof.h"
#include "tension.h"
scalar f[], * interfaces = {f};
#endif

#define DIAMETER 0.4
#define MU sqrt(DIAMETER/LAPLACE)
#define TMAX (sq(DIAMETER)/MU)

int LEVEL;
double LAPLACE = 120.;
double DC = 0.;
FILE * fp = NULL;

int main() {
  
  TOLERANCE = 1e-6 [*];
  stokes = true;
#if CLSVOF || HF2D
  d.sigmaf = sigmaf;
#else
  f.sigma = 1;
#endif

  periodic (right);
  LEVEL = 6;
  N = 1 << LEVEL;
  for (LAPLACE = 120; LAPLACE <= 12000; LAPLACE *= 10)
    run();
}

event init (i = 0) {

  /**
  We set the constant viscosity field... */

  const face vector muc[] = {MU,MU};
  mu = muc;

  /**
  ... open a new file to store the evolution of the amplitude of
  spurious currents for the various LAPLACE... */

  char name[80];
  sprintf (name, "La-%g", LAPLACE);
  if (fp)
    fclose (fp);
  fp = fopen (name, "w");
  
  /**
  ... and initialise the shape of the interface and the initial volume
  fraction field. */

#if CLSVOF || HF2D
  foreach(){
    d[] = - sqrt (sq(x - L0/2.) + sq(y - L0/2.)) + DIAMETER/2.;
    sigmaf[] = 1.;
    u.x[] = 1.;
  }  
#else
  fraction (f, sq(DIAMETER/2) - sq(x - L0/2.) - sq(y - L0/2.));
  foreach()
    u.x[] = 1.;
#endif
}

event logfile (i++; t <= TMAX)
{
  scalar un[];
  foreach()
    un[] = norm(u) - 1.;
  fprintf (fp, "%g %10.10f %10.10f\n",MU*t/sq(DIAMETER), normf(un).max, normf(un).rms);
}

event error (t = end) {
  
  scalar un[];
  foreach() {
    un[] = norm(u) - 1.;
  }
  
  fprintf (stderr, "%d %10.10f %10.10f %10.10f\n", 
           LEVEL, LAPLACE,normf(un).max*sqrt(DIAMETER),normf(un).rms*sqrt(DIAMETER));
  
}

/**
We use an adaptive mesh with a constant (maximum) resolution along the
interface. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){0}, maxlevel = LEVEL, minlevel = 0);
}
#endif

/**
## Results

~~~pythonplot Spurious currents

import numpy
import sys 
import string
import math
import glob
from pylab import *
from matplotlib.colors import LogNorm, Normalize,ListedColormap
from scipy.interpolate import CubicSpline
from matplotlib.ticker import LogLocator

params = {
    'axes.labelsize': 14,
    'font.size': 12,
    'mathtext.fontset' : 'cm',
    'font.family' : 'sans-serif',
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'xtick.major.size' : 2.,
    'ytick.major.size' : 2.,
    'xtick.minor.size' : 1.,
    'ytick.minor.size' : 1.,
    'text.usetex': False,
#    'figure.figsize': [2.15,1.85],
    'axes.linewidth' : 0.75,
    'figure.subplot.left' : 0.1,
    'figure.subplot.bottom' : 0.18,
    'figure.subplot.right' : 0.96,
    'figure.subplot.top' : 0.98,
    'savefig.dpi' : 300,
}
rcParams.update(params)

La = ['120','1200','12000']

fig, ax = plt.subplots(1, 3, figsize=(24/2.54, 5.5/2.54), gridspec_kw={'wspace': 0.5,'width_ratios': [1., 1, 1.]},sharey=False)

folders=["./","../spuriousMov-clsvof","../spuriousMov-hf2d"]
cmap=["#440154","#3b528b","#21918c"]
figno=np.array([0,1,2])

ax[0].set_title("CSF", loc='center')
ax[1].set_title("CLSVOF", loc='center')
ax[2].set_title("HF2D", loc='center')

for i in figno:
    ax[i].set_yscale('log')
    ax[i].set_xlim(0.,1)
    ax[i].set_xlabel(r'$t/t_\mu$')
    ax[i].set_ylabel(r'$\epsilon_{RMS} (u)$')
    for j in range(len(La)):
        datafile = "%s/La-%s" % (folders[i],La[j])
        data = np.loadtxt(datafile)
        ax[i].plot(data[::,0],abs(data[::,2]),color=cmap[j],label=La[j])
    ax[i].legend(loc=0,frameon=False,fontsize=9)
    
fig.text(0.02, 0.98, '$(a)$', fontsize=14, fontweight='bold')
fig.text(0.34, 0.98, '$(b)$', fontsize=14, fontweight='bold')
fig.text(0.66, 0.98, '$(c)$', fontsize=14, fontweight='bold')

savefig("spurious.pdf",bbox_inches='tight', pad_inches=0.3/2.54)

~~~
## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/spurious.html)
*/
