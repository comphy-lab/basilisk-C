/**
# Capillary wave testcase with integral Height functions method
*/

# include "grid/multigrid.h"
# include "navier-stokes/centered.h"
#if CLSVOF
#if HF2D
# include "two-phase-HF.h"
#else
# include "two-phase-clsvof.h"
#endif
# include "integral.h"
# include "curvature.h"
#else
# include "vof.h"
scalar f[], * interfaces = {f};
#if INTHF
#include "integral-HF.h"
scalar sigma[];
#else
# include "tension.h"
#endif
#endif

#include "test/prosperetti.h"

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

double se = 0; int ne = 0;

int main() {

  size (2. [1]);
  Y0 = -L0/2.;
#if CLSVOF
  const scalar sigma[] = 1.;
  d.sigmaf = sigma;
#else
#if INTHF
  f.sigmaf = sigma;
#else
  f.sigma = 1.;
#endif
#endif
  TOLERANCE = 1e-6 [*];
  const face vector muc[] = {0.0182571749236, 0.0182571749236};
  mu = muc;

  for (N = 16; N <= 128; N *= 2) {
    se = 0, ne = 0;
    run();
  }
}

event init (t = 0) {
  double k = 2., a = 0.01;
#if CLSVOF
  foreach()
    d[] = y - a*cos (k*pi*x);
#else
  fraction (f, y - a*cos (k*pi*x));
#endif

#if INTHF
  foreach()
    sigma[] = 1.;
#endif
}

event vof (i++, first);

event amplitude (t += 3.04290519077e-3; t <= 2.2426211256) {

  scalar pos[];
  position (f, pos, {0,1 [0]});
  double max = statsf(pos).max;
  
  char name[80];
  sprintf (name, "wave-%d", N);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t*11.1366559937, max);
  fflush (fp);

  se += sq(max - prosperetti[ne][1]); ne++;
}

event error (t = end)
  fprintf (stderr, "%g %g\n", N/L0, sqrt(se/ne)/0.01);

#if 0
event gfsview (i += 1) {
  static FILE * fp = popen ("gfsview2D -s ../capwave.gfv", "w");
  output_gfs (fp);
}
#endif

/**
## Results

~~~gnuplot Evolution of the amplitude of the capillary wave as a function of non-dimensional time $\tau=\omega_0 t$
set xlabel '{/Symbol t}'
set ylabel 'Relative amplitude'
plot '../prosperetti.h' u 2:4 w l t "Prosperetti", \
     '../capclsvof/wave-64' every 10 w p t "CLSVOF", \
     '../caphf/wave-64' every 10 w p t "HF", \
     '../caphf2d/wave-64' every 10 w p t "HF2D",\
     'wave-64' every 10 w p t "VOF"
~~~

~~~gnuplot Convergence of the RMS error as a function of resolution (number of grid points per wavelength)
set xlabel 'Number of grid points'
set ylabel 'Relative RMS error'
set logscale y
set logscale x 2
set grid
plot [5:200][1e-4:1]\
    2./x**2 t "Second order",\
     '../capclsvof/log' t "CLSVOF" w lp, \
     '../caphf/log' t "HF" w lp, \
     '../caphf2d/log' t "HF2D" w lp,\
     'log' t "VOF" w lp
     
~~~

~~~pythonplot Amplitude of perturbation and Convergence
import numpy
import sys 
import string
import math
import glob
from pylab import *
from matplotlib.colors import LogNorm, Normalize,ListedColormap
from scipy.interpolate import CubicSpline

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

fig, ax = plt.subplots(1, 2, figsize=(6.29921, 2.5), gridspec_kw={'wspace': 0.35})

datafile1 = "wave-64"
data1=numpy.loadtxt(datafile1)

datafile2 = "../capclsvof/wave-64"
data2=numpy.loadtxt(datafile2)

datafile3 = "../caphf/wave-64"
data3=numpy.loadtxt(datafile3)

datafile4 = "../prosperetti.dat"
data4=numpy.loadtxt(datafile4)

datafile5 = "../caphf2d/wave-64"
data5=numpy.loadtxt(datafile5)

datafile11 = "log"
data11=numpy.loadtxt(datafile11)

datafile21 = "../capclsvof/log"
data21=numpy.loadtxt(datafile21)

datafile31 = "../caphf/log"
data31=numpy.loadtxt(datafile31)

datafile51 = "../caphf2d/log"
data51=numpy.loadtxt(datafile51)

NUM = 400
    
ax[0].set_ylabel(''r'$a$',labelpad=0)
ax[0].set_xlabel(r'$\omega_0 t$',labelpad=2)

ax[0].set_xticks([0,5,10,15,20,25])
ax[0].set_yticks([0,0.002,0.004,0.006,0.008,0.01])

ax[0].set_xlim(0,25)
ax[0].set_ylim(0,0.01)
    
#plot(data1[:,0],data1[:,1],'-',markersize=1.5,color = cmap[0])
ax[0].plot(data2[:,0],data2[:,1],'-',linewidth=4.5,color = "#ce2d4f",label="CLSVOF")
ax[0].plot(data3[:,0],data3[:,1],'-',linewidth=3.5,color = "#5aae61",label="HF")
ax[0].plot(data5[:,0],data5[:,1],'-',linewidth=1.25,color = "#0000FF",label="HF2D")

ax[0].plot(data4[0:140:2,0],data4[0:140:2,1],'o',markersize=2,color = "#000000",label="Prosperetti")
ax[0].plot(data4[141::5,0],data4[141::5,1],'o',markersize=2,color = "#000000")

ax[0].legend(ncol = 1, loc=1,frameon=False,fontsize=10)

x = np.logspace(0,2,num=NUM)

ax[1].plot(data21[:,0],data21[:,1]/2.,'-+',markersize=4,linewidth=2,color = "#ce2d4f",label="CLSVOF")
ax[1].plot(data31[:,0],data31[:,1]/2.,'-*',markersize=4,linewidth=2,color = "#5aae61",label="HF")
ax[1].plot(data51[:,0],data51[:,1]/2.,'-o',markersize=4,linewidth=2,color = "#0000FF",label="HF2D")
ax[1].plot(x,30/x/x/2,'-',linewidth=2,color = "#000000",label="Order 2")

ax[1].set_xlim(3,100)
ax[1].set_ylim(0.001,1)

ax[1].set_xscale("log")
ax[1].set_yscale("log")

ax[1].legend(ncol = 1, loc=1,frameon=False,fontsize=10)

ax[1].set_ylabel(r'$\epsilon_{L2}$',labelpad=0)
ax[1].set_xlabel(r'$N/\lambda$',labelpad=2)

fig.text(-0.035, 0.98, '$(a)$', fontsize=14, fontweight='bold')
fig.text(0.49, 0.98, '$(b)$', fontsize=14, fontweight='bold')

savefig("capwave.svg",bbox_inches='tight', pad_inches=1/2.54)

~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/capwave.html)
*/
