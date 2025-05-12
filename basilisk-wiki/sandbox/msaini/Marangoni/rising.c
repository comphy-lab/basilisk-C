//#include "axi.h" // fixme: does not run with -catch
#include "navier-stokes/centered.h"


#if CLSVOF
# include "two-phase-clsvof.h"
# include "integral.h"
#elif HF2D
# include "two-phase-HF.h"
#include "integral.h"
#elif INTHF
#include "two-phase.h"
#include "integral-HF.h"
scalar sigma[];
#else
#include "two-phase.h"
#include "tension.h"
#endif

#ifndef LEVEL
# define LEVEL 8
#endif

u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
uf.n[bottom] = 0.;
uf.n[top] = 0.;

int main() {

  size (2 [1]);
  DT = 1. [0,1];
  init_grid (1 << LEVEL);
  
  rho1 = 1000.[0], mu1 = 10.;  // works also with rho1 = [-3,0,1]
  rho2 = 100., mu2 = 1.;
  
#if CLSVOF || HF2D
  const scalar sigma[] = 24.5;
  d.sigmaf = sigma;
#elif INTHF
  f.sigmaf = sigma;
#else
  f.sigma = 24.5;  
#endif
  
  TOLERANCE = 1e-4 [*];

  run();
}

event init (t = 0) {

  mask (y > 0.5 ? top : none);

#if HF2D || CLSVOF
  foreach()
    d[] = sqrt (sq(x - 0.5) + sq(y)) - 0.25;
#else
  fraction (f, sq(x - 0.5) + sq(y) - sq(0.25));
#endif

#if INTHF
  foreach(){
    sigma[] = 24.5;
  }
#endif

}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}

event logfile (i++) {
  double xb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vb += u.x[]*dv;
    xb += x*dv;
    sb += dv;
  }
  static double sb0 = 0.;
  if (i == 0) {
    sb0 = sb;
  }
  printf ("%g %g %g %g %g %g %g %g \n", 
	  t, (sb - sb0)/sb0, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);
  fflush (stdout);
}

event interface (t = 3.) {

  output_facets (f, stderr);
}
#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-4,1e-3,1e-3}, LEVEL);
}
#endif

/**
## Results

~~~gnuplot Shapes at tend
reset

unset border

set style line 1 linecolor rgb '#ce2d4f' dt 1 linewidth 6 ps 0.2 pt 7
set style line 2 linecolor rgb '#5aae61' dt 1 linewidth 4 ps 0.2 pt 7
set style line 3 linecolor rgb '#2b6cbb' dt 1 linewidth 2 ps 0.2 pt 7
set style line 4 linecolor rgb '#000000' dt 2 linewidth 2 ps 0.2 pt 7

set key top outside horizontal font "sans-serif,10"
set size ratio -1

unset xtics
unset ytics

set lmargin 1
set bmargin 2.5
set tmargin 3
set rmargin 2

p "../rising-clsvof/log" u 2:1 w l ls 1 t "CLSVOF",'' u (-$2):1 w l ls 1 notitle,\
  "../rising-hf/log" u 2:1 w l ls 2 t "HF",'' u (-$2):1 w l ls 2 notitle,\
  "../rising-hf2d/log" u 2:1 w l ls 3 t "HF2D",'' u (-$2):1 w l ls 3 notitle,\
  "log" u 2:1 every 10 w p ls 4 t "VOF",'' u (-$2):1 every 10 w p ls 4 notitle
     
~~~


~~~pythonplot Rise velocity
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

datafile1 = "../rising-hf/out"
data1=numpy.loadtxt(datafile1)

datafile2 = "../rising-hf2d/out"
data2=numpy.loadtxt(datafile2)

datafile3 = "../rising-clsvof/out"
data3=numpy.loadtxt(datafile3)

datafile4 = "out"
data4=numpy.loadtxt(datafile4)

fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))

ax.set_ylabel(r'$u_y$',labelpad=1)
ax.set_xlabel(r'$t$',labelpad=1)

ax.set_xlim(0.0,3)
ax.set_ylim(0.0,0.25)

ax.plot(data3[::,0],data3[::,4],'-',linewidth=4.5,color = "#ce2d4f",label="CLSVOF")
ax.plot(data1[::,0],data1[::,4],'-',linewidth=3.5,color = "#5aae61",label="HF")
ax.plot(data2[::,0],data2[::,4],'-',linewidth=1.25,color = "#0000FF",label="HF2D")
ax.plot(data4[::30,0],data4[::30,4],'o',markersize=2,color = "#000000",label="CSF")

ax.legend(ncol = 1, loc=0,frameon=False)

savefig("./rising.svg",bbox_inches='tight', pad_inches=0.3/2.54)

~~~

*/
