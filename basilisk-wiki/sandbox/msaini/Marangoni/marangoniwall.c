/**
In this case, we study the effect of wall on the drop driven by Marangoni flows. The test is discussed in detail in our new paper. Note that the results here are obtained lower resolution as compared to the paper because of limited computational time on the Basilisk.fr servers. Therefore, small differences persist as compared to the figures the paper.
 */

#include "axi.h"
#include "navier-stokes/centered.h"
#if CLSVOF
# include "two-phase-clsvof.h"
# include "integral.h"
#elif HF2D
# include "two-phase-HF.h"
#include "integral.h"
#else
#include "two-phase.h"
#include "integral-HF.h"
#endif

#include "advdiff.h"
#include "view.h"
int LEVEL = 8;

u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);

PhiC[left] = dirichlet(L0);
PhiC[right] = dirichlet(0);

const double R = 1. [1], NablaT = 1., Mu = 1., Rho = 1. [0];
const double Re =  0.08;
#if CASE2
const double Ca = 0.15;
#else
const double Ca = 4.55;
#endif
const double Gamma_T = Re*sq(Mu)/(Rho*sq(R)*NablaT);
const double Gamma_0 = (Gamma_T*R*NablaT)/Ca;
const double t0 = Mu/(Gamma_T*NablaT);
const double Cdrop = 1., Cbulk = 1.;
double U_drop;

scalar sigmaf[];

int main()
{
  size (16*R);
  rho1 = rho2 = Rho;
  mu1 = Mu; mu2 = Mu*0.01;

  Diff_C1 = 100000.;
  Diff_C2 = 100000.*0.01;

#if CLSVOF || HF2D
  d.sigmaf = sigmaf;
#else
  f.sigmaf = sigmaf;
#endif
  
  TOLERANCE = 1e-4 [*];
  
  U_drop = - 2./((2. + 3.*mu2/mu1)*(2. + Diff_C2/Diff_C1))*Gamma_T*R*NablaT/mu1;

  N = 1 << LEVEL;  
  run();
}

event init (t = 0)
{
  vertex scalar dist[];
  foreach_vertex()
    dist[] = sqrt (sq(x - 4.*R) + sq(y)) - R;
  
#if CLSVOF || HF2D
  foreach() {
    d[] = sqrt (sq(x - 4.*R) + sq(y)) - R;
  }
#else
  fractions(dist,f);
#endif    
  
  foreach() {
    PhiC[] = 0;   
    sigmaf[] = Gamma_0 - Gamma_T*(PhiC[] - L0);
  }
}

event properties(i++){
  foreach()
    sigmaf[] = Gamma_0 - Gamma_T*(PhiC[] - L0);
}

double u_drop = 0.;

event logfile (i += 2)
{
  double xb = 0., vb = 0., sb = 0.;
  static double xb0 = 0., previous = 0.;
  if (t == 0.)
    previous = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vb += u.x[]*dv;
    xb += x*dv;
    sb += dv;
  }
  static double sb0 = 0.;
  if (i == 0) {
    sb0 = sb;
    fprintf (ferr, "#dsb xb vb/U_drop ta u_drop/U_drop dt\n");
  }
  if(sb0 == 0)
    sb0 = 1;
  u_drop = t > previous ? (xb/sb - xb0)/(t - previous) : 0.;
  fprintf (ferr, "%g %g %g %g %g %g %g\n",
	   t/t0, (sb - sb0)/sb0, xb/sb, vb/sb/U_drop,
	   (t + previous)/2./t0, u_drop/U_drop,
	   dt);
  xb0 = xb/sb, previous = t;
}

event movie (t += 7.5*t0/120.){

#if !HF2D
  clear();
  view(fov = 10., ty = -0.02, quat = {0,0,-cos(pi/4.),cos(pi/4.)}, width = 1980, height = 1980);
  draw_vof("f",lw=4);
  squares("PhiC", linear=true, map = cool_warm);
  vectors (u = "u", scale = 0.1);
  mirror({0,1}) {
    draw_vof("f",lw=4);
    squares("PhiC", linear=true, map = cool_warm);
    vectors (u = "u", scale = 0.05);
  }
  char fname[80];
  sprintf(fname,"movie-%2.2f.mp4",Ca);
  save(fname);
#endif
}

event end (t = 7.5*t0)
{
  double U = - 2./((2. + 3.*mu2/mu1)*(2. + Diff_C2/Diff_C1))*Gamma_T*R*NablaT/mu1;
  fprintf (stderr, "#%d %g %g %g\n", N/16, u_drop/U, U_drop/U,Ca);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-2, 1e-6, 1e-6}, LEVEL);
}
#endif

/**
Here are some of the results for HF method: 

<video width="512" height="512" controls>
<source src="marangoniwall/movie-4.55.mp4" type="video/mp4">
</video> 

<video width="512" height="512" controls>
<source src="marangoniwall2/movie-0.15.mp4" type="video/mp4">
</video> 

Note how the oscillations are higher at lower Capillary numbers (Right video) as compared to the high capillary number (Left video). Additionally, the oscillations are concentrated around the points where $n_x = n_y$ where normal to interface $\mathbf{n} = [n_x,n_y]$. This indicates that the origin of the oscillations is change in the orientation of height functions. This problem does not exist in the CLSVOF method and the velocity fields are relatively smoother as seen from the follwing videos.

<video width="512" height="512" controls>
<source src="marangoniwall-clsvof/movie-4.55.mp4" type="video/mp4">
</video> 

<video width="512" height="512" controls>
<source src="marangoniwall-clsvof2/movie-0.15.mp4" type="video/mp4">
</video> 

~~~pythonplot Velocity as function of drop position and its comparison with baseline solutions.

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

datafile1 = "../marangoniwall-clsvof2/log"
data1=numpy.loadtxt(datafile1)

datafile11 = "../marangoniwall2/log"
data11=numpy.loadtxt(datafile11)

datafile111 = "../marangoniwall-hf2d2/log"
data111=numpy.loadtxt(datafile111)

datafile2 = "../marangoniwall-clsvof/log"
data2=numpy.loadtxt(datafile2)

datafile22 = "./log"
data22=numpy.loadtxt(datafile22)

datafile222 = "../marangoniwall-hf2d/log"
data222=numpy.loadtxt(datafile222)

fig, ax = plt.subplots(1, 2, figsize=(6.29921, 2.5), gridspec_kw={'wspace': 0.42})

ax[0].plot(data1[:4800:10,2],data1[:4800:10,3],'--',linewidth=2,color='#FF0000',label="CLSVOF")
ax[0].plot(data11[::10,2],data11[::10,3],'--',linewidth=2,color='#00FF00',label="HF")
ax[0].plot(data111[::10,2],data111[::10,3],'--',linewidth=2,color='#0000FF',label="HF2D")

ax[0].set_xlim(0,3)
ax[0].set_ylim(0,1.1)

ax[0].set_xlabel(r"$D/R$")
ax[0].set_ylabel(r"$u/u_{YGB}$")

datafile4 = "../theoretical.dat"
data4=numpy.loadtxt(datafile4)

ax[0].plot(data4[::,0],data4[::,1],'x',markersize=3,linewidth=2,color='#000000',label="Meyappan et al.")

ax[0].legend(loc='lower right',frameon=False,fontsize=9,borderpad=-0.275)

ax[1].plot(data2[:4800:10,2],data2[:4800:10,3],'--',linewidth=2,color='#FF0000',label="CLSVOF")
ax[1].plot(data22[::10,2],data22[::10,3],'--',linewidth=2,color='#00FF00',label="HF")
ax[1].plot(data222[::10,2],data222[::10,3],'--',linewidth=2,color='#0000FF',label="HF2D")
#last plot is to check convergence

ax[1].set_xlim(0,3)
ax[1].set_ylim(0,1.1)

ax[1].set_xlabel(r"$D/R$")
ax[1].set_ylabel(r"$u/u_{YGB}$")

datafile5 = "../numerical.dat"
data5=numpy.loadtxt(datafile5)

ax[1].plot(data5[::,0],data5[::,1],'o',markersize=3,linewidth=2,color='#000000',label="Ascoli and Leal")
ax[1].legend(loc=0,frameon=False,fontsize=9,borderpad=-0.275)

fig.text(0.0, 0.98, '$(a)$', fontsize=14, fontweight='bold')
fig.text(0.48, 0.98, '$(b)$', fontsize=14, fontweight='bold')

savefig("velwall.svg",bbox_inches='tight', pad_inches=0.25/2.54)
~~~

*/
