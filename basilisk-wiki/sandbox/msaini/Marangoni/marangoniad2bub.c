/**
In this case, we study the effect of wall on the drop driven by Marangoni flows. The test is discussed in detail in our new paper. Note that the results here are obtained lower resolution as compared to the paper because of limited computational time on the Basilisk.fr servers. As you refine the grid, you will obtain better results with lower oscillations.
 */

#include "axi.h"
#include "navier-stokes/centered.h"
#if CLSVOF
# include "two-phase-clsvof.h"
# include "integral.h"
#else
# include "two-phase-HF.h"
#include "integral.h"
#endif

#include "advdiff.h"
#include "view.h"

#include "tag.h"

int LEVEL = 8;

PhiC[left] = dirichlet(L0);
PhiC[right] = dirichlet(0);

const double R = 1. [1], NablaT = 1., Mu = 1., Rho = 1. [0];
const double Re =  0.066, Ca = 0.66;
const double Gamma_T = Re*sq(Mu)/(Rho*sq(R)*NablaT);
const double Gamma_0 = (Gamma_T*R*NablaT)/Ca;
const double t0 = Mu/(Gamma_T*NablaT);
const double Cdrop = 1., Cbulk = 1.;
double U_drop;

scalar sigmaf[];
static FILE * fp;

int main()
{
  size (12*R);
  rho1 = rho2 = Rho;
  mu1 = Mu; mu2 = Mu*0.0005;

  Diff_C1 = 1.e+4;
  Diff_C2 = 1.e+4*0.0005;

  d.sigmaf = sigmaf;
    
  TOLERANCE = 1e-4 [*];
  U_drop = - 2./((2. + 3.*mu2/mu1)*(2. + Diff_C2/Diff_C1))*Gamma_T*R*NablaT/mu1;

  N = 1 << (LEVEL - 3);
  run();
}

/**
We initialize the signed distance *d* and the surface tension gradient. */

event init (t = 0)
{
  fp = fopen ("bubbles.dat", "w");
  
  if (!restore (file = "restart")) {
    for(int ii = 0; ii < 3; ii++)
      refine(y < (L0/2 + 4.*Delta) && level < (LEVEL-(2. - ii)));
    
    vertex scalar dist[];
    foreach_vertex()
      dist[] = min(sqrt (sq(x - L0/2. - R) + sq(y)) - R,sqrt (sq(x - L0/2. + R/3. + R*0.5) + sq(y)) - R/3.);
    
    fractions(dist,f);
    
    foreach() {
      d[] = min(sqrt (sq(x - L0/2. - R) + sq(y)) - R,sqrt (sq(x - L0/2. + R/3. + R*0.5) + sq(y)) - R/3.);
    }
    
    foreach(){
      PhiC[] = 0;
      sigmaf[] = Gamma_0 - Gamma_T*(PhiC[] - L0);
    }
  }
}

event properties(i++){
  foreach()
    sigmaf[] = Gamma_0 - Gamma_T*(PhiC[] - L0);
}

double u_drop = 0.;

event statistics(t+=2.*t0/2000.){
  double U_drops[2] = {- 2./((2. + 3.*mu2/mu1)*(2. + Diff_C2/Diff_C1))*Gamma_T*R/3.*NablaT/mu1,- 2./((2. + 3.*mu2/mu1)*(2. + Diff_C2/Diff_C1))*Gamma_T*R*NablaT/mu1};

  scalar m[];
  foreach()
    m[] = (1. - f[]) > 1e-3;
  int n = tag (m);
  double v[n];
  coord b[n];
  coord ub[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = ub[j].x = ub[j].y =0.;
  foreach (serial)
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*(1. - f[]);
      coord p = {x,y,z};
      foreach_dimension(){
	b[j].x += dv()*(1. - f[])*p.x;
	ub[j].x += dv()*(1. - f[])*u.x[];
      }
    }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ub, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  if(pid() ==0 && n > 1){
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g ", t/t0,b[j].x/v[j],ub[j].x/v[j]/U_drops[j]);
  }
  fprintf(fp,"\n");
}

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
    fprintf (stdout, "#dsb xb vb/U_drop ta u_drop/U_drop dt perf.t perf.speed\n");
  }
  u_drop = t > previous ? (xb/sb - xb0)/(t - previous) : 0.;
  fprintf (ferr, "%g %g %g %g %g %g %g\n",
	   t/t0, (sb - sb0)/sb0, xb/sb, vb/sb/U_drop,
	   (t + previous)/2./t0, u_drop/U_drop,
	   dt);
  xb0 = xb/sb, previous = t;
}

event movie (t += 1.68*t0/120.){
  clear();
  view(fov = 9., ty = -0.5, quat = {0,0,-cos(pi/4.),cos(pi/4.)}, width = 1980, height = 1980);
  draw_vof("f",lw=4);
  squares("PhiC", linear=true, map = cool_warm);
  vectors (u = "u", scale = 0.1);
  mirror({0,1}) {
    draw_vof("f",lw=4);
    squares("PhiC", linear=true, map = cool_warm);
    vectors (u = "u", scale = 0.05);
  }
  save("movie.mp4");
}

event graphics (t = 1.5*t0)
{
  double U = -2./((2. + 3.*mu2/mu1)*(2. + Diff_C2/Diff_C1))*Gamma_T*R*NablaT/mu1;;
  fprintf (stderr, "#%d %g %g %g\n", N/16, u_drop/U, U_drop/U,Diff_C2/Diff_C1);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-2, 1e-5, 1e-5}, LEVEL);
}
#endif

/**
<center>
<video width="426" height="660" controls>
<source src="marangoniad2bub/movie.mp4" type="video/mp4">
</video> 
</center>
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

colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#FFC107', 
          '#6A3D9A', '#00CED1', '#4D4D4D', '#FB8072']

datafile1 = "./bubbles.dat"
data1=numpy.loadtxt(datafile1)

datafile = "../marangoniad2bub-clsvof/bubbles.dat"
data=numpy.loadtxt(datafile)

datafile2 = "../thSmall.dat"
data2=numpy.loadtxt(datafile2)

datafile3 = "../thBig.dat"
data3=numpy.loadtxt(datafile3)

fig, ax = plt.subplots(figsize=(6/2.54, 5.5/2.54))

xdata = np.linspace(0.15,1,num=200)

ax.plot((-data[::5,1] + data[::5,4] - 1. - 1./3.),data[::5,2],'--',linewidth=3,color=colors[0],label=r"CLSVOF, $u_S$")
ax.plot((-data1[::4,1] + data1[::4,4] - 1. - 1./3.),data1[::4,2],'--',linewidth=2.5,color=colors[1],label=r"HF2D, $u_S$")
ax.plot(xdata,1.+2.*(3./(3.*xdata+4.))**3.,'-.',linewidth=2,color=colors[2],label="Equation (50)")
ax.plot(data2[:,0],data2[:,1],'--',linewidth=2,color=colors[3],label=r"Meyyappan et al. ($u_S$)")

ax.set_xlim(0,0.6)
ax.set_ylim(0.9,2.1)

ax.set_xlabel(r"$H/R$")
ax.set_ylabel(r"$u/u_{YGB}$")

ax.plot((-data[::5,1] + data[::5,4] - 1. - 1./3.),data[::5,5],'--',linewidth=2.5,color=colors[4],label=r"CLSVOF, $u_L$")
ax.plot((-data1[::4,1] + data1[::4,4] - 1. - 1./3.),data1[::4,5],'--',linewidth=3,color=colors[5],label=r"HF2D, $u_L$")
ax.plot(xdata,1.-2./3.*(1/(3.*xdata+4.))**3.,'-.',linewidth=2,color=colors[6],label="Equation (51)")
ax.plot(data3[:,0],data3[:,1],'--',linewidth=2,color=colors[7],label=r"Meyyappan et al. ($u_L$)")

ax.legend(loc=0,frameon=False,fontsize=9,bbox_to_anchor=(1, 0.15))

savefig("vel2bub.svg",bbox_inches='tight', pad_inches=0.3/2.54)

~~~

*/
