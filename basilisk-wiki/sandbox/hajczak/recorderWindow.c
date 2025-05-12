 /**
# Quasi 2D incompressible flow in a recorder window
In flute-like instruments, an air jet flows out of a channel and impinges on an edge (called labium) after a distance $W$. A hard-wall pipe is located under the edge (the resonator). Above the edge is the free-field. This defines the so-called window of the flute.

We reproduce the bidimensional test case from [Auvray *et al.* (2013)](http://www.lam.jussieu.fr/Publications/Conferences/SMAC_2013_Auvray.pdf). The pipe diameter is taken as $4W$ and the channel is ignored in the geometry. It is interesting to determine the potential flow in the recorder window since the production of acoustic waves is related to the rate at which vortices cross the streamlines of the potential flow (cf. [Howe 2002, Lectures on the Theory of Vortex sound](https://link.springer.com/chapter/10.1007/3-540-45880-8_2))

The 2D potential flow is solved using the Poisson-Helmholtz solver with embedded boundaries. The geometry is a recorder labium of angle 15°.

Reference data has been (roughly) digitized from the publication pdf using [WebPlotDigitizer](https://apps.automeris.io/wpd/) 
*/

#include "embed.h"
#include "poisson.h"
#include "view.h"

double W = 1.;
double alpha_b = 15*M_PI/180.;

/**
This function is necessary to impose Dirichlet BC. 
*/

static
double dirichlet_homogeneous_bc (Point point, Point neighbor,
				 scalar s, void * data) {
  return 0.;
}
/** 
*/
int main()
{
  N=512; // Fixme : the code crashes if N=1024 but still works when N<=512. Why ?
  size (10. [0]); 
  origin (-1, -4.); // Set the origin at the labium edge
  init_grid (N);					
  solid (cs, fs,union (y - x*tan(alpha_b), -y));
  /* Fixme : when using the foreach_vertex phi[] = ...
     fractions(cs,fs,phi)... we obtain a geometry with
     cs equal either to 0 or 1 (fractions not computed ???)
  */
  cm = cs;
  fm = fs;

  scalar a[], b[]; // Define the potential field

/** The boundary conditions are as follows:

* Neumann ($\partial \phi / \partial n = 0$) on the bottom, left and embedded boundaries

* Dirichlet $\phi = -x$ at outflow and $y < 0$

* Dirichlet $\phi = ln(r)$ at top and outflow $y>(L_0 / 2)\text{tan}(\alpha)$

* NB : this is different from Auvray who solves the streamfunction problem and thus imposes Dirichlet BC everywhere on $\psi$

* The geometry is also slightly different since there is a portion of the left boundary which is not a hard wall in Auvray.

* NB : I am not sure if the Delta/2 shifts are important in the BC

*/
  
  a[bottom] = neumann(0.);
  a[left] = neumann(0.);
  a[top] = dirichlet(log(sqrt(sq(x)+sq(y))));
  a[right] = (y < Delta/2.) ? dirichlet(-9.) : (y > tan(alpha_b)*(x))-Delta/2. ? dirichlet(log(sqrt(sq(x)+sq(y)))) : 0.;
  //a.boundary_homogeneous[right] = dirichlet_homogeneous_bc; 
  a[embed] = neumann(0.); 
  /**
Use "third-order interpolation" as in the example.
  */
  
  a.third = true;

  /**
  */
  
  struct Poisson p;
  p.alpha = fs;
  p.lambda = zeroc;
  p.embed_flux = embed_flux;
  scalar res[]; 
  double maxp = residual ({a}, {b}, {res}, &p), maxf = 0.;  
  foreach()  
    if (cs[] == 1. && fabs(res[]) > maxf)  // NB : large residual values on boundaries
      maxf = fabs(res[]);  
  fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxf, maxp);  
	  
  // FIXME: need to set minlevel to 4
  timer t = timer_start();
  mgstats s = poisson (a, b, alpha = fs, tolerance = 1e-4, minlevel = 4);
  double dt = timer_elapsed (t);
  printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
    
  /**
     Compute velocity field from potential and integrate back to get the streamfunction
     
     Fixme : compute directly $\psi$ using Dirichlet BC as Auvray...
  */

  vector u[];
  gradients ({a}, {u});
  
  scalar psi[];
  psi[bottom] = 0.;  
  foreach()
    {
      psi[] = -u.y[]*Delta + psi[-1,0];
    }

  FILE * fp = fopen ("streamlines.dat","w");
  foreach() 
    fprintf(fp,"%g %g %g %g %g %g \n", x, y, a[], u.x[], u.y[], psi[]);
  fflush (fp);

}

/**
Red marks are streamfunction isocontours obtained from Auvray et al.
~~~pythonplot
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

a=np.loadtxt('streamlines.dat')
x=a[:,0]
y=a[:,1]
phi=a[:,2]
ux=a[:,3]
uy=a[:,4]
psi=a[:,5]

xi = np.linspace(x.min()-1, x.max()+1, 100)
yi = np.linspace(y.min()-1, y.max()+1, 100)
zi = griddata((x, y), psi, (xi[None,:], yi[:,None]), method='linear')
b=np.loadtxt('../refAuvray.dat')
plt.scatter(b[:,0],b[:,1],marker='+',c='red',s=5)
plt.contour(xi,yi,zi,25,linewidths=1.,colors='k',linestyles='solid')
#plt.quiver(x,y,ux,uy,units='inches')
plt.fill([0.,1.,1.],[0.,0.,np.tan(15.*np.pi/180.)],'k')
plt.axis('square')
plt.xlim(left=-1, right=1)
plt.ylim(bottom=-1, top=1)
plt.savefig('streamlines.png')
~~~

* Some streamlines cross the labium edge due to insufficient resolution in that region. Fixme : manage to refine without crashing the code


*/

/**
## References

* Howe, M. S. (2002). Lectures on the theory of vortex-sound. In Sound-flow interactions (pp. 31-111). Berlin, Heidelberg: Springer Berlin Heidelberg.

* Auvray, R., Ernoult, A., Fabre, B., & Lagrée, P. Y. (2014). Time-domain simulation of flute-like instruments: Comparison of jet-drive and discrete-vortex models. The Journal of the Acoustical Society of America, 136(1), 389-400.

*/
