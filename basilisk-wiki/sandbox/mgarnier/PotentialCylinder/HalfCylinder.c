/**
# Potential flow around a half cylinder

## Case 1 : Dirichlet conditions

In this case we are studying the case of a half cylinder. We take the code of [Potential flow around a cylinder](http://basilisk.fr/sandbox/mgarnier/PotentialCylinder/cyclindre.c) and we just shift the domain of computation of the solution.
*/


/**
As a starting point, i used my code [Poisson equation on complex domains](http://basilisk.fr/src/test/neumann.c). 

*/

#include "embed.h"
#include "poisson.h"
#include "view.h"

static double Rayon_cyl (void) {
  return 1;
}
static double Vitesse (void) {
  return 1;
}

static double exact (double x, double y) {
  double theta = atan2 (y, x), r2 = x*x + y*y;
  return Vitesse()*sqrt(r2)*cos(theta)*(1+(Rayon_cyl()*Rayon_cyl())/r2);
}

double exact_gradient (Point point, double theta, double r)
{
  coord n = facet_normal (point, cs, fs);
  normalize (&n);
  double dsdtheta = -Vitesse()*sin(theta)*(1+(Rayon_cyl()*Rayon_cyl())/(r*r));
  double dsdr = Vitesse()*cos(theta)*(1-(Rayon_cyl()*Rayon_cyl())/(r*r));
  return (n.x*(dsdr*cos(theta) - dsdtheta*sin(theta)) +
	  n.y*(dsdr*sin(theta) + dsdtheta*cos(theta)));  
}

static
double dirichlet_homogeneous_bc (Point point, Point neighbor,
				 scalar s, void * data) {
  return 0.;
}

/**
We just need to modify the origin to the middle of the left side of the box of calcul. 
*/

int main()
{
  for (N = 32; N <= 512; N *= 2) {
    size (6. [0]); // dimensionless
    origin (-L0, -L0/2);
    init_grid (N);

    
    vertex scalar phi[];
    foreach_vertex() {
     double r = sqrt (sq(x) + sq(y));

    phi[] = r - Rayon_cyl() ;

    }
    fractions (phi, cs, fs);  
#if TREE
    cs.refine = cs.prolongation = fraction_refine;
#endif
    restriction ({cs,fs});

    cm = cs;
    fm = fs;
    
    scalar adir[], b[], sol[];

    adir[left]   = exact (x - Delta/2., y);
    adir.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
    adir[right]  = exact (x + Delta/2., y);
    adir.boundary_homogeneous[right] = dirichlet_homogeneous_bc;
    adir[top]    = exact (x, y + Delta/2.);
    adir.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
    adir[bottom] = exact (x, y - Delta/2.);
    adir.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;
 
    /**
    The boundary conditions on the embedded boundary are Dirichlet  */


    adir[embed] = dirichlet(exact(x,y));

    adir.third = true;

    /**
    We need to solve this equation
    $$
    \Delta\phi = 0
    $$
    */
    
    foreach() {
      adir[] = cs[] > 0. ? exact (x, y) : nodata;
      
      
      double xc = x, yc = y;
      if (cs[] > 0. && cs[] < 1.) {
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	line_center (n, alpha, cs[], &p);
	xc += p.x*Delta, yc += p.y*Delta;
      }
      //      fprintf (stderr, "xc %g %g\n", xc, yc);
      b[] = 0;
    }

#if 0
    output_cells (stdout);
    output_facets (cs, stdout, fs);

    scalar e[];
    foreach() {
      if (cs[] > 0. && cs[] < 1.) {
	scalar s = a;
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	double length = line_length_center (n, alpha, &p);
	x += p.x*Delta, y += p.y*Delta;
	double theta = atan2(y,x), r = sqrt(x*x + y*y);
	
	double dsdtheta = - 3.*cube(r)*sin (3.*theta);
	double dsdr = 4.*cube(r)*cos (3.*theta);
	double nn = sqrt (sq(n.x) + sq(n.y));
	n.x /= nn, n.y /= nn;
	double dsdn = (n.x*(dsdr*cos(theta) - dsdtheta*sin(theta)) +
		       n.y*(dsdr*sin(theta) + dsdtheta*cos(theta)));

	e[] = dsdn - dirichlet_gradient (point, s, cs, n, p, exact (x, y));
#if 1
       fprintf (stderr, "g %g %g %g %g\n",
		x, y, dsdn,
		dirichlet_gradient (point, s, cs, n, p, exact (x, y)));
#endif
      }
      else
	e[] = nodata;
    }

    norm n = normf (e);
    fprintf (stderr, "%d %g %g\n",
	     N, n.rms, n.max);
#else

    /**
    The Poisson equation is solved. */
    
    struct Poisson p;
    p.alpha = fs;
    p.lambda = zeroc;
    p.embed_flux = embed_flux;
    scalar resdir[];
    double maxpdir = residual ({adir}, {b}, {resdir}, &p), maxfdir = 0.;
    foreach()
      if (cs[] == 1. && fabs(resdir[]) > maxfdir)
	maxfdir = fabs(resdir[]);
    fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxfdir, maxpdir);

    // FIXME: need to set minlevel to 4
    timer t = timer_start();
    mgstats s = poisson (adir, b, alpha = fs, tolerance = 1e-6, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

    /**
    The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
    fields and their norms are computed. */
    
    scalar edir[], epdir[], efdir[];
    foreach() {
      if (cs[] == 0.)
	epdir[] = efdir[] = edir[] = nodata;
      else {
	edir[] = fabs(adir[] - exact (x, y));
	epdir[] = cs[] < 1. ? edir[] : nodata;//Erreur sur les frontières immmergées
	efdir[] = cs[] >= 1. ? edir[] : nodata;//Erreur sur le domaine de résolution
      }
    }
    norm ndir = normf (edir), npdir = normf (epdir), nfdir = normf (efdir);
    fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %d %d\n",
	     N, ndir.avg, ndir.max, npdir.avg, npdir.max, nfdir.avg, nfdir.max, s.i, s.nrelax);
      
  FILE * fp2 = fopen ("errors_dir.dat","a");
  fprintf(fp2,"%g %g %g %g %g %g %d \n", ndir.avg, ndir.max, npdir.avg, npdir.max, nfdir.avg, nfdir.max, N);
  fflush (fp2);

    FILE * fp1 = fopen ("potentialdir.dat","w");
    foreach() 
      fprintf(fp1,"%g %g %g %g \n", x, y, adir[], sol[]);
    fflush (fp1);




    /**
## Case 2 : Neumann conditions
 */

    scalar aneu[];

    aneu[left]   = exact (x - Delta/2., y);
    aneu.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
    aneu[right]  = exact (x + Delta/2., y);
    aneu.boundary_homogeneous[right] = dirichlet_homogeneous_bc;
    aneu[top]    = exact (x, y + Delta/2.);
    aneu.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
    aneu[bottom] = exact (x, y - Delta/2.);
    aneu.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;

/**
    The boundary conditions on the embedded boundary are Neumann  */

    aneu[embed] = neumann (exact_gradient (point, atan2(y, x), sqrt(x*x + y*y)));


    aneu.third = true;

    foreach() {
      aneu[] = cs[] > 0. ? exact (x, y) : nodata;
      
      
      double xc = x, yc = y;
      if (cs[] > 0. && cs[] < 1.) {
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	line_center (n, alpha, cs[], &p);
	xc += p.x*Delta, yc += p.y*Delta;
      }
      //      fprintf (stderr, "xc %g %g\n", xc, yc);
      b[] = 0;
    }

    struct Poisson p1;
    p1.alpha = fs;
    p1.lambda = zeroc;
    p1.embed_flux = embed_flux;
    scalar resneu[];
    double maxpneu = residual ({aneu}, {b}, {resneu}, &p1), maxfneu = 0.;
    foreach()
      if (cs[] == 1. && fabs(resneu[]) > maxfneu)
	maxfneu = fabs(resneu[]);
    fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxfneu, maxpneu);

    // FIXME: need to set minlevel to 4
    timer t1 = timer_start();
    mgstats s1 = poisson (aneu, b, alpha = fs, tolerance = 1e-6, minlevel = 4);
    double dt1 = timer_elapsed (t1);
    printf ("%d %g %d %d\n", N, dt1, s1.i, s1.nrelax);


    /**
    The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
    fields and their norms are computed. */
    
    scalar eneu[], epneu[], efneu[];
    foreach() {
      if (cs[] == 0.)
	epneu[] = efneu[] = eneu[] = nodata;
      else {
	eneu[] = fabs(aneu[] - exact (x, y));
	epneu[] = cs[] < 1. ? eneu[] : nodata;
	efneu[] = cs[] >= 1. ? eneu[] : nodata;
      }
    }
    foreach() {
      if (cs[] == 0.)
      sol[] = nodata;
    else
      sol[] = exact(x,y);  
    }


    norm nneu = normf (eneu), npneu = normf (epneu), nfneu = normf (efneu);
    fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %d %d\n",
	     N, nneu.avg, nneu.max, npneu.avg, npneu.max, nfneu.avg, nfneu.max, s.i, s.nrelax);

  FILE * fp4 = fopen ("errors_neu.dat","a");
  fprintf(fp4,"%g %g %g %g %g %g %d \n", nneu.avg, nneu.max, npneu.avg, npneu.max, nfneu.avg, nfneu.max, N);
  fflush (fp4);

    FILE * fp = fopen ("potentialneu.dat","w");
    foreach() 
      fprintf(fp,"%g %g %g %g \n", x, y, aneu[], sol[]);
    fflush (fp);
#endif
    dump ("dump");
  }
} 

/** 

~~~pythonplot Legend

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

L = 6
l = L/2


a_neu = np.loadtxt('potentialneu.dat')

x_neu = a_neu[:,0]
y_neu = a_neu[:,1]
phi_neu = a_neu[:,2]
th_neu = a_neu[:,3]


a_dir = np.loadtxt('potentialdir.dat')

x_dir = a_dir[:,0]
y_dir = a_dir[:,1]
phi_dir = a_dir[:,2]


xi_neu = np.linspace(x_neu.min() - 1, x_neu.max() + 1, 5000)
yi_neu = np.linspace(y_neu.min() - 1, y_neu.max() + 1, 5000)
xi_neu, yi_neu = np.meshgrid(xi_neu, yi_neu)


zi1_neu = griddata((x_neu/l, y_neu/l), th_neu, (xi_neu, yi_neu), method='linear')
zi2_neu = griddata((x_neu/l, y_neu/l), phi_neu, (xi_neu, yi_neu), method='linear')


xi_dir = np.linspace(x_dir.min() - 1, x_dir.max() + 1, 5000)
yi_dir = np.linspace(y_dir.min() - 1, y_dir.max() + 1, 5000)
xi_dir, yi_dir = np.meshgrid(xi_dir, yi_dir)

zi1_dir = griddata((x_dir/l, y_dir/l), th_neu, (xi_dir, yi_dir), method='linear')
zi3_dir = griddata((x_dir/l, y_dir/l), phi_dir, (xi_dir, yi_dir), method='linear')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

contour1_neu = ax1.contour(xi_neu, yi_neu, zi1_neu, levels=np.arange(-6.16, 0, 0.664), linewidths=2., colors='r', linestyles='solid')
contour2_neu = ax1.contour(xi_neu, yi_neu, zi2_neu, levels=np.arange(-6.16, 0, 0.664), linewidths=2., colors='g', linestyles='dashed')

theta = np.linspace(0, 2 * np.pi, 100)
x_circ = np.cos(theta)
y_circ = np.sin(theta)
ax1.fill(x_circ/l, y_circ/l, color='black')

ax1.axis('square')
ax1.set_xlim(left=-2, right=0)
ax1.set_ylim(bottom=-1, top=1)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Neumann boundary conditions')

contour1_dir = ax2.contour(xi_dir, yi_dir, zi1_dir, levels=np.arange(-6.16, 0, 0.664), linewidths=2., colors='r', linestyles='solid')
contour2_dir = ax2.contour(xi_dir, yi_dir, zi3_dir, levels=np.arange(-6.16, 0, 0.664), linewidths=2., colors='g', linestyles='dashed')

ax2.fill(x_circ/l, y_circ/l, color='black')

ax2.axis('square')
ax2.set_xlim(left=-2, right=0)
ax2.set_ylim(bottom=-1, top=1)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('Dirichlet boundary conditions')

plt.tight_layout()
plt.savefig('isoline_halfcylinder.png')
plt.show()
~~~
*/
