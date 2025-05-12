/**
# Potential flow around a cylinder

## Case 1 : Dirichlet conditions
 */
#include "embed.h"
#include "poisson.h"
#include "view.h"





/**
As a starting point, i used the code [Poisson equation on complex domains](http://basilisk.fr/src/test/neumann.c). 


First of all, we know the exact solution of this problem. The solution is written in the following form : 
$$
\phi(r,\theta) = Urcos(\theta) \left( 1+\frac{R^2}{r^2} \right)
$$
with R the radius of the cylinder and U the velocity of the input flow. 

In the beginning of the code, we create two functions wich will allow us to easily modifiy the values of the two parameters $U$ and $R$. 
*/
static double Rayon_cyl (void) {
  return 1;
}
static double Vitesse (void) {
  return 1;
}

/**
The tree following functions will be useful in order to define the different boundary conditions in our problem. 
*/

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

int main()
{
  for (N = 32; N <= 512; N *= 2) {
    size (6. [0]); // dimensionless
    origin (-L0/2, -L0/2);
    init_grid (N);

    /**
    The shape of the domain is given by : \
    $$
    \Omega = {(r,\theta): r \geq R}
    $$
    */
    
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
    
        /**
    In our case, the boundary condition on the box boundaries defined here. We also define $a$ the numerical solution and $sol$ the theorical solution, so that we can compare them. */
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

    /**
    We use "third-order" [face flux interpolation](/src/embed.h). */

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

    scalar respdir[], resfdir[];
    foreach() {
      if (cs[] == 0.)
	respdir[] = resfdir[] = resdir[] = nodata;
      else {
	respdir[] = cs[] < 1. ? resdir[] : nodata;
	resfdir[] = cs[] >= 1. ? resdir[] : nodata;
      }
    }

  norm nresdir = normf(resdir), nrespdir = normf(respdir), nresfdir = normf(resfdir);
  FILE * fp6 = fopen ("res_dir.dat","a");
  fprintf(fp6,"%g %g %g %g %g %g %d \n", nresdir.avg, nresdir.max, nrespdir.avg, nrespdir.max, nresfdir.avg, nresfdir.max, N);
  fflush (fp6);

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

  scalar respneu[], resfneu[];
    foreach() {
      if (cs[] == 0.)
	respneu[] = resfneu[] = resneu[] = nodata;
      else {
	respneu[] = cs[] < 1. ? resneu[] : nodata;
	resfneu[] = cs[] >= 1. ? resneu[] : nodata;
      }
    }

  norm nresneu = normf(resneu), nrespneu = normf(respneu), nresfneu = normf(resfneu);
  FILE * fp5 = fopen ("res_neu.dat","a");
  fprintf(fp5,"%g %g %g %g %g %g %d \n", nresneu.avg, nresneu.max, nrespneu.avg, nrespneu.max, nresfneu.avg, nresfneu.max, N);
  fflush (fp5);

    /**
    The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
    fields and their norms are computed. */
    
    scalar eneu[], epneu[], efneu[];
    foreach() {
      if (cs[] == 0.)
	epneu[] = efneu[] = eneu[] = nodata;
      else {
	eneu[] = fabs(aneu[] - exact (x, y));
	epneu[] = cs[] < 1. ? eneu[] : nodata;//Erreur sur les frontières immmergées
	efneu[] = cs[] >= 1. ? eneu[] : nodata;//Erreur sur le domaine de résolution
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
Here are the results. For each case (numerical solution with neumann boundary conditions and numerical solution with dirichlet boundary conditions). I choose to plot the isoline of potential beacause when i plotted teh solution it was impossible to see any difference with the naked eye. Here are the results. 

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


contour1_neu = ax1.contour(xi_neu, yi_neu, zi1_neu, levels=np.arange(-3.32, 3.32, 0.664), linewidths=2., colors='r', linestyles='solid')
contour2_neu = ax1.contour(xi_neu, yi_neu, zi2_neu, levels=np.arange(-3.32, 3.32, 0.664), linewidths=2., colors='g', linestyles='dashed')

theta = np.linspace(0, 2 * np.pi, 100)
x_circ = np.cos(theta)
y_circ = np.sin(theta)
ax1.fill(x_circ/l, y_circ/l, color='black')

ax1.axis('square')
ax1.set_xlim(left=-1, right=1)
ax1.set_ylim(bottom=-1, top=1)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Neumann boundary conditions')


contour1_dir = ax2.contour(xi_dir, yi_dir, zi1_dir, levels=np.arange(-3.32, 3.32, 0.664), linewidths=2., colors='r', linestyles='solid')
contour2_dir = ax2.contour(xi_dir, yi_dir, zi3_dir, levels=np.arange(-3.32, 3.32, 0.664), linewidths=2., colors='g', linestyles='dashed')

ax2.fill(x_circ/l, y_circ/l, color='black')

ax2.axis('square')
ax2.set_xlim(left=-1, right=1)
ax2.set_ylim(bottom=-1, top=1)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('Dirichlet boundary conditions')

plt.tight_layout()
plt.savefig('contour_plots.png')

~~~

## Mesh convergence with Neumann conditions

To finsh, we are looking at the influence of the mesh size on the computation of the solution. 

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_law(x, a, b):
    return a * x**b

data = np.loadtxt('errors_neu.dat')

nneu_avg = data[:, 0]
nneu_max = data[:, 1]
npneu_avg = data[:, 2]
npneu_max = data[:, 3]
nfneu_avg = data[:, 4]
nfneu_max = data[:, 5]
N = data[:, 6]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

popt_nneu_avg, _ = curve_fit(power_law, N, nneu_avg)
popt_npneu_avg, _ = curve_fit(power_law, N, npneu_avg)
popt_nfneu_avg, _ = curve_fit(power_law, N, nfneu_avg)

popt_nneu_max, _ = curve_fit(power_law, N, nneu_max)
popt_npneu_max, _ = curve_fit(power_law, N, npneu_max)
popt_nfneu_max, _ = curve_fit(power_law, N, nfneu_max)

ax1.plot(N, nneu_avg, 'o', label='all cells')
ax1.plot(N, npneu_avg, 's', label='partial cells')
ax1.plot(N, nfneu_avg, '^', label='full cells')

ax1.plot(N, power_law(N, *popt_nneu_avg), '-', label='{:.2f}x^({:.2f})'.format(popt_nneu_avg[0], popt_nneu_avg[1]))
ax1.plot(N, power_law(N, *popt_npneu_avg), '--', label='{:.2f}x^({:.2f})'.format(popt_npneu_avg[0], popt_npneu_avg[1]))

ax1.set_xscale('log', base=2)
ax1.set_yscale('log', base=10)
ax1.set_xlabel('N (log2)')
ax1.set_ylabel('Error (log10)')
ax1.set_title('Errors (avg) with Neumann conditions')
ax1.legend()
ax1.grid(True, which="both", ls="--")

ax2.plot(N, nneu_max, 'o', label='all cells')
ax2.plot(N, npneu_max, 's', label='partial cells')
ax2.plot(N, nfneu_max, '^', label='full cells')

ax2.plot(N, power_law(N, *popt_npneu_max), '--', label='{:.2f}x^({:.2f})'.format(popt_npneu_max[0], popt_npneu_max[1]))

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=10)
ax2.set_xlabel('N (log2)')
ax2.set_ylabel('Error (log10)')
ax2.set_title('Errors (max) with Neumann conditions')
ax2.legend()
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.savefig('mesh_neumann.png')
~~~

We can also plot the residues. 
  
~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def power_law(x, a, b):
    return a * x**b


data = np.loadtxt('res_neu.dat')


residue_avg_all_cells = data[:, 0]
residue_max_all_cells = data[:, 1]
residue_avg_partial_cells = data[:, 2]
residue_max_partial_cells = data[:, 3]
residue_avg_full_cells = data[:, 4]
residue_max_full_cells = data[:, 5]
N = data[:, 6]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

popt_residue_avg_all_cells, _ = curve_fit(power_law, N, residue_avg_all_cells)
popt_residue_avg_partial_cells, _ = curve_fit(power_law, N, residue_avg_partial_cells)
popt_residue_avg_full_cells, _ = curve_fit(power_law, N, residue_avg_full_cells)

popt_residue_max_all_cells, _ = curve_fit(power_law, N, residue_max_all_cells)
popt_residue_max_partial_cells, _ = curve_fit(power_law, N, residue_max_partial_cells)
popt_residue_max_full_cells, _ = curve_fit(power_law, N, residue_max_full_cells)

ax1.plot(N, residue_avg_all_cells, 'o', label='all cells')
ax1.plot(N, residue_avg_partial_cells, 's', label='partial cells')
ax1.plot(N, residue_avg_full_cells, '^', label='full cells')

ax1.plot(N, power_law(N, *popt_residue_avg_all_cells), '-', label='{:.2f}x^({:.2f})'.format(popt_residue_avg_all_cells[0], popt_residue_avg_all_cells[1]))
ax1.plot(N, power_law(N, *popt_residue_avg_partial_cells), '--', label='{:.2f}x^({:.2f})'.format(popt_residue_avg_partial_cells[0], popt_residue_avg_partial_cells[1]))
ax1.plot(N, power_law(N, *popt_residue_avg_full_cells), '-.', label='{:.2f}x^({:.2f})'.format(popt_residue_avg_full_cells[0], popt_residue_avg_full_cells[1]))

ax1.set_xscale('log', base=2)
ax1.set_yscale('log', base=10)
ax1.set_xlabel('N (log2)')
ax1.set_ylabel('Residue (log10)')
ax1.set_title('Residues (avg) with Neumann conditions')
ax1.legend()
ax1.grid(True, which="both", ls="--")

ax2.plot(N, residue_max_all_cells, 's', label='all cells')
ax2.plot(N, residue_max_partial_cells, 'o', label='partial cells')
ax2.plot(N, residue_max_full_cells, '^', label='full cells')

ax2.plot(N, power_law(N, *popt_residue_max_partial_cells), '--', label='{:.2f}x^({:.2f})'.format(popt_residue_max_partial_cells[0], popt_residue_max_partial_cells[1]))
ax2.plot(N, power_law(N, *popt_residue_max_full_cells), '-.', label='{:.2f}x^({:.2f})'.format(popt_residue_max_full_cells[0], popt_residue_max_full_cells[1]))

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=10)
ax2.set_xlabel('N (log2)')
ax2.set_ylabel('Residue (log10)')
ax2.set_title('Residues (max) with Neumann conditions')
ax2.legend()
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.savefig('res_neumann.png')

~~~

## Mesh convergence with Dirichlet conditions

We can do the same with Dirichlet boundary conditions : 

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_law(x, a, b):
    return a * x**b
    
data = np.loadtxt('errors_dir.dat')

ndir_avg = data[:, 0]
ndir_max = data[:, 1]
npdir_avg = data[:, 2]
npdir_max = data[:, 3]
nfdir_avg = data[:, 4]
nfdir_max = data[:, 5]
N = data[:, 6]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

popt_ndir_avg, _ = curve_fit(power_law, N, ndir_avg)
popt_npdir_avg, _ = curve_fit(power_law, N, npdir_avg)
popt_nfdir_avg, _ = curve_fit(power_law, N, nfdir_avg)

popt_ndir_max, _ = curve_fit(power_law, N, ndir_max)
popt_npdir_max, _ = curve_fit(power_law, N, npdir_max)
popt_nfdir_max, _ = curve_fit(power_law, N, nfdir_max)

ax1.plot(N, ndir_avg, 'o', label='all cells')
ax1.plot(N, npdir_avg, 's', label='partial cells')
ax1.plot(N, nfdir_avg, '^', label='full cells')

ax1.plot(N, power_law(N, *popt_ndir_avg), '-', label='{:.2f}x^({:.2f})'.format(popt_ndir_avg[0], popt_ndir_avg[1]))
ax1.plot(N, power_law(N, *popt_npdir_avg), '--', label='{:.2f}x^({:.2f})'.format(popt_npdir_avg[0], popt_npdir_avg[1]))

ax1.set_xscale('log', base=2)
ax1.set_yscale('log', base=10)
ax1.set_xlabel('N (log2)')
ax1.set_ylabel('Error (log10)')
ax1.set_title('Errors (avg) with Dirichlet conditions')
ax1.legend()
ax1.grid(True, which="both", ls="--")

ax2.plot(N, ndir_max, 'o', label='all cells')
ax2.plot(N, npdir_max, 's', label='partial cells')
ax2.plot(N, nfdir_max, '^', label='full cells')

ax2.plot(N, power_law(N, *popt_npdir_max), '--', label='{:.2f}x^({:.2f})'.format(popt_npdir_max[0], popt_npdir_max[1]))

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=10)
ax2.set_xlabel('N (log2)')
ax2.set_ylabel('Error (log10)')
ax2.set_title('Errors (max) with Dirichlet conditions')
ax2.legend()
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.savefig('mesh_dirichlet.png')

~~~

Same for residue : 

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_law(x, a, b):
    return a * x**b

data = np.loadtxt('res_dir.dat')

residue_avg_all_cells = data[:, 0]
residue_max_all_cells = data[:, 1]
residue_avg_partial_cells = data[:, 2]
residue_max_partial_cells = data[:, 3]
residue_avg_full_cells = data[:, 4]
residue_max_full_cells = data[:, 5]
N = data[:, 6]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

popt_residue_avg_all_cells, _ = curve_fit(power_law, N, residue_avg_all_cells)
popt_residue_avg_partial_cells, _ = curve_fit(power_law, N, residue_avg_partial_cells)
popt_residue_avg_full_cells, _ = curve_fit(power_law, N, residue_avg_full_cells)

popt_residue_max_all_cells, _ = curve_fit(power_law, N, residue_max_all_cells)
popt_residue_max_partial_cells, _ = curve_fit(power_law, N, residue_max_partial_cells)
popt_residue_max_full_cells, _ = curve_fit(power_law, N, residue_max_full_cells)

ax1.plot(N, residue_avg_all_cells, 'o', label='all cells')
ax1.plot(N, residue_avg_partial_cells, 's', label='partial cells')
ax1.plot(N, residue_avg_full_cells, '^', label='full cells')

ax1.plot(N, power_law(N, *popt_residue_avg_all_cells), '-', label='{:.2f}x^({:.2f})'.format(popt_residue_avg_all_cells[0], popt_residue_avg_all_cells[1]))
ax1.plot(N, power_law(N, *popt_residue_avg_partial_cells), '--', label='{:.2f}x^({:.2f})'.format(popt_residue_avg_partial_cells[0], popt_residue_avg_partial_cells[1]))
ax1.plot(N, power_law(N, *popt_residue_avg_full_cells), '-.', label='{:.2f}x^({:.2f})'.format(popt_residue_avg_full_cells[0], popt_residue_avg_full_cells[1]))

ax1.set_xscale('log', base=2)
ax1.set_yscale('log', base=10)
ax1.set_xlabel('N (log2)')
ax1.set_ylabel('Residue (log10)')
ax1.set_title('Residues (avg) with Dirichlet conditions')
ax1.legend()
ax1.grid(True, which="both", ls="--")

ax2.plot(N, residue_max_all_cells, 's', label='all cells')
ax2.plot(N, residue_max_partial_cells, 'o', label='partial cells')
ax2.plot(N, residue_max_full_cells, '^', label='full cells')

ax2.plot(N, power_law(N, *popt_residue_max_partial_cells), '--', label='{:.2f}x^({:.2f})'.format(popt_residue_max_partial_cells[0], popt_residue_max_partial_cells[1]))
ax2.plot(N, power_law(N, *popt_residue_max_full_cells), '-.', label='{:.2f}x^({:.2f})'.format(popt_residue_max_full_cells[0], popt_residue_max_full_cells[1]))

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=10)
ax2.set_xlabel('N (log2)')
ax2.set_ylabel('Residue (log10)')
ax2.set_title('Residues (max) with Dirichlet conditions')
ax2.legend()
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.savefig('residue_dir.png')

~~~

*/