/**
# Curvature of a circular/spherical interface

This test evaluates the accuracy of the generalised height-function
curvature calculation. It is similar to the case presented in
[Popinet, 2009](/src/references.bib#popinet2009) (Figure 5). The
curvatures of circles/spheres with a randomised position and varying
radii are computed and statistics on the error are gathered and
displayed on the graph below. 

We use the
[Vofi](http://www.ida.upmc.fr/~zaleski/paris/Vofi-1.0.tar.gz) library
for accurate initialisation of interfacial shapes [Bna et al,
2015](/src/references.bib#bna2015). */

#include <vofi.h>
#pragma autolink -L$HOME/local/lib -lvofi

#include "fractions.h"
#include "curvature.h"

/**
We pass arguments to *Vofi* through these global variables. */

static double xc, yc, zc, Rc;

static double sphere (creal p[dimension])
{
#if dimension == 3
  return sq(p[0] - xc) + sq(p[1] - yc) + sq(p[2] - zc) - sq(Rc);
#else
  return sq(p[0] - xc) + sq(p[1] - yc) - sq(Rc);
#endif
}

static void vofi (scalar c, int levelmax)
{
  double fh = Get_fh (sphere, NULL, 1./(1 << levelmax), dimension, 0);
  foreach() {
    creal p[3] = {x - Delta/2., y - Delta/2., z - Delta/2.};
    c[] = Get_cc (sphere, p, Delta, fh, dimension);
  }
}

/**
The function below is called with different arguments for "coarse" and
"fine" interface resolution in order to minimize the computation
times. The number of randomised positions is *nr*, the radius of the
circle is *R*, the maximum refinement level to consider is *levelmax*,
and the statistics for each level are stored in the array *n*, while
the statistics on which method is used are stored in *sc*. */

void sample_circles (int nr, double R, int levelmax, norm * n, cstats * sc)
{
  while (nr--) {

    /**
    We refine the grid down to *levelmax* but only around the interface. */

    scalar c[], kappa[];
    c.refine = c.prolongation = fraction_refine;

    init_grid (1 << 5);
    xc = noise()/8., yc = noise()/8., zc = noise()/8., Rc = R;
    vofi (c, 5);
    for (int l = 6; l <= levelmax; l++) {
      refine (c[] > 0. && c[] < 1. && level < l);
      vofi (c, l);
    }
    
    /**
    We then successively coarsen this fine initial grid to compute the
    curvature on coarser and coarser grids (thus saving on the
    expensive initial condition). */

    restriction ({c});
    for (int l = levelmax; l >= 3; l--) {
      unrefine (level >= l || c[] <= 0. || c[] >= 1.);
      
      cstats s = curvature (c, kappa, sigma = 2.);
      
      /**
      We store statistics on the methods used for curvature
      computation... */

      sc[l].h += s.h; sc[l].f += s.f; sc[l].a += s.a; sc[l].c += s.c;
      foreach()
	if (c[] > 0. && c[] < 1.) {

          /**
	  ...and error statistics (for a given level of refinement *l*). */
           
	  double e = fabs(kappa[]/2. - (dimension - 1)/R)*R/(dimension - 1);
	  n[l].volume += dv();
	  n[l].avg += dv()*e;
	  n[l].rms += dv()*e*e;
	  if (e > n[l].max)
	    n[l].max = e;
	}
    }
  }  
}

int main()
{
  size (1. [0]); // dimensionless
  origin (-0.5, -0.5, -0.5);
  init_grid (N);
  
  /**
  We try a wide enough range of radii. */

  for (double R = 0.1; R <= 0.2; R *= 1.2) {

    /**
    We initialize the arrays required to store the statistics for each
    level of refinement. */

    int levelmax = 8;
    norm n[levelmax + 1];
    cstats sc[levelmax + 1];
    for (int i = 0; i <= levelmax; i++) {
      n[i].volume = n[i].avg = n[i].rms = n[i].max = 0;
      sc[i].h = sc[i].f = sc[i].a = sc[i].c = 0.;
    }

    /**
    We can limit randomisation for the higher resolutions (since we
    expect less "special cases" on fine meshes). We thus limit the
    total runtime by sampling many (100) locations on coarse meshes
    but only few (1) location on the finest mesh. */

#if dimension == 2
    sample_circles (1000, R, 4, n, sc);
    sample_circles (100, R, 6, n, sc);
    sample_circles (10, R, levelmax, n, sc);
#else // dimension == 3
    sample_circles (100, R, 4, n, sc);
    sample_circles (10, R, 6, n, sc);
    sample_circles (1, R, levelmax, n, sc);
#endif
    
    /**
    Finally we output the statistics for this particular radius and for
    each level of refinement. */

    for (int l = levelmax; l >= 3; l--)
      if (n[l].volume) {
	n[l].avg /= n[l].volume;
	n[l].rms = sqrt(n[l].rms/n[l].volume);
	double t = sc[l].h + sc[l].f + sc[l].a + sc[l].c;
	fprintf (stderr, "%g %g %g %g %g %g %g %g\n",
		 2.*R*(1 << l),
		 n[l].avg, n[l].rms, n[l].max,
		 sc[l].h/t, sc[l].f/t, sc[l].a/t, sc[l].c/t);
      }
  }

  /**
  At the end of the run, we sort the data by increasing order of
  diameter (in grid points). */

  fflush (stderr);
  system ("sort -k1,2 -n log > log.s; mv -f log.s log");
}

/**
The results are summarised in the figure below. There are two sets of
points: error norms (bottom three curves) and percentages of
curvatures computed with each method (top four curves). As expected
the second-order convergence of the max and RMS norms is recovered for
the pure HF method, for a diameter greater than about 15 grid
points. The RMS norm shows consistent second-order convergence across
the whole range of diameters. There is a significant degradation of
the max norm between 7 and 15 grid points, corresponding with the
introduction of the "averaging" and "(mixed) HF fit" methods. This
degradation is significantly less pronounced for the method
implemented in [Popinet, 2009](/src/references.bib#popinet2009).

The grading of the methods used for curvature calculation follows what
is expected: 

* exclusively HF for $D > 15$, 
* a combination of HF, nearest-neighbor average and "mixed HF fit" for
$3 < D < 15$
* and exclusively "centroids fit" for $D < 3$.

~~~gnuplot Relative curvature error as a function of resolution
set logscale
set grid
set key top right
set xlabel 'Diameter (grid points)'
set ylabel 'Relative curvature error / percentage' 
f(x)=(x > 0. ? 100.*x : 1e1000)
set yrange [:100]
plot 2./(x*x) t '2/x^{2}', 'log' u 1:4 w lp t 'Max', '' u 1:3 w lp t 'RMS', \
  '../popinet.csv' u ($1*2):2 w lp t 'Popinet (2009)', \
  'log' u 1:(f($5)) w lp t 'HF', '' u 1:(f($6)) w lp t 'HF fit', \
  '' u 1:(f($7)) w lp t 'Average', '' u 1:(f($8)) w lp t 'Centroids'
~~~

Similar results are obtained in three dimensions, but with a much
wider intermediate range of low-order convergence for the max-norm.

~~~gnuplot Relative curvature error as a function of resolution in three dimensions
set key bottom left
plot 2./(x*x) t '2/x^{2}', '../curvature.3D/log' u 1:4 w lp t 'Max',	\
  '' u 1:3 w lp t 'RMS', '../curvature.3D/log' u 1:(f($5)) w lp t 'HF',	\
  '' u 1:(f($6)) w lp t 'HF fit',					\
  '' u 1:(f($7)) w lp t 'Average', '' u 1:(f($8)) w lp t 'Centroids'
~~~
*/
