/**
# Embedded boundary conditions

We make sure here that the proper boundary conditions are used. */

#include "../myembed.h"
#include "utils.h"
#include "view.h"

/**
## Exact solution

We define here the exact solution, evaluated at the center of each
cell. */

static double exact (double x, double y)
{
  return x - 2.*y;
}

double exact_gradient (Point point, double xc, double yc)
{
  coord p, n;
  embed_geometry (point, &p, &n);
  
  return n.x*(x) - n.y*(2.*y);
}

/**
We also define the shape of the domain. */

#if GEOM == 1 // concave
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = r - (0.25 + 0.05*cos(6.*theta));
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}
#elif GEOM == 2 // convex
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = (0.25 + 0.05*cos(6.*theta)) - r;
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);  
}
#elif GEOM == 3 // concave
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = r - (0.3 + 0.15*cos(6.*theta));
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);  
}
#else // convex
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = (0.3 + 0.15*cos(6.*theta)) - r;
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);  
}
#endif // GEOM

#if TREE
/**
When using *TREE*, we try to increase the accuracy of the restriction
operation in pathological cases by defining the gradient of *a* at the
center of the cell. */

void a_embed_gradient (Point point, scalar s, coord * g)
{
  g->x = 1.;
  g->y = -2.;
}
#endif // TREE

/**
## Setup */

int lvl;

int main()
{
  /**
  The domain is $1\times 1$. */
  
  size (1 [0]);
  origin (-L0/2., -L0/2.);

  for (lvl = 6; lvl <= 11; lvl++) { // minlevel = 4

    /**
    We initialize the grid. */

    N = 1 << (lvl);
    init_grid (N);
    
    /**
    We initialize the embedded boundary. */

#if TREE
    /**
    When using *TREE*, it is necessary to modify the refinement and
    prolongation of the volume and face fractions to account for
    embedded boundaries. In this case, using *embed_fraction_refine*
    for prolongation (as well as refine) seems to overall improve the
    accuracy of the results. */

    cs.refine = embed_fraction_refine;
    cs.prolongation = fraction_refine;

    foreach_dimension()
      fs.x.prolongation = embed_face_fraction_refine_x;
#endif // TREE
    
    cm   = cs;
    fm   = fs;
    csm1 = cs;

#if TREE
    /**
    When using *TREE*, we refine the mesh around the embedded
    boundary. */

    astats ss;
    int ic = 0;
    do {
      ic++;
      p_shape (cs, fs);
      ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			  maxlevel = (lvl), minlevel = (lvl) - 2);
    } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
    
    p_shape (cs, fs);

    /**
    We define *a* and set its boundary conditions. */
    
    scalar a[];
    foreach()
      a[] = cs[] > 0. ? exact (x, y) : nodata;

#if GEOM == 1 || GEOM == 3
    a[left]   = exact (x - Delta/2., y);
    a[right]  = exact (x + Delta/2., y);
    a[top]    = exact (x, y + Delta/2.);
    a[bottom] = exact (x, y - Delta/2.);
#endif // GEOM == 1 || GEOM == 3

#if TREE
    /**
    On *TREE*, we also modify the prolongation an restriction
    operators for *a*. */
    
    a.refine = a.prolongation = refine_embed_linear;
    a.restriction = restriction_embed_linear;
#endif // TREE
    
    /**
    We use "third-order" [face flux interpolation](/src/embed.h). */

#if ORDER2
    a.third = false;
#else
    a.third = true;
#endif // ORDER2

#if TREE
    /**
    We also define the gradient of *a* at the full cell center of
    cut-cells. */

    a.embed_gradient = a_embed_gradient;
#endif // TREE

    a[embed] = dirichlet (exact (x, y));

    foreach() {
      if (cs[] > 0. && cs[] < 1.) {
	coord b, n;
        embed_geometry (point, &b, &n);
	bool dirichlet;
        double ab = (a.boundary[embed] (point, point, a, &dirichlet));
	assert (dirichlet);
	/* fprintf (stderr, "%g %.10f %.10f\n", cs[], ab, exact (x + Delta*b.x, y + Delta*b.y)), fflush (stderr); */
	assert (ab == exact (x + Delta*b.x, y + Delta*b.y));
      }
    }

    a[embed] = neumann (exact_gradient (point, x, y));

    foreach() {
      if (cs[] > 0. && cs[] < 1.) {
	coord b, n;
        embed_geometry (point, &b, &n);
	bool dirichlet;
        double ab = (a.boundary[embed] (point, point, a, &dirichlet));
	assert (!dirichlet);
	/* fprintf (stderr, "%g %.10f %.10f\n", cs[], ab, exact_gradient (point, x + Delta*b.x, y + Delta*b.y)), fflush (stderr); */
	assert (ab == exact_gradient (point, x + Delta*b.x, y + Delta*b.y));
      }
    }
  }
}
