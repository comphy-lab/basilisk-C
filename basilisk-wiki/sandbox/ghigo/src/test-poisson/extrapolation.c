/**
# 2D extrapolation 

We assess here the order of the extrapolation function
*embed_extrapolate()*. */

#include "../myembed.h"
#include "utils.h"
#include "view.h"

/**
## Exact solution

We define here the exact solution, evaluated at the center of each
cell. The solution is quadratic so the extrapolation should be
second-order when using *embed_extrapolate()* and exact with
*embed_extrapolate_ls()*. When using mesh adaptation, the
extrapolation should be second-order since the
restriction/prolongation functions are second-order. */

static double exact (double x, double y)
{
  return pi + x - 2.*y + 3.*x*x - 4.*x*y + 5.*y*y;
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
  g->x = 1. + 6.*x - 4.*y;
  g->y = -2. - 4.*x + 10.*y;
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

    a[embed] = dirichlet (exact (x, y));

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

    /**
    We now compute *a* in each cut-cell using the *embed_extrapolate*
    function. We perform several iterations as the cut-cells are
    initialized with the exact solution. */

    int ni = 0;
    while (ni < 10) {
      ni ++;
      foreach() {
	if (cs[] > 0. && cs[] < 1.) {
	  // Cell centroid, barycenter and normal of the embedded fragment
	  coord c = {0., 0.}, b, n;
	  embed_geometry (point, &b, &n);
#if LS
	  a[] = embed_extrapolate_ls (point, a, cs, c, false);
#else
	  bool dirichlet;
	  double ab = (a.boundary[embed] (point, point, a, &dirichlet));
	  assert (dirichlet);
	  a[] = embed_extrapolate (point, a, cs, n, c, ab);
#endif // LS
	}
      }
    }
        
    /**
    The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
    fields and their norms are computed. */
    
    scalar e[], ep[], ef[];
    foreach() {
      if (cs[] == 0.)
	ep[] = ef[] = e[] = nodata;
      else {
	e[] = fabs (a[] - exact (x, y));
	ep[] = cs[] < 1. ? e[] : nodata;
	ef[] = cs[] >= 1. ? e[] : nodata;
      }
    }
    norm n = normf (e), np = normf (ep), nf = normf (ef);
    fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g\n",
	     N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max);
    fflush (stderr);
    
    /**
    The solution and the error are displayed using bview. */

    if (lvl == 7) {
      clear ();
      view (fov = 18);
      
      squares ("a", spread = -1);
      draw_vof ("cs", "fs");
      save ("a.png");
      
      squares ("e", spread = -1);
      draw_vof ("cs", "fs");
      save ("e.png");

      cells ();
      draw_vof ("cs", "fs");
      save ("mesh.png");
    }
  }
}

/**
## Results for for the default domain $\Omega_4$

![Solution for *l=7*](extrapolation/a.png)

![Error for *l=7*](extrapolation/e.png)

![Mesh for *l=7*](extrapolation/mesh.png)

#### Errors

~~~gnuplot Average error convergence
set terminal svg font ",16"
set key top right spacing 1.1
set xtics 32,4,4096
set grid ytics
set ytics format "%.0e" 1.e-20,10000,1.e2 
set xlabel 'N'
set ylabel '||error||_{1}'
set xrange [32:4096]
set yrange [1.e-20:1.e2]
set logscale
plot 'log' u 1:4 w p ps 1.25 pt 7 lc rgb "black" t 'cut-cells', \
     ''    u 1:6 w p ps 1.25 pt 5 lc rgb "blue" t 'full cells', \
     ''    u 1:2 w p ps 1.25 pt 2 lc rgb "red" t 'all cells'
~~~

~~~gnuplot Maximum error convergence
set ylabel '||error||_{inf}'
plot '' u 1:5 w p ps 1.25 pt 7 lc rgb "black" t 'cut-cells', \
     '' u 1:7 w p ps 1.25 pt 5 lc rgb "blue" t 'full cells', \
     '' u 1:3 w p ps 1.25 pt 2 lc rgb "red" t 'all cells'
~~~

#### Order of convergence

~~~gnuplot Order of convergence of the  average error
reset
set terminal svg font ",16"
set key bottom right spacing 1.1
set xtics 32,4,4096
set ytics -4,1,4
set grid ytics
set xlabel 'N'
set ylabel 'Order of ||error||_{1}'
set xrange [32:4096]
set yrange [-4:4]
set logscale x

# Asymptotic order of convergence

ftitle(b) = sprintf(", asymptotic order n=%4.2f", -b);
Nmin = log(1024);

f1(x) = a1 + b1*x; # cut-cells
f3(x) = a3 + b3*x; # all cells

fit [Nmin:*][*:*] f1(x) '< sort -k 1,1n log' u (log($1)):(log($4)) via a1,b1; # cut-cells 
fit [Nmin:*][*:*] f3(x) '< sort -k 1,1n log' u (log($1)):(log($2)) via a3,b3; # all cells

plot '< sort -k 1,1n log | awk -f ../data/order-extrapolation.awk' u 1:4 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '< sort -k 1,1n log | awk -f ../data/order-extrapolation.awk' u 1:2 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3)
~~~

~~~gnuplot Order of convergence of the maximum error
set ylabel 'Order of ||error||_{inf}'

# Asymptotic order of convergence

fit [Nmin:*][*:*] f1(x) '< sort -k 1,1n log' u (log($1)):(log($5)) via a1,b1; # cut-cells 
fit [Nmin:*][*:*] f3(x) '< sort -k 1,1n log' u (log($1)):(log($3)) via a3,b3; # all cells

plot '< sort -k 1,1n log | awk -f ../data/order-extrapolation.awk' u 1:5 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '< sort -k 1,1n log | awk -f ../data/order-extrapolation.awk' u 1:3 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3)
~~~

## References

~~~bib
@article{johansen1998,
  title={A Cartesian grid embedded boundary method for Poisson's
  equation on irregular domains},
  author={Johansen, Hans and Colella, Phillip},
  journal={Journal of Computational Physics},
  volume={147},
  number={1},
  pages={60--85},
  year={1998},
  publisher={Elsevier},
  url={https://pdfs.semanticscholar.org/17cd/babecd054d58da05c2ba009cccb3c687f58f.pdf}
}
~~~
*/
