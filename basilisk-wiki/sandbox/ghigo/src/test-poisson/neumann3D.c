/**
# Poisson equation on 3D complex domains with Dirichlet and Neumann boundary conditions

We reproduce a test case proposed by [Schwarz et al.,
2006](#schwarz2006). The following Poisson equation is solved:
$$
\Delta a = -14 f, \quad f(x,y,z) = \sin(x)\sin(2y)\sin(3z)
$$
on the domains:
$$
\Omega_1 = \left\{(r): r < 0.392\right\},
$$
$$
\Omega_2 = \left\{(r): r > 0.392\right\}.
$$
The exact solution is:
$$
a(x,y,z) = f.
$$
*/

#include "../myembed.h"
#include "../mypoisson.h"
#include "view.h"

/**
## Exact solution

We define here the exact solution. */

static double exact (double x, double y, double z)
{
  return sin(x)*sin(2.*y)*sin(3.*z);
}

/**
Next, we define the gradient of the exact solution. */

double exact_gradient (Point point, double xc, double yc, double zc)
{
  coord p, n;
  embed_geometry (point, &p, &n);
  
  return (n.x*cos(xc)*sin(2.*yc)*sin(3.*zc) +
	  n.y*2.*sin(xc)*cos(2.*yc)*sin(3.*zc) +
	  n.z*3.*sin(xc)*sin(2.*yc)*cos(3.*zc));
}

/**
We also define the shape of the domain. */

#if GEOM // concave
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = sq(x) + sq(y) + sq(z) - sq(0.392);
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
    phi[] = sq(0.392) - sq(x) - sq(y) - sq(z);
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}
#endif // GEOM

#if GEOM
/**
In the case of the domain $\Omega_2$, the poisson problem requires the
homogeneous equivalent of the boundary conditions imposed on the
non-embedded boundaries of the domain. We therefore define a function
for homogeneous Dirichlet condition as we are not using the
*dirichlet* function to impose the boundary conditions for *a* on the
domain boundaries. */

static double dirichlet_homogeneous_bc (Point point, Point neighbor,
					scalar s, void * data)
{
  return 0.;
}
#endif // GEOM

#if TREE
/**
When using *TREE*, we try to increase the accuracy of the restriction
operation in pathological cases by defining the gradient of *a* at the
center of the cell. */

void a_embed_gradient (Point point, scalar s, coord * g)
{
  g->x = cos(x)*sin(2.*y)*sin(3.*z);
  g->y = 2.*sin(x)*cos(2.*y)*sin(3.*z);
  g->z = 3.*sin(x)*sin(2.*y)*cos(3.*z);
}
#endif // TREE

/**
## Setup */

int lvl;

int main()
{
  /**
  The domain is $1\times 1\times 1$. */

  size (1 [0]);
  origin (-L0/2., -L0/2., -L0/2.);

  for (lvl = 4; lvl <= 8; lvl++) { // minlevel = 2 (3.1pt/D)

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
      a[] = cs[] > 0. ? exact (x, y, z) : nodata;

#if GEOM
    a[left]   = exact (x - Delta/2., y, z);
    a.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
    a[right]  = exact (x + Delta/2., y, z);
    a.boundary_homogeneous[right] = dirichlet_homogeneous_bc;
    a[top]    = exact (x, y + Delta/2., z);
    a.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
    a[bottom] = exact (x, y - Delta/2., z);
    a.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;
    a[front]  = exact (x, y, z + Delta/2.);
    a.boundary_homogeneous[front] = dirichlet_homogeneous_bc;
    a[back]   = exact (x, y, z - Delta/2.);
    a.boundary_homogeneous[back] = dirichlet_homogeneous_bc;
#endif // GEOM    

#if DIRICHLET
    a[embed] = dirichlet (exact (x, y, z));
#else
    a[embed] = neumann (exact_gradient (point, x, y, z));
#endif

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

#if TREE && !NGRAD
    /**
    We also define the gradient of *a* at the full cell center of
    cut-cells. */

    a.embed_gradient = a_embed_gradient;
#endif // TREE
    
    /**
    We define the right hand-side *b*:
    $$
    \Delta a = -14 f
    $$
    using the coordinates of the barycenter of the cut cell (xc,yc),
    which is calculated from the cell and surface fractions. */
    
    scalar b[];

#if TREE
    /**
    In theory for *TREE*, we should also use embed-specific
    prolongation and restriction operators for *b*. However, *b* is
    only used to compute the residual in leaf cells and its stencil is
    never accessed. Furtheremore, and more importantly, if a *res*
    variable is not given as a parameter of the *poisson* function
    (which is the case here), than the prolongation and restriction
    operators for the *res* variable used in *mg_solve* are copied
    from those of *b*. However, for the multigrid solver to work
    properly, the prolongation and restriction operators for the *res*
    variable should be similar to the multigrid prolongation and
    restriction operators. We therefore use the default operators for
    *b* here. */
    
    /* b.refine = b.prolongation = refine_embed_linear; */
    /* b.restriction = restriction_embed_linear; */
#endif // TREE
    
    foreach() {
      double xc = x, yc = y, zc = z;
      if (cs[] > 0. && cs[] < 1.) {
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	plane_center (n, alpha, cs[], &p); // Different from line_area_center
	xc += p.x*Delta, yc += p.y*Delta, zc += p.z*Delta;
      }
      b[] = - 14.*exact (xc, yc, zc)*cs[];
    }
    
    /**
    We compute the trunction error, or residual. */
      
    struct Poisson p;
    p.alpha = fs;
    p.lambda = zeroc;
    scalar res[];
    residual ({a}, {b}, {res}, &p);
    double maxp = 0., maxc = 0., maxf = 0.;
    foreach(reduction(max:maxf) reduction(max:maxc) reduction(max:maxp)) {
      if (cs[] == 1. && fabs(res[]) > maxf)
	maxf = fabs(res[]);
      if (cs[] > 0. && cs[] < 1. && fabs(res[]) > maxp) {
	maxp = fabs(res[]);
	maxc = cs[];
      }
    }
    fprintf (stdout, "%d %.3g %.3g %.3g\n", N, maxf, maxp, maxp*(maxc + SEPS));
    fflush (stdout);
    
    /**
    We solve the Poisson equation. */

    mgstats s = poisson (a, b, alpha = fs, tolerance = 1e-6, minlevel = 1);
    
    /**
    The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
    fields and their norms are computed. */
    
    scalar e[], ep[], ef[];
    foreach() {
      if (cs[] == 0.)
	ep[] = ef[] = e[] = nodata;
      else {
	e[] = fabs (a[] - exact (x, y, z));
	ep[] = cs[] < 1. ? e[] : nodata;
	ef[] = cs[] >= 1. ? e[] : nodata;
      }
    }
    norm n = normf (e), np = normf (ep), nf = normf (ef);
    fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %d %d\n",
	     N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max, s.i, s.nrelax);
    fflush (stderr);
    
    /**
    The solution and the error are displayed using bview. */

    if (lvl == 6) {
      clear ();
      view (fov = 32.2073, quat = {-0.309062,0.243301,0.0992085,0.914026},
	    tx = 0.0122768, ty = 0.0604286, bg = {1,1,1},
	    width = 800, height = 800);

      box ();
      draw_vof("cs", "fs", color = "e");
      save ("e.png");

      cells (n = {0, 0, 1}, alpha = -0.392/2.);
      squares ("e", n = {0, 1, 0}, alpha = -0.392/2., spread = -1);
      save ("e-slice.png");

      box ();
      draw_vof("cs", "fs", color = "a");
      save ("a.png");

      cells (n = {0, 0, 1}, alpha = -0.392/2.);
      squares ("a", n = {0, 1, 0}, alpha = -0.392/2., spread = -1);
      save ("a-slice.png");
    }
  }
}

/**
## Results for the default domain $\Omega_1$

![Solution for *l=6*](neumann3D/a.png)

![Solution for *l=6*](neumann3D/a-slice.png)

![Error for *l=6*](neumann3D/e.png)

![Error for *l=6*](neumann3D/e-slice.png)

#### Errors

~~~gnuplot Maximum truncation error (residual)
set terminal svg font ",16"
set key bottom left spacing 1.1
set xtics 8,4,512
set grid ytics
set ytics format "%.0e" 1.e-12,100,1.e2
set xlabel 'N'
set ylabel '||res||_{inf}'
set xrange [8:512]
set yrange [1.e-6:1.e1]
set logscale
plot 'out' u 1:3 w p ps 1.25 pt 7 lc rgb "black" t 'cut-cells', \
     ''    u 1:2 w p ps 1.25 pt 5 lc rgb "blue" t 'full cells'
~~~

~~~gnuplot Average error convergence
set ylabel '||error||_{1}'
set yrange [1.e-9:1.e-2]
plot 'log' u 1:4 w p ps 1.25 pt 7 lc rgb "black" t 'cut-cells', \
     '' u 1:6 w p ps 1.25 pt 5 lc rgb "blue" t 'full cells', \
     '' u 1:2 w p ps 1.25 pt 2 lc rgb "red" t 'all cells'
~~~

~~~gnuplot Maximum error convergence
set ylabel '||error||_{inf}'
plot '' u 1:5 w p ps 1.25 pt 7 lc rgb "black" t 'cut-cells', \
     '' u 1:7 w p ps 1.25 pt 5 lc rgb "blue" t 'full cells', \
     '' u 1:3 w p ps 1.25 pt 2 lc rgb "red" t 'all cells'
~~~

#### Order of convergence

Following [Johansen and Collela, 1998](#johansen1998), [Schwarz et al.,
2006](#schwarz2006) determined the following
estimates for the residual, or truncation error:

* for full cells, even though gradients at the faces of the cell are
estimated with an $O\left(\Delta^2\right)$ accuracy, the standard
cancellation of the error when using cell-centered values to
discretize the Laplacian gives:
$$
\mathrm{res}_f \sim \Delta^2
$$

* for cut-cells, because the gradients on the faces of the cell and
the embedded boundary are estimated with an $O\left(\Delta^2\right)$
accuracy, but the error doesn't cancel out in the Laplacian, we have:
$$
\mathrm{res}_p \sim \frac{\Delta}{c_s}
$$

We recover here the expected order of convergence for full cells only,
but have an order of convergence lower than 1 for cut-cells. The same
behaviour is observed in the equivalent 2D test case when we do not
account for the neighboring face fractions when computing the fluxes
on the faces of each cut-cell.

~~~gnuplot Order of convergence of the maximum residual
reset
set terminal svg font ",16"
set key bottom right spacing 1.1
set xtics 8,4,512
set ytics -4,1,4
set grid ytics
set xlabel 'N'
set ylabel 'Order of ||res||_{inf}'
set xrange [8:512]
set yrange [-4:4.5]
set logscale x

# Asymptotic order of convergence

ftitle(b) = sprintf(", asymptotic order n=%4.2f", -b);
Nmin = log(128);

f1(x) = a1 + b1*x; # cut-cells
f2(x) = a2 + b2*x; # full cells
f3(x) = a3 + b3*x; # all cells

fit [Nmin:*][*:*] f1(x) '< sort -k 1,1n out' u (log($1)):(log($3)) via a1,b1; # cut-cells 
fit [Nmin:*][*:*] f2(x) '< sort -k 1,1n out' u (log($1)):(log($2)) via a2,b2; # full cells

plot '< sort -k 1,1n out | awk -f ../data/order-res.awk' u 1:3 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '< sort -k 1,1n out | awk -f ../data/order-res.awk' u 1:2 w lp ps 1.25 pt 5 lc rgb "blue" t 'full cells'.ftitle(b2)
~~~

[Schwarz et al., 2006](#schwarz2006) also determined the following
estimates for the error of the solution of the Poisson problem:

* for full cells:
$$
\mathrm{err}_f \sim \Delta^2
$$

* for cut-cells:
$$
\mathrm{err}_p \sim \Delta^2
$$

We recover here the expected orders of convergence. Note that when
using mesh adaptation, the order of convergence of both the residual
and the error are determined by the order of the
prolongataion/refinement functions, which are second-order accurate in
Basilisk. However, in the complex "convex" geometry presented here,
the restriction operation requires the user to input an
*embed_gradient* to maintain second-order accuracy.

~~~gnuplot Order of convergence of the average error
set ylabel 'Order of ||error||_{1}'

# Asymptotic order of convergence

fit [Nmin:*][*:*] f1(x) '< sort -k 1,1n log' u (log($1)):(log($4)) via a1,b1; # cut-cells 
fit [Nmin:*][*:*] f2(x) '< sort -k 1,1n log' u (log($1)):(log($6)) via a2,b2; # full-cells
fit [Nmin:*][*:*] f3(x) '< sort -k 1,1n log' u (log($1)):(log($2)) via a3,b3; # all cells

plot '< sort -k 1,1n log | awk -f ../data/order.awk' u 1:4 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '< sort -k 1,1n log | awk -f ../data/order.awk' u 1:6 w lp ps 1.25 pt 5 lc rgb "blue" t 'full cells'.ftitle(b2), \
     '< sort -k 1,1n log | awk -f ../data/order.awk' u 1:2 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3)
~~~

~~~gnuplot Order of convergence of the maximum error
set ylabel 'Order of ||error||_{inf}'

# Asymptotic order of convergence

fit [Nmin:*][*:*] f1(x) '< sort -k 1,1n log' u (log($1)):(log($5)) via a1,b1; # cut-cells 
fit [Nmin:*][*:*] f2(x) '< sort -k 1,1n log' u (log($1)):(log($7)) via a2,b2; # full-cells
fit [Nmin:*][*:*] f3(x) '< sort -k 1,1n log' u (log($1)):(log($3)) via a3,b3; # all cells

plot '< sort -k 1,1n log | awk -f ../data/order.awk' u 1:5 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '< sort -k 1,1n log | awk -f ../data/order.awk' u 1:7 w lp ps 1.25 pt 5 lc rgb "blue" t 'full cells'.ftitle(b2), \
     '< sort -k 1,1n log | awk -f ../data/order.awk' u 1:3 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3)
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

@article{schwartz2006,
  title={A Cartesian grid embedded boundary method for the heat equation 
  and Poisson's equation in three dimensions},
  author={Schwartz, Peter and Barad, Michael and Colella, Phillip and Ligocki, 
  Terry},
  journal={Journal of Computational Physics},
  volume={211},
  number={2},
  pages={531--550},
  year={2006},
  publisher={Elsevier},
  url={https://cloudfront.escholarship.org/dist/prd/content/qt0fp606kk/qt0fp606kk.pdf}
}
~~~
*/
