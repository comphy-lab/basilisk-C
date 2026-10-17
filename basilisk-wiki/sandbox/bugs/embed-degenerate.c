/**
# Poisson solver with degenerate embedded Dirichlet cells

For some cut cells, the third-order Dirichlet gradient of
[dirichlet_gradient()](/src/embed.h#dirichlet_gradient) cannot be
computed using the neighboring cells (they are not in the fluid). The
gradient then depends on the value of the cell itself and
[embed_flux()](/src/embed.h#embed_flux) returns a non-zero coefficient
$e$ such that the embedded flux is $c + e\,a[]$.

This coefficient must be taken into account consistently by both the
relaxation and the residual functions of the
[Poisson solver](/src/poisson.h). This test checks this using a
"star" with sharp branches, which contains degenerate cells.

The solution is trivial: with $a=1$ on the embedded boundary and on the
domain boundaries and a zero right-hand-side, the exact solution is
$a=1$ everywhere. The residual of the exact solution must then be zero
(to round-off) and the multigrid solver must converge to it.

<div class="message">
The patch that fixes this issue was published in [this thread](https://groups.google.com/g/basilisk-fr/c/R_vdpFltUl8/m/Grp8A-JHAQAJ) on the user forum. </div>
*/

#include "embed.h"
#include "poisson.h"
#include "view.h"

int main()
{
  double maxres = 0., maxresa = 0.;
  for (N = 32; N <= 256; N *= 2) {
    size (1. [0]); // dimensionless
    origin (-L0/2., -L0/2.);
    init_grid (N);

    vertex scalar phi[];
    foreach_vertex() {
      double theta = atan2 (y, x), r = sqrt (sq(x) + sq(y));
      phi[] = r - (0.25 + 0.15*cos(6.*theta));
    }
    fractions (phi, cs, fs);
#if TREE
    cs.refine = cs.prolongation = fraction_refine;
#endif
    restriction ({cs, fs});

    cm = cs;
    fm = fs;

    scalar a[], b[];
    a[embed] = dirichlet (1.);
    a[left]   = dirichlet (1.);
    a[right]  = dirichlet (1.);
    a[top]    = dirichlet (1.);
    a[bottom] = dirichlet (1.);
    a.third = true;

    /** 
    We count the number of degenerate cells i.e. cells for which
    *embed_flux()* returns a non-zero coefficient. */

    int ndeg = 0;
    scalar deg[];
    foreach (reduction(+:ndeg)) {
      double c, e = embed_flux (point, a, fs, &c);
      deg[] = e ? 1. : cs[] > 0. ? 0. : nodata;
      if (e) 
        ndeg++;
    }   

    /**
    We compute the residual of the exact solution $a=1$. */

    foreach()
      a[] = 1., b[] = 0.;
    scalar res[];
    struct Poisson p;
    p.a = a, p.b = b, p.alpha = fs, p.lambda = zeroc;
    p.embed_flux = embed_flux;
    maxres = max (maxres, residual ({a}, {b}, {res}, &p));

    /**
    We solve the Poisson equation starting from $a=0$ and compute the
    maximum error. */

    foreach()
      a[] = 0.;
    mgstats s = poisson (a, b, alpha = fs, tolerance = 1e-6);
    maxresa = max (maxresa, s.resa);
    double maxerr = 0.;
    scalar err[];
    foreach (reduction(max:maxerr)) {
      err[] = cs[] > 0. ? a[] - 1. : nodata;
      if (cs[] > 0. && fabs (a[] - 1.) > maxerr)
        maxerr = fabs (a[] - 1.);
    }

    fprintf (stderr, "%d %d %d %.2g %.2g\n",
	     N, ndeg, s.i, s.resa, maxerr);

    /**
    We display the degenerate cells and the error for $N=64$, the
    resolution with the largest number of degenerate cells. */

    if (N == 64) {
      view (fov = 18);
      squares ("deg", spread = -1);
      draw_vof ("cs", "fs", lw = 2);
      cells();
      save ("deg.png");

      clear();
      squares ("err", spread = -1);
      draw_vof ("cs", "fs", lw = 2);
      save ("err.png");
    }
  }

  /**
  For all resolutions, the residual of the exact solution must vanish
  (to round-off) and the solver must converge. These checks are done
  after the loop so that the figures are also generated when the
  solver diverges. */

  assert (maxres < 1e-9);
  assert (maxresa < 1e-6);
}

/**
The columns are: the resolution, the number of degenerate cells, the
number of multigrid cycles, the final residual and the maximum error.

If the residual is not consistent with the relaxation function (i.e.
if $e$ is not accounted for with the same sign in both), the solver
diverges. 

The degenerate cells are located at the tips of the narrow fluid
wedges between the branches of the star, where the stencil of
*dirichlet_gradient()* does not fit within the fluid.

![Degenerate cells (red) in the fluid (blue) for $N=64$](embed-degenerate/deg.png)

With a consistent residual, the error is of the order of the tolerance
of the solver. Otherwise it is concentrated in the degenerate cells
and grows without bound.

![Error $a-1$ for $N=64$](embed-degenerate/err.png)
*/
