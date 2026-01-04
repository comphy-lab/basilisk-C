/**
# Poisson equation on complex domains in 3D

We reproduce the test cases initially proposed by [Schwarz et al.,
2006](#schwarz2006), section 4.1, with Dirichlet or Neumann boundary
conditions. */

#include "grid/multigrid3D.h"
#include "embed.h"
#include "poisson.h"
#include "view.h"

/**
The exact solution is given by
$$
a(x,y,z) = \sin(x) \sin(2y) \sin(3z)
$$
*/

static double exact (double x, double y, double z) {
  return sin(x)*sin(2.*y)*sin(3.*z);
}

double exact_gradient (Point point, double xc, double yc, double zc)
{
  coord n = facet_normal (point, cs, fs);
  normalize (&n);
  return (n.x*cos(xc)*sin(2.*yc)*sin(3.*zc) +
	  n.y*2.*sin(xc)*cos(2.*yc)*sin(3.*zc) +
	  n.z*3.*sin(xc)*sin(2.*yc)*cos(3.*zc));
}

int main()
{
  for (N = 32; N <= 256; N *= 2) {
    size (1. [0]); // dimensionless space
    origin (-L0/2., -L0/2., -L0/2.);
    init_grid (N);

    /**
    The domain is a sphere of radius 0.392 centered at the origin. */

    solid (cs, fs, sq(0.392) - sq(x) - sq(y) - sq(z));
#if TREE
    cs.refine = cs.prolongation = fraction_refine;
#endif
    restriction ({cs,fs});

    cm = cs;
    fm = fs;
    
    scalar a[], b[];
    
    /**
    The boundary conditions on the embedded boundary are Dirichlet and
    Neumann, respectively. */

#if DIRICHLET
    a[embed] = dirichlet (exact (x, y, z));
#else
    a[embed] = neumann (exact_gradient (point, x, y, z));
#endif

    /**
    We use "third-order" [face flux interpolation](/src/embed.h). */
    
    a.third = true;
    
    /**
    The right-hand-side
    $$
    \Delta a = -14 a
    $$
    is defined using the coordinates of the
    barycenter of the cut cell (xc,yc), which is calculated from the
    cell and surface fractions. */
    
    foreach() {
      a[] = cs[] > 0. ? exact (x, y, z) : nodata;
      
      double xc = x, yc = y, zc = z;
      if (cs[] > 0. && cs[] < 1.) {
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	plane_center (n, alpha, cs[], &p);
	xc += p.x*Delta, yc += p.y*Delta, zc += p.z*Delta;
      }
      // fprintf (stderr, "xc %g %g %g\n", xc, yc, zc);
      b[] = - 14.*exact (xc, yc, zc)*cs[];
    }

#if 0    
    foreach_face (z)
      if (z == 0. && fs.z[] > 0. && fs.z[] < 1.) {
	// fprintf (stderr, "fs %g %g %g %g\n", x, y, z, fs.z[]);
	embed_face_gradient_z (point, a, 0);
      }
    exit (0);
#endif
    
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
    scalar res[];
    double maxp = residual ({a}, {b}, {res}, &p), maxf = 0.;
    foreach()
      if (cs[] == 1. && fabs(res[]) > maxf)
	maxf = fabs(res[]);
    fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxf, maxp);

    timer t = timer_start();
    mgstats s = poisson (a, b, alpha = fs,
			 flux = a.boundary[embed] != symmetry ? embed_flux : NULL,
			 tolerance = 1e-6);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

    /**
    The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
    fields and their norms are computed. */
    
    scalar e[], ep[], ef[];
    foreach() {
      if (cs[] == 0.)
	ep[] = ef[] = e[] = nodata;
      else {
	e[] = a[] - exact (x, y, z);
	ep[] = cs[] < 1. ? e[] : nodata;
	ef[] = cs[] >= 1. ? e[] : nodata;
      }
    }
    norm n = normf (e), np = normf (ep), nf = normf (ef);
    fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %d %d\n",
	     N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max, s.i, s.nrelax);

    /**
    The solution error is displayed using bview. */

    view (fov = 16.1659, quat = {-0.270921,0.342698,0.106093,0.893256},
	  tx = -0.00535896, ty = 0.000132663, bg = {1,1,1},
	  width = 600, height = 600, samples = 4);
    draw_vof("cs", "fs", color = "e");
    save ("e.png");
#endif
    
    //    dump ("dump"); // too big
  }
}

/**
## Results

### Dirichlet boundary condition

![Error on the finest mesh](dirichlet3D/e.png)

For Dirichlet boundary conditions, the residual converges at first order 
in partial cells.

~~~gnuplot Maximum residual convergence
set xrange [*:*]
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) '< grep maxres ../dirichlet3D/log' u (log($2)):(log($3)) via a,b
f2(x) = a2 + b2*x
fit f2(x) '' u (log($2)):(log($4)) via a2,b2
set xlabel 'Resolution'
set logscale
set cbrange [1:2]
set xtics 16,2,1024
set grid ytics
set ytics format "% .0e"
set xrange [16:512]
set ylabel 'Maximum residual'
set yrange [1e-6:]
set key bottom left
plot '' u 2:3 t 'full cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 2:4 t 'partial cells', exp(f2(log(x))) t ftitle(a2,b2)
~~~

This leads to third-order overall convergence.

~~~gnuplot Error convergence
set xrange [*:*]
fit f(x) '../dirichlet3D/log' u (log($1)):(log($3)) via a,b
fit f2(x) '' u (log($1)):(log($4)) via a2,b2
f3(x) = a3 + b3*x
fit f3(x) '' u (log($1)):(log($6)) via a3,b3
set xrange [16:512]
set ylabel 'Error'
set yrange [*:*]
plot '' u 1:3 pt 6 t 'max all cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 1:4 t 'avg partial cells', exp(f2(log(x))) t ftitle(a2,b2), \
     '' u 1:6 t 'avg full cells', exp(f3(log(x))) t ftitle(a3,b3)
~~~

### Neumann boundary condition

![Error on the finest mesh](neumann3D/e.png)

For Neumann boundary conditions, the residual converges at less than
first order in partial cells (which may be worth looking into).

~~~gnuplot Maximum residual convergence
set xrange [*:*]
fit f(x) '< grep maxres log' u (log($2)):(log($3)) via a,b
fit f2(x) '' u (log($2)):(log($4)) via a2,b2
set ylabel 'Maximum residual'
set xrange [16:512]
set yrange [1e-6:]
set key bottom left
plot '' u 2:3 t 'full cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 2:4 t 'partial cells', exp(f2(log(x))) t ftitle(a2,b2)
~~~

But this now leads to second-order overall convergence.

~~~gnuplot Maximum error convergence
set xrange [*:*]
fit f(x) 'log' u (log($1)):(log($7)) via a,b
fit f2(x) '' u (log($1)):(log($5)) via a2,b2
set ylabel 'Maximum error'
set xrange [16:512]
set yrange [*:*]
plot '' u 1:7 pt 6 t 'full cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 1:5 t 'partial cells', exp(f2(log(x))) t ftitle(a2,b2)
~~~

## References

~~~bib
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
