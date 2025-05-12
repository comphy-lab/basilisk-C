/**
# Poisson equation on complex domains

We reproduce the test cases initially proposed by [Johansen and
Collela, 1998](#johansen1998), Problem 1, p. 14, with Dirichlet
boundary conditions and Problem 3, p. 19, with Neumann boundary
conditions. But now we solve it on both side of the interface ! */

#include "embed.h"
#include "poisson.h"
#include "view.h"

/**
The exact solution is given by
$$
\phi(r,\theta) = r^4\cos 3\theta
$$
*/

static double exact (double x, double y) {
  double theta = atan2 (y, x), r2 = x*x + y*y;
  return sq(r2)*cos (3.*theta);
}

double exact_gradient (Point point, double theta, double r)
{
  coord n = facet_normal (point, cs, fs);
  normalize (&n);
  double dsdtheta = - 3.*cube(r)*sin (3.*theta);
  double dsdr = 4.*cube(r)*cos (3.*theta);
  return (n.x*(dsdr*cos(theta) - dsdtheta*sin(theta)) +
	  n.y*(dsdr*sin(theta) + dsdtheta*cos(theta)));  
}

/**
We also need a function for homogeneous Dirichlet condition. */

static
double dirichlet_homogeneous_bc (Point point, Point neighbor,
				 scalar s, void * data) {
  return 0.;
}

int main()
{
  for (N = 128; N <= 512; N *= 2) {
    origin (-L0/2., -L0/2.);
    init_grid (N);

    /**
    The shape of the domain is given by
    $$
    \Omega = {(r,\theta): r \leq 0.30 + 0.15\cos 6\theta}
    $$
    for Problem 1 and
    $$
    \Omega = {(r,\theta): r \geq 0.25 + 0.05\cos 6\theta}
    $$
    for Problem 2.
    */
    
    vertex scalar phi[];
    foreach_vertex() {
      double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
#if DIRICHLET
      phi[] = 0.30 + 0.15*cos(6.*theta) - r;
#else
      phi[] = r - (0.25 + 0.05*cos(6.*theta));
#endif
    }
    boundary ({phi});
    fractions (phi, cs, fs);  
#if TREE
    cs.refine = cs.prolongation = fraction_refine;
#endif
    boundary ({cs,fs});
    restriction ({cs,fs});

    cm = cs;
    fm = fs;
    
    /**
    Conditions on the box boundaries are set (only relevant for Problem 3). */
    
    scalar a[], b[];
    scalar a2[], b2[];
    
    a[left]   = exact (x - Delta/2., y);
    a.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
    a[right]  = exact (x + Delta/2., y);
    a.boundary_homogeneous[right] = dirichlet_homogeneous_bc;
    a[top]    = exact (x, y + Delta/2.);
    a.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
    a[bottom] = exact (x, y - Delta/2.);
    a.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;

    /**
    The boundary conditions on the embedded boundary are Dirichlet and
    Neumann for Problem 1 and 3, respectively. */

#if DIRICHLET
    a[embed] = dirichlet (exact (x,y));
#else
    a[embed] = neumann (exact_gradient (point, atan2(y, x), sqrt(x*x + y*y)));
#endif

    /**
    We use "third-order" [face flux interpolation](/src/embed.h). */

    a.third = true;

    /**
    The right-hand-side
    $$
    \Delta\phi = 7r^2\cos 3\theta
    $$
    is defined using the coordinates of the
    barycenter of the cut cell (xc,yc), which is calculated from the
    cell and surface fractions. */
    
    int j = 0;
    scalar e[], ep[], ef[];
    scalar e2[];
    mgstats s;
    norm n, np , nf ;
    
/**
    Now we're ready to solve our equations on both side of the interface.
*/

    for (j = 0; j<=1; j++){

      foreach() {
        a[] = cs[] > 0. ? exact (x, y) : nodata;
        
        double xc = x, yc = y;
        if (cs[] > 0. && cs[] < 1.) {
    coord n = facet_normal (point, cs, fs), p;
    double alpha = plane_alpha (cs[], n);
    line_center (n, alpha, cs[], &p);
    xc += p.x*Delta, yc += p.y*Delta;
        }
        //      fprintf (stderr, "xc %g %g\n", xc, yc);
        double theta = atan2(yc, xc), r2 = sq(xc) + sq(yc);
        b[] = 7.*r2*cos (3.*theta)*cs[];
      }
      boundary ({a,b});

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
      if(j==0)
        fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxf, maxp);
      else
        fprintf (stderr, "maxres2 %d %.3g %.3g\n", N, maxf, maxp);


      // FIXME: need to set minlevel to 4
      timer t = timer_start();
      s = poisson (a, b, alpha = fs,
         embed_flux =
         a.boundary[embed] != symmetry ? embed_flux : NULL,
         tolerance = 1e-6, minlevel = 4);
      double dt = timer_elapsed (t);
      printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

      /**
      The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
      fields and their norms are computed. */
      
      foreach() {
        if (cs[] == 0.)
    ep[] = ef[] = e[] = nodata;
        else {
    e[] = a[] - exact (x, y);
    ep[] = cs[] <  1. ? e[] : nodata;
    ef[] = cs[] >= 1. ? e[] : nodata;
        }
      }

      n = normf (e), np = normf (ep), nf = normf (ef);
      fprintf (stderr, "err%d %d %.3g %.3g %.3g %.3g %.3g %.3g %d %d\n",
     j, N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max, s.i, s.nrelax);
  
      if(j==0){
          foreach(){
          a2[]  = a[];
          b2[]  = b[];
          e2[]  = e[];
          cs[]  = 1. - cs[];

        }

        foreach_face(){
          fs.x[] = 1.-fs.x[];
        }
        boundary({a2,b2,cs,fs});
        restriction({a2,b2,cs,fs});
      }


    }
    /**
    The solution is displayed using bview. */
    view (fov = 18);
    squares ("a", spread = -1);
    squares ("a2", spread = -1);
    draw_vof ("cs", "fs");
    save ("a.png");

    clear();
    squares ("e", spread = -1);
    squares ("e2", spread = -1);
    draw_vof ("cs", "fs");
    save ("e.png");
#endif
    
    dump ("dump");
  }
}

/**
## Results

### Problem 1

![Solution on the finest mesh](dirichlet/a.png)

![Error on the finest mesh](dirichlet/e.png)

For Dirichlet boundary conditions, the residual converges at first order 
in partial cells.

~~~gnuplot Maximum residual convergence
set xrange [*:*]
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) '< grep maxres ../dirichlet/log' u (log($2)):(log($3)) via a,b
f2(x) = a2 + b2*x
fit f2(x) '' u (log($2)):(log($4)) via a2,b2
f3(x) = a3 + b3*x
fit f3(x) '< grep maxres2 ../dirichlet/log' u (log($2)):(log($3)) via a3,b3
f4(x) = a4 + b4*x
fit f4(x) '' u (log($2)):(log($4)) via a4,b4
set xlabel 'Resolution'
set logscale
set cbrange [1:2]
set xtics 32,2,1024
set grid ytics
set ytics format "% .0e"
set xrange [32:1024]
set ylabel 'Maximum residual'
set key top right
plot '< grep maxres ../dirichlet/log' u 2:3 t 'full cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 2:4 t 'partial cells', exp(f2(log(x))) t ftitle(a2,b2), \
     '< grep maxres2 ../dirichlet/log' u 2:3 t 'full cells2', exp(f3(log(x))) t ftitle(a3,b3), \
     '' u 2:4 t 'partial cells2', exp(f4(log(x))) t ftitle(a4,b4), \
     'neumann.table1' u 1:4 w lp t 'full cells (J and C, 1998)', \
     'neumann.table1' u 1:2 w lp t 'partial cells (J and C, 1998)'
~~~

This leads to third-order overall convergence.

~~~gnuplot Error convergence (phase 1)
set xrange [*:*]
fit f(x) '< grep err0 ../dirichlet/log' u (log($2)):(log($4)) via a,b
fit f2(x) '' u (log($2)):(log($5)) via a2,b2
fit f3(x) '' u (log($2)):(log($7)) via a3,b3

set xrange [32:1024]
set ylabel 'Error'
set yrange [*:*]
plot '' u 2:4 pt 6 t 'max all cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 2:5 t 'avg partial cells', exp(f2(log(x))) t ftitle(a2,b2), \
     '' u 2:7 t 'avg full cells', exp(f3(log(x))) t ftitle(a3,b3), \
     'neumann.table2' u 1:2 w lp t 'max all cells (J and C, 1998)',\
     '' u 1:3 w lp t 'avg partial cells (J and C, 1998)', \
     '' u 1:4 w lp t 'avg full cells (J and C, 1998)'
~~~

~~~gnuplot Error convergence (phase 2)
set xrange [*:*]
fit f(x) '< grep err1 ../dirichlet/log' u (log($2)):(log($4)) via a,b
fit f2(x) '' u (log($2)):(log($5)) via a2,b2
fit f3(x) '' u (log($2)):(log($7)) via a3,b3

set xrange [32:1024]
set ylabel 'Error'
set yrange [*:*]
plot '' u 2:4 pt 6 t 'max all cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 2:5 t 'avg partial cells', exp(f2(log(x))) t ftitle(a2,b2), \
     '' u 2:7 t 'avg full cells', exp(f3(log(x))) t ftitle(a3,b3), \
     'neumann.table2' u 1:2 w lp t 'max all cells (J and C, 1998)',\
     '' u 1:3 w lp t 'avg partial cells (J and C, 1998)', \
     '' u 1:4 w lp t 'avg full cells (J and C, 1998)'
~~~


### Problem 3

![Solution on the finest mesh](neumann_double/a.png)

![Error on the finest mesh](neumann_double/e.png)

For Neumann boundary conditions, the residual also converges at first order 
in partial cells.

~~~gnuplot Maximum residual convergence
set xrange [*:*]
fit f(x) '< grep maxres log' u (log($2)):(log($3)) via a,b
fit f2(x) '' u (log($2)):(log($4)) via a2,b2
fit f3(x) '< grep maxres2 log' u (log($2)):(log($3)) via a3,b3
fit f4(x) '' u (log($2)):(log($4)) via a4,b4
set ylabel 'Maximum residual'
set xrange [32:1024]
set yrange [1e-6:]
set key top right
plot '< grep maxres log' u 2:3 t 'full cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 2:4 t 'partial cells', exp(f2(log(x))) t ftitle(a2,b2), \
     '< grep maxres2 log' u 2:3 t 'full cells', exp(f3(log(x))) t ftitle(a3,b3), \
     '' u 2:4 t 'partial cells', exp(f4(log(x))) t ftitle(a4,b4), \
     'neumann.table5' u 1:2 w lp t 'full cells (J and C, 1998)', \
     'neumann.table5' u 1:3 w lp t 'partial cells (J and C, 1998)'
~~~

But this now leads to second-order overall convergence.

~~~gnuplot Maximum error convergence (phase 1)
set xrange [*:*]
fit f(x) '< grep err0 log' u (log($2)):(log($4)) via a,b
fit f2(x) '' u (log($2)):(log($6)) via a2,b2
fit f3(x) '< grep err1 log' u (log($2)):(log($4)) via a3,b3
fit f4(x) '' u (log($2)):(log($6)) via a4,b4

set ylabel 'Maximum error'
set xrange [32:1024]
set yrange [*:*]
plot '< grep err0 log' u 2:4 pt 6 t 'all cells (1)', exp(f(log(x))) t ftitle(a,b), \
     '' u 2:6 t 'partial cells (1)', exp(f2(log(x))) t ftitle(a2,b2), \
     '< grep err1 log' u 2:4 pt 6 t 'all cells (2)', exp(f3(log(x))) t ftitle(a3,b3), \
     '' u 2:6 t 'partial cells (2)', exp(f4(log(x))) t ftitle(a4,b4), \
     'neumann.table5' u 1:5 w lp t 'all cells (J and C, 1998)', \
     'neumann.table5' u 1:4 w lp t 'partial cells (J and C, 1998)'
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
