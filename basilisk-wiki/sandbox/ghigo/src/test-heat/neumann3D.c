/**
# Heat equation in 3D fixed, expanding and shrinking complex domains with Dirichlet and Neumann boundary conditions

We reproduce a test case proposed by [Schwarz et al.,
2006](#schwarz2006). The following heat equation is solved:
$$
\left\{
\begin{aligned}
&
a_t = \Delta a + f \\
&
f(x,y,z) = 4\frac{x^2 + y^2 + z^2 + 5(t + 1)}{125 \pi
(t+1)^3}\exp\left( -\frac{x^2 + y^2 + z^2}{5(t + 1)}\right)
\end{aligned}
\right.
$$
on three different domains (fixed, shrinking and expanding):
$$
\Omega_1 = \left\{(x,y,z): r < 0.392\right\},
$$
$$
\Omega_2 = \left\{(x,y,z): r < 0.392 - t\right\},
$$
$$
\Omega_3 = \left\{(x,y,z): r < 0.392 + t\right\}.
$$
The exact solution in all cases is:
$$
a(x,y,z) = \frac{4}{5\pi (t + 1)}\exp\left( -\frac{x^2 + y^2 +
z^2}{5(t + 1)} \right).
$$
*/

#include "../myembed.h"
#include "run.h"
#include "../mydiffusion.h"
mgstats mgp, mgpf, mgu;
#include "../myperfs.h"
#include "view.h"

/**
## Exact solution

We define here the exact solution. */

static double exact (double t, double x, double y, double z)
{
  return 4./(5.*pi*(t + 1))*exp (-(sq (x) + sq (y) + sq (z))/(5*(t + 1)));
}

/**
Next, we define the gradient of the exact solution. */

double exact_gradient (Point point, double t, double xc, double yc, double zc)
{
  coord p, n;
  embed_geometry (point, &p, &n);

  return (n.x*(-2.*xc)/(5.*(t + 1))*exact (t, xc, yc, zc) +
	  n.y*(-2.*yc)/(5.*(t + 1))*exact (t, xc, yc, zc) +
	  n.z*(-2.*zc)/(5.*(t + 1))*exact (t, xc, yc, zc));
}

/**
We then define the forcing term, evaluated at the geometric center of
each cell (which is different from the center of a cell for
cut-cells). */

double f (double t, double x, double y, double z)
{
  return 4.*(sq (x) + sq (y) + sq (z) + 5.*(t + 1))/(125.*pi*cube (t + 1))*exp (-(sq (x) + sq (y) + sq (z))/(5*(t + 1)));
}

/**
We also define the shape of the domain. */

#if GEOM // concave
void p_shape (scalar c, face vector f, double r)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = sq(x) + sq(y) + sq(z) - sq(r);
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}
#else // convex
void p_shape (scalar c, face vector f, double r)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = sq(r) - sq(x) - sq(y) - sq(z);
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}
#endif // GEOM

#if GEOM
/**
In this case, the poisson problem requires the homogeneous equivalent
of the boundary conditions imposed on the non-embedded boundaries of
the domain. We therefore define a function for homogeneous Dirichlet
condition as we are not using the *dirichlet* function to impose the
boundary conditions for *a* on the domain boundaries. */

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
  g->x = (-2.*x)/(5.*(t + 1))*exact (t, x, y, z);
  g->y = (-2.*y)/(5.*(t + 1))*exact (t, x, y, z);
  g->z = (-2.*z)/(5.*(t + 1))*exact (t, x, y, z); 
}
#endif // TREE

/**
## Setup 

We define the scalar *a* for the temperature and the variable *radius*
to track the evolution of the radius of the sphere, if necessary. We
also record the statistics of the Poisson solver. */

scalar a[];
double radius;
mgstats s;

int lvl;

int main()
{
  /**
  The domain is $1\times 1\times 1$. */
  
  size (1 [0]);
  origin (-L0/2., -L0/2., -L0/2.);

  /**
  We set the time step to respect the diffusive time step for the
  smallest mesh size considered and to minimize splitting errors. */
  
  DT = 1.e-5;

  /**
  We set the tolerance and the Poisson solver. */
    
  TOLERANCE = 1.e-6;
  NITERMAX  = 250; // For the initial iterations of cases with mesh adaptation

  for (lvl = 4; lvl <= 7; lvl++) { // minlevel = 2 (3.1pt/D, 3.93pt/D and 2.33pt/D)

    /**
    We initialize the grid. */

    N = 1 << (lvl);
    init_grid (N);
    
    /**
    Finally, we initialize the value of radius. */

    radius = 0.392;

    run ();
  } 
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */

#if ORDER2
  a.third = false;
#else
  a.third = true;
#endif // ORDER2

#if TREE    
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  also modify the prolongation an restriction operators for *a*. The
  prolongation and restriction operators for the volume and face
  fractions *cs* and *fs* are modified in the [*event
  metric*](/src/embed.h). */
    
  a.refine = a.prolongation = refine_embed_linear;
  a.restriction = restriction_embed_linear;

  /**
  We also define the gradient of *a* at the full cell center of
  cut-cells. */
    
  a.embed_gradient = a_embed_gradient;
#endif // TREE
  
  /**
  We initialize the embedded boundary. */

#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, radius);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			maxlevel = (lvl), minlevel = (lvl) - 2);
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, radius);

  /**
  We also initialize the volume fraction at the previous timestep
  *csm1*. */

  trash ({csm1});
  foreach()
    csm1[] = cs[];
  restriction ({csm1}); // Since restriction/prolongation depend on csm1
  
  /**
  We initialize *a*. */
  
  foreach()
    a[] = cs[] > 0. ? exact (0., x, y, z) : nodata;

  /**
  We define the boundary conditions for *a* on the domain boundary, if
  necessary. */

#if GEOM
  a[left]   = exact (t + dt, x - Delta/2., y, z);
  a.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
  a[right]  = exact (t + dt, x + Delta/2., y, z);
  a.boundary_homogeneous[right] = dirichlet_homogeneous_bc;
  a[top]    = exact (t + dt, x, y + Delta/2., z);
  a.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
  a[bottom] = exact (t + dt, x, y - Delta/2., z);
  a.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;
  a[front]  = exact (t + dt, x, y, z + Delta/2.);
  a.boundary_homogeneous[front] = dirichlet_homogeneous_bc;
  a[back]   = exact (t + dt, x, y, z - Delta/2.);
  a.boundary_homogeneous[back] = dirichlet_homogeneous_bc;
#endif // GEOM    

#if DIRICHLET
  a[embed] = dirichlet (exact (t + dt, x, y, z));
#else
  a[embed] = neumann (exact_gradient (point, t + dt, x, y, z));
#endif
}

/**
## Outputs */

scalar e[], ep[], ef[];
double eavg  = 0., emax  = 0.;
double epavg = 0., epmax = 0.;
double efavg = 0., efmax = 0.;

event init (i = 0)
{
  eavg  = 0., emax  = 0.;
  epavg = 0., epmax = 0.;
  efavg = 0., efmax = 0.;
}

event error (i = 1; i++)
{
  /**
  The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
  fields and their norms are computed. */
    
  foreach() {
    if (cs[] == 0.)
      ep[] = ef[] = e[] = nodata;
    else {
      e[]  = fabs (a[] - exact (t, x, y, z));
      ep[] = cs[] < 1. ? e[] : nodata;
      ef[] = cs[] >= 1. ? e[] : nodata;
    }
  }
  norm n = normf (e), np = normf (ep), nf = normf (ef);

  // All cells
  eavg  += dt*n.avg;
  emax  += dt*n.max;
  // Partial cells
  epavg += dt*np.avg;
  epmax += dt*np.max;
  // Full cells
  efavg += dt*nf.avg;
  efmax += dt*nf.max;
}

/**
We set the final time small enough to make sure that when using the
expanding domain $\Omega_3$, the embedded boundaries do not overflow
outside of the computational domain. */

event logfile (t = 0.1)
{
  /**
  The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
  are saved. */
    
  fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %d %d %d %g %g\n",
	   N,
	   (eavg/(t + SEPS)),  (emax/(t + SEPS)),
	   (epavg/(t + SEPS)), (epmax/(t + SEPS)),
	   (efavg/(t + SEPS)), (efmax/(t + SEPS)),
	   s.i, s.nrelax,
	   i, t, dt);
  fflush (stderr);

  /**
  The solution and the error are displayed using bview. */

  if (lvl == 6) {
    clear ();
    view (fov = 32.2073, quat = {-0.309062,0.243301,0.0992085,0.914026},
	  tx = 0.0122768, ty = 0.0604286, bg = {1,1,1},
	  width = 800, height = 800);

    box ();
    draw_vof ("cs", "fs");
    save ("vof.png");

    box ();
    draw_vof ("cs", "fs", color = "e");
    save ("e.png");

    cells (n = {0, 0, 1}, alpha = -radius/2.);
#if GEOM
    draw_vof ("cs", "fs", color = "e");
#endif // GEOM
    squares ("e", n = {0, 1, 0}, alpha = -radius/2., spread = -1);
    save ("e-slice.png");

    box ();
    draw_vof ("cs", "fs", color = "a");
    save ("a.png");

    cells (n = {0, 0, 1}, alpha = -radius/2.);
#if GEOM
    draw_vof ("cs", "fs", color = "a");
#endif // GEOM
    squares ("a", n = {0, 1, 0}, alpha = -radius/2., spread = -1);
    save ("a-slice.png");    
  }
}

/**
## Time integration */

event integration (i++)
{
  dt = dtnext (DT);

#if EXPANDING || SHRINKING
  /**
  We make sure that the time step verifies the CFL condition for the
  displacement of the embedded boundary. */
  
  foreach(reduction(min:dt)) {
    if (cs[] > 0. && cs[] < 1.) {

      // Local maximum velocity
      double umax = 1; // u = ((radius +/- dt) - radius)/dt; 
      
      // Non-restrictive timestep (independent of *cs* and *fs*)
      double dte = Delta/(fabs (umax) + SEPS);
      if (dte < dt)
	dt = dte;
    }
  }

  /**
  We then store the previous volume fraction in *csm1*. */

  trash ({csm1});
  foreach()
    csm1[] = cs[];
  restriction ({csm1}); // Since restriction/prolongation depend on csm1

  /**
  We finally modify the shape of the embedded boundary. */

#if EXPANDING
  radius += dt;
#elif SHRINKING
  radius -= dt;
#endif // EXPANDING && SHRINKING
  p_shape (cs, fs, radius);
  
  /**
  When considering the expanding domain $\Omega_3$, as the radius of
  the sphere increases, emerged cells appear that have no time
  history. We therefore initialize the values of *a* in these newly
  emerged cells by extrapolating their values in the direction of the
  normal to the embedded boundary, while making sure not to use any
  values from other newly emerged cells by setting the flag *emerged*
  to false. This is a necessary step as these values are used in the
  rhs of the diffusion equation (i.e. it is not just an initial
  guess). */

  emerged = false;

  foreach() {
    if (cs[] <= 0.) // Solid cells
      a[] = nodata;
    else if (cs[] > 0. && csm1[] <= 0.) { // Emerged cells

      assert (cs[] < 1.);
      
      // Cell centroid, barycenter and normal of the embedded fragment
      coord c = {0., 0., 0.}, b, n;
      embed_geometry (point, &b, &n);

      // Boundary condition on the embedded boundary
      bool dirichlet = true;
      double ab = (a.boundary[embed] (point, point, a, &dirichlet));
      if (!dirichlet) {
	double coef = 0.;
	ab = neumann_scalar (point, a, cs, n, b, ab, &coef);
	// a[] is undefined here so we use the average over the neighbors
	// ab += coef*a[];
	int navg = 0;
	double aavg = 0.;
	foreach_neighbor(1)
	  if (cs[] > 0. && (emerged || csm1[] > 0.)) {
	    navg += 1;
	    aavg += a[];
	  }
	ab += coef*(aavg/(navg + SEPS));
      }
 
      // Emerged cell value at time t
#if EXTRAPOLATE == 1
      a[] = exact (t, x, y, z);
#else // 0
      a[] = embed_extrapolate (point, a, cs, n, c, ab);
#endif // EXTRAPOLATE
    }
  }
  
  /**
  Before using the *boundary* function, we set the *emerged* flag to
  true to indicate that all emerged cells have been updated and can
  now be used. */
  
  emerged = true;
#endif // EXPANDING || SHRINKING
    
  /**
  We define the right hand-side *r* at time *t* on the new geometry,
  using the coordinates of the barycenter of the cut cell (xc,yc,zc),
  which is calculated from the cell and surface fractions.

  Note here that we don't multiply *r* by *cs* as this will done in
  [mydiffusion.h](src/mydiffusion.h). */
    
  scalar r[];

#if TREE
  /**
  In theory for *TREE*, we should also use embed-specific prolongation
  and restriction operators for *r*. However, *r* is only used to
  compute the residual in leaf cells and its stencil is never
  accessed. Furtheremore, and more importantly, if a *res* variable is
  not given as a parameter of the *poisson* function (which is the
  case here), than the prolongation and restriction operators for the
  *res* variable used in *mg_solve* are copied from those of
  *r*. However, for the multigrid solver to work properly, the
  prolongation and restriction operators for the *res* variable should
  be similar to the multigrid prolongation and restriction
  operators. We therefore use the default operators for *r* here. */
    
  /* r.refine = r.prolongation = refine_embed_linear; */
  /* r.restriction = restriction_embed_linear; */
#endif // TREE
  
  foreach() {
    double xc = x, yc = y, zc = z;
    if (cs[] > 0. && cs[] < 1.) {
      coord n = facet_normal (point, cs, fs), p;
      double alpha = plane_alpha (cs[], n);
      plane_center (n, alpha, cs[], &p); // Different from line_area_center
      xc += p.x*Delta, yc += p.y*Delta, zc += p.z*Delta;
    }
    r[] = f (t, xc, yc, zc);
  }

  /**
  We define a unity diffusion coefficient, accounting for the
  metric. */
  
  face vector D = fm;

  /**
  We then use the diffusion solver to advance the solution to time
  *t+dt*. */

  s = diffusion (a, dt, D, r = r);
  mgp = mgpf = mgu = s;
}

/**
## Mesh adaptation */

#if TREE
event adapt (i++)
{
  /**
  We refine the mesh around the embedded boundary. */

  adapt_wavelet ({cs}, (double[]) {1.e-30},
		 maxlevel = (lvl), minlevel = (lvl) - 2);
  fractions_cleanup (cs, fs,
		     smin = 1.e-14, cmin = 1.e-14);
}
#endif // TREE
  
/**
## Results

![Embedded boundaries for *l=6*](neumann3D/vof.png)

![Solution for *l=6*](neumann3D/a.png)

![Solution for *l=6*](neumann3D/a-slice.png)

![Error for *l=6*](neumann3D/e.png)

![Error for *l=6*](neumann3D/e-slice.png)

#### Errors

~~~gnuplot Average error convergence
set terminal svg font ",16"
set key bottom left spacing 1.1
set xtics 4,4,512
set grid ytics
set ytics format "%.0e" 1.e-12,100,1.e2
set xlabel 'N'
set ylabel '||error||_{1}'
set xrange [8:256]
set yrange [1.e-10:1.e-2]
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

[Schwarz et al., 2006](#schwarz2006) determined the following
estimates for the error of the solution of the Poisson problem:

* for full cells:
$$
\mathrm{err}_f \sim \Delta^2
$$

* for cut-cells:
$$
\mathrm{err}_p \sim \Delta^2
$$

By choosing a timestep *DT* small enough to minimize splitting errors,
we recover the expected orders of convergence. Note that when using
mesh adaptation, the order of convergence of the error is determined
by the order of the prolongataion/refinement functions, which are
second-order accurate in Basilisk. However, in the complex "convex"
geometry presented here, the restriction operation requires the user
to input an *embed_gradient* to maintain second-order
accuracy. Furthermore, for the case where the isotropic expansion of
the sphere uncovers emerged cells, neighboring emerged cells are
likely to appear simultaneously. In this case, the number of available
cells in the stencil of these emerged cells is too small to use a
second-order extrapolation and a first-order extrapolation is used
instead.

~~~gnuplot Order of convergence of the average error
reset
set terminal svg font ",16"
set key bottom right spacing 1.1
set xtics 4,4,256
set ytics -10,1,10
set grid ytics
set xlabel 'N'
set ylabel 'Order of ||error||_{1}'
set xrange [8:256]
set yrange [-4:6]
set logscale x

# Asymptotic order of convergence

ftitle(b) = sprintf(", asymptotic order n=%4.2f", -b);
Nmin = log(64);

f1(x) = a1 + b1*x; # cut-cells
f2(x) = a2 + b2*x; # full cells
f3(x) = a3 + b3*x; # all cells

fit [Nmin:*][*:*] f1(x) '< sort -k 1,1n log | awk "!/#/{print }"' u (log($1)):(log($4)) via a1,b1; # cut-cells 
fit [Nmin:*][*:*] f2(x) '< sort -k 1,1n log | awk "!/#/{print }"' u (log($1)):(log($6)) via a2,b2; # full-cells
fit [Nmin:*][*:*] f3(x) '< sort -k 1,1n log | awk "!/#/{print }"' u (log($1)):(log($2)) via a3,b3; # all cells

plot '< sort -k 1,1n log | awk "!/#/{print }" | awk -f ../data/order.awk' u 1:4 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '< sort -k 1,1n log | awk "!/#/{print }" | awk -f ../data/order.awk' u 1:6 w lp ps 1.25 pt 5 lc rgb "blue" t 'full cells'.ftitle(b2), \
     '< sort -k 1,1n log | awk "!/#/{print }" | awk -f ../data/order.awk' u 1:2 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3)
~~~

~~~gnuplot Order of convergence of the maximum error
set ylabel 'Order of ||error||_{inf}'

# Asymptotic order of convergence

fit [Nmin:*][*:*] f1(x) '< sort -k 1,1n log | awk "!/#/{print }"' u (log($1)):(log($5)) via a1,b1; # cut-cells 
fit [Nmin:*][*:*] f2(x) '< sort -k 1,1n log | awk "!/#/{print }"' u (log($1)):(log($7)) via a2,b2; # full-cells
fit [Nmin:*][*:*] f3(x) '< sort -k 1,1n log | awk "!/#/{print }"' u (log($1)):(log($3)) via a3,b3; # all cells

plot '< sort -k 1,1n log | awk "!/#/{print }" | awk -f ../data/order.awk' u 1:5 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '< sort -k 1,1n log | awk "!/#/{print }" | awk -f ../data/order.awk' u 1:7 w lp ps 1.25 pt 5 lc rgb "blue" t 'full cells'.ftitle(b2), \
     '< sort -k 1,1n log | awk "!/#/{print }" | awk -f ../data/order.awk' u 1:3 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3)
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
