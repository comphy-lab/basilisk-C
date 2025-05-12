/**
# Dewetting

Implemented from [Mahady, Afkhami, Kondic, Phys. Fluids 28, 062002
(2016)](http://doi.org/10.1063/1.4949522). 

Several things need to be improved: adaptivity, more robust height estimation 
(as done by [VariablePosition](http://gerris.dalembert.upmc.fr/gfsvariableposition.html)).
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

double R = 0.4, hs = 0.01, thetai = pi/2., Oh = 1./sqrt(10.);

int main() {
  N = 128;
  f.sigma = 1;
  mu1 = mu2 = sq(Oh)*f.sigma*2.*R;
  u.t[bottom] = dirichlet(0);
  f[bottom] = 1.;
  run();
}

event init (t = 0) {
#if 0
  fraction (f, 5.*hs*(1. + 0.01*sin(4.*M_PI*x)) - y);
#elif 0
  fraction (f, max(sq(R) - sq(x) - sq(y + R*cos(thetai) - hs),
  		   hs - y));
#else
  fraction (f, 5.*hs*(1. + 0.01*noise()) - y);
#endif
}

event acceleration (i++)
{
  scalar hy[];

  foreach()
    f[] = clamp (f[], 0, 1);
  
  position (f, hy, {0,1});

#if TREE
  f.prolongation = p.prolongation;
  f.dirty = true;
#endif

  // Eq. 10 of Mahady et al, 2016
  double theta_eq = pi/2., m = 3, n = 2;
  double xi = (1. - cos(theta_eq))/hs*(m - 1.)*(n - 1.)/(m - n);
  
  face vector st = a;
  foreach_face()
    if (f[] != f[-1]) {
      double h = 
	(hy[] < nodata && hy[-1] < nodata) ?
	(hy[] + hy[-1])/2. :
	hy[] < nodata ? hy[] :
	hy[-1] < nodata ? hy[-1] :
	nodata;
      
      if (h != nodata)
	st.x[] -= alpha.x[]/fm.x[]*f.sigma*xi*(pow(hs/h, m) - pow(hs/h, n))
	  *(f[] - f[-1])/Delta;
    }
  
#if TREE
  f.prolongation = fraction_refine;
  f.dirty = true;
#endif
}

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D -s dewetting.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

event logfile (i++) {
  fprintf (stderr, "%g %g\n", t, statsf (u.x).max);
}

event interface (i += 200; t <= 4) {
  output_facets (f);
}

/**
# Results

~~~gnuplot Evolution of the interface
plot 'out' w l t '', 0.01
~~~
*/
