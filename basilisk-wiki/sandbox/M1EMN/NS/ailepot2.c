/**
# Poisson equation on complex domains

We reproduce the test cases initially proposed by [Johansen and
Collela, 1998](#johansen1998), Problem 1, p. 14, with Dirichlet
boundary conditions and Problem 3, p. 19, with Neumann boundary
conditions. */

#include "embed.h"
#include "poisson.h"
 

/**
The exact solution is given by
$$
\phi(r,\theta) = r^4\cos 3\theta
$$
*/

static double exact (double x, double y) {
  return  1.;
}

static double bord (double x, double y) {
  return  0;
}


static double  extra(double x){
  return 0.17735*sqrt(x)-0.075597*x 
       - 0.212836*(x*x)+0.17363*(x*x*x)-0.06254*(x*x*x*x);}
static double intra(double x){
  return   -(0.17735*sqrt(x)-0.075597*x 
        -0.212836*(x*x)+0.17363*(x*x*x)-0.06254*(x*x*x*x));} 

/**
We also need a function for homogeneous Dirichlet condition. */

static
double dirichlet_homogeneous_bc (Point point, Point neighbor,
				 scalar s, void * data) {
  return 0.;
}

int main()
{
  N=32;
  L0=10;
   {
    origin (-L0/2., -L0/2.);
    init_grid (N);

    /**
    The shape of the domain is given by
 
    */
    
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = ((sq(x) + sq(y)) - (2*2));
     //phi[] = (y>0? ((x>0) && (x<1 ) ? extra(x)  - y : 1) :((x>0) && (x<1 ) ? - intra (x)  + y : 1)) ;
    }
    boundary ({phi});


    

    fractions (phi, cs, fs);  

    boundary ({cs,fs});
  //  restriction ({cs,fs});

/*foreach() {
     printf ("%g %g %g \n", x, y,cs[]);}
*/

    cm = cs;
    fm = fs;
    
    /**
    Conditions on the box boundaries are set (only relevant for Problem 3). */
    
    scalar a[], b[];
    
    a[left]   = 0;
    a.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
    a[right]  = 0;
    a.boundary_homogeneous[right] = dirichlet_homogeneous_bc;
    a[top]    = 0;
    a.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
    a[bottom] = 0;
    a.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;

    /**
    The boundary conditions on the embedded boundary are Dirichlet   */

 
    a[embed] = dirichlet (exact(x,y));
 
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
    
    foreach() {
      a[] = cs[] > 0. ? exact (x, y) : nodata;
      
      double xc = x, yc = y;
      if (cs[] > 0. && cs[] < 1.) {
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	line_center (n, alpha, cs[], &p);
	xc += p.x*Delta, yc += p.y*Delta;
      }
      b[] = 0. ;
    }
    boundary ({a,b});

    
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

    // FIXME: need to set minlevel to 4
    timer t = timer_start();
    mgstats s = poisson (a, b, alpha = fs,
			 embed_flux =
			 a.boundary[embed] != symmetry ? embed_flux : NULL,
			 tolerance = 1e-4, minlevel = 3);
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
	e[] = a[] - exact (x, y);
	ep[] = cs[] < 1. ? e[] : nodata;
	ef[] = cs[] >= 1. ? e[] : nodata;
      }
    }
   

    /**
    The solution is displayed   */
foreach() {
    printf ("%g %g %g %g\n", x, y, a[], cs[]);
}
   
    
  }
}

/**
## Results


http://www.lmm.jussieu.fr/~lagree/SOURCES/GERRIS/AILE/steady/index.html

*/
