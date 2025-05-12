/**
# Couette flow between rotating cylinders

We test embedded boundaries by solving the (Stokes) Couette flow
between two rotating cylinders. If SLIP or NAVIER is defined, 
the outer cylinder is either free of viscous stress or 
the viscous tresses is equal to the tangential velocity at the solid boundary.
This [test](/sandbox/lopez/slip.c) has been proposed by Lopez and Ghigo to validate slip boundary condition */

#include "grid/multigrid.h"
//#include "embed.h"
#include "tavares/embed-navier.h"
//#include "../embed_C.h" // Lopez version
#include "navier-stokes/centered.h"
#include "view.h"
#define SLIP 0
#define NAVIER 1

int main()
{
  size (1. [0]);
  DT = 1. [0];
  
  origin (-L0/2., -L0/2.);
  
  stokes = true;
  //DT = 1e-4;
  //TOLERANCE = 1e-4;

  int NF = 64;
  for (N = 16; N <= NF; N *= 2)
    run();
}

double lambda = 0; //equivalent to impose Dirichlet BC
scalar un[];
#define WIDTH 0.5

event init (t = 0) {

  /**
  The viscosity is unity. */
  
  mu = fm;

  /**
  The geometry of the embedded boundary is defined as two cylinders of
  radii 0.5 and 0.25. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = difference (sq(0.5) - sq(x) - sq(y),
			sq(0.25) - sq(x) - sq(y));
  boundary ({phi});
  fractions (phi, cs, fs);

  /**
  The outer cylinder is fixed and the inner cylinder is rotating with
  an angular velocity unity. */
 
#if SLIP // OR NAVIER
  u.n[embed] = dirichlet (x*x + y*y > 0.14 ? nodata : -y);
  u.t[embed] = dirichlet (x*x + y*y > 0.14 ? nodata :  x);
#elif NAVIER
  u.n[embed] = dirichlet (x*x + y*y > 0.14 ? nodata : -y);
  u.t[embed] = dirichlet (x*x + y*y > 0.14 ? nodata :  x);
   //u.n[embed] = (x*x + y*y > 0.14 ? navier (0.25) : dirichlet (-y)); Lopez version
   //u.t[embed] = (x*x + y*y > 0.14 ? navier (0.25) : dirichlet (x));
#else// DIRICHLET 
  u.n[embed] = dirichlet (x*x + y*y > 0.14 ? 0. : -y);
  u.t[embed] = dirichlet (x*x + y*y > 0.14 ? 0. :   x);
#endif
  /**
  We initialize the reference velocity field. */
  
  foreach()
    un[] = u.y[];
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

/**
We compute error norms and display the angular velocity, pressure and
error fields using bview. */

#if SLIP
  #define powerlaw(r) (sq(0.25)/(sq(0.25) + sq(0.5))*(sq(0.5)/r + r))
#elif NAVIER  
  #define powerlaw(r) (r*(0.5+lambda) + sq(0.5)*(lambda-0.5)/r)*(sq(0.25)/(sq(0.25)*(0.5+lambda)+sq(0.5)*(lambda-0.5)))
#else
  #define powerlaw(r) (r*(sq(0.5/r) - 1.)/(sq(0.5/0.25) - 1.))
#endif

event profile (t = end)
{
  scalar utheta[], e[];
  foreach() {
    double theta = atan2(y, x), r = sqrt(x*x + y*y);
    if (cs[] > 0.) {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      e[] = utheta[] - powerlaw (r);
    }
    else
      e[] = p[] = utheta[] = nodata;
  }

  norm n = normf (e);
  #if NAVIER
    char name[50];
    sprintf (name, "file-Navier-lambda%g.dat", lambda);
    static FILE * fp = fopen (name,"w");
  #elif SLIP
    char name[50];
    sprintf (name, "file-slip.dat");
    static FILE * fp = fopen (name,"w");
  #else
    char name[50];
    sprintf (name, "file-dirichlet.dat");
    static FILE * fp = fopen (name,"w");
  #endif  

  fprintf (fp, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  dump();
  
  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("utheta", spread = -1);
  save ("utheta-dirichlet.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("p", spread = -1);
  save ("p-dirichlet.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("e", spread = -1);
  save ("e-dirichlet.png");

  if (N == 32){
    #if SLIP
      char name1[50];
      sprintf (name1, "out%d-slip.dat", N);
      static FILE * fp1 = fopen (name1,"w");
    #elif NAVIER
      char name1[50];
      sprintf (name1, "out", N, lambda);
      static FILE * fp1 = fopen (name1,"w");
    #else
      char name1[50];
      sprintf (name1, "out%d-dirichlet.dat", N);
      static FILE * fp1 = fopen (name1,"w");
    #endif

    foreach() {
      double theta = atan2(y, x), r = sqrt(x*x + y*y);
      fprintf (fp1, "%g %g %g %g %g %g %g %g %g\n",
	       r, theta, u.x[], u.y[], p[], utheta[], cs[], e[], dt);
    }
  }  
}

/**
## Results

![Angular velocity](slip/utheta-dirichlet.png)

![Pressure field](slip/p-dirichlet.png)

![Error field](slip/e-dirichlet.png)

~~~gnuplot Velocity profile (N = 32)
set xlabel 'r'
set ylabel 'u_theta'
powerlawnavier(r,z)= (r*(0.5+lambda) + (0.5**2)*(lambda-0.5)/r)*((0.25**2)/((0.25**2)*(0.5+lambda)+(0.5**2)*(lambda-0.5)))
set grid
set arrow from 0.25, graph 0 to 0.25, graph 1 nohead
set arrow from 0.5, graph 0 to 0.5, graph 1 nohead

plot [0.2:0.55][-0.05:0.5] 'out' u 1:6 t 'numerics' w p pt 2 ps 2 lt -1 lw 2,\
powerlawnavier(x, lambda=0) t 'theory' lt -1 dt 2 lc rgb "black" lw 2.5


~~~

*/
