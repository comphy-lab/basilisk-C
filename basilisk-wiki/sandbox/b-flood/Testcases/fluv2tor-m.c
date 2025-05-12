/**
# Transcritical test case with Manning friction

## Declarations

We call the Saint-Venant solver on a 1D grid and we add the Manning
friction term. */

#include "grid/cartesian1D.h"
#include "b-flood/saint-venant-topo.h"
#include "b-flood/manning.h"

int LEVEL;

scalar e[];
norm nerror;
double tmax = 1000, q0 = 2, z0;
double dx;

// Analytical solution for h(x) and dh/dx
double hex (double x) {
  if (x <= 500)
    return pow(4/G,1/3.)*(1 - 1/3.*tanh(3*(x/1000.-1/2.)));
  else
    return pow(4/G,1/3.)*(1 - 1/6.*tanh(6*(x/1000.-1/2.)));
}

double dhex (double x) {
  if (x <= 500)
    return -pow(4/G,1/3.)/(1000*sq(cosh(3/2.-3*x/1000.)));
  else
    return -pow(4/G,1/3.)/(1000*sq(cosh(3-3*x/500.)));
}

// Manning's friction term
double sfm (double x) {
  return -sq(n)*sq(q0)/pow(hex(x),10/3.);
}

// Z and dz/dx
// We use Runge-Kuta 4 algo to fix the topography

double dzex (double x) {
  return (sq(q0)/(G*cube(hex(x)))-1)*dhex(x) + sfm(x);
}

double zex (double x, double z) { 
  return z + dx/4.*(dzex(x-dx)+2*dzex(x-0.5*dx)+dzex(x));
}

/**
## Parameters

Definition of parameters and calling of the saint-venant-topo subroutine run().*/

int main()
{
  n = 0.0218;
  L0 = 1000.;
  X0 = 0;
  G = 9.81;
  tmax = 1000;
  for (LEVEL = 4; LEVEL <= 9; LEVEL++) {
    N = 1 << LEVEL;
    dx = L0/N;
    run();
    fprintf (stderr, "%d %g %g\n", N, nerror.avg, nerror.rms);
  }
}

/**
## Boundary condition

We fix h and u at the left boundary (fluvial). */

h[left] = dirichlet(max(hex(0),0));
eta[left] =  dirichlet(max(hex(0)+zb[],zb[]));
u.n[left] = dirichlet(max(q0/hex(0),0));

/**
We fix a free exit condition on the right boundary (torrential). */

u.n[right] = neumann(0);
h[right] =  neumann(0);
eta[right] = neumann(0);

/**
## Initial conditions */

event init (i = 0)
{
  // Because the slope is initially dry, we set a maximum time-step. 
  DT = 1e-2;
  // the topography start at the altitude z = 0 at the left of the domain
  z0 = 0;
  foreach() {
    zb[] = zex(x,z0);
    z0 = zb[];
    u.x[] = 0;
    h[] = 0;
  }
  boundary (all);
}

/**
## Error norms

We compute the differents error norms. */

event error (i++; t <= tmax)
{
  foreach()
    e[] = h[] - hex(x);  
  nerror = normf (e);
}

/**
## Gnuplot output

We print the water profile along the channel at the final time for each
resolution. */

event printprofile (t = tmax)
{
  char name[100];
  FILE * fp;
  sprintf (name, "profil-%i.dat", N);
  fp = fopen (name, "w");
  foreach()
    fprintf (fp,"%g\t%g\t%g\t%g\t%g\n", x, h[], zb[], hex(x), u.x[]);
  fclose (fp);
}

/**
## References 

~~~bib
@article{macdonald1997analytic,
  title={Analytic benchmark solutions for open-channel flows},
  author={MacDonald, I and Baines, MJ and Nichols, NK and Samuels, PG},
  journal={Journal of Hydraulic Engineering},
  volume={123},
  number={11},
  pages={1041--1045},
  year={1997},
  publisher={American Society of Civil Engineers},
  doi = {10.1061/(ASCE)0733-9429(1997)123:11(1041)}
}
~~~

## Results

~~~gnuplot Free surface and topography
set xlabel 'L (m)' 
set ylabel 'Height (m)' 
set xtics 
set ytics
plot [][] './profil-512.dat' u 1:($3+$4) w l lw 0.5 \
                axes x1y1 t 'exact solution :Zb + he', \
	  './profil-512.dat' u 1:($2+$3) w l lt 0 lw 7 \
                axes x1y1 t 'N=512 :Zb + h', \
	  './profil-32.dat' u 1:($2+$3) axes x1y1 t 'N=32: Zb + h', \
	  './profil-512.dat' u 1:3 w l axes x1y1 t 'topo: Zb'
~~~

~~~gnuplot Error convergence
set logscale
set xlabel 'Number of cells N' 
set ylabel 'Error norm (m)' 
set xtics 
set ytics
set cbrange [1:2]
ftitle(a,b) = sprintf("order %4.2f", -b)
f1(x) = a1 + b1*x
fit f1(x) 'log' u (log($1)):(log($2)) via a1,b1
f2(x) = a2 + b2*x
fit f2(x) 'log' u (log($1)):(log($3)) via a2,b2
plot exp (f1(log(x))) t ftitle(a1,b1), './log' u 1:2 t 'average error', \
     exp (f2(log(x))) t ftitle(a2,b2), './log' u 1:3 t 'rms error'
~~~
  		  
## Link to the homepage

* [Homepage](/sandbox/b-flood/Readme)
*/
