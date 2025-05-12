/**
# Transonic shock with Darcy friction

## Declarations

We use the Saint-Venant solver on a 1D grid and we add the Darcy
friction term. */

#include "grid/cartesian1D.h"
#include "b-flood/saint-venant-topo.h"
#include "b-flood/darcy.h"

int LEVEL;

scalar e[];
norm nerror;
int ne = 0;
double tmax = 200, q0 = 2, z0;
double dx ;
 
double a1 = 0.674202, a2 = 21.7112, a3 = 14.492, a4 = 1.4305;

// Analytical solution for h(x) and dh/dx
double hex (double x) {
  if (x <= 200/3.)
    return pow(4/G,1/3.)*((4/3.-x/100.)-9*x/1000.*(x/100.-2/3.));
  else
    return pow(4/G,1/3.)*(a1*pow(x/100.-2/3.,4) +
			  a1*cube(x/100.-2/3.) -
			  a2*sq(x/100.-2/3.) +
			  a3*(x/100.-2/3.) + a4);
}

double dhex (double x) {
  if (x <= 200/3.)
    return pow(4/G,1/3.)*(-9*x - 200)/50000.;
  else
    return pow(4/G,1/3.)*(a1*cube(x)/25000000. - a1*sq(x)/200000. +
			  a1*x/7500. + a1/675. - a2*x/5000. + a2/75. + a3/100.);
}

// Darcy friction term in kinematic formulation
double sfd(double x){
  return  -f/(8*G)*sq(q0)/cube(hex(x));
}

// Z and dz/dx
// We use RK4 to solve the topography
double dzex(double x) {
  return (sq(q0)/(G*cube(hex(x))) - 1.)*dhex(x) + sfd(x);
}

double zex(double x, double z) {
  return z + dx/4.*(dzex(x - dx) + 2.*dzex(x - 0.5*dx) + dzex(x));
}

/**
## Parameters

Definition of parameters and calling of the Saint-Venant subroutine
run(). */

int main()
{
  f = 0.093;
  L0 = 100.;
  X0 = 0;
  G = 9.81;
  for (LEVEL = 4; LEVEL <= 9; LEVEL++) {  
    N = 1 << LEVEL;
    dx = L0/N;
    run();
    fprintf (stderr, "%d %g %g\n", N, nerror.avg, nerror.rms);
  }
}

/**
## Boundary condition

We set h and q (u) at both boundaries (torrential upstream and fluvial
downstream). */

h[left] = dirichlet(max(hex(0),0));
eta[left] =  dirichlet(max(hex(0)+zb[],zb[]));
u.n[left] = dirichlet(max(q0/hex(0),0));

h[right] = dirichlet(max(hex(100),0));
eta[right] =  dirichlet(max(hex(100)+zb[],zb[]));
u.n[right] = dirichlet(max(q0/hex(100),0));

/**
## Initial conditions */

event init (i = 0)
{
  // Because the slope is initially dry, we set a maximum time-step. 
  DT = 1e-2;
  z0 = 0;
  foreach(){
    zb[] = zex(x,z0);
    z0 = zb[];
    u.x[] = 0;
  }
  boundary(all);
}

/**
## Error norms

We compute the different error norms */

event error (i++; t <= tmax) {
  foreach()
    e[] = h[] - hex(x);  
  nerror = normf (e);
}

/**
## Gnuplot output

We print the water profile along the channel at final time.*/

event printprofile (t = tmax)
{
  char name[100];
  FILE * fp;
  sprintf (name, "profil-%i.dat", N);
  fp = fopen(name, "w");
  foreach()
    fprintf (fp, "%g\t%g\t%g\t%g\t%g\n",
	     x, h[], zb[], hex(x), u.x[]);
  fclose(fp);
}

/**
## References 

~~~bib
@article{macdonald1997,
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

@Article{popinet2011,
  author = 	 {S. Popinet},
  title = 	 {Quadtree-adaptive tsunami modelling},
  journal = 	 {Ocean Dynamics},
  year = 	 {2011},
  url =          {http://gerris.dalembert.upmc.fr/papers/tsunami.pdf},
  volume =       {61},
  number =       {9},
  pages =        {1261-1285}
}
~~~

## Results

~~~gnuplot Free surface and topography
set xlabel 'L (m)' 
set ylabel 'Height (m)' 
set xtics 
set ytics
set y2label 'Error (m)'
set y2tics
set key l b
plot [][] './profil-512.dat' u 1:($3+$4) w l lw 0.5 \
                axes x1y1 t 'exact solution :Zb + he', \
	  './profil-512.dat' u 1:($2+$3) w l lt 0 lw 7 \
                axes x1y1 t 'N=512 :Zb + h', \
	  './profil-32.dat' u 1:($2+$3) axes x1y1 t 'N=32: Zb + h', \
	  './profil-512.dat' u 1:3 w l axes x1y1 t 'topo: Zb', \
	  './profil-512.dat' u 1:($2-$4) w l \
	        axes x1y2 t'error N=512: h - he (right axis)'
~~~

~~~gnuplot Error convergence
reset
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
