/**
# Making a radial profile
Using the *interfaceaverage.h* file, we can now make a profile of a field in the radial direction (i.e. among other things). Almost all required geometric functions were already in the source code.
*/

# include "grid/octree.h"
# include "fractions.h"
# include "utils.h"
# include "interfaceaverage.h"

/**
This is a generic function of a sphere that will be usefull to tell the program that we want radial averages.
*/
#define interface   (sq(x)+sq(y)+sq(z)-sq(R))
#define func(x,y,z) (sq(x)+sq(y)+sq(z)-sq(R))
double R;
scalar f[],g[];
/** 
We use a simple all-or-nothing boundary conditions for the volume fraction field
*/
f[left] =   sign(func(x-Delta,y,z));
f[right] =  sign(func(x+Delta,y,z));
f[top] =    sign(func(x,y+Delta,z));
f[bottom] = sign(func(x,y-Delta,z));
f[front] =  sign(func(x,y,z-Delta));
f[back] =   sign(func(x,y,z+Delta));
int main(){
  FILE * fp = fopen("radprof.dat","w");
  L0=1.;
  X0=Y0=Z0=-L0/2;
  init_grid(64);
  /**
  Here the centred scalar field (*g*) that we will diagnose is defined. 
  */
  foreach()
    g[]=exp(-(sq(x)+sq(y)+sq(z)));
  boundary({g});
  /**
  A loop is used to cycle over different radii for the radial profiling.
  */
  for (R=0.01;R<1.3;R+=0.01){
    fraction(f,interface);
    boundary({f}); 
    fprintf(fp,"%g\t%g\t%g\n",R,interface_average(g,f),exp(-sq(R)));
  }
  fclose(fp);
}
/**
## Results
We can display the diagnosed profile and compare it with the analytical answer.


~~~gnuplot
set yr [0.0 : 1.01]

set xlabel 'r / (L0)'
set ylabel '<g>'
plot "radprof.dat" using 1:2 title "Diagnosed profile",\
     exp(-x*x) w lines title "Analytical profile"
~~~

When the interface is entirely within the domain (i.e. $r < L0/2$) the algorithm does reasonably well. The algorithm does not make any sense if the interface is entirely outside the domain (i.e for $r>L0\sqrt{3/4}$), but is does not crash. Also note that we have not set proper boundary conditions for *g* (and *f*), this influences the results for interfaces that lay near the domain's boundaries.

## Warning
This page only shows the general consistency of the algorithm. but note:

- It does not work with MPI (yet). 
- Interfaces that cross the domain boundaries are an issue at the moment. This is due to inproper boundary conditions for the used volume-fraction field *f*.  
- The algorithm does not have proper convergence characteristics, even if the entire inface if away from the domain boundaries. I think this is enherited from the second-order accurate interface reconstruction, i.e. *per face*.   
*/