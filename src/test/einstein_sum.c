/**
# Tests for the `einstein_sum()` macro */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

/**
We create a second-order tensor structure. */

typedef struct{
  coord x;
  coord y;
  coord z;
} coord2;

/**
We create a third-order tensor structure. */

typedef struct{
  coord2 x;
  coord2 y;
  coord2 z;
} coord3;

int main ()
{
  init_grid (64);

  // define time and length scale
  double gamma = 1 [0,-1];
  double length = 1 [1,0];

  // define arbitrary second and third order tensors
  coord2 E = {{gamma,0.1*gamma}, {0.2*gamma, - gamma}};
  coord2 E2 = {{gamma/length, 0.2* gamma/length},
                {0.7*gamma/length, - 0.3*gamma/length}};
  coord3 K = {E2, E2};
  vector u2[];
  
  /** 
  We now initialize a quadratic velocity field around the origin.
  Formally, we have
  $$ \textbf{u}(\textbf{x}) = \textbf{E}\cdot \textbf{x}+\textbf{E2}: \textbf{xx}$$
  This can be encoded in indices notation as follows: */
  
  foreach(){
    coord X = {x,y,z};
    einstein_sum (i,k,j) {
      // factorized expression
      u.i[] = (E.i.j + K.i.j.k*X.k)*X.j;
      // developed expression
      u2.i[] = E.i.j*X.j + K.i.j.k*X.k*X.j;
    }
  }
  coord2 R =  (coord2) {{0}} ;
  double kinetic_energy = 0; 
  coord2 R2 =  (coord2) {{0}} ;
  double kinetic_energy2 = 0; 
  foreach()
    einstein_sum (i,j) {
      R.i.j += u.i[]*u.j[];
      kinetic_energy += u.i[]*u.i[]/2.;
      R2.i.j += u2.i[]*u2.j[];
      kinetic_energy2 += u2.i[]*u2.i[]/2.;
    }

  fprintf (stderr, "%g %g\n", kinetic_energy, kinetic_energy2);
  
  einstein_sum (i,j) {
    fprintf (stderr, "%g %g\n", R.i.j, R2.i.j);
  }
}
