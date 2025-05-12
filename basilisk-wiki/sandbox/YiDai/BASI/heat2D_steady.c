/** 
#Steady 2D heat equation 

Steady 2D heat equation is basically laplace equation
$$\frac{\partial T}{\partial t}=\alpha\left(\frac{\partial^{2} T}{\partial x^{2}}+\frac{\partial^{2} T}{\partial y^{2}}\right)$$

$$0=\left(\frac{\partial^{2} T}{\partial x^{2}}+\frac{\partial^{2} T}{\partial y^{2}}\right)$$
*/

#include "grid/multigrid.h"
#include "utils.h"
#include "poisson.h"

scalar T[];
scalar f[];


T[left] = sin(y);
T[right] = dirichlet (0);
T[top] = sin(x);
T[bottom] = dirichlet (0);



int main(){
  L0 = 10;
  X0 = -L0/2.;
  Y0 = -L0/2.;
  init_grid(1<<8);
  foreach() {
      T[] = 0.;
      f[] = 0.;
  }
  boundary({T});
  TOLERANCE = 1e-4;
  poisson(T, f);
  foreach(){
    printf("%g %g\n", x, T[]);
  }
  output_ppm (T, file = "T.png", spread = 3, linear = true);
}

/**
![solution](heat2D_steady/T.png)

*/


