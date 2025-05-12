/**
# Potential flow continuum equation solver

The steady continuum equations reads,
$$\nabla \cdot \rho \mathbf{v} = S, $$
with density $\rho$, flowfield $\mathbf{v}$ and source term $S$. Assuming $\mathbf{v} = -\nabla P$, one can rewite, 
$$-\rho \nabla^2 P - \nabla P \cdot  \nabla\rho = S $$

Which we will solve for $P$, with given $S$ and $\rho$.

## Results
Example problem:

![$\rho$](continuum/rho.png)

![$S$](continuum/S.png)

![$P$](continuum/P.png)

## Method

We use a Miltigrid-accelerated Poisson-equation solver
*/
#include "grid/multigrid.h"
#include "poisson.h"
#include "utils.h"

// Declare fields and set some boundary conditions
scalar P[], rho[], S[];
P[left] = dirichlet (0);
P[top] = dirichlet (0);
P[top] = neumann (0);

/**
We iteratively solve the corresponding Poisson problem and update the gradient field iteratively, bringing it to the right-hand side. 

You can see the adjustment stages:
![Each field solves the intermediate Poisson problem](continuum/P.mp4)
*/
mgstats solve_continuum (scalar P, scalar rho, scalar S) {
  // Compute the gradient of rho
  vector drho[];
  gradients ({rho}, {drho});
  // Log convergence
  mgstats MG = {0};
  MG.resb = 1;
  // A new right-hand-side scalar field
  scalar St[];
  while (MG.resb > TOLERANCE) {
    foreach() {
      St[] = S[];
      foreach_dimension()
	St[] -= drho.x[]*(P[1] - P[-1])/(2*Delta);
      St[] /= -rho[];
    }
    MG = poisson (P, St);
    output_ppm (P, file = "P.mp4", n = 300,
                min = -1, max = 1, opt = "-r 3");
    if (MG.resa > 1.1*MG.resb) {
      fprintf (stderr, "# Unstable iterations\n No solution was found\n");
      return MG;
    }
  }
  return MG;
}

int main() {
  L0 = 15;
  X0 = Y0 = -L0/2;
  init_grid (N);
  double R = 8;
  /**
  Initialize the fields
  */
  foreach() {
    P[] = 0;
    S[] = exp(-sq(x - 1) - sq(y - 1)) - 1.5*exp(-sq(x + 2) - sq(y + 2));
    rho[] = exp (-(sq(x/2) + sq(y))/sq(R));
  }
  /**
  Solve and output images
  */
  solve_continuum (P, rho, S);
  output_ppm (rho, file = "rho.png", n = 300);
  output_ppm (S, file = "S.png", n = 300);
  output_ppm (P, file = "P.png", n = 300);
  
}

