/**
## A diffusion-equation solver

For a species-field concentration $c$, the diffusion equations reads
(or a version thereof),

$$\frac{\partial c}{\partial t} = \nabla \cdot \left( \kappa \nabla c
\right).$$

With $\kappa$ the diffusivity of the medium. Because of reasons, we
wish to use a finite-volume discretization to solve for the temporal
evolution of $c(\mathbf{x},t)$. We recognize the flux vector ($\mathbf{F}$):

$$\frac{\partial c}{\partial t} + \nabla \cdot \mathbf{F} = 0,$$

as,

$$\mathbf{F} = -\kappa \nabla c.$$

Given the alligment of the faces, we can make a function that writes
the fluxes ($F$) of $c$ in a face vector field. We also request its
user to provide the diffusivity on faces, which could be `const`ant:
*/

void flux_diffusion (scalar c, (const) face vector kappa, face vector F) {
  /**
 Neighboring cells are possibly behind a boundary of some
     sort. As such, we need to compute the relevant solution values. 
*/
  boundary ({c});
  /**
By convention, a face separates a cell from its left, bottom, front
neighbor in the `x,y` and `z` direction, respectively. The
`foreach_face()` function will automatically rotate over all
dimensions.
   */
  foreach_face()
    F.x[] = -kappa.x[]*(c[] - c[-1])/Delta; //2nd order accurate gradient *estimation*
}

/**
   Using this function we can easily compute the tendency for $c$ in a
   cell ($j$) with $N$ faces:

$$   \left(\frac{\partial c}{\partial t}\right)_j = \frac{1}{V_j} \sum_{i=1}^N \mathbf{F_i}
   \cdot \mathbf{A_i}.$$

Assuming that the ratio $\frac{V_j}{A_i}=\Delta$. Which is true for
one-dimensional, square or cubic cells.
*/

void tendency_from_flux (face vector F, scalar dc) {
  boundary_flux ({F});  //flux on levels
  foreach() {
    dc[] = 0;
    foreach_dimension() //Rotates over the dimensions
      dc[] +=  (F.x[] - F.x[1])/Delta; 
  }
}
/**
   With the tendency, the solution can be advanced in time with
   timestep `dt`.
*/

void advance (scalar c, scalar dc, double dt) {
  foreach()
    c[] += dt*dc[];
}

/**
Now we can construct a user-interface function, `diffusion()` that
employs a fordward-Euler scheme.
 */

void diffusion (scalar c, double dt, face vector kappa) {
  face vector F[];
  scalar dc[];
  flux_diffusion (c, kappa, F);
  tendency_from_flux (F, dc);
  advance (c, dc, dt);
}
/**
   At the cost of more computational effort and memory, we could also
   choose to use the mid-point rule to update the solution. It is
   second-order accurate in time!
 */

void diffusion_midpoint (scalar c, double dt, face vector kappa) {
  face vector F[];
  scalar dc[], c_temp[];
  foreach()
    c_temp[] = c[];                 //create a scratch
  flux_diffusion(c_temp, kappa, F); 
  tendency_from_flux (F, dc);
  advance (c_temp, dc, dt/2.);      //advance to the mid point.
  flux_diffusion(c_temp, kappa, F); //Re-estimate the fluxes at the mid point,
  tendency_from_flux (F, dc);       //and the tendency there.
  advance (c, dc, dt);              //update the solution
}

/**
# Further, 

We could use the [`Runke-Kutta`](/src/runge-kutta.h) schemes (upto 4th
order) or a [predictor corrector](/src/predictor-corrector.h) scheme.
 */



