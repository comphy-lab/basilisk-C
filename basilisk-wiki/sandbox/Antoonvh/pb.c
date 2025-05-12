/**
# A dirichlet BC trick for the Poisson solver

We solve $\nabla^2 f = s$ for $f$ and want to prescribe some value at
the top boundary. *Normally* it would look like this,
*/
#include "poisson.h"
#include "utils.h"

scalar f[], s[];

int main() {
  periodic (left);
  L0 = 10;
  X0 = -5.;
  Y0 = -2.5;
  init_grid (128);
  foreach()
    s[] = exp(-sq(x) - sq(y));
  output_ppm (s, file = "s.png", n = 256, linear = true);
  /**
![The source term looks like this](pb/s.png)

We set boundary conditions and solve the problem
  */
  f[top] =  dirichlet (1.);
  f[bottom] = dirichlet (0.);
  poisson (f, s);
  output_ppm (f, file = "f.png", n = 256,
	      linear = true, min = -2, max = 1.5);
  /**
     ![The solution for $f$ on the square domain](pb/f.png)
   
## Boundary within the domain 
     
     We want to prescribe the value of $f$ halfway the domain (i.e. $y =$
     `Y0 + L0/2`) and for some reason we do not use mask or embedded
     boundaries. We must then first set a homogeneous dirichlet
     condition at the top, so that we solve for a modified field $f_{\mathrm{mod}} =
     f-1$.
*/
  f[top] =  dirichlet (0.);
  f[bottom] = dirichlet (-1.);
  /**
The dirichlet condition is extented into the domain via a face vector
field and a scalar field:
   */
  face vector as[];
  scalar ls[];
# define OUTSIDE_DOMAIN (y >= (Y0 + L0/2.)) 
  foreach() {
    if (OUTSIDE_DOMAIN) { // Outide the internal domain:
      s[] = 0;            // Source term is zero
      ls[] = -10;         // Linear term 
    }else
      ls[] = 0;           // No linear term
  }
  foreach_face()          // Set the default face weights 
    as.x[] = 1;
  foreach_face(x)         // Except for the boundary-perpendicular component,
    if (OUTSIDE_DOMAIN)   // outside the internal domain.
      as.x[] = 0.;
  boundary((scalar*) all);
  /**
     We now solve the modified problem:
  */
  poisson (f, s, as, ls);
  output_ppm (f, file = "f2.png", n = 256,
	      linear = true, min = -3, max = 0.5);
  /**
![The solution looks as if it could be true](pb/f2.png)

we compare it against the proper mask implementation,
  */
  mask (OUTSIDE_DOMAIN ? top : none);
  boundary(all);
  poisson (f, s);
  output_ppm (f, file = "f3.png", n = 256,
	      linear = true, min = -3, max = 0.5);
  /**
  ![Its OK-ish](pb/f3.png)
  */  
  
}
