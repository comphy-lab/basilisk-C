/**
# Velocity potential solver for the solid phase velocity field
The calculation of the solid phase velocity field 'ubf' is performed through the
solution of a Poisson equation for the velocity potential 'psi'.
The source term of the Poisson equation is given by the product of the
reaction rate 'omega', the volume fraction field 'f' and the
shrinkage factor 'zeta', divided by the solid phase density 'rhoS'.
*/

extern scalar omega;
extern scalar zeta;
extern double rhoS;

#include "diffusion.h"

bool shift_prod = false;
bool diffuse_prod = false;

scalar psi[];
face vector ubf[];
mgstats mgpsf;
double TOLERANCE_SOLID = 1e-5;
extern face vector fsS;

void shift_field (scalar fts, scalar f, int dir);

scalar prod[];
extern scalar fS;

/**
## Projection method for the solid phase velocity field
We solve the Poisson equation for the velocity potential 'psi' and
then compute the solid phase velocity field 'ubf' as the negative gradient
of 'psi'.
*/

trace
mgstats project_sv (face vector ubf, scalar psi,
    (const) face vector alpha = unityf,
    int nrelax = 4)
{
  /**
  We compute the source term for the Poisson equation
  */

  foreach()
    prod[] = omega[]*f[]*zeta[]*cm[]/rhoS;
 
  /**
  We optionally shift and/or diffuse the source term to improve stability
  when the heat exchange at the interface is high and the reaction rate
  is very localized.
  */

  if (shift_prod) 
    shift_field (prod, f, 1);
  
  if (diffuse_prod) {
    face vector D[];
    scalar theta[];

    foreach()
      theta[] = cm[]*max(fS[], F_ERR);

    foreach_face()
      D.x[] = fm.x[]*fsS.x[]*1e-5;

    diffusion (prod, dt, D=D, theta=theta);
  }

  mgstats mgp = poisson (psi, prod, alpha,
      tolerance = TOLERANCE_SOLID, nrelax = nrelax);

  foreach_face()
    ubf.x[] = -alpha.x[]*face_gradient_x (psi, 0);

  return mgp;
}

/**
## Boundary conditions
*/

psi[right] = neumann (neumann_pressure(ghost));
psi[left]  = neumann (- neumann_pressure(0));

#if AXI
ubf.n[bottom] = 0.;
ubf.t[bottom] = dirichlet(0);
psi[top]    = neumann (neumann_pressure(0));
#else // !AXI
#  if dimension > 1
psi[top]    = neumann (0);
psi[bottom] = neumann (0);
#  endif
#  if dimension > 2
psi[front]  = neumann (0);
psi[back]   = neumann (0);
#  endif
#endif // !AXI

event defaults (i=0) {
  psi.nodump = true;
  #if TREE
  ubf.x.refine = refine_face_solenoidal;
  #endif
}
