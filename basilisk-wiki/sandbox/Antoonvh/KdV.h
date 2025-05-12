/**
# A solver for the Korteweg-De Vries Equation

The Korteweg-De Vries(KdV) equation reads:

$$\frac{\partial c}{\partial t} = 6c\frac{\partial c}{\partial x} -
\frac{\partial^3 c}{\partial x^3} $$

When we integrate over a *finite volume* ($V$) the evolution equation
for the cell averaged value of $c$ (C) reads $\frac{\partial
C}{\partial t} = \frac{1}{V}\Sigma F_i\cdot A_i$, where the right-hand
side concerns a summation of fluxes over the cell's face areas
($A_i$). For $F$ we employ the divergence theorem and write:

$$F_i = 3c_i^2-\left(\frac{\partial^2 c}{\partial x^2}\right)_i$$

These fluxes can be *aproximated* from the spatial structure of $C(x)$
in the neighborhood of a cell's face employing a discretized approach
and the numerical solution of $C$ can be advanced in time, using a
discrete timestep (`dt`) according to an fancy implicit, symplectic
and eight-order accurate Gauss-Legendre time-integration method.
 */
#define BGHOSTS 2
#define RKORDER 8
#include "GLrk.h"
/**
   This function computes the tendencies with second order spatial
   accuracy. This requires two layers of ghost cells:
*/
double six = 6;
static void KdV2 (scalar * s, scalar *kl){
  //boundary(s);
  for (int j = 0; j < list_len(s); j++){
    scalar m = s[j];
    scalar k = kl[j];
    face vector flx[];
    foreach_face()
      flx.x[] = 0.75*(six/6)*sq((m[] + m[-1])) - (m[-2] - m[-1] - m[] + m[1])/(2.*sq(Delta));
    //boundary_flux({flx});
    foreach(){
      k[] = 0;
      foreach_dimension()
	k[] += (flx.x[1] - flx.x[])/Delta;
    }
  }
}
/**
   The forward-in-time integrator is prone to instabilities and since
   the user may not be aware of any relevant details we take some care
   here: The discretized timestep `Dt` can be compared against the
   grid-cell size ($\Delta$) and the tendency from each term. For an
   equation of the form:
   
   $$\partial_tc = ac\partial_xc$$ 
   
   with $a$ a constant with units $[a] = LT^{-1}[c]^{-1}$. thus we
   formulate a Courant number:
   
   $$Co = \frac{a\mathrm{Dt}c}{\Delta}$$
   
   Similarly for an equation of the form,
   
   $$\partial_tc = b\partial_x^3c,$$
  
   with b a constant with units $[b] = L^3T^{-1}$ we define the KdV
   number ($KdVn$):
   
   $$KdVn = \frac{\mathrm{Dt}b}{\Delta^3}.$$
   
   $\mathrm{Dt}$ should be chosen such that both $Co$ and $KdVn$ remain within
   limited bounds. By default we use:
*/
double Co = .8;
double KdVn = .5;
/**
If the user-requested timestep `dt`$>\mathrm{Dt}$, the time
integration is split into `it` equal parts to ensure a limited value
of $Co$ and $KdVn$. The function below provides the described user
interface and returns the number of iterations (`it`). The input is a
list of scalar fields and the time step `dt`.
*/
int KdV(scalar * s, double dt){
  double delt = L0/(double)(1 << depth());
  double maxs = 1E-30;
  for (scalar m in s){
    double max = fabs(statsf(m).max);
    if (max > maxs)
      max = maxs;
  }
  double Dtmin = min(Co*delt/(six*maxs), KdVn*cube(delt));
  int it = 1;
  if (Dtmin < dt)
    it = (int)((dt/(Dtmin))+0.99);
  double Dt = dt/((double)(it));
  for (int iter = 0; iter < it; iter++)
    A_Time_Step(s, Dt, KdV2, 1e-5); 
  return it;
}
