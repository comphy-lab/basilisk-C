/** Test case case following the problem outlined by Howell \& Stone 2005 about the lack of marginal pinching for viscous (stokes flow) 2D films attached to a plateau border 
*/
/**##Setup*/
#include "grid/multigrid1D.h"
#include "../hydro-tension.h"
#include "../nh.h"
#include "../viscous_surface.h"
#include "layered/remap.h"



#define lay 40
#define LEVEL 8

double T0 = 2E3;
//Length of domain and Plateau boder curvature (linearised)
double L = 20E0 [1];
double detadx2 = 1. [-1];

//Transverse integral for vertical velocity to compute the vertical velocity at the domain boundary 
//from the continuity equation and the horizontal velocity BC here a constant flux is hard coded 
double w_integral(Point point,int l,scalar s) {
  if (l>0)
    return (s[0,0,l]+w_integral(point,l-1,s));
  else
    return s[0,0,l];
  
}

// Macro for the definition of the 2nd layer of ghost cell functions
@define ghost_layer (neighbor.i < GHOSTS ? (GHOSTS - neighbor.i) : neighbor.i - (1 << level) - GHOSTS + 1) // gives 1 for the first layer, 2 for the second, deepest, layer



int main() {
  
  N = 1 << LEVEL;
  L0 = L;
  G = 0. [1,-2]; //no gravity
  nl = lay;
  TOLERANCE = 1e-12 [*];
  CFL_H = 0.1;
  //CFL_H = 0.04;
  nu= 0.01 [2,-1];
  
  h_diffusion= true; //horizontal diffusion is activated
  symmetric_bathymetry = true; //bathyletry is a flat symmetry plane
  max_slope = 10; // max slope in increased to 10 -> 84° ugly bodge not really a solution for a skewed mesh
  
  run();
  
  
  
}

event init (i = 0)
{
  foreach() {
    zb[] = 0. [1];
    foreach_layer(){
      h[] = (0.5 [1]+ (x<(L0/64.)? sq((L0/64-x))*1.[-1]:0))/nl;
    }
  }
  
  //Film is thickness is constant to the right
  eta[right]=0.5 [1];
  h[right]=0.5 [1]/nl;
  
  //Film thickness is parabolic to the left (second layer of ghost cells does nothing)
  eta[left]=(2.-ghost_layer)*(sq(Delta)*detadx2+2.*eta[]-eta[1])+(ghost_layer-1.)*(3.*sq(Delta)*detadx2+3.*eta[-1]-2.*eta[]);
  h[left]=(2.-ghost_layer)*(sq(Delta)*detadx2/nl+2.*h[]-h[1])+(ghost_layer-1.)*(3.*sq(Delta)*detadx2/nl+3.*h[-1]-2.*h[]);
  
  //Zero gradient condition on horizontal velocity to the right forces a zero condition on vertical velocity
  u.n[right]=neumann(0.);
  w[right]=dirichlet(0.);
  
  //Constant horizontal flux  condition on the left
  u.n[left]=u.n[]*h[]/(sq(Delta)*detadx2/nl+2.*h[]-h[1]);
  w[left]=-((u.n[]/2.+w_integral(point,point.l-1,u.n))/Delta)*(sq(Delta)*detadx2/nl+h[]-h[1])*(sq(Delta)*detadx2/nl+3.*h[]-h[1])/(sq(Delta)*detadx2/nl+2.*h[]-h[1])-w[];
}

/**
Boundary conditions for face acceleration and flux are in principle not necessary since these fields are computed at the domain boundaries from the boundary conditions on h and u. In practice though imposing the zero acceleration condition by setting a constant curvature condition using the second layer of ghost cells is too complicated. 
*/
event half_advection (i++)
{
  
  ha.n[left] = 0.;
  ha.n[right] = 0.;
  //The zero neumann conditon is a direct result of the set u and h BCs so this is commented 
  //hu.n[left]=neumann(0.);
  //hu.n[right]=neumann(0.);
  
  //boundary ((scalar *){ha, hu});
}

/**
##Outputs
*/

event eta_profile (t+=0.01) {
  
  char name[80];
  sprintf (name, "profile_eta_Oh_%g_L_%g_N_%d_nl_%d",nu,L0,N,nl);
  static FILE * fpp = fopen (name, "w");
  if (i ==0){
    fprintf (fpp, "%g",t);
    foreach()
      fprintf (fpp, " %g",x);
    fprintf (fpp, "\n");
  }
  fprintf (fpp, "%g",t);
  foreach(){
    //fprintf (fpp, "%g",x);
    //foreach_layer()
    fprintf (fpp, " %g",eta[]);
    //fprintf (fpp, "\n");
  }
  fprintf (fpp, "\n");
  fflush (fpp);

}

event end (t = T0);