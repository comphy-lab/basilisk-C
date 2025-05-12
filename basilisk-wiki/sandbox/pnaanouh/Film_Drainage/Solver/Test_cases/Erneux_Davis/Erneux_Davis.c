/**
Test case following the problem presented by Erneux \& Davis 1993 of the relaxation of an infinite symmetric 2D film under the effect of surface tension and Van der Waals attraction
*/

#include "grid/multigrid1D.h"

//Due to how the macros are loaded for the surface forces van-der-waals.h includes hydro-tension.h that in turn includes hydro.h
//#include "../hydro.h"
#include "../van-der-waals.h"
//#include "../hydro-tension.h"
#include "../nh.h"
#include "../viscous_surface.h"
#include "layered/remap.h"



#define lay 2
#define delta 0.05 [1]
double T0 = 80.E4;
double k = 1.E-2 [-1];
double sg= 8E0  ;
double diff =0.0000;

int main() {
  
  N = 64;
  L0 = 2./k;
  G = 0. [1,-2];
  nl = lay;
  TOLERANCE = 1e-8 [*];
  CFL_H = 0.1;
  nu= 1. [2,-1];
  HAM =diff*sq(pi*k)/6.*1 [7,-2];
  
  h_diffusion= true;
  symmetric_bathymetry=true;
  
  //Simulation time is set to 12 characteristic periods
  T0=-24[0,1]/(-4 [2]*sq(k*pi)+sqrt(16[4]*sq(sq(k*pi))-6[2]*sq(k*pi)*(1 [6,0]*sq(k*pi)*sg/(3*sq(nu))*(1-diff))))/(nu*1[-2,1]);

  //Test  matrix: Outer loop Oh*2 every time and inner loop increases the effect of the Van der Waals force. 
  //Van der Waals force is taken to be proportional to the surface tension such that diff=0 no van der Waals and diff=1 film is stationary
  double inv;
  for (int j=0;j<4;j++){
    inv=1;
    for (int H=0;H<8;H++){
      diff=1.-inv;
      HAM = diff*sq(pi*k)*sg/6.*1 [7,0];
      T0=-24[0,1]/(-4 [2]*sq(k*pi)+sqrt(16[4]*sq(sq(k*pi))-6[2]*sq(k*pi)*(1 [6,0]*sq(k*pi)*sg/(3*sq(nu)))*(1-diff)))/(nu*1[-2,1]);
      run();
      inv /=2.;
    }
    sg/=2.;
  } 
  
  //run();
  
  
  
}
/**
#Initialization 
We the initialize the surface perturbation to a cosine and the velocity and pressure to the solution of the linear attenuation mode outlined by Erneux \& Davis 
*/

event init (i = 0)
{
  const scalar sig[]=sg;
  sigma=sig;
  periodic(right);
  double u_init=0;
  double w_init=0;
  foreach() {
    zb[] = 0. [1];
    u_init = -2. [-1,-1]*((-4 [2]*sq(k*pi)+sqrt(16[4]*sq(sq(k*pi))-6[2]*sq(k*pi)*(1 [6,0]*sq(k*pi)*sg/(3*sq(nu)))*(1-diff)))/(nu*1[-2,1]))/2.*sin(k*x*pi)/(k*pi)*delta;
    w_init = +2. [-1,-1]*((-4 [2]*sq(k*pi)+sqrt(16[4]*sq(sq(k*pi))-6[2]*sq(k*pi)*(1 [6,0]*sq(k*pi)*sg/(3*sq(nu)))*(1-diff)))/(nu*1[-2,1]))/2.*cos(k*x*pi)*delta;
    foreach_layer(){
      h[] = (0.5 [1] +delta*cos(k*x*pi))/nl;
      u.x[]=u_init;
      w[]=w_init*(0.5 [1] +delta*cos(k*x*pi))*(point.l+0.5)/nl;
      phi[]=nu*w_init;
    }
  }
}

/**
#Outputs
We output the temporal evolution of the vertical velocity, perturbation amplitude, and the time derivative of the log of the perturbation amplitude
*/

event v_out (t+=T0/1000)
{
  char name[80];
  sprintf (name, "profile_w_nu_%g_k_%g_sig_%g_diff_%g_N_%d_L_%d", nu,k,sg,diff,N,nl);
  static FILE * fpp1 = fopen (name, "w");
  foreach(){
    fprintf (fpp1, "%g %g",x,t);
    //foreach_layer()
    fprintf (fpp1, " %g",w[]);
    fprintf (fpp1, "\n");
  }
  fflush (fpp1);
}

event amplitude (t += T0/1000) {
  double max = statsf(eta).max;
  double mean = normf(eta).avg;
  
  char name[80];
  sprintf (name, "ampSimu1_nu_%g_k_%g_sig_%g_diff_%g_N_%d_L_%d", nu,k,sg,diff,N,nl);
  static FILE * fpw = fopen (name, "w");
  fprintf (fpw, "%g %g %g\n", t, (max-0.5)/0.5,(mean-0.5)/0.5);
  fflush (fpw);

  }


double pt=0;
double pmeta;

event slope (i += 10) {
  double max = statsf(eta).max;
  if (pt !=0) { 
    char name[80];
    sprintf (name, "slopeSimu4_nu_%g_k_%g_sig_%g_diff_%g_N_%d_L_%d", nu,k,sg,diff,N,nl);
    static FILE * fpsm = fopen (name, "w");
    fprintf (fpsm, "%g %g %g\n", t, t-pt,log((max-0.5)/0.5)-log((pmeta-0.5)/0.5));
    fflush (fpsm);
  }
  pt=t;
  pmeta=max;
}

event end (t = T0);

