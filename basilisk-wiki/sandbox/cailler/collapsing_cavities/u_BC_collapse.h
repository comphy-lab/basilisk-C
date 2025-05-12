/** 
This file is a copy of the code
[`u_BC_dipolar_flow.h`](http://basilisk.fr/sandbox/cailler/sierou_lister/u_BC_dipolar_flow.h). 
Please refer to the link provided for a detailed documentation on this file.
*/


#include "elliptic.h"

#define PHI_0_A phi_0_A(THETA_0, MU_0)
#define PHI_0_OQ phi_0_OQ(THETA_0, MU_0)


                            /* LIQUID PHASE */

double uxA (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*cos(th)*legendreP_half(th);
  double B = sin(th)*legendreP_prime_half(th);

  return ( PHI_0_A*(A - B)/sqrt(R) );
}

double uyA (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*sin(th)*legendreP_half(th);
  double B = cos(th)*legendreP_prime_half(th);

  return ( PHI_0_A*(A + B)/sqrt(R) );
}

// ∂uxA/∂x
double duxA_dx (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 8.*sq(R)*legendreP_half(th);
  double B = 24.*x*R*legendreP_three_half(th);
  double C = 15.*sq(R)*legendreP_five_half(th);

  return ( 0.25*PHI_0_A*( A - B + C )/pow(R, 7./2.) );
}

// ∂uxA/∂y
double duxA_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 4.*x*sq(R)*legendreP_half(th);
  double B = R*( 9.*sq(x) + sq(y) )*legendreP_three_half(th);
  double C = 5.*x*sq(R)*legendreP_five_half(th);

  return ( -0.75*PHI_0_A*( A - B + C )/( y*pow(R, 7./2.)) );
}

// ∂uyA/∂x
double duyA_dx (double x, double y){
  return ( duxA_dy(x,y) );
}

// ∂uyA/∂y
double duyA_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = ( 15.*pow(x, 4.) + 14.*sq(x)*sq(y) - pow(y, 4.) )*legendreP_half(th);
  double B = -2.*R*( 5.*sq(x) + sq(y) )*legendreP_three_half(th);
  double C = 5.*x*sq(R)*legendreP_five_half(th);

  return ( 0.25*PHI_0_A*( A + 3.*x*( B + C ) )/( sq(y)* pow(R, 7./2.)) );
}


                              /* GAS PHASE */

double uxO (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*cos(th)*legendreQ_half(th);
  double B = sin(th)*legendreQ_prime_half(th);

  return ( PHI_0_OQ*(A - B)/sqrt(R) );
}

double uyO (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*sin(th)*legendreQ_half(th);
  double B = cos(th)*legendreQ_prime_half(th);

  return ( PHI_0_OQ*(A + B)/sqrt(R) );
}

// ∂uxO/∂x
double duxO_dx (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 8.*sq(R)*legendreQ_half(th);
  double B = 24.*x*R*legendreQ_three_half(th);
  double C = 15.*sq(R)*legendreQ_five_half(th);

  return ( 0.25*PHI_0_OQ*( A - B + C )/pow(R, 7./2.) );
}

// ∂uxO/∂y
double duxO_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 4.*x*sq(R)*legendreQ_half(th);
  double B = R*( 9.*sq(x) + sq(y) )*legendreQ_three_half(th);
  double C = 5.*x*sq(R)*legendreQ_five_half(th);

  return ( -0.75*PHI_0_OQ*( A - B + C )/( y*pow(R, 7./2.)) );
}

// ∂uyO/∂x
double duyO_dx (double x, double y){
  return ( duxO_dy(x,y) );
}

// ∂uyO/∂y
double duyO_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = ( 15.*pow(x, 4.) + 14.*sq(x)*sq(y) - pow(y, 4.) )*legendreQ_half(th);
  double B = -2.*R*( 5.*sq(x) + sq(y) )*legendreQ_three_half(th);
  double C = 5.*x*sq(R)*legendreQ_five_half(th);

  return ( 0.25*PHI_0_OQ*( A + 3.*x*( B + C ) )/( sq(y)* pow(R, 7./2.)) );
}

