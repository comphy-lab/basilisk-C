#ifndef Stokes5_H
#define Stokes5_H

#define eta Stokes5_eta // rename functions
#define u Stokes5_u
#define v Stokes5_v

// #include<stdio.h>
// #include<math.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif



typedef struct{
/* Coefficients depending on wave and environmental properties */
  double A11, A22, A31, A33, A42, A44, A51, A53, A55,
  B22, B31, B42, B44, B53, B55,
  C0, C2, C4,
  D2, D4,
  E2, E4,
  S, eps, wave_length, wave_height, current, gravity, depth, k, kd, kH, c, cs, u_mean, Q, R,
  x0, y0;
} Stokes5;


/* Stokes 5th order waves based on Fenton's paper */

void set_stokes5_properties(Stokes5 *stokes5, double wave_length, double wave_height, double current, double depth, double gravity, double x0, double y0)

/* Compute coefficients based on wave properties and depth */
{

  double eps, S, k, kd;


  // compute basic properties, dimensional numbers, and set attributes of stokes5  
  k = 2.*PI/wave_length;
  kd = k*depth;
  eps = k*wave_height/2.;
  S = 1./cosh(2.*kd);

  stokes5->k = k;
  stokes5->kd = kd;
  stokes5->wave_height = wave_height;
  stokes5->wave_length = wave_length;
  stokes5->current = current;
  stokes5->depth = depth;
  stokes5->gravity = gravity;
  stokes5->kH = k*wave_height;
  stokes5->eps = eps;
  stokes5->S = S;
  stokes5->x0 = x0;
  stokes5->y0 = y0;

  // compute coefficients
  stokes5->A11 = 1/sinh(kd);
  stokes5->A22 = 3*pow(S,2)/(2*pow(1-S,2));
  stokes5->A31 = (-4 -20*S + 10*pow(S,2) - 13*pow(S,3))/(8*sinh(kd)*pow(1-S,3));
  stokes5->A33 = (-2*pow(S,2) + 11*pow(S,3))/(8*sinh(kd)*pow(1-S,3));
  stokes5->A42 = (12*S - 14*pow(S,2) - 264*pow(S,3) - 45*pow(S,4) - 13*pow(S,5))/(24*pow(1-S,5));
  stokes5->A44 = (10*pow(S,3) - 174*pow(S,4) + 291*pow(S,5) + 278*pow(S,6))/(48*(3+2*S)*pow(1-S,5));
  stokes5->A51 = (-1184 + 32*S + 13232*pow(S,2) + 21712*pow(S,3) + 20940*pow(S,4) + 12554*pow(S,5) - 500*pow(S,6) - 3341*pow(S,7) - 670*pow(S,8))/(64*sinh(kd)*(3+2*S)*(4+S)*pow(1-S,6));
  stokes5->A53 = (4*S + 105*pow(S,2) + 198*pow(S,3) - 1376*pow(S,4) - 1302*pow(S,5) - 117*pow(S,6) + 58*pow(S,7))/(32*sinh(kd)*(3+2*S)*(pow(1-S,6)));
  stokes5->A55 = (-6*pow(S,3) + 272*pow(S,4) - 1552*pow(S,5) + 852*pow(S,6) + 2029*pow(S,7) + 430*pow(S,8))/(64*sinh(kd)*(3+2*S)*(4+S)*pow(1-S,6));

  stokes5->B22 = (cosh(kd)/sinh(kd))*(1+2*S)/(2*(1-S));
  stokes5->B31 = -3*(1 + 3*S + 3*pow(S,2) + 2*pow(S,3))/(8*pow(1-S,3));
  stokes5->B42 = (cosh(kd)/sinh(kd))*(6 - 26*S - 182*pow(S,2) - 204*pow(S,3) - 25*pow(S,4) + 26*pow(S,5))/(6*(3+2*S)*pow(1-S,4));
  stokes5->B44 = (cosh(kd)/sinh(kd))*(24 + 92*S + 122*pow(S,2) + 66*pow(S,3) + 67*pow(S,4) + 34*pow(S,5))/(24*(3+2*S)*pow(1-S,4));
  stokes5->B53 = 9*(132 + 17*S - 2216*pow(S,2) - 5897*pow(S,3) - 6292*pow(S,4) - 2687*pow(S,5) + 194*pow(S,6) + 467*pow(S,7) + 82*pow(S,8))/(128*(3+2*S)*(4+S)*pow(1-S,6));
  stokes5->B55 = 5*(300 + 1579*S + 3176*pow(S,2) + 2949*pow(S,3) + 1188*pow(S,4) + 675*pow(S,5) + 1326*pow(S,6) + 827*pow(S,7) + 130*pow(S,8))/(384*(3+2*S)*(4+S)*pow(1-S,6));
  stokes5->C0 = pow(tanh(kd),0.5);
  stokes5->C2 = pow(tanh(kd),0.5)*(2+7*pow(S,2))/(4*pow(1-S,2));
  stokes5->C4 = pow(tanh(kd),0.5)*(4 + 32*S - 116*pow(S,2) - 400*pow(S,3) - 71*pow(S,4) + 146*pow(S,5))/(32*pow(1-S,5));
  stokes5->D2 = -pow(cosh(kd)/sinh(kd),0.5)/2;
  stokes5->D4 = pow(cosh(kd)/sinh(kd),0.5)*(2 + 4*S + pow(S,2) + 2*pow(S,3))/(8*pow(1-S, 3));
  stokes5->E2 = tanh(kd)*(2 + 2*S + 5*pow(S,2))/(4*pow(1-S,2));
  stokes5->E4 = tanh(kd)*(8 + 12*S - 152*pow(S,2) - 308*pow(S,3) - 42*pow(S,4) + 77*pow(S,5))/(32*pow(1-S,5));

  stokes5->u_mean = (stokes5->C0 + pow(eps,2)*stokes5->C2 + pow(eps,4)*stokes5->C4)/sqrt(k/gravity);
  stokes5->c = stokes5->u_mean + stokes5->current; // wave celerity
  stokes5->Q = (stokes5->u_mean)*depth + (stokes5->D2*pow(eps,2) + stokes5->D4*pow(eps,4))/sqrt(pow(k,3)/gravity);
  stokes5->cs = stokes5->c - (stokes5->Q)/depth; // cs = Stokes drift (set the current to Q/d (Stokes drift=0) for closed wave tanks
  stokes5->R = 0.5*(gravity/k)*pow(stokes5->C0,2) + kd + (stokes5->E2)*pow(eps,2) + (stokes5->E4)*pow(eps,4); // Bernoulli constant

}


double eta(Stokes5 *stokes5, double t, double X)

/* Given an initialized stokes wave (stokes5), compute the surface elevation eta at time t and coordinate X */
{
  double keta, kx, eps;
  eps = stokes5->eps;
  kx = (stokes5->k)*(X - (stokes5->c)*t - (stokes5->x0));
  keta = stokes5->kd + eps*cos(kx) + 
              pow(eps,2)*(stokes5->B22)*cos(2*kx) + 
              pow(eps,3)*(stokes5->B31)*(cos(kx) - cos(3*kx)) + 
              pow(eps,4)*((stokes5->B42)*cos(2*kx) + (stokes5->B44)*cos(4*kx)) + 
              pow(eps,5)*(-((stokes5->B53) + (stokes5->B55))*cos(kx) + (stokes5->B53)*cos(3*kx) + (stokes5->B55)*cos(5*kx));

  return stokes5->y0 + keta/(stokes5->k) - stokes5->depth;
}



double u(Stokes5 *stokes5, double t, double X, double Y)

/* Given an initialized stokes wave (stokes5), compute the horizontal velocity at time t and coordinate X, Y.  Values above the surface are computed uncritically. */
{

  double k, kx, ky, eps;
  eps = stokes5->eps;
  k = stokes5->k;
  kx = (stokes5->k)*(X - (stokes5->c)*t - (stokes5->x0));
  ky = (stokes5->k)*(Y + stokes5->depth - (stokes5->y0));

  return stokes5->c - stokes5->u_mean + 
         (stokes5->C0)*sqrt((stokes5->gravity)/pow(k,3))*(pow(eps,1)*(1*k*(stokes5->A11)*cosh(1*ky)*cos(1*kx)) +
         pow(eps,2)*(2*k*(stokes5->A22)*cosh(2*ky)*cos(2*kx)) + 
         pow(eps,3)*(1*k*(stokes5->A31)*cosh(1*ky)*cos(1*kx) + 3*k*(stokes5->A33)*cosh(3*ky)*cos(3*kx)) + 
         pow(eps,4)*(2*k*(stokes5->A42)*cosh(2*ky)*cos(2*kx) + 4*k*(stokes5->A44)*cosh(4*ky)*cos(4*kx)) + 
         pow(eps,5)*(1*k*(stokes5->A51)*cosh(1*ky)*cos(1*kx) + 3*k*(stokes5->A53)*cosh(3*ky)*cos(3*kx) + 
         5*k*(stokes5->A55)*cosh(5*ky)*cos(5*kx)));
}




double v(Stokes5 *stokes5, double t, double X, double Y)

/* Given an initialized stokes wave (stokes5), compute the vertical velocity at time t and coordinate X, Y.  Values above the surface are computed, but should later be set to zero. */
{

  double k, kx, ky, eps;
  eps = stokes5->eps;
  k = stokes5->k;
  kx = (stokes5->k)*(X - (stokes5->c)*t - (stokes5->x0));
  ky = (stokes5->k)*(Y + stokes5->depth - (stokes5->y0));

  return (stokes5->C0)*sqrt((stokes5->gravity)/pow(k,3))*(pow(eps,1)*(1*k*(stokes5->A11)*sinh(1*ky)*sin(1*kx)) +
         pow(eps,2)*(2*k*(stokes5->A22)*sinh(2*ky)*sin(2*kx)) + 
         pow(eps,3)*(1*k*(stokes5->A31)*sinh(1*ky)*sin(1*kx) + 3*k*(stokes5->A33)*sinh(3*ky)*sin(3*kx)) + 
         pow(eps,4)*(2*k*(stokes5->A42)*sinh(2*ky)*sin(2*kx) + 4*k*(stokes5->A44)*sinh(4*ky)*sin(4*kx)) + 
         pow(eps,5)*(1*k*(stokes5->A51)*sinh(1*ky)*sin(1*kx) + 3*k*(stokes5->A53)*sinh(3*ky)*sin(3*kx) + 
         5*k*(stokes5->A55)*sinh(5*ky)*sin(5*kx)));
}




#endif
