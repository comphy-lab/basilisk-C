/**
# Shock tube with Riemann solver

We solve the Euler equations for a compressible gas. We also need to
compute volume fractions for the initial condition. */

#include "grid/multigrid.h"
#include "compressible.h"

double rhoL = 10., rhoR = 0.125;
double pL = 10., pR = 0.1;
double gamma1 = 1.4, gamma2 = 1.4;
double tend = 1.;

int main() {

  /**
  We make boundary conditions free outflow. */

  foreach_dimension() {
    w.n[right] = neumann(0);
    w.n[left]  = neumann(0);
  }
  
  L0 = 10;
  X0 = Y0 = -L0/2.;
  N = 256;
  run(); 
}

event init (t = 0)
{ 
  scalar f[];
  
  /**
  Left and right initial states for $\rho$, $\mathbf{w}$ and energy
  $E = \rho \mathbf{u}^2/2 + p/(\gamma-1)$. */
  
  foreach() {
    f[] = (x < 0.);
    rho[] = rhoL*f[] + rhoR*(1. - f[]);
    foreach_dimension()
      w.x[] = 0.;
    E[] = (pL*f[] + pR*(1. - f[]))/(gammao - 1.);
  }
}

event print (t = tend)
{

  foreach() 
    printf("%g %g %g %g \n", x, (gammao - 1.)*(E[] - 1./2.*rho[]*pow(w.x[]/rho[],2)) ,rho[], w.x[]/rho[]);
}

/* Theoretical solution */

double get_u4(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return (p1 - p2)*sqrt((1.-m)/(r*(p1 + m*p2))); 
}

double get_du4dp(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return sqrt((1.-m)/(r*(p1 + m*p2)))*(1. - 1./2.*(p1-p2)/(p1 + m*p2)); 
}


double get_u2(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return (pow(p1,(g-1)/(2*g)) - pow(p2,(g-1)/(2*g)))*sqrt((1.-m*m)*pow(p1,1/g)/(r*m*m)); 
}

double get_du2dp(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return - (g-1)/(2*g)*pow(p2,-(g+1)/(2*g))*sqrt((1.-m*m)*pow(p1,1/g)/(r*m*m)); 
}


double get_p3 ()
{

  double a,b;
  a = pR;
  b = (pR+pL)/2;
  double p3 = (a + b)/2;

  int Imax = 1000, i = 1;
  double error = 1.e10, tol = 1.e-3;

  //Newton-Raphson algorithm
  while (error > tol && i < Imax) {
    double u4 = get_u4(p3,pR,gamma2,rhoR);
    double u2 = get_u2(pL,p3,gamma1, rhoL);
    error = fabs(u4 - u2);
    p3 -= (u4 - u2)/(get_du4dp(p3,pR,gamma2,rhoR)- get_du2dp(pL,p3,gamma1, rhoL));
    if (p3 < a)
      p3 = a;
    else if (p3 > b)
      p3 = b;

    i++;
  }

  return p3;
}

event theoreticalsolution (t= tend) 
{
    FILE * fp1;
    fp1    = fopen ("theory.dat", "w");

    double p3 = get_p3();
    double p4 = p3;
    double u4 = get_u4(p4,pR,gamma2,rhoR);
    double u3 = u4;
    double m2 = (gamma2 - 1)/(gamma2 + 1);
    double m1 = (gamma1 - 1)/(gamma1 + 1);
    double rho4 = rhoR*(p4 + m2*pR)/(pR + m2*p4);
    double rho3 = rhoL*pow(p3/pL,1./gamma1);
    double ushock =  u4*rho4/(rho4-rhoR);

    double cson1 = sqrt(gamma1*pL/rhoL);
    double chi1 = -cson1;
    double chi2 = (u3/(1. - m1) - cson1);

    double u2, rho2, p2;
    foreach () { 
      double chi = x/t;
      if (chi < chi1) 
        fprintf(fp1, "%g %g %g %g \n", chi, pL, rhoL, 0.);
      else if (chi < chi2) {
        u2 = (1. - m1)*(chi + cson1);
        rho2 = pow(pow(rhoL,gamma1)/(gamma1*pL)*pow(u2 - chi,2),1./(gamma1-1.));
        p2 = pL*pow(rho2/rhoL, gamma1);
        fprintf(fp1, "%g %g %g %g \n", chi, p2, rho2, u2);
      }
    }
    u2 = (1. - m1)*(chi2 + cson1);
    rho2 = pow(pow(rhoL,gamma1)/(gamma1*pL)*pow(u2 - chi2,2),1./(gamma1-1.));
    p2 = pL*pow(rho2/rhoL, gamma1);
    fprintf(fp1, "%g %g %g %g \n", chi2, p2, rho2, u2);
    fprintf(fp1, "%g %g %g %g \n", u4, p3, rho3, u3);
    fprintf(fp1, "%g %g %g %g \n", u4, p4, rho4, u4);
    fprintf(fp1, "%g %g %g %g \n", ushock, p4, rho4, u4);
    fprintf(fp1, "%g %g %g %g \n", ushock, pR, rhoR, 0.);
    fprintf(fp1, "%g %g %g %g \n", 1.5*ushock, pR, rhoR, 0.);

    fclose(fp1);
}

/**
 *
~~~gnuplot Pressure profile
set output 'p.png'
set xrange[-5:5]
set xlabel 'x/t'
set ylabel 'p'
set cblabel 't'
p "theory.dat" u 1:2 not w l lc 0 lw 3,  "out" u 1:2 not w p 
~~~ 
 
~~~gnuplot Density profile
set output 'r.png'
set xrange[-5:5]
set xlabel 'x/t'
set ylabel '{/Symbol r}'
set cblabel 't'
p "theory.dat" u 1:3 not w l lc 0 lw 3,  "out" u 1:3 not w p
~~~ 

~~~gnuplot Velocity profile
set output 'u.png'
set xrange[-5:5]
set xlabel 'x/t'
set ylabel 'u'
set cblabel 't'
p "theory.dat" u 1:4 not w l lc 0 lw 3,  "out" u 1:4 not w p
~~~ 

*/