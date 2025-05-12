/**
# Shock tube problem for a single ideal gas (strong shock wave): Modified formulation

*/

#include "grid/multigrid.h"
#include "all-mach.h"

////////////////////////////////////////////////////////////

#include "vof.h"
#include "tension.h"

scalar c[], rhov[], rho1[], rho2[], Ei1[], Ei2[], * interfaces = {c}, * interfaces1 = {c};

face vector alphav[];

scalar rhoc2v[];

double CFLacous = 0.5;

double rhoL = 10., rhoR = 0.125;
double pL = 10., pR = 0.1;
double gammaL = 1.4, gammaR = 1.4;
double tend = 1.;


event vof (i++) {
  vector q1 = q, q2[];
  foreach()
    foreach_dimension() {
      double u = q.x[]/rho[];
      q1.x[] = rho1[]*u;
      q2.x[] = rho2[]*u;
    }
  boundary ((scalar *){q1,q2});
  theta = 1.;
  foreach_dimension() {
    q2.x.inverse = true;
    q1.x.gradient = q2.x.gradient = minmod2;
  }
  rho2.inverse = true;
  rho1.gradient = rho2.gradient = minmod2;
  Ei2.inverse = true;
  Ei1.gradient = Ei2.gradient = minmod2;
  c.tracers = {rho1,rho2,q1,q2,Ei1,Ei2};
  vof_advection ({c}, i);
  foreach() 
    foreach_dimension() 
      q.x[] = q1.x[] + q2.x[];

  boundary ((scalar *){q});

  vector u[];
  foreach() 
    foreach_dimension() 
      u.x[] = q.x[]/(rho1[] + rho2[]);
  boundary ((scalar *){u});

  foreach() {
    double div = 0.;
    foreach_dimension() 
      div += u.x[1] - u.x[];
    Ei1[] -= p[]*div/Delta*c[]*dt;
    Ei2[] -= p[]*div/Delta*(1. - c[])*dt;
 }

  interfaces = NULL;
}

event acceleration (i++) {
  interfaces = interfaces1;
}

double cstate (double rho, double rhoref, double pref,
         double B, double gamma) {
  return (1. + B)*pref/rhoref*gamma*pow(rho/rhoref, gamma - 1.);
}

event stability (i++) {

  double dt;
  foreach() {
    if (c[] > 0.5) {
      dt = Delta/sqrt(cstate (rho1[]/c[], rhoL, pL, 0., gammaL));
    } else {
      dt = Delta/sqrt(cstate (rho2[]/(1. - c[]), rhoR, pR, 0., gammaR));
    }
    dt *= CFLacous;
    if (dt < dtmax)
      dtmax = dt;
  }

}

event properties (i++) {
  alpha = alphav;
  rhoc2 = rhoc2v;
  rho = rhov;
  
  foreach() {
    rhov[] = max(rho1[] + rho2[], 1e-6);
    ps[] = (gammaR - 1.)*(Ei1[]+Ei2[]);
    rhoc2v[] = gammaR*ps[];
  }
  boundary ({rhov, rho1,rho2});
  
  foreach_face()
    alphav.x[] = 2./(rho[] + rho[-1]);
}

FILE * fp;

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2.;
  N = 256;
  fp    = fopen ("numerical.dat", "w");
  run();
  fclose(fp);
}

event defaults (i = 0) {
  foreach() {
    rho1[] = rho2[] = c[] = 1.;
    Ei1[] = Ei2[] = 0.;
  }
  boundary ({rho1,rho2,c, Ei1, Ei2});
}

event init (i = 0) {
  foreach() {
    c[] = (x < 0.);
    rho1[] = c[]*rhoL;
    rho2[] = (1. - c[])*rhoR;
    Ei1[] = c[]*pL/(gammaL-1.);
    Ei2[] = (1. - c[])*pR/(gammaR-1.);
  }
  boundary ({rho1,rho2,Ei1,Ei2});
}

event outputdata (t= tend)
//event outputdata (i += 1; t= tend)
{
  foreach () 
      fprintf(fp, "%g %g %g %g %g %g \n", x, t, p[], rho1[]+rho2[], q.x[]/(rho1[]+rho2[]), Ei1[]+Ei2[]);

  fprintf(fp, " \n");
}

/** Theoretical solution */

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
    double u4 = get_u4(p3,pR,gammaR,rhoR);
    double u2 = get_u2(pL,p3,gammaL, rhoL);
    error = fabs(u4 - u2);
    p3 -= (u4 - u2)/(get_du4dp(p3,pR,gammaR,rhoR)- get_du2dp(pL,p3,gammaL, rhoL));
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
    double u4 = get_u4(p4,pR,gammaR,rhoR);
    double u3 = u4;
    double mR = (gammaR - 1)/(gammaR + 1);
    double mL = (gammaL - 1)/(gammaL + 1);
    double rho4 = rhoR*(p4 + mR*pR)/(pR + mR*p4);
    double rho3 = rhoL*pow(p3/pL,1./gammaL);
    double ushock =  u4*rho4/(rho4-rhoR);

    double csonL = sqrt(gammaL*pL/rhoL);
    double chi1 = -csonL;
    double chi2 = (u3/(1. - mL) - csonL);

    double u2, rho2, p2;
    foreach () { 
      double chi = x/t;
      if (chi < chi1) 
        fprintf(fp1, "%g %g %g %g \n", chi, pL, rhoL, 0.);
      else if (chi < chi2) {
        u2 = (1. - mL)*(chi + csonL);
        rho2 = pow(pow(rhoL,gammaL)/(gammaL*pL)*pow(u2 - chi,2),1./(gammaL-1.));
        p2 = pL*pow(rho2/rhoL, gammaL);
        fprintf(fp1, "%g %g %g %g \n", chi, p2, rho2, u2);
      }
    }
    u2 = (1. - mL)*(chi2 + csonL);
    rho2 = pow(pow(rhoL,gammaL)/(gammaL*pL)*pow(u2 - chi2,2),1./(gammaL-1.));
    p2 = pL*pow(rho2/rhoL, gammaL);
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
p "theory.dat" u 1:2 not w l lc 0 lw 3,  "./numerical.dat" u ($1/$2):3 not w l lw 3
~~~ 
 
~~~gnuplot Density profile
set output 'r.png'
set xrange[-5:5]
set xlabel 'x/t'
set ylabel '{/Symbol r}'
set cblabel 't'
p "theory.dat" u 1:3 not w l lc 0 lw 3,  "./numerical.dat" u ($1/$2):4 not w l lw 3
~~~ 

~~~gnuplot Velocity profile
set output 'u.png'
set xrange[-5:5]
set xlabel 'x/t'
set ylabel 'u'
set cblabel 't'
p "theory.dat" u 1:4 not w l lc 0 lw 3,  "./numerical.dat" u ($1/$2):5 not w l lw 3
~~~ 

 */  
