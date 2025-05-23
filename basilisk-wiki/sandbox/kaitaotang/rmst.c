/**
# Compressible Richtmyer-Meshkov instability with surface tension 

We use the compressible solver with VOF interface tracking and surface tension developed in Fuster and Popinet (2018). */

#include "grid/multigrid.h"
#include "two-phase-compressible.h"
#include "compressible-tension.h"

/**
   Set length scale and resolution for the problem. */
#define D 1.0
#define LEVEL 9
/**
   Define Mach number and Weber numbers. */
#define M 1.2
#define We 6

/**
   Some useful reference quantities. */
double Ms = M;
double rhoR = 0.1;
double rhoW = 1;
double pR = 1.;
double tend = 4;
double gam = 1.4;

/**
   Calculate the postshock conditions based on the Mach number. */
double prat(double Ms, double gam){
  return 2.*gam*sq(Ms)/(gam+1.) - (gam-1.)/(gam+1.);
}

double rhorat(double Ms, double gam){
  return (gam + 1.)*sq(Ms)/((gam - 1.)*sq(Ms) + 2.);
}

double urat(double Ms, double gam){
  return 1./rhorat(Ms, gam);
}

double us(double Ms, double gam, double rhoR, double pR){
  return Ms*sqrt(gam*pR/rhoR);
}

/**
   Some useful fields. */
FILE * fp = NULL;

double rho1, rho2, A, dv;
double p0R = 1.;

double rhoL, uL, pL, ushock, qL, EL, eta0, l, x0;
scalar denslap[];
vector densgrad[];
scalar adensgrad[];
scalar schl[];
scalar totdens[];
scalar vort[];

/**
   CFL number needs to be reduced to 0.25 for obtaining physical results at low Mach numbers. */
double CFLac = 0.25;

/**
   Make sure we meet the CFL condition by correctly estimating maximum signal speeds (|u|+|c|). */
event stability (i++) {
  stats pstats = statsf(p);
  stats fr1stats = statsf(frho1);
  stats qstats = statsf(q.x);
  
  double cson = fmax(sqrt(gamma1*(pstats.max)/rhoL), sqrt(gamma2*(pstats.max + PI2)/rhoW));
  double u_c = fmax(fabs(qstats.max/rhoL), fabs(qstats.max/rhoW));

  foreach () {
      DT    = L0*CFLac/(cson + u_c)/pow(2,LEVEL);
      dtmax = L0*CFLac/(cson + u_c)/pow(2,LEVEL);
  }
}

/**
   Boundary conditions.*/
p[left] = dirichlet(pL);
frho1[left] = dirichlet(rhoL);
q.n[left] = dirichlet(qL);
q.t[left] = dirichlet(0.);
fE1[left] = dirichlet(EL);
q.n[bottom] = neumann(0.);
q.n[top] = neumann(0.);
p[top] = neumann(0.);
p[bottom] = neumann(0.);
frho1[top] = neumann(0.);
frho1[bottom] = neumann(0.);

/**
   Main routine, set all the relevant parameters and domain size, etc. */

int main() {
  gamma1 = gam;
  gamma2 = gam;
  PI2 = 0.;
  l=D;

  eta0= 0.01*l;
  rhoL = rhorat(Ms, gam)*rhoR;
  pL = pR*prat(Ms, gam);
  ushock = us(Ms, gam, rhoR, pR);
  uL = ushock*(1.-urat(Ms, gam));
  qL = uL*rhoL;
  EL = pL/(gam-1.) + 0.5*sq(qL)/rhoL;
  fprintf(ferr, "%g %g %g %g\n", pL, rhoL, qL, EL);
  periodic(top);

  /** Size of the domain: */
  size (11.*D);
  x0=-eta0-L0/324; origin (x0, -L0/11.);

  /**
  The density is variable. */

  rho1 = rhoR;
  rho2 = rhoW;
  A=(rho2-rho1)/(rho2+rho1);

  /**
     Set viscosity (mu1, mu2) and surface tension (f.sigma) via pre-shock Weber numbers. */
  dv = Ms*sqrt(gam*pR/rhoR)-sqrt(gam*pR*prat(Ms,gam)/(rhoR*rhorat(Ms,gam)))*sqrt(((gam-1)*sq(Ms)+2)/(2*gam*sq(Ms)-(gam-1))); 
  fprintf(ferr, "dv=%g\n", dv);
  f.sigma = (rho1+rho2)*sq(A*dv)*eta0/We;

  mu1=0; mu2=0;
  f.gradient = frho1.gradient = frho2.gradient = p.gradient = q.x.gradient = q.y.gradient = minmod;

  N = (1 << LEVEL)*11;
  fprintf(ferr, "uL=%g %g\n", uL, ushock);

  /**
     We open a file indexed by the level to store the time evolution of
     the kinetic energy. */

  char name[80];
  sprintf (name, "k-%d", LEVEL);
  fp = fopen (name, "w");
  run();
  fclose (fp);

  /**
  We use *grep* to filter the lines generated by gnuplot containing
  the results of the fits (see below). */
  system ("grep ^fit out >> log");
}

event init (i = 0) {
  fprintf(ferr, "start\n");   
  fprintf(ferr, "l=%g\n", l); fprintf(ferr, "L0=%g\n", L0); fprintf(ferr, "sigma=%g\n", f.sigma);
  fprintf(ferr, "mu=%g\n", mu1);

  /**
  We initialise the shape of the interface. */
  fraction (f, -x+eta0*cos(2*pi*y/l) );
  boundary({f});

  foreach() {
    double pLap = p0R;
    double m = (x < -7.);
    p[] = f[]*(m*pL + (1. - m)*pR) + (1.-f[])*pLap;
    frho1[] = f[]*(m*rhoL + (1. - m)*rhoR);
    frho2[] = (1. - f[])*rhoW;
    q.x[] = m*frho1[]*uL;
    q.y[] = 0.;
    fE1[] = f[]*(p[]/(gam - 1.) + 0.5*pow(q.x[],2)/rhoL);
    fE2[] = (1.-f[])*(pLap + PI2*gamma2)/(gamma2-1.);
  }
  boundary ((scalar *){q.x, q.y, frho1, frho2, p, fE1, fE2});
}

/**
At each timestep we output the kinetic energy. */

//event logfile (i++; t <= 3.0) {
event logfile (i++; t <= 12.0) {
  stats pstats = statsf(p);
  double ke = 0.;
  foreach (reduction(+:ke))
    ke += sq(Delta)*(sq(q.x[]) + sq(q.y[]))/(frho1[]+frho2[]);
  fprintf (fp, "%g %g %g %d\n", t, ke, mgp.i, pstats.max);
  fflush (fp);
}

event outputInterface(t += 0.02) {
  char resultname[100];
  sprintf( resultname, "results%4.2f_%d.dat", t, pid() );
  FILE * fp = fopen(resultname, "w");
  scalar xpos[];
  scalar ypos[];
  position (f, xpos, {1, 0});
  position (f, ypos, {0, 1});
  foreach()
      if (xpos[] != nodata)
        fprintf (fp, "%g %g\n", xpos[], ypos[]);
  fclose (fp);
}

event output (t = 0.00; t += 0.01){
  srand(time(NULL));
  scalar vort[];
  vector uu[];
  dump("dump");
  fprintf(ferr, "t=%g\n", t);

  
  foreach()
    {
      totdens[] = frho1[] + frho2[];
      vort[] = (q.y[1] - q.y[-1])/(2.*Delta) - (q.x[0,1]-q.x[0,-1])/(2.*Delta);
      foreach_dimension() {
	    adensgrad[] += sq((totdens[1] - totdens[-1])/(2.*Delta));
      }
      adensgrad[] = sqrt(adensgrad[]);
      denslap[] = -fabs(totdens[1] + totdens[0,1] + totdens[-1] + totdens[0,-1] - 4.*totdens[])/sq(Delta);
    }
  
  stats sadensg = statsf(adensgrad);
  foreach()
    {
      double kscale = 40.;
      schl[] = exp( -kscale * adensgrad[] / sadensg.max);
    }
  
  static FILE * fp = fopen("p.ppm", "w");
  static FILE * fuf = fopen("f.ppm", "w");
  static FILE * fr1 = fopen("frho1.ppm", "w");
  static FILE * fr2 = fopen("frho2.ppm", "w");
  static FILE * fqx = fopen("fqx.ppm", "w");
  static FILE * fe1 = fopen("fe1.ppm", "w");
  static FILE * fe2 = fopen("fe2.ppm", "w"); 
  static FILE * fdl = fopen("fdl.ppm", "w");
  static FILE * fsc = fopen("fsc.ppm", "w");
  static FILE * fom = fopen("fom.ppm", "w");
  output_ppm (p, fp, n=512*11, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (f, fuf, n=512*11, min=0., max=1., box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (frho1, fr1, n=512*11, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (frho2, fr2, n=512*11, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (q.x, fqx, n=512*11, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (fE1, fe1, n=512*11, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (fE2, fe2, n=512*11, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (denslap, fdl, n=512*11, map=gray, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (schl, fsc, n=512*11, map=gray, box = {{x0,-L0/11},{x0+L0,0}});
  output_ppm (vort, fom, n=512*11, box = {{x0,-L0/11},{x0+L0,0}});
}

event end (t==tend) {
}

#if TREE
event adapt (i++) {
  scalar ffp[];
  foreach()
    ffp[] = 2.*f[] - 1.;
  adapt_wavelet ({ffp, frho1, p}, (double[]){0.01, 0.01, 0.01}, LEVEL, 5);
}

#endif //TREE