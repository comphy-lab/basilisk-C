/**
   This is a standard file to reproduce the simulations presented in: [Mostert and Deike, Inertial energy dissipation in shallow-water breaking waves. Journal of Fluid Mechanics, 890:A12, 2020](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/inertial-energy-dissipation-in-shallowwater-breaking-waves/8B818C95BBB1C8BF2723DB611AFB9DCD).

   We use Navier-Stokes in 2D with surface tension, and with the momentum conserving VOF scheme. */

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"
#include "tag.h"

/**
   Include profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"

/**
   Set maximum and minimum refinement levels. The "minimum" level is
   only really used in initialization: the adaptation during runtime
   uses level 9 as the minimum. */

#define LEVEL 14
#define MINLEVEL 5

/**
   Switch on movie creation if desired. */
#define MOVIES 0
/**
   Switch on beach definition (required) - when IB methods come, we
   will work on this. */
#define BEACH 1
/**
   Define useful physical parameters. */
#define ALPA 3.0*pi/180. //Beach slope
#define BO 1000.0 //Bond number
#define RE 40000.0 //Reynolds number
#define RATIO 1.0/850.0 //density ratio, water to air
#define MURATIO 17.4e-6/8.9e-4 //dynamic viscosity ratio, water to air
#define AMPLITUDE 0.3 // Amplitude
#define MAXTIME 70.0 // Maximum runtime.
/**
   First, we set the scales.
   We use the initial depth h_ as the length scale,
   and gravity g_ as the velocity scale. This should
   completely characterize the system. */
double h_ = 1.0;
double g_ = 1.0;
double a_ = AMPLITUDE;
/**
   Keep the beach off the bottom of the domain. Initialize the centre of the wave at distance xWave_ from the left boundary. The beach begins to slope at distance xChange_. The size of the domain is xextent_. */
double yOffset_ = 0.5;
double xWave_ = 5.0;
double xextent_ = 40.0;
double xChange_ = 10.0;
//double xextent_ = 30.0;
double patm_ = 1.0;
scalar beach[];
scalar dissrate[];
/**
   Define a secondary tracer which is the fluid kept within the beach. */
scalar fbb[];
scalar * tracers = {fbb};
//----------------------MAIN-----------------------//
//Begin main function
int main()
{
  size (xextent_);
  origin (0.0, -h_-yOffset_);
  rho1 = 1.0;
  rho2 = RATIO;
  /**
     For calculating viscosities, interpret Reynolds number as defined at depth of unity. */
  mu1 = 1.0/RE;
  mu2 = 1.0/RE*MURATIO;
  /**
     Use depth length scale for Bond number as well. */
  f.sigma = 1.0/BO;
  init_grid(1 << (LEVEL-5));
  /**
     Acceleration using reduced gravity. */
  G.y = -g_;
  run();
}
//----------------------HELPER FUNCTIONS---------------------//
/** 
    We use the moving cylinder example to inform the implementation
    of the beach.*/
#ifdef BEACH
double beachProfile(double x, double y) {
//Define beach profile H(x)
  double H = 0.0;
  if (x < xChange_)
    H = -h_;
  else
    H = -h_ + ALPA*(x-xChange_);
  return H;
}

/**
  Set so that the "beach" variable is 1.0 for beachProfile<0,
  and 0.0 for beachProfile > 1. */
event set_beach(i=0;i++) {
  fraction (beach, beachProfile(x,y)-y);
  /** 
      And set the corresponding momentum B.C.s within the beach region:*/
  foreach(){
    foreach_dimension()
      u.x[] = (1.0 - beach[])*u.x[];
  }
  boundary ((scalar *){u});
}

/**
   At a time sufficiently close to zero, update the beach to include water. */
event updatef(i=20) {
  foreach(){
    fbb[] = beach[];
    f[] = f[]+ fbb[];
  }
}
#endif //BEACH
/**
   Necessary user-defined functions to get this show off the road!*/
double sech(double qval) {
  return 1.0/cosh(qval);
}

double maxv (scalar a) {
  double maxi = - 1.0e100;
  /**
     Calculate the maximum quantity in the water. */
  foreach (reduction(max:maxi)){
    if (fabs(a[]*f[]) > maxi)
      maxi = fabs(a[]*f[]);
  }
  return maxi;
}

//---------------------INITIALIZATION------------------------//
double waveGN( double x, double y )
{
  // Try "analytical solution of G-N equations":
  double k = sqrt(3.*a_)/(2.*h_*sqrt(h_*(1. + a_/h_)));
  return a_*sq(sech(k*x)) - y;
}

double detax( double x )
{
  double tmp1 = -sqrt(3.0)*a_*sqrt(a_)/sqrt(h_*(1.0 + a_/h_));
  double tmp2 = sqrt(3.0*a_)/(2.0*h_*sqrt(h_*(1.0 + a_/h_)));
  double tmp3 = sq(sech(tmp2*x));
  double tmp4 = tanh(tmp2*x);
  double deta = tmp1 * tmp3 * tmp4;
  return deta;
}

/*Write initialization event*/
event init (i=0)
{
  if (!restore("restart")){
    double femax = 2e-4;
    double uemax = 2e-2;
    double bemax = 8e-3;
    int jdx = 0;
    do {
      jdx += 1;
      fprintf(ferr, "Initializing, %d\n", jdx );
      fraction ( f, waveGN(x - xWave_,y) ); //Define interface
      /* //Define the wave speed and initialize all fluid to move with this speed */
      double cspeed = sqrt(g_*h_*(1.0 + a_/h_));
      foreach() { //Green-Naghdi analytical soliton solution
    	double eta = waveGN (x-xWave_,y) + y;
    	double deta = detax (x-xWave_);
    	u.x[] = cspeed*eta / (h_ + eta) * f[];
    	u.y[] = -(y + h_)*cspeed*deta/(eta + h_) * (1.0 - eta/(h_+eta) ) * f[];
	f[] = f[] + beach[];
      }
      boundary ((scalar *){u});
    }
    while (adapt_wavelet ({f,u,beach}, (double[]){femax,uemax, uemax, uemax, bemax}, LEVEL, MINLEVEL).nf);
  }
}

//-------------------ADAPTIVITY---------------------//
/*Adapt once on error in volume fraction, velocity field, and beach fraction*/
event adapt(i++) {
  //double uemax = 1e-5;
  double femax = 1e-8;
  double uemax = 2e-3;
  double bemax = 1e-2;
  adapt_wavelet ({f,u,beach}, (double[]){femax,uemax,uemax,uemax,bemax}, LEVEL, 9);
}

//------------------DIAGNOSTICS---------------------//
/*Define functions for determining kinetic and potential energy*/
int dissipation_rate (vector u, double* rates, scalar drate)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1] - u.x[-1])/(2.0*Delta);
    double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]   - u.y[-1]  )/(2.*Delta);
    double dvdy = (u.y[0,1] - u.y[0,-1])/(2.0*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.0*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.0*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			      sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			      sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz));
    double localRate = mu1/rho[]*f[]*sqterm*(1.0-beach[]);
    drate[] = localRate;
    rateWater += localRate; //water
    rateAir += mu2/rho[]*(1.0-f[])*sqterm*(1.0-beach[]); //air
  }
  //fprintf (fix, "\n");
  //fprintf (fiy, "\n");
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

/**
We log the evolution of the kinetic and potential energies and
dissipation rate as functions of the non-dimensional time. */

event graphs (i++) {
  scalar umag[];
  static FILE * fpwater = fopen("budgetWater.dat", "w");
  static FILE * fpair = fopen("budgetAir.dat", "w");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv()*(1.0-beach[]);
    keAir += rho[]*norm2*(1.0-f[])*dv()*(1.0-beach[]);
    gpe += rho[]*g_*y*f[]*dv()*(1.0-beach[]);
    gpeAir += rho[]*g_*y*(1.0-f[])*dv()*(1.0-beach[]);
    umag[] = sqrt(sq(u.x[]) + sq(u.y[]));
  }
  double rates[2];
  dissipation_rate(u, rates, dissrate);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
  }
  double maxs = maxv(umag);
  fprintf (fpwater, "%g %g %g %g %g\n",
	   t, ke/2., gpe, dissWater, maxs);
  fprintf (fpair, "%g %g %g %g %g\n",
	   t, keAir/2., gpeAir, dissAir, maxs);
  fprintf (ferr, "%g %g %g %g %g\n",
	   t, ke/2., gpe, dissWater, maxs);
}


//------------------------IMAGING--------------------------
#define POPEN(name, mode) fopen (name ".ppm", mode)

#ifdef MOVIES
event movie (t += 0.05) {
  double c = sqrt(g_*h_*(1.0+a_/h_));
  double xbeach = pow(ALPA, -1);
  double xbreak = xChange_ + xbeach;
  double prefac = 5.0;
  scalar m[];
  scalar pid[];
  foreach(){
#ifdef BEACH
    m[] = 0.5 - beach[];
#else
    m[] = 0.0;
#endif
    pid[] = pid();
  }
  static FILE * fp = POPEN ("f", "w");
  static FILE * fpg = POPEN ("fg", "w");
  static FILE * fpid = POPEN ("fpid", "w");
  static FILE * fpb = POPEN ("fpb", "w");
  output_ppm(f,fpg, mask=m, min=0.0,max=1.0, n=1024, 
             box= {{c*t, -1.0 }, {c*t + 2.0*xWave_,2.0}});
  output_ppm (f, fp, mask=m, min =0.0, max=1.0, n=1024, box= {{0.0,-1.0}, {L0,2.0}});
  output_ppm (f, fpb, mask=m, min=0.0, max=1.0, n=2048, box= {{xbreak-prefac*xWave_,-1.0}, {xbreak+prefac*xWave_,2.0}});
  output_ppm (pid, fpid, mask=m, n=1024, box= {{0.0,-1.0}, {L0,2.0}});
  scalar omega[];
  scalar l[];
  vorticity (u, omega);
  static FILE *fp1 = POPEN ("omega", "w");
  output_ppm (omega, fp1, n=1024, box= {{0.0,-1.0}, {L0,2.0}});
  static FILE *fp1b = POPEN ("omegab", "w");
  output_ppm (omega, fp1b, mask=m, n=2048, box= {{xbreak-prefac*xWave_,-1.0}, {xbreak+prefac*xWave_,2.0}});
  static FILE * fp2 = POPEN ("ux", "w");
  output_ppm (u.x, fp2, mask=m, n=1024, box= {{0.0,-1.0}, {L0,2.0}});
  static FILE * fp2b = POPEN ("uxb", "w");
  output_ppm (u.x, fp2b, mask=m, n=2048, box= {{xbreak-prefac*xWave_,-1.0}, {xbreak+prefac*xWave_,2.0}});
  static FILE * fp3 = POPEN ("uy", "w");
  output_ppm (u.y, fp3, mask=m, n=1024, box= {{0.0,-1.0}, {L0,2.0}});
  static FILE * fp3b = POPEN ("uyb", "w");
  output_ppm (u.y, fp3b, mask=m, n=2048, box= {{xbreak-prefac*xWave_,-1.0}, {xbreak+prefac*xWave_,2.0}});
#if TREE
  foreach()
    l[] = level;
  static FILE * fp4 = POPEN ("l", "w");
  output_ppm (l, fp4, box= {{0.0,-1.0}, {L0,2.0}}, min=1, max=LEVEL);
  static FILE * fp4c = POPEN ("lg", "w");
  output_ppm(l,fp4c,  n=1024, 
             box= {{c*t, -1.0 }, {c*t + 2.0*xWave_,2.0}}, min=1, max=LEVEL);
#endif
  static FILE * fp5 = POPEN ("p", "w");
  output_ppm (p, fp5, mask=m, n=1024, box= {{0.0,-1.0}, {L0,2.0}});
}
#endif //MOVIES

event dumpstep(t += 0.05) {
  char dname[100];
  sprintf(dname,"dump%g.ppm",t);  
  dump(dname);
}

event outputInterface(t += 0.05) {
  char resultname[100];
  sprintf( resultname, "results%4.2f_%d.dat", t, pid() );
  FILE * fp = fopen(resultname, "w");
  scalar xpos[];
  scalar ypos[];
  position (f, xpos, {1, 0});
  position (f, ypos, {0, 1});
  foreach()
    {
      if (xpos[] != nodata){
	fprintf (fp, "%g %g %g %g\n", x, y, xpos[], ypos[]);
      }
    }
  fclose (fp);
}

event end (t=MAXTIME) { 
  printf ("i = %d t=%g\n", i, t);
  dump ("end");
}
