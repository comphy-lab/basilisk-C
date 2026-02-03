/**
This simulation studies the effect of surfactants on gravity-capillary wave based on breaking wave example [wave.c](/sandbox/popinet/wave.c) and Marangoni flow package by Palas Kumar Farsoiya et al. [farsoiya/marangoni_surfactant/](/sandbox/farsoiya/marangoni_surfactant/).
*/

/*
Requirement: redistance2.h, surfactant-transport.h
*/

/**
# Breaking wave

We solve the two-phase Navier--Stokes equations with surface tension
and using a momentum-conserving transport of each phase. Gravity is
taken into account using the "reduced gravity approach" and the
results are visualised using Basilisk view. */



#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "curvature.h"
#include "view.h"
#include "tag.h"

//Surfactants
#include "surfactant-transport.h"

/**
We log some profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"
#include <sys/stat.h>

//Surfactants
double beta = 1.0;
double dsigma0 = 0.5;
double gamma0 = 1., gamma_inf = 1.;
scalar sigmaf[];
vector iforce[];

/**
The primary parameters are the wave steepness $ak$, the Bond and
Reynolds numbers. */

double ak = 0.3;
double BO = 10.;
double RE = 40000.;
double TEND = 1.0;    // end of the simulation

/**
The default maximum level of refinement depends on the dimension. */

int LEVEL = dimension == 2 ? 8 : 6;

/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.005;

/**
The density and viscosity ratios are those of air and water. */

#define RATIO (1./850.)
#define MURATIO (17.4e-6/8.9e-4)

/**
Define if we want to use a Dirac viscous layer initialization. */

int DIRAC = 0;

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. */

#define k_  (2.*pi)
#define h_   0.5
#define g_   1.
double MA = 1;

/**
The program takes optional arguments which are the level of
refinement, steepness, Bond and Reynolds numbers, and optional Dirac
initialisation. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    beta = atof (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    RE = atof(argv[4]);    
  if (argc > 5)
    DIRAC = atof(argv[5]);
  if (argc > 6)
    LEVEL = atof(argv[6]);
  if (argc > 7)
    TEND = atof(argv[7]);
  if (argc > 8)
    dsigma0 = atof(argv[8]);
  MA = beta / (BO * k_ * sqrt(k_)/RE);
  /**
  The domain is a cubic box centered on the origin and of length
  $L0=1$, periodic in the x- and z-directions. */
   
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
#if dimension > 2
  periodic (front);
#endif

  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;

  TOLERANCE = 1e-4;
  D_s = 0.01; //Surface Diffusivity of surfactant D_s
  surfactant = 1;

  d.sigmaf = sigmaf;
  tracers = list_concat(tracers, {c1, gamma2, d2, pfield});
  reinit_skip_steps = 100; //this needs consideration if there are topological changes
  /**
  When we use adaptive refinement, we start with a coarse mesh which
  will be refined as required when initialising the wave. */
  
#if TREE  
  N = 1 << LEVEL;
#else
  N = 1 << LEVEL;
#endif


  struct stat st = {0};
  char name[80];
  sprintf (name, "dat");
  char newname[80];
  sprintf (name, "dat");

  if (stat(name, &st) == -1)
    {
      mkdir(name, 0700);
    }

  sprintf (name, "eta");
  if (stat(name, &st) == -1)
    {
      mkdir(name, 0700);
    }

    /**
      Preserve old stats restoring the simulation with a timestamp suffixed. */		
  sprintf (name, "stats.dat");
	
  if (stat(name, &st) == 0)
    {
      sprintf (newname, "stats-%ld.dat",time(0));
	rename(name, newname);
    }
  sprintf (name, "perfs");
	
  if (stat(name, &st) == 0)
    {
      sprintf (newname, "perfs-%ld",time(0));
	rename(name, newname);
    }

  run();
}

//Surfactants (Gravity)
event acceleration (i++) {

  face vector av = a;
  
  foreach_face(y)
    av.y[] += (-1. + alphav.y[]*rho1);
  
  boundary ((scalar *){av});
}

//Surfactants (isotherm)
event stability (i++){

  foreach(){
       
      double deltas = (pfield[]*(1. - pfield[]))/EPSILON;
      if (deltas > 1.e-1 && surfactant == 1){
	  double gamma = c1[]*4.*EPSILON;
	  sigmaf[] =  1.0/(BO*sq(k_))*(1.   - dsigma0*tanh(MA/dsigma0 * BO * k_ * sqrt(k_)/RE * gamma/gamma0));
      }
      else
        sigmaf[] = 1.0/(BO*sq(k_));
      
  }
  boundary({sigmaf});

 /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;
    }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  The maximum timestep is set using the sum of surface tension
  coefficients. */

  for (scalar c in interfaces)
  
    foreach(serial){
      if (sigmaf[] > 0) {
        double dt = sqrt (rhom*cube(dmin)/(pi*sigmaf[]))/10; //some prefactors maybe needed for stability
        if (dt < dtmax)
          dtmax = dt;
      }
   }

}

event test (i++){
  if (i == 0)
  {
      foreach ()
      {
          double deltas = (pfield[] * (1. - pfield[])) / EPSILON;
          // initially uniform surfactant concentration
          if (deltas > 1.e-1 && surfactant == 1){
            c1[] = gamma0 * deltas;
            double gamma = c1[]*4.*EPSILON;
	    sigmaf[] =  1.0/(BO*sq(k_))*(1.   - dsigma0*tanh(MA/dsigma0 * BO * k_ * sqrt(k_)/RE * gamma/gamma0));
          }
          else
            sigmaf[] = 1./(BO*sq(k_));
      }
  }
}

/**
## Initial conditions

These functions return the shape of a third-order Stokes wave with the
wavenumber and steepness given by the parameters above ($ak$ and
$_k_$). */

double wave (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

double eta (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}

/**
We also calculate an approximation to a Dirac distribution on the wave
surface.  This allows us to calculate a vortex sheet on the surface to
provide a boundary layer in the air above the water surface. */

double gaus (double y, double yc, double T){
  double deltaw = sqrt(2.0/RE)/k_;
  double deltaa = sqrt(2.0/RE*MURATIO/RATIO)/k_;
  double r = y - yc;
  return 2.0/(sqrt(2.0*pi*sq(deltaa)) + sqrt(2.0*pi*sq(deltaw))) *
    (T*exp(-sq(r)/(2.0*sq(deltaw))) + (1.0 - T)*exp(-sq(r)/(2.0*sq(deltaa))));
}

/**
We either restart (if a "restart" file exists), or initialise the wave
using the third-order Stokes wave solution. */

# define POPEN(name, mode) fopen (name ".ppm", mode)

event init (i = 0)
{
  if (!restore ("restart")) {
    do {

      foreach(){
        d[] = wave(x,y);
      }
      fraction (f, wave(x,y));

      /**
      To initialise the velocity field, we first define the potential. */
      
      scalar Phi[];
      foreach() {
	double alpa = 1./tanh(k_*h_);
	double a_ = ak/k_;
	double sgma = sqrt(g_*k_*tanh(k_*h_)*
			   (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					      (sq(alpa) - 1.) + sq(alpa))));
	double A_ = a_*g_/sgma;
	double phi1 = A_*cosh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x);
	double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
	  cosh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_);
	double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
	  (9.*sq(alpa) - 13.)*
	  cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
	Phi[] = phi1 + ak*phi2 + ak*ak*phi3;
      }
      boundary ({Phi});
      
      if (DIRAC) {

        /**
        We calculate the vorticity in the Dirac layer. We need a separate
        foreach here because we need the derivative of the potential phi.*/

        scalar vort2[];
        scalar psi[];
        foreach() {
          vort2[] = -2.0*gaus(y,wave(x,y)+y,f[])*(Phi[1,0]-Phi[-1,0])/(2.*Delta);
          psi[] = 0.0;
        }
        boundary ({vort2,psi});
        psi[top] = dirichlet(0.);
        psi[bottom] = dirichlet(0.);
        
        /**
        Solve the Poisson problem for the streamfunction psi given the
        vorticity field.*/
        
        poisson (psi, vort2);
        
        /**
        And then define the velocity field using centered-differencing
        of the streamfunction. */
            
        foreach() {
          u.x[] = (psi[0,1] - psi[0,-1])/(2.*Delta);
          u.y[] = -(psi[1] - psi[-1])/(2.*Delta);
        }
      }
      else {

        /**
        If we choose not to use the Dirac layer, instead initialize
        in the water only according to the potential already calculated.*/
        
        foreach()
          foreach_dimension()
            u.x[] = (Phi[1] - Phi[-1])/(2.0*Delta) * f[];
            //u.x[] = f[]*10;
      }
      boundary ((scalar *){u});

  //Surfactants
      event("properties2"); //call this so that pfield can be populated and then c1 can be initialized	
    foreach(){
      double deltas = (pfield[]*(1. - pfield[]))/EPSILON; //

    // double gamma0 = gamma_inf*0.01;
      c1[] = gamma0*deltas; //convert volume to surface
      if (deltas > 1.e-2){
	  double gamma = c1[]*4.*EPSILON;
	  sigmaf[] =  1.0/(BO*sq(k_))*(1.   - dsigma0*tanh(MA/dsigma0 * BO * k_ * sqrt(k_)/RE * gamma/gamma0));
      }
      else
        sigmaf[] =  1./(BO*sq(k_));
      
    }
    boundary({sigmaf});

    }

    /**
    On trees, we repeat this initialisation until mesh adaptation does
    not refine the mesh anymore. */

#if TREE  
    while (adapt_wavelet ({f,u,pfield,c1},
			  (double[]){0.01,uemax,uemax,uemax,0.001,0.001}, LEVEL, 5).nf);
#else
    while (0);
#endif
  }
}

/**
## Outputs

We are interested in the viscous dissipation rate in both water and air. */

int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    rateWater += mu1/rho[]*f[]*sqterm; //water
    rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

/**
We also want to count the drops and bubbles in the flow. */

event countDropsBubble(i++)
{
  scalar m1[]; //droplets
  scalar m2[]; //bubbles
  foreach(){
    m1[] = f[] > 0.5; //i.e. set m true if f[] is close to unity (droplets)
    m2[] = f[] < 0.5; //m true if f[] close to zero (bubbles)
  }
  int n1 = tag(m1);
  int n2 = tag(m2);
  
  /**
  Having counted the bubbles, now we find their size. This example
  is similar to the jet atomization problem. We are interested in
  the volumes and positions of each droplet/bubble.*/

  double v1[n1]; //droplet
  coord b1[n1];  //droplet
  double v2[n2]; //bubble
  coord b2[n2];  //bubble
  
  /**
  We initialize: */
  
  for (int j=0; j<n1; j++)
      v1[j] = b1[j].x = b1[j].y = b1[j].z = 0.0;
  for (int j=0; j<n2; j++)
      v2[j] = b2[j].x = b2[j].y = b2[j].z = 0.0;
  
  /**
  We proceed with calculation. */

  foreach_leaf() {
    // droplets
    if (m1[] > 0) {
      int j = m1[] - 1;
      v1[j] += dv()*f[]; //increment the volume of the droplet
      coord p = {x,y,z};
      foreach_dimension()
	b1[j].x += dv()*f[]*p.x;
    }
    // bubbles
    if (m2[] > 0) {
      int j = m2[] - 1;
      v2[j] += dv()*(1.0-f[]);
      coord p = {x,y,z};
      foreach_dimension()
	b2[j].x += dv()*(1.0-f[])*p.x;
    }
  }
  
  /**
  Reduce for MPI. */
  
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b1, 3*n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, v2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /**
  Output the volume and position of each droplet to file. */
  
  static FILE * fdrop = fopen("droplets.dat","w");
  static FILE * fbubb = fopen("bubbles.dat","w");
  for (int j=0; j<n1; j++)
    fprintf (fdrop, "%d %g %d %g %g %g\n", i, t,
	     j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j]);
  for (int j=0; j<n2; j++)
    fprintf (fbubb, "%d %g %d %g %g %g\n", i, t,
	     j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j]);
}

/**
We log the evolution of the kinetic and potential energies and
dissipation rate as functions of the non-dimensional time. */

event graphs (i++) {
  static FILE * fpwater = fopen("budgetWater.dat", "w");
  static FILE * fpair = fopen("budgetAir.dat", "w");
  double ke = 0., gpe = 0.;
  double spe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) reduction(+:spe)
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv();
    keAir += rho[]*norm2*(1.0-f[])*dv();
    gpe += rho1*g_*y*f[]*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();
    if (interfacial (point, f)) {
      coord n      = interface_normal(point, f), pp;
      double alpha = plane_alpha (f[], n);
      double ar    = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
      spe += ar*sigmaf[];
    }

  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
    if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation spe\n");
    fprintf (fpair, "t ke gpe dissipation spe\n");
    }
  fprintf (fpwater, "%g %g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater, spe);
  fprintf (fpair, "%g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), keAir/2., gpeAir + 0.125, dissAir);
  fprintf (ferr, "%g %g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater, spe);
}

/**
## Visualisation
*/

event movies (t += 0.01) {

  scalar omega[];
  vorticity (u, omega);
  foreach(){
      gamma2[] = c1[]*4.*EPSILON/gamma0;
  }
#if dimension == 2
  view (width = 800, height = 600, fov = 18.8);
  clear();

  /**
  We repeat the drawing periodically in the x-direction. */
  
  for (double x = -L0; x <= L0; x += L0)
    translate (x) {
      draw_vof ("f");
      squares ("omega", linear = true, map = cool_warm , min=-10.0, max=10.0);
    }

  {
    static FILE * fp = POPEN ("vor", "a");
    save (fp = fp);
    fflush (fp);
  }

  /**
  This gives the following movie.
  ![](movie.mp4)
  */
  
#else // dimension == 3

  view (width = 1600, height = 1200, theta = pi/4, phi = -pi/6, fov = 20);
  clear();
  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
      for (double z = -3*L0; z <= L0; z += L0)
	translate (z = z)
	  draw_vof ("f");
    }
  save ("below.mp4");

  view (width = 1600, height = 1200, theta = pi/4, phi = pi/6, fov = 20);
  clear();

  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
      for (double z = -3*L0; z <= L0; z += L0)
	translate (z = z)
	  draw_vof ("f");
    }
#endif // dimension == 3
save ("movie.mp4");
}

/**
## Dump/restore

To be able to restart, we dump the entire simulation at regular
intervals. */
scalar stress2[];
scalar stx[];
scalar sty[];
event snapshot (t += 0.5) {
    foreach(){
      gamma2[] = c1[]*4.*EPSILON;
    }
  char filename[100];
  sprintf(filename, "dat/Ma%g-Bo%g-ak%g-%d-%g",MA,BO,ak,(1 << LEVEL),t);
  foreach_dimension () {
    iforce.x.nodump = false;
  }
  
  dump(filename);
  foreach_dimension () {
    iforce.x.nodump = true;
  }
}

/**
## End 

The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
(alternatively 4) periods. */

event end (t = TEND*k_/sqrt(g_*k_)) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

/**
## Mesh adaptation

On trees, we adapt the mesh according to the error on volume fraction
and velocity. */

#if TREE
event adapt (i++) {

 double femax = 0.001;
 double uxemax = 0.005;
 double uyemax = 0.005;
 double uzemax = 0.005;
 double c1emax = 0.001;
 double pfemax = 0.001;
 
 

  adapt_wavelet ({f,u,pfield,c1}, (double[]){femax,uxemax,uyemax,uzemax,pfemax,c1emax}, LEVEL); 
}
#endif

/**
## Functions to extract surface quantities
*/

static int compare (const void * a, const void * b)
{
  const double * A = (const double *) a;
  const double * B = (const double *) b;
  // sort by first column (x)
  if (A[0] < B[0]) return -1;
  if (A[0] > B[0]) return  1;
  return 0;
}

/** We want to compute some quantities at the interface. */
void output_int_qtn (char * fname, int istep,
                     scalar c, scalar sca,
                     double stp_eta)
{
  // We first loop over all the interfacial points
  // and we count them (per processor)
  int int_pt = 0;

  foreach (serial) {
    if (interfacial (point, c)) {
      if (point.level == LEVEL) {
        coord n = interface_normal (point, c), pp;
        double alpha1 = plane_alpha (c[], n);
        plane_area_center (n, alpha1, &pp);

        // we compute centroid coordinates (used later)
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        (void) xc; (void) yc; // silence warnings in case you later edit outputs

#if dimension > 2
        double zc = z + Delta*pp.z;
        (void) zc;
#endif

        // Stay on the current interfacial cell (no locate/Point reassignment)
        int_pt++;
      }
    }
  }

  int tot_column = 4;
  double t_mat[int_pt][tot_column];

  for (int j = 0; j < tot_column; j++)
    for (int i = 0; i < int_pt; i++)
      t_mat[i][j] = 0.;

  fprintf (stderr, "First pass over interfacial points\n");

  scalar pos[];
#if dimension > 2
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#else
  coord G = {0.,1.},   Z = {0.,0.};
#endif
  position (c, pos, G, Z);

  int count = 0;

  foreach (serial) {
    if (interfacial (point, c)) {
      if (point.level == LEVEL) {
        coord n = interface_normal (point, c), pp;
        double alpha1 = plane_alpha (c[], n);
        plane_area_center (n, alpha1, &pp);

        double eta = pos[];

        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        (void) yc; // yc not otherwise used unless you re-enable interpolation

#if dimension > 2
        double zc = z + Delta*pp.z;
#endif

        t_mat[count][0] = xc;
#if dimension > 2
        t_mat[count][1] = zc;
#else
        t_mat[count][1] = 0.;
#endif
        t_mat[count][2] = sca[];
        t_mat[count][3] = eta; // already defined at the interface
        count++;
      }
    }
  }

  fprintf (stderr, "Second pass over interfacial points\n");

  // We sort locally t_mat by the x coordinate (the first, i.e. 0-th, column)
  qsort (t_mat, int_pt, tot_column*sizeof(double), compare);
  fprintf (stderr, "First sort\n");

  double tot_row;

#if _MPI
  // On multiple cores, we gather all the int_pt to the root pid
  // and, then, we broadcast this information to all the processes
  int nproc;
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  int counts_it[nproc];

  if (pid() == 0) {
    MPI_Gather (&int_pt, 1, MPI_INT, counts_it, 1, MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    MPI_Gather (&int_pt, 1, MPI_INT, NULL,      1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  MPI_Bcast (counts_it, nproc, MPI_INT, 0, MPI_COMM_WORLD);

  // Each processor knows the int_pt of the others. So we can compute
  // the displacement, the total number of interfacial points and
  // the number of elements owned by each processor
  int tot_int_pt = 0;
  int tot_el_p[nproc], disp_r[nproc], disp[nproc];

  for (int i = 0; i < nproc; i++) {
    disp_r[i] = tot_int_pt;
    disp[i]   = disp_r[i]*tot_column;
    tot_int_pt += counts_it[i];
    tot_el_p[i] = counts_it[i]*tot_column;
  }

  tot_row = tot_int_pt;

  // --> Gather to the root pid
  // --> Sort by first column
  double t_mat_tot[tot_int_pt][tot_column];

  if (pid() == 0) {
    int tot_el_p0[nproc], disp0[nproc];
    for (int i = 0; i < nproc; i++) {
      tot_el_p0[i] = tot_el_p[i];
      disp0[i]     = disp[i];
    }

    MPI_Gatherv (&t_mat, int_pt*tot_column, MPI_DOUBLE,
                 t_mat_tot, tot_el_p0, disp0, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    fprintf (stderr, "GatherV 0\n");

    qsort (t_mat_tot, tot_int_pt, tot_column*sizeof(double), compare);
    fprintf (stderr, "Second sort 0\n");
  }
  else {
    int tot_el_p1[nproc], disp1[nproc];
    for (int i = 0; i < nproc; i++) {
      tot_el_p1[i] = tot_el_p[i];
      disp1[i]     = disp[i];
    }

    MPI_Gatherv (&t_mat, int_pt*tot_column, MPI_DOUBLE,
                 NULL, tot_el_p1, disp1, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
  }

#else
  tot_row = int_pt;

  double t_mat_tot[int_pt][tot_column];
  for (int j = 0; j < tot_column; j++)
    for (int i = 0; i < int_pt; i++)
      t_mat_tot[i][j] = t_mat[i][j];
#endif

  if (pid() == 0) {
    fflush (stderr);
    FILE * eta_loc = fopen (fname, "w");
    // Print in ASCII format
    for (int i = 0; i < tot_row; i++) {
      fprintf (eta_loc, "%8E %8E %8E %8E\n",
               t_mat_tot[i][0], t_mat_tot[i][1],
               t_mat_tot[i][2], t_mat_tot[i][3]);
    }
    fclose (eta_loc);
    fflush (eta_loc);
  }

  fprintf (stderr, "Print\n");
}

event eta_loc (t += 0.01)
{
  fflush (stderr);

  char etaname_in[100];
  sprintf (etaname_in, "./eta/eta_loc_t%09d.out", i);

  double stp_eta = 0.0*(L0/pow(2.0,LEVEL)); // it corresponds to 0*Delta

  boundary ({gamma2});
  output_int_qtn (etaname_in, i, f, gamma2, stp_eta);

  fflush (stderr);

  if (pid() == 0) {
    char name_1[80];
    sprintf (name_1,"./log_eta.out");
    FILE * log_sim = fopen (name_1,"a");
    fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
    fclose (log_sim);
  }
}