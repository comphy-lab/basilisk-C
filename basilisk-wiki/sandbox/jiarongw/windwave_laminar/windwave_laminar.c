/** 
This is the wave wind interaction simulation with linear wind profile and adaptive grid. Pressure is output on the run by `output_matrix_mpi` function and also dumped into dump file. For velocity initialization, velocity at the top is initialized with slope*height and kept the same across x position. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "navier-stokes/perfs.h"

/**
Functions to output pressure matrix on the fly (from Antoon). */

void output_matrix_mpi (struct OutputMatrix p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));

  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
        field[i][j] = interpolate (p.f, xp, yp);
      }
      else {
        Point point = locate (xp, yp);
        field[i][j] = point.level >= 0 ? val(p.f) : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif


    fwrite (&fn, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < p.n; i++){
      float xp = Delta*i + X0 + Delta/2.;
      fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif

  matrix_free (field);
}


/**
   We log some profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"

/**
   The primary parameters are the wave steepness $ak$, the Bond and
   Reynolds numbers. */

double ak = 0.05;
double BO = 200.;
double RE = 40000.;

/**
   The default maximum level of refinement depends on the dimension. */

int LEVEL = dimension == 2 ? 10 : 6;

/**
   The error on the components of the velocity field used for adaptive
   refinement. */

double uemax = 0.01;
double femax = 0.00001;

/**
   The density and viscosity ratios are those of air and water. */

double RATIO = (1./850.);
double MURATIO = (17.4e-6/8.9e-4);

/**
   Define if we want to use a Dirac viscous layer initialization. */
int DIRAC = 0;

/**
   Define the velocity output file number. */
int countvelocity = 0;

/**
   The wave number, fluid depth and acceleration of gravity are set to
   these values. */
double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 1.;
double T0_ = 1;

/** Define the wind profile related parameters. */

double m = 5.;  // vary between 5 and 8
double B = 0.;
double Karman = 0.41;   // Karman universal turbulence constant
double UstarRATIO = 1;   // Ratio between Ustar and c
double Ustar;
double Utop;
double y_1 = 0.;
double Udrift = 0.;   

/**
   The program takes optional arguments which are: the level of
   refinement **LEVEL**, steepness **ak**, Bond number **BO**, Reynolds
   numbers **RE**, the inverse wave age ($u_* / c$) **UstarRATIO**. The **m**, 
   **B** and **DIRAC** parameters are here but were not used.  */

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    RE = atof(argv[4]);   
  if (argc > 5)
    m = atof(argv[5]);
  if (argc > 6)
    B = atof(argv[6]);
  if (argc > 7)
    UstarRATIO = atof(argv[7]);
  if (argc > 8)
    DIRAC = atoi(argv[8]);

  /**
     Here we set the densities and viscosities corresponding to the
     parameters above. Note that these variables are defined in two-phase.h already.*/
 
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; // using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;  
  T0_ = 2.*pi/sqrt(g_*k_ + f.sigma*k_*k_*k_/rho1);
  
  /**
     The domain is a cubic box centered on the origin and of length
     $L_0=1$, periodic in the x- and z-directions. The lin-log profile with m 
     and B parameters was not used. */
  
  Ustar = sqrt(g_/k_+f.sigma*k_)*UstarRATIO;
  Utop = sq(Ustar)/(mu2/rho2)*L0/2.;
  
  // y_1 = m*mu2/rho2/Ustar; 
  // Udrift = B*Ustar;  
  
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
  u.n[top] = dirichlet(0);
  u.t[top] = dirichlet(Utop);
#if dimension > 2
  periodic (front);
#endif

  /**
     When we use adaptive refinement, we start with a coarse mesh which
     will be refined as required when initialising the wave. */
  
#if TREE  
  N = 128; // If coarsened or not
#else
  N = 1 << LEVEL;
#endif
  run();
}

/**
## Initial conditions

   These functions return the shape of a third-order Stokes wave with the
   wavenumber and steepness given by the parameters above ($ak$ and
   $k_$). */

double wave (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 - y;
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
  return eta1 + eta2 + eta3;
}
/**

   We also calculate an approximation to a Dirac distribution on the wave surface.
   This allows us to calculate a vortex sheet on the surface to provide a boundary
   layer in the air above the water surface. (NOTICE: this is not actually used
   because we already have prescribed air velocity.) */

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
event init (i = 0)
{
  fprintf(stderr, "UstarRATIO=%g B=%g\n m=%g ak=%g", UstarRATIO, B, m, ak);

  if (!restore ("restart")) {
    do {
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
      	  (9.*sq(alpa) - 13.)*cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
      	Phi[] = phi1 + phi2 + phi3;
      } 
      boundary ({Phi});
      if (DIRAC){
	/** 
	    We calculate the vorticity in the Dirac layer. We need a separate
	    foreach here because we need the derivative of the potential phi.*/
	scalar vort2[];
	scalar psi[];
	foreach() {
	  vort2[] = -2.0*gaus(y,wave(x,y)+y,f[])*(Phi[1,0]-Phi[-1,0])/(2.*Delta);
	  psi[] = 0.0;
	}  
     
	boundary({vort2,psi});
	psi[top] = dirichlet(0.);
	psi[bottom] = dirichlet(0.);
  
	/**
	   Solve the Poisson problem for the streamfunction psi given the 
           vorticity field.*/
        
	poisson(psi, vort2);
        
	/**
	   And then define the velocity field using centered-differencing
	   of the streamfunction. */

	foreach(){
	  u.x[] = (psi[0,1] - psi[0,-1])/(2.*Delta);
	  u.y[] = -(psi[1] - psi[-1])/(2.*Delta);
	}
      }
      else{
	/**
	   If we choose not to use the Dirac layer, instead initialize
	   in the water only according to the potential already calculated.*/
	foreach(){
	  foreach_dimension()
	    u.x[] = (Phi[1] - Phi[-1])/(2.0*Delta) * f[]; 
	}
	/**
	   Superpose the air velocity on top (linear profile). */
	foreach(){
	  u.x[] += Udrift + Utop / (L0/2.-eta(x,y)) * (y-eta(x,y)) * (1-f[]);
	}
      }
      boundary ((scalar *){u});    
    }

    /**
       On trees, we repeat this initialisation until mesh adaptation does
       not refine the mesh anymore. */

#if TREE  
    while (adapt_wavelet ({f, u},
			  (double[]){femax,100.*uemax,100.*uemax,100.*uemax}, LEVEL, 5).nf); // if not adapting anymore, return zero
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
   We log the evolution of the kinetic and potential energies and
   dissipation rate as functions of the non-dimensional time. */

event graphs (i++) {
  static FILE * fpwater = fopen("budgetWaterwind.dat", "a");
  static FILE * fpair = fopen("budgetAirwind.dat", "a");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv();
    keAir += rho[]*norm2*(1.0-f[])*dv();
    gpe += rho1*g_*y*f[]*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();
  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
  }
  fprintf (fpwater, "%g %g %g %g\n",
	   t, ke/2., gpe + 0.125*L0*L0*L0, dissWater);
  fprintf (fpair, "%g %g %g %g\n",
	   t, keAir/2., gpeAir + 0.125*L0*L0*L0, dissAir);
  fprintf (ferr, "%g %g %g %g\n",
	   t, ke/2., gpe + 0.125*L0*L0*L0, dissWater);
}

/**
## Visualisation

   We use Basilisk view (and output_ppm()) to display animations of the
   results.

   On some parallel systems, pipes tend to cause problems, so we switch
   to simple uncompressed PPM outputs when running with MPI. Otherwise,
   we use MPEG-4 file compression (which requires working pipes and
   ffmpeg). */

#define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t += 0.1) {

  /**
     We do simple movies of the volume fraction, level of
     refinement fields. */

  scalar omega[];
  vorticity (u, omega);
  view (width = 800, height = 600, fov = 18.8);
  clear();

  /**
     We repeat the drawing periodically in the x-direction. */
  
  for (double x = -L0; x <= L0; x += L0)
    translate (x) {
      draw_vof ("f");
      squares ("omega", linear = true, max=50, min=-50);
    }
  {
    static FILE * fp = POPEN ("omega", "a");
    save (fp = fp);
  }
}

/**
## End 

   The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
   (alternatively 4) periods. */

event end (t = 10.*T0_) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}


/** 
   Dump every 1/32 period.  */

event dumpstep (t += T0_/32.) {
  char dname[100];
  sprintf (dname, "dump%g", t/T0_);
  p.nodump = false;
  /* scalar p_air[],p_air_round[]; */
  /* foreach(){ */
  /*   p_air[] = p[]*(1-f[]); */
  /*   if (f[]<0.5){ */
  /*     p_air_round[] = p[]; */
  /*   } */
  /*   else{ */
  /*     p_air_round[] = 0; */
  /*   } */
  /* } */
  dump (dname);
  // Output pressure
  /* char pressurename1[100], pressurename2[100], pressurename3[100]; */
  /* sprintf (pressurename1, "./pressure/p_matrix%g.dat", t/(k_/sqrt(g_*k_))); */
  /* sprintf (pressurename2, "./pressure/p_air_matrix%g.dat", t/(k_/sqrt(g_*k_))); */
  /* sprintf (pressurename3, "./pressure/p_air_round_matrix%g.dat", t/(k_/sqrt(g_*k_))); */
  /* FILE * fpressure1 = fopen (pressurename1, "w");   */
  /* output_matrix_mpi (p, fpressure1, n = 512); */
  /* FILE * fpressure2 = fopen (pressurename2, "w");   */
  /* output_matrix_mpi (p_air, fpressure2, n = 512); */
  /* FILE * fpressure3 = fopen (pressurename3, "w");   */
  /* output_matrix_mpi (p_air_round, fpressure3, n = 512);   */
}

/**
## Mesh adaptation

   On trees, we adapt the mesh according to the error on volume fraction
   and velocity. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
}
#endif


/**
## Running in parallel

   This file will work in 2D or 3D, either with parallel multigrid
   (without adaptivity), using for example:

   ~~~bash
   qcc -source -D_MPI=1 -grid=multigrid3D wave.c
   scp _wave.c occigen.cines.fr:
   ~~~

   and then following a recipe similar to that of the
   [atomisation](/src/examples/atomisation.c#on-occigen) example.

   To use adaptivity, just do something like:

   ~~~bash
   qcc -source -D_MPI=1 -grid=octree wave.c
   scp _wave.c occigen.cines.fr:
   ~~~
*/
