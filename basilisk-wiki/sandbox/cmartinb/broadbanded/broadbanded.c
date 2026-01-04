/** 
   This is the simulation file for the wind wave interaction problem. */

//#include "navier-stokes/centered.h"
#include "sandbox/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "view.h"
#include "tag.h"
#include "sandbox/profile6.h"   // From Antoon
#include "sandbox/maxruntime.h"   
#include "sandbox/my_funcs.h"   
#include "sandbox/colormap.h"   

/** 
   Input parameters of the simulations. */

double Re_ast = 720;                   // Friction Reynolds number
double BO = 200;                       // Bond number
double Re_wave = 2.55e4;               // Wave Reynolds number
double UstarRATIO = 0.5;               // Friction velocity over wave speed
double ak = 0.3;                       // initial steepness
double r_L0lam = 4.0;                  // ratio between the box size and one wavelength
double rho_r = 1.225/1000.0;           // reference density (air)
double mu_r  = 2.2471881948940954e-06; // reference dynamic viscosity (air)
int MAXLEVEL = 10;                     // max level of the simulation
int MINLEVEL = 6;                      // min level of the simulation
double RELEASETIME = 1000;             // physical time before which the precursor simulation should be run
double uemaxRATIO = 0.3;               // Ratio to control the maximum error on the air velocity
double end_sim = 1000;                 // end of the simulation (physical time unit)
int do_eta_loc = 1;                    // output or not eta_loc
int do_profile = 1;                    // output or not profiles
int do_fields  = 1;                    // output or not fields
int do_tagging = 1;                    // output or not tagging
int prt_res = 9;                       // printing resolution
int from_pr = 1;                       // from precursor simulation
int st_wave = 1;                       // initialize with a Stokes Wave (1) or import a prescribed distribution
int RANDOM = 2;                        // random number (employed in spectrum.h)

/**
   We define these values: the wave number, fluid depth, wave period, gravity acceleration,
   the wave speed, the friction velocity, the wavelength, alter_MU. */

double k_  = 1.0; // we change later 
double h_  = 1.0; // we change later 
double T0_ = 1.0; // we change later
double g_  = 1.0; // we change later
double c_  = 1.0; // we change later
double Ustar = 1.0; // we change later
double lam_ = 1.0; // we change later
double alter_MU = 1.0; // we change later

/**
   For some statistics. */

bool is_first = true;
int i_sp = 0;
double t_sp = 0.0;
#include "sandbox/perfs.h"  

/**
   Clean the VoF from bubbles/droplets. */

bool do_remove   = true;   // do it or not
int min_size     = 3;      // size of the bubble/droplet in grid cells
double threshold = 1.0e-4; // VoF and 1.0-VoF threshold
double cirp_th   = 0.20;   // minimum threshold for position
double eta_m0    = 1.0;    // initial mean eta0

/**
   Define some output frequencies. */

double tout_glo_my = 64.00; // output frequency of global observables (don't go beyond otherwise it does not print)
double tout_tag_my = 64.00; // output frequency of tagging
double tout_eta_my = 32.00; // output frequency of interfacial quantity
double tout_slc_my = 32.00; // output frequency of slices
double tout_pro_my = 32.00; // output frequency of 1d profile
double tout_mov_my = 4.00;  // output frequency of the movie
double tout_cut_my = 4.00;  // output frequency of some 2D cuts
double tout_res_my = 2.00;  // output frequency of restart
double tout_rbk_my = 2.00;  // output frequency of the backup restarting files

/**
   For the restarting step. */

int counter_max = 2;
int counter = 0;

/**
   The error on the components of the velocity field used for adaptive
   refinement. */

double uemax  = 0.01;   // we are going to change later
double uwemax = 0.001;  // we are going to change later
double femax  = 0.0001; // we keep this value

/**
   The density and viscosity ratios are those of air and water. 
   note: the ratio is air properties over water properties */

double RHO_RATIO = 1.225/1000.;
double MU_RATIO  = 18.31e-6/10.0e-4;

/** 
   We define a water velocity field for the refinement criteria */

vector u_water[]; 

/**
   We need to store the variable forcing field av[]. */

face vector av[];

/**
   Input paramters are passed in from the command line. */

int main(int argc, char *argv[]) {

  /**
     Provide the inputs */

  maxruntime (&argc, argv);

  if (argc > 1)
    Re_ast  = atof(argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    Re_wave = atof(argv[3]);
  if (argc > 4)
    UstarRATIO = atof(argv[4]);
  if (argc > 5)
    ak = atof(argv[5]);
  if (argc > 6)
    r_L0lam = atof(argv[6]);
  if (argc > 7)
    rho_r = atof(argv[7]);
  if (argc > 8)
    mu_r = atof(argv[8]);
  if (argc > 9)
    MAXLEVEL = atoi(argv[9]);
  if (argc > 10)
    MINLEVEL = atoi(argv[10]);
  if (argc > 11)
    RELEASETIME = atof(argv[11]);
  if (argc > 12)
    uemaxRATIO = atof(argv[12]);
  if (argc > 13)
    end_sim  = atof(argv[13]);
  if (argc > 14)
    do_eta_loc = atoi(argv[14]); // do_eta_loc = 1 --> true, else --> false
  if (argc > 15)
    do_profile = atoi(argv[15]); // do_profile = 1 --> true, else --> false
  if (argc > 16)
    do_fields  = atoi(argv[16]); // do_fields  = 1 --> true, else --> false
  if (argc > 17)
    do_tagging = atoi(argv[17]); // do_tagging = 1 --> true, else --> false
  if (argc > 18)
    prt_res    = atoi(argv[18]);
  if (argc > 19)
    from_pr    = atoi(argv[19]);
  if (argc > 20)
    st_wave    = atoi(argv[20]);
  if (argc > 21)
    RANDOM     = atoi(argv[21]);
  
  if (argc < 22) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 22-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "maximum runtime = %.10e seconds\n", _maxruntime), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " Re_ast = %.10e\n ", Re_ast), fflush (stderr);
  fprintf(stderr, " Bo     = %.10e\n ", BO), fflush (stderr);
  fprintf(stderr, " Re_wave = %.10e\n ", Re_wave), fflush (stderr);
  fprintf(stderr, " UstarRATIO = %.10e\n ", UstarRATIO), fflush (stderr);
  fprintf(stderr, " a_0k  = %.10e\n ", ak), fflush (stderr);
  fprintf(stderr, " r_L0lam = %.10e\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAX LEVEL  = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL  = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " RELEASETIME = %.10e\n ", RELEASETIME), fflush (stderr);
  fprintf(stderr, " uemaxRatio = %.10e\n ", uemaxRATIO), fflush (stderr);
  fprintf(stderr, " tend  = %.10e\n ", end_sim), fflush (stderr);
  fprintf(stderr, " do_eta_loc = %d\n ", do_eta_loc), fflush (stderr);
  fprintf(stderr, " do_profile = %d\n ", do_profile), fflush (stderr);
  fprintf(stderr, " do_fields  = %d\n ", do_fields), fflush (stderr);
  fprintf(stderr, " do_tagging = %d\n ", do_tagging), fflush (stderr);
  fprintf(stderr, " Printing resolution = %d\n ", prt_res), fflush (stderr);
  fprintf(stderr, " Restart from precursor = %d\n ", from_pr), fflush (stderr);
  fprintf(stderr, " Initialize or not with a 3rd order Stokes wave = %d\n ", st_wave), fflush (stderr);
  fprintf(stderr, " RANDOM = %d\n ", RANDOM), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     The domain is a cubic box centered on the origin and of length
     $L_0=2*\pi$, periodic in the x- and z-directions. */

  L0   = 2.0*pi;
  h_   = L0/(2.0*pi); // Water depth set L0/(2*pi) following Wu, Popinet & Deike JFM2022
  lam_ = L0/r_L0lam;
  k_   = 2.0*pi/lam_;
  origin (-L0/2., 0, -L0/2.);

  /**
     Set the boundary conditions. */
  
  u.r[top]    = neumann(0);
  u.r[bottom] = neumann(0);
  u.n[top]    = dirichlet(0);
  u.n[bottom] = dirichlet(0);
  u.t[top]    = neumann(0);
  u.t[bottom] = neumann(0);
  p[top]      = neumann(0);
  p[bottom]   = neumann(0);

  periodic (right);
  periodic (front);

  /**
     Here we set the densities and viscosities corresponding to the
     parameters above. Note that these variables are defined in two-phase.h already. 
     To set the dimensional parameters we can use three strategies:
       1. $u_\ast$ to a fixed value and we modify $c$. This strategy allows to keep the same $Re_\ast$ and (1a) modify either the wave Reynolds number and keep fixed the viscosity ratio or (2a) modify the viscosity ratio and keep fixed the wave Reynolds number (strategy of Wu, Popinet and Deike, JFM2022);
       2. $c$ to a fixed value and we modify $u_\ast$. This strategy allows to keep the same viscosity ratio and wave Reynolds number, but we change $Re_\ast$.
     */

  rho2     = rho_r;
  mu2      = mu_r;
  rho1     = rho2/RHO_RATIO;
  Ustar    = (mu2*Re_ast)/(rho2*(L0-h_));
  c_       = Ustar/UstarRATIO;
  g_       = k_*sq(c_)/( 1.0+(rho1-rho2)/(rho1*BO) ); // tune g_ 
  f.sigma  = g_*(rho1-rho2)/(BO*sq(k_));
  G.y      = -g_;
  mu1      = rho1*c_*(2.*pi/k_)/Re_wave;
  alter_MU = (mu2/mu1)/MU_RATIO;
  T0_      = 2.*pi/sqrt(g_*k_ + f.sigma*k_*k_*k_/rho1);

  /**
     Give the address of av[] to a[] so that acceleration can be changed */

  a = av;
  
  /**
     We set the refinement level */

  uemax  = uemaxRATIO*Ustar;
  uwemax = 0.001*c_;
  fprintf (stderr, "RELEASETIME = %.10e, uemax = %.10e, uwemax = %.10e, femax = %.10e\n", 
		   RELEASETIME, uemax, uwemax, femax);
  
  /**
     We run! */

  run();

}

/** 
   Specify the interface shape using a third-order Stokes wave. */

double WaveProfile (double x, double y) {

  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 + h_;

}

/**
   Set profile either with a third-order Stokes wave (set_profile function)
   or wave spectrum (import_profile function). */

void set_profile (scalar f) {
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = -(y-WaveProfile(x,z));
  }
  fractions (phi,f);
  delete ({phi});
  stats stat1 = statsf(f); // log some stats
  fprintf (stderr, "f_min, f_sum and f_max, %8E %8E %8E \n",
                   stat1.min, stat1.sum, stat1.max), fflush (stderr);
}

#include "./spectrum.h" // Used for new input method (the spectrum info)
void import_profile (scalar f) {
  power_input();
  dkx_ = kx_[1] - kx_[0];
  dky_ = ky_[1] - ky_[0];
  fprintf (stderr, "dkx = %g, dky = %g\n", dkx_, dky_), fflush (stderr);
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = -(y-(wave(x,z)+h_));
  }
  fractions (phi,f);
  delete ({phi});
  stats stat1 = statsf(f); // log some stats
  fprintf (stderr, "f_min, f_sum and f_max, %8E %8E %8E \n",
                   stat1.min, stat1.sum, stat1.max), fflush (stderr);
}

# define POPEN(name, mode) fopen (name ".ppm", mode)

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " Re_ast = %.10e\n ", rho2*Ustar*(L0-h_)/mu2), fflush (stderr);
  fprintf(stderr, " Bo    = %.10e\n ", g_*(rho1-rho2)/(f.sigma*sq(k_))), fflush (stderr);
  fprintf(stderr, " Re_wave = %.10e\n ", rho1*c_*(2.0*pi/k_)/mu1), fflush (stderr);
  fprintf(stderr, " UstarRATIO = %.10e\n ", Ustar/c_), fflush (stderr);
  fprintf(stderr, " a_0k  = %.10e\n ", ak), fflush (stderr);
  fprintf(stderr, " r_L0lam = %.10e\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAX LEVEL  = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL  = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " RELEASETIME  = %.10e\n ", RELEASETIME), fflush (stderr);
  fprintf(stderr, " uemaxRATIO  = %.10e\n ", uemaxRATIO), fflush (stderr);
  fprintf(stderr, " end_sim  = %.10e\n ", end_sim), fflush (stderr);
  fprintf(stderr, " do_eta_loc = %d\n ", do_eta_loc), fflush (stderr);
  fprintf(stderr, " do_profile = %d\n ", do_profile), fflush (stderr);
  fprintf(stderr, " do_fields  = %d\n ", do_fields), fflush (stderr);
  fprintf(stderr, " do_tagging  = %d\n ", do_tagging), fflush (stderr);
  fprintf(stderr, " gravity = %.10e\n ", g_), fflush (stderr);
  fprintf(stderr, " Re_wave = %.10e\n ", rho1*c_*(2*pi/k_)/mu1), fflush (stderr);
  fprintf(stderr, " rho_r = %.10e\n ", (rho1/rho2)), fflush (stderr);
  fprintf(stderr, " mu_r  = %.10e\n ", (mu1/mu2)), fflush (stderr);
  fprintf(stderr, " alter_MU = %.10e\n ", alter_MU), fflush (stderr);
  fprintf(stderr, " wave speed = %.10e\n ", c_), fflush (stderr);
  fprintf(stderr, " T0 = %.10e\n ", T0_), fflush (stderr);

  fprintf(stderr, " We_w = %.10e\n ", rho1*sq(c_)*(2.0*pi/k_)/f.sigma ), fflush (stderr);
  fprintf(stderr, " Fr_w = %.10e\n ", sq(c_)/(g_*2.0*pi/k_) ), fflush (stderr);
  fprintf(stderr, " We_a = %.10e\n ", rho2*sq(Ustar)*(L0-h_)/f.sigma ), fflush (stderr);
  fprintf(stderr, " Fr_a = %.10e\n ", sq(Ustar)/(g_*(L0-h_)) ), fflush (stderr);
  fprintf(stderr, " Surface tension = %.10e\n ", f.sigma), fflush (stderr);

  fprintf(stderr, " Printing resolution = %d\n ", prt_res), fflush (stderr);
  fprintf(stderr, " Restart from precursor = %d\n ", from_pr), fflush (stderr);
  fprintf(stderr, " Initialize or not with a 3rd order Stokes wave = %d\n ", st_wave), fflush (stderr);
  fprintf(stderr, " RANDOM = %d\n ", RANDOM), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     Create directories for better organization of the output */
  
  fprintf(stderr, "Create directories if needed\n"), fflush (stderr);

  char comm[80];
  sprintf (comm, "mkdir -p profiles");
  system(comm);
  sprintf (comm, "mkdir -p eta");
  system(comm);
  sprintf (comm, "mkdir -p eta/eta_loc");
  system(comm);
  sprintf (comm, "mkdir -p field");
  system(comm);
  sprintf (comm, "mkdir -p restart_bk");
  system(comm);
  sprintf (comm, "mkdir -p tagging");
  system(comm);
  sprintf (comm, "mkdir -p tagging/tagging_f1");
  system(comm);
  sprintf (comm, "mkdir -p tagging/tagging_f2");
  system(comm);
  sprintf (comm, "mkdir -p budgets");
  system(comm);
  sprintf (comm, "mkdir -p slices/");
  system(comm);

  /**
     Restart from the latest dump files or initialize a new simulation */

  if (from_pr == 1) {

    // Restore the simulation
    scalar cs = new scalar; // temporary file just for restore 
    fprintf(stderr, "Scalar created\n"), fflush (stderr);
    restore ("restart_sp.bin");
    delete ({cs}); // deleted since we use f
    fprintf(stderr, "Simulation restarts from a dumped file (SINGLE-PHASE)\n"), fflush (stderr);

    // Initalize the VoF profile
    if(st_wave == 1) {
      fprintf(stderr, "We use a 3rd-order Stokes wave!\n"), fflush (stderr);
      set_profile (f);
    }
    else {
      fprintf(stderr, "We import a profile previously generated!\n"), fflush (stderr);
      import_profile (f);
    }

    // Adjust the density based on the VoF profile, rho1 and rho2 
    foreach() {
      rhov[] = cm[]*rho(sf[]); 
    }
    
    // Check air velocity
    vector uA_tmp = new vector;
    foreach() {
      foreach_dimension() {
	uA_tmp.x[] = u.x[]*(1.0-f[]);
      }
    }
    stats uxA_tmp = statsf (uA_tmp.x);
    stats uyA_tmp = statsf (uA_tmp.y);
    stats uzA_tmp = statsf (uA_tmp.z);
    delete({uA_tmp.x,uA_tmp.y,uA_tmp.z});
    fprintf(stderr, "Air Max ux = %.10e\n", uxA_tmp.max), fflush (stderr);
    fprintf(stderr, "Air Min ux = %.10e\n", uxA_tmp.min), fflush (stderr);
    fprintf(stderr, "Air Max uy = %.10e\n", uyA_tmp.max), fflush (stderr);
    fprintf(stderr, "Air Min uy = %.10e\n", uyA_tmp.min), fflush (stderr);
    fprintf(stderr, "Air Max uz = %.10e\n", uzA_tmp.max), fflush (stderr);
    fprintf(stderr, "Air Min uz = %.10e\n", uzA_tmp.min), fflush (stderr);

    // Output an image for checking 
    char stg[80];
    char file[99];
    clear();
    view (fov = 27.5,  
          theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
          width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
    box(false, lc = {1,1,1}, lw = 0.1);
    draw_vof ("f", color = "u.x");
    //squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
    //squares ("u.z", linear = true, n = {1,0,0}, alpha = -L0/2.0);
    squares ("u.x", linear = true, map = rain_cm   , n = {0,0,1}, alpha = -L0/2.0);
    squares ("u.z", linear = true, map = balance_cm, n = {1,0,0}, alpha = -L0/2.0);
    cells (n = {1,0,0}, alpha = -L0/2.0);
    sprintf (stg, "t = %0.3f", 2.0*pi*t/T0_);
    draw_string (stg, size = 30); 
    sprintf (file, "./3D_sp_%09d.ppm", i);
    save(file);
    /* 
    u_water.x.nodump = true;
    u_water.y.nodump = true;
    u_water.z.nodump = true;
    dump ("test_tw.bin");
    //return 1;
    */

    /**
       Add a check on RELEASETIME. */

    if (t > RELEASETIME) {
      fprintf(stderr, "The restarting time t is larger than the RELEASETIME\n"), fflush (stderr);
      fprintf(stderr, "Please increase RELEASETIME otherwise the water velocity is not initialized\n"), fflush (stderr);
      return 1;
    }

  }
  else {
    if (restore ("restart.bin")) {
      fprintf(stderr, "Simulation restarts from a dumped file (TWO-PHASE)\n"), fflush (stderr);
    }
    else {
      fprintf(stderr, "Restarting files absent! Please provide them.\n"), fflush (stderr);
      return 1;
    }
  }
   
}

/**
   We log some useful information. */

event log_simulation (i += 10) {
 
  int min_level = +100, max_level = -100;
  foreach(reduction(min:min_level),reduction(max:max_level)) {
    if (level > max_level) max_level = level;
    if (level < min_level) min_level = level;
  }

  if (pid() == 0) { 
    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_simulation.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %d %d\n", t, 1.0*i, dt, max_level, min_level);
    fclose(log_sim);
  }

}

/**
  Keep the wave profile and set the wave velocity to 0, until the RELEASETIME. */
  
event set_wave(i=0; i++; t<RELEASETIME) {

  if(st_wave == 1) {
    fprintf(stderr, "We use a 3rd-order Stokes wave!\n"), fflush (stderr);
    set_profile (f);
  }
  else {
    fprintf(stderr, "We import a profile previously generated!\n"), fflush (stderr);
    import_profile (f);
  }

  foreach() {
    foreach_dimension() {
      u.x[] = (1.0 - f[])*u.x[];
    }
  }

}

/**
   Release the wave at RELEASETIME. We don't do any adaptation at this step. 
   */

// TO-DO: maybe include g correction by Bond number. High Bond should be fine.

double u_x (double x, double y) {

  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*3.*k_*cos(3.*k_*x);

}

double u_y (double x, double y) {

  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*sin(3.*k_*x);

}

event start(t = RELEASETIME) {
  // A slightly changed version of the test/stokes.h wave function as y = 0 at the bottom now so y+h -> y
 
  /*	
  fprintf(stderr, "t = RELEASETIME, we initialize the liquid velocity\n"), fflush (stderr);
  fraction (f, WaveProfile(x,z)-y);
  foreach () {
    u.x[] += u_x(x, y-h_)*f[];
    u.y[] += u_y(x, y-h_)*f[];
  }
  //boundary ((scalar *){u});
  */
  
  if(st_wave == 1) {
    fprintf(stderr, "We initialize the liquid velocity with u_x_st and u_y_st\n"), fflush (stderr);
    set_profile (f);
    foreach () {
      u.x[] += u_x(x, y-h_)*f[]; // stokes wave profile
      u.y[] += u_y(x, y-h_)*f[]; // stokes wave profile
    }
  }
  else {
    fprintf(stderr, "We initialize the liquid velocity with u_x_ge, u_y_ge and u_z_ge\n"), fflush (stderr);
    import_profile (f);
    foreach () {
      u.x[] += u_x_ge(x, z, y-h_)*f[]; // profile coherent with the wave spectrum
      //u.y[] += u_y_ge(x, z, y-h_)*f[]; // profile coherent with the wave spectrum
      u.y[] += u_z_ge(x, z, y-h_)*f[]; // profile coherent with the wave spectrum
      u.z[] += u_y_ge(x, z, y-h_)*f[]; // profile coherent with the wave spectrum
    }
  }

}

/**
   Forcing term to sustain the flow equivalent to the pressure gradient in x direction. */

event acceleration (i++) {

  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x) {
    double ff = (f[] + f[-1])/2.;
    av.x[] += ampl*(1.0-ff);
  }

}

/** 
   Output video and field in uncompressed format 
   (we can compress later using convert .ppm to .mpg and open with mplayer) */

//# define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_mov_my) {
//event movies (i += 2) {

  //fprintf(stderr, "I output the movies every t=%.10e\n", T0_/tout_mov); fflush (stderr);
  fprintf(stderr, "I output the movies every t=%.10e\n", T0_/tout_mov_my); fflush (stderr);
  
  scalar l2[];
  //lambda2 (u, l2);
  scalar omega[];
  vorticity (u, omega);

  char stg[80];

  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "u.x");
  //squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  //squares ("omega", linear = true, n = {1,0,0}, alpha = -L0/2.0);
  squares ("u.x", linear = true, map = rain_cm   , n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.z", linear = true, map = balance_cm, n = {1,0,0}, alpha = -L0/2.0);
  //squares ("omega", linear = true, map = jet       , n = {1,0,0}, alpha = -L0/2.0);
  //cells (n = {1,0,0}, alpha = -L0/2.0);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30); 
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
    fflush (fp);
  }

  clear();
  view (fov = 22.5, camera = "front", ty = -0.50,
        //theta = 0.0, phi = -pi/240.0, psi = 0.0, ty = -0.5,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("f", color = "omega");
  //squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.x", linear = true, map = rain_cm, n = {0,0,1}, alpha = -L0/2.0);
  sprintf (stg, "t = %0.3f", 2.0*pi*(t-RELEASETIME)/T0_);
  draw_string (stg, size = 30, pos = 3); 
  {
    static FILE * fp = POPEN ("int", "a");
    save (fp = fp);
    fflush (fp);
  }

}

event out_pro (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_pro_my) {
//event out_pro (i += 2) {

  if (do_profile == 1) {

    //fprintf(stderr, "I output the vertical profiles file every t=%.10e\n", T0_/tout_pro);
    fprintf(stderr, "I output the vertical profiles file every t=%.10e\n", T0_/tout_pro_my);
    profile_output (u, i, MAXLEVEL);
    fflush (stderr);
   
    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./profiles/log_pro.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%.10e %.10e\n", t, 1.0*i);
      fclose(log_sim);
    
    }

  }

}

event slice_stat (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_slc_my) {
//event slice_stat (i += 2) {

  if (do_fields == 1) {

    //fprintf(stderr, "I output the slices every t=%.10e\n", T0_/tout_slc), fflush (stderr);
    fprintf(stderr, "I output the slices every t=%.10e\n", T0_/tout_slc_my), fflush (stderr);

    char filename[100];
    int res        = prt_res;
    bool do_linear = true;
    bool print_bin = true;

    sprintf (filename, "./field/ux_2d_avg_%09d.bin", i); // x-vel
    output_2d_span_avg (filename,u.x  ,res, do_linear, print_bin);
    sprintf (filename, "./field/uy_2d_avg_%09d.bin", i); // y-vel
    output_2d_span_avg (filename,u.y  ,res, do_linear, print_bin);
    sprintf (filename, "./field/uz_2d_avg_%09d.bin", i); // z-vel
    output_2d_span_avg (filename,u.z  ,res, do_linear, print_bin);
    sprintf (filename, "./field/fv_2d_avg_%09d.bin", i); // VoF
    output_2d_span_avg (filename,f    ,res, do_linear, print_bin);
    sprintf (filename, "./field/pr_2d_avg_%09d.bin", i); // pressure
    output_2d_span_avg (filename,p    ,res, do_linear, print_bin);

    // Compute some quantities (we defined as new scalar so we can delete as soon as they are used)
    scalar omega = new scalar;
    vorticity (u, omega);
    sprintf (filename, "./field/ow_2d_avg_%09d.bin", i); // omega
    output_2d_span_avg (filename,omega,res, do_linear, print_bin); delete ({omega});
    
    scalar uv = new scalar; scalar duy = new scalar; 
    scalar diss = new scalar; scalar ke = new scalar;
    scalar vke = new scalar;
    foreach () {
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
      double sqterm = 2.*( sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
          		   sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
        		   sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz) );
      uv[]   = u.x[]*u.y[];
      duy[]  = dudy; 
      diss[] = sqterm;
      ke[]   = sq(u.x[])+sq(u.y[])+sq(u.z[]);
      vke[]  = u.y[]*ke[];
    }

    sprintf (filename, "./field/uv_2d_avg_%09d.bin", i); // uv
    output_2d_span_avg (filename,uv   ,res, do_linear, print_bin); delete({uv});
    sprintf (filename, "./field/du_2d_avg_%09d.bin", i); // dudy
    output_2d_span_avg (filename,duy  ,res, do_linear, print_bin); delete({duy});
    sprintf (filename, "./field/ke_2d_avg_%09d.bin", i); // ke
    output_2d_span_avg (filename,ke   ,res, do_linear, print_bin); delete({ke});
    sprintf (filename, "./field/vk_2d_avg_%09d.bin", i); // vke
    output_2d_span_avg (filename,vke  ,res, do_linear, print_bin); delete({vke});
    sprintf (filename, "./field/di_2d_avg_%09d.bin", i); // dissipation
    output_2d_span_avg (filename,diss ,res, do_linear, print_bin); delete({diss});
    
    scalar dkdy = new scalar;
    foreach() {
      dkdy[] = (ke[0,1]-ke[0,-1])/(2.*Delta);
    }
    sprintf (filename, "./field/dk_2d_avg_%09d.bin", i); // dkdy
    output_2d_span_avg (filename,dkdy ,res, do_linear, print_bin); delete({dkdy});

    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./field/log_field.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%.10e %.10e\n", t, 1.0*i);
      fclose(log_sim);
    
    }

  }
  //return 1;

}

event slice_spec (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_cut_my) {
//event slice_spec (i += 2) {

  if (do_fields == 1) {

    //fprintf(stderr, "I output the slices every t=%.10e\n", T0_/tout_slc), fflush (stderr);
    fprintf(stderr, "I output the cuts every t=%.10e\n", T0_/tout_slc_my), fflush (stderr);

    char filename[100];
    int res        = prt_res;
    bool do_linear = true;

    int len = 2;
    double pos[len]; 
    for (int i = 0; i <= len; i++) {
      double stp = L0/4.0;
      pos[i] = -L0/4.0+i*stp;
    }

    scalar uv[];
    foreach () {
      uv[] = u.x[]*u.y[];
    }

    scalar omega[];
    vorticity (u, omega);

    for (int s = 0; s <= len; s++) {
      double stp = pos[s];
      sprintf (filename, "./slices/ux_2d_%09d_%03d.bin", i, s); // x-vel
      sliceXY (filename, u.x, stp, res, do_linear);
      sprintf (filename, "./slices/uy_2d_%09d_%03d.bin", i, s); // y-vel
      sliceXY (filename, u.y, stp, res, do_linear);
      sprintf (filename, "./slices/uz_2d_%09d_%03d.bin", i, s); // z-vel
      sliceXY (filename, u.z, stp, res, do_linear);
      sprintf (filename, "./slices/fv_2d_%09d_%03d.bin", i, s); // VoF
      sliceXY (filename, f  , stp, res, do_linear);
      sprintf (filename, "./slices/uv_2d_%09d_%03d.bin", i, s); // uv
      sliceXY (filename, uv , stp, res, do_linear);
      sprintf (filename, "./slices/om_2d_%09d_%03d.bin", i, s); // vorticity
      sliceXY (filename, omega, stp, res, do_linear);
    }

    if (pid() == 0) {
    
      char name_1[80];
      sprintf (name_1,"./slices/log_cut.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%.10e %.10e\n", t, 1.0*i);
      fclose(log_sim);
    
    }

  }
  //return 1;

}

/**
   We want to compute some global observables */

void output_global_obs_1 (char * fname, int istep, scalar c, scalar p_a, double stp_eta) {

  // Preliminary calculations
  
  // --> position
  scalar pos[];
#if dimension > 2
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#else
  coord G = {0.,1.}, Z = {0.,0.};
#endif
  position (c, pos, G, Z);

  // --> my diagnosis of eta_my, area_my, amp_my (without offsets)
  double area_my = 0; double eta_my = 0; 
  foreach(reduction(+:area_my), reduction(+:eta_my)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
      //if (abs(pos[]-eta_m0) < cirp_th) {
        coord n      = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar    = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
	area_my += ar;
	eta_my  += ar*pos[]; // already defined at the interface
      }
    }
  }
  area_my += 1.0e-12; // to prevent arithmetic failure if area = 0
  eta_my   = eta_my/area_my;
  fprintf(stderr, "area %.10e, eta %.10e\n", area_my, eta_my);

  // --> my diagnosis of eta_my, area_my, amp_my (we want to just do it at the interface)
  double amp2_my = 0;
  foreach(reduction(+:amp2_my)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
      //if (abs(pos[]-eta_m0) < cirp_th) {
        coord n = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
        amp2_my += ar*2.0*sq(pos[]-eta_my);
      }
    }
  }
  amp2_my = amp2_my/area_my;

  // --> compute some mean quantities to be used later (with offset)
  double area = 0; double eta_m = 0; double pr_m = 0;
  double u_xi = 0; double u_yi  = 0; double u_zi = 0;
  foreach(reduction(+:area), reduction(+:eta_m), reduction(+:pr_m),
	  reduction(+:u_xi), reduction(+:u_yi) , reduction(+:u_zi)) {
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
      //if (abs(pos[]-eta_m0) < cirp_th) {
	coord n      = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar    = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
	double eta   = pos[]; // must be here since here it is defined
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        double zc = z + Delta*pp.z;
#if dimension > 2
	Point point1 = locate (xc, yc, zc);
#else
	Point point1 = locate (xc, yc);
#endif
	if (point1.level > 0) {
        
	  POINT_VARIABLES;
	
	  area  += ar;
	  eta_m += ar*eta; // already defined at the interface

          pr_m += ar*my_interpolation(point, p_a, xc, yc, zc);
          u_xi += ar*my_interpolation(point, u.x, xc, yc, zc);
          u_yi += ar*my_interpolation(point, u.y, xc, yc, zc);
#if dimension > 2
          u_zi += ar*my_interpolation(point, u.z, xc, yc, zc);
#else
	  u_zi = 0.;
#endif

	}
      }
    }
  }
  area  += 1.0e-12; // to prevent arithmetic failure if area = 0
  eta_m  = eta_m/area; // we need the mean surface elevation
  pr_m   = pr_m/area;  // we need the mean pressure
  u_xi   = u_xi/area; u_yi = u_yi/area; u_zi = u_zi/area;
  //fprintf(stderr, "area %.10e, eta %.10e, pr_m %.10e\n", area, eta_m, pr_m);
  //fprintf(stderr, "u.x %.10e, u.y %.10e, u.z %.10e\n", u_xi, u_yi, u_zi);

  // --> quantities to be probed at the interface
  // n.b.: for the stress we can use just u since it contains derivatives
  scalar Sxx[]; scalar Syy[]; scalar Szz[];
  scalar Sxy[]; scalar Sxz[]; scalar Syz[];
  foreach() {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    Sxx[] = dudx+dudx; 
    Syy[] = dvdy+dvdy;
    Szz[] = dwdz+dwdz;
    Sxy[] = dudy+dvdx;
    Sxz[] = dudz+dwdx;
    Syz[] = dvdz+dwdy;
  }
  //boundary({Sxx,Syy,Szz,Sxy,Sxz,Syz}); // must be kept
  boundary({Sxx,Syy,Szz,Sxy,Sxz,Syz,p_a,u.x,u.y,u.z}); // must be kept

  // --> compute the amplitude, stress, velocity at the interface, energy fluxes
  double amp2  = 0.0;
  double mf_px = 0.0; double mf_py = 0.0; double mf_pz = 0.0; // mom. flux pressure
  double mp_px = 0.0; double mp_py = 0.0; double mp_pz = 0.0; // mean pressure contribution
  double mf_vx = 0.0; double mf_vy = 0.0; double mf_vz = 0.0; // mom. flux viscous stress
  double ef_pn = 0.0; double ef_pp = 0.0; // en. flux pressure (without-with subtraction)
  double ef_fn = 0.0; double ef_fp = 0.0; // en. flux pressure (no mean pressure)
  double ef_vn = 0.0; double ef_vp = 0.0; // en. flux viscous stress
  double Sxx_gp1 = 0.0; double Syy_gp1 = 0.0; double Szz_gp1 = 0.0; // stress component 11-22-33 (projected along the normal)
  double Sxy_gp1 = 0.0; double Sxz_gp1 = 0.0; double Syz_gp1 = 0.0; // stress component 11-13-23 (projected along the normal)
  double Sxy_gp2 = 0.0; double Sxz_gp2 = 0.0; double Syz_gp2 = 0.0; // stress component 11-13-23 (projected along the normal)
  double Sxx_gu1 = 0.0; double Syy_gu1 = 0.0; double Szz_gu1 = 0.0; // stress component 11-22-33 (unprojected along the normal)
  double Sxy_gu1 = 0.0; double Sxz_gu1 = 0.0; double Syz_gu1 = 0.0; // stress component 11-13-23 (unprojected along the normal)
  double Sxy_gu2 = 0.0; double Sxz_gu2 = 0.0; double Syz_gu2 = 0.0; // stress component 11-13-23 (unprojected along the normal)
  double ep_gs = 0.0;
  foreach(reduction(+:amp2), // amplitude 
	  reduction(+:mf_px), reduction(+:mf_py), reduction(+:mf_pz), // pressure - momentum flux (x,y,z)
	  reduction(+:mp_px), reduction(+:mp_py), reduction(+:mp_pz), // pressure - momentum flux (x,y,z) - (mean pressure)
	  reduction(+:mf_vx), reduction(+:mf_vy), reduction(+:mf_vz), // viscous  - momentum flux (x,y,z)
	  reduction(+:ef_pn), reduction(+:ef_pp), // pressure en flux (without-with subtraction)
	  reduction(+:ef_fn), reduction(+:ef_fp), // pressure en flux (without the mean pressure)
	  reduction(+:ef_vn), reduction(+:ef_vp), // viscous en flux (without-with subtraction)
	  reduction(+:Sxx_gp1), reduction(+:Syy_gp1), reduction(+:Szz_gp1), // stress component 11-22-33 (projected)
	  reduction(+:Sxy_gp1), reduction(+:Sxz_gp1), reduction(+:Syz_gp1), // stress component 11-13-23 (projected)
	  reduction(+:Sxy_gp2), reduction(+:Sxz_gp2), reduction(+:Syz_gp2), // stress component 11-13-23 (projected)
	  reduction(+:Sxx_gu1), reduction(+:Syy_gu1), reduction(+:Szz_gu1), // stress component 11-22-33 (unprojected)
	  reduction(+:Sxy_gu1), reduction(+:Sxz_gu1), reduction(+:Syz_gu1), // stress component 11-13-23 (unprojected)
	  reduction(+:Sxy_gu2), reduction(+:Sxz_gu2), reduction(+:Syz_gu2), // stress component 11-13-23 (unprojected)
	  reduction(+:ep_gs)) { // gravitational energy due to surface tension
    if (interfacial (point, c)) {
      //if (point.level == MAXLEVEL) {
      if (point.level == MAXLEVEL && abs(pos[]-eta_m0) < cirp_th) {
      //if (abs(pos[]-eta_m0) < cirp_th) {
	coord n      = interface_normal(point, c), pp;
        double alpha = plane_alpha (c[], n);
        double ar    = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &pp);
	double eta   = pos[]; // must be here since here it is defined
	normalize (&n); // |n| = 1, we should normalize since we use the normals for online calculations
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        double zc = z + Delta*pp.z;
#if dimension > 2
	Point point1 = locate (xc, yc, zc);
#else
	Point point1 = locate (xc, yc);
#endif
	if (point1.level > 0) {
        
	  POINT_VARIABLES;
	
	  amp2 += ar*2.0*sq(eta-eta_m); // amplitude
 
          double pr_int = 0;   
	  double Sxx_int = 0; double Syy_int = 0; double Szz_int = 0;
	  double Sxy_int = 0; double Sxz_int = 0; double Syz_int = 0;
	  double ux_int  = 0; double uy_int  = 0; double uz_int  = 0;

          pr_int  = my_interpolation(point, p_a, xc, yc, zc);
	  Sxx_int = my_interpolation(point, Sxx, xc, yc, zc);
	  Syy_int = my_interpolation(point, Syy, xc, yc, zc);
	  Szz_int = my_interpolation(point, Szz, xc, yc, zc);
	  Sxy_int = my_interpolation(point, Sxy, xc, yc, zc);
	  Sxz_int = my_interpolation(point, Sxz, xc, yc, zc);
	  Syz_int = my_interpolation(point, Syz, xc, yc, zc);
	  ux_int  = my_interpolation(point, u.x, xc, yc, zc);
	  uy_int  = my_interpolation(point, u.y, xc, yc, zc);
#if dimension > 2
          uz_int  = my_interpolation(point, u.z, xc, yc, zc);
#else
	  uz_int  = 0.;
#endif

	  // stress components (projected)
	  Sxx_gp1 += ar*mu2*Sxx_int*n.x;
	  Syy_gp1 += ar*mu2*Syy_int*n.y;
	  Szz_gp1 += ar*mu2*Szz_int*n.z;
	  Sxy_gp1 += ar*mu2*Sxy_int*n.y;
	  Sxz_gp1 += ar*mu2*Sxz_int*n.z;
	  Syz_gp1 += ar*mu2*Syz_int*n.z;
	  Sxy_gp2 += ar*mu2*Sxy_int*n.x;
	  Sxz_gp2 += ar*mu2*Sxz_int*n.x;
	  Syz_gp2 += ar*mu2*Syz_int*n.y;

	  // stress components (unprojected)
	  Sxx_gu1 += ar*mu2*Sxx_int;
	  Syy_gu1 += ar*mu2*Syy_int;
	  Szz_gu1 += ar*mu2*Szz_int;
	  Sxy_gu1 += ar*mu2*Sxy_int;
	  Sxz_gu1 += ar*mu2*Sxz_int;
	  Syz_gu1 += ar*mu2*Syz_int;
	  Sxy_gu2 += ar*mu2*Sxy_int;
	  Sxz_gu2 += ar*mu2*Sxz_int;
	  Syz_gu2 += ar*mu2*Syz_int;
	  
          // momentum flux --> pressure	  
	  mf_px += ar*( -(pr_int)*n.x ); // x-Mom. flux due to pressure work
	  mf_py += ar*( -(pr_int)*n.y ); // y-Mom. flux due to pressure work
	  mf_pz += ar*( -(pr_int)*n.z ); // z-Mom. flux due to pressure work
	  
	  mp_px += ar*( +(pr_m)*n.x ); // mean pressure contribution - x
	  mp_py += ar*( +(pr_m)*n.y ); // mean pressure contribution - y
	  mp_pz += ar*( +(pr_m)*n.z ); // mean pressure contribution - z
	  
	  // momentum flux --> viscous stress
	  mf_vx += ar*mu2*(n.x*Sxx_int + n.y*Sxy_int + n.z*Sxz_int); // x-Mom. flux due to viscous dissipation
          mf_vy += ar*mu2*(n.x*Sxy_int + n.y*Syy_int + n.z*Syz_int); // y-Mom. flux due to viscous dissipation
          mf_vz += ar*mu2*(n.x*Sxz_int + n.y*Syz_int + n.z*Szz_int); // z-Mom. flux due to viscous dissipation
          
	  // energy flux --> pressure
	  ef_pn += ar*( -(pr_int-pr_m)*( (ux_int)*n.x+(uy_int)*n.y+(uz_int)*n.z) ); // En. flux due to pressure (no vel. subtraction)
	  ef_pp += ar*( -(pr_int-pr_m)*( (ux_int-u_xi)*n.x+(uy_int-u_yi)*n.y+(uz_int-u_zi)*n.z) ); // En. flux due to pressure
	  ef_fn += ar*( -(pr_int)*( (ux_int)*n.x+(uy_int)*n.y+(uz_int)*n.z) ); // En. flux due to pressure (no vel. subtraction)
	  ef_fp += ar*( -(pr_int)*( (ux_int-u_xi)*n.x+(uy_int-u_yi)*n.y+(uz_int-u_zi)*n.z) ); // En. flux due to pressure

	  // energy flux --> viscous stress
	  ef_vn += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*(ux_int) +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*(uy_int) +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*(uz_int) ); // En. flux due to viscous dissipation (no vel. subtraction)
          ef_vp += ar*mu2*( ( Sxx_int*n.x + Sxy_int*n.y + Sxz_int*n.z )*(ux_int-u_xi) +
	                    ( Sxy_int*n.x + Syy_int*n.y + Syz_int*n.z )*(uy_int-u_yi) +
	                    ( Sxz_int*n.x + Syz_int*n.y + Szz_int*n.z )*(uz_int-u_zi) ); // En. flux due to viscous dissipation
	  
	  ep_gs += ar*f.sigma*(sqrt(1.0+sq(n.x)+sq(n.z))-1.0); // potential energy due to surface tension
 
	}
      }
    }
  }
  amp2 = amp2/area;

  // --> compute quantities at the top/bottom boundary (for y-mom in air and water)
  /*
  double pre_a_m = 0.0, pre_w_m = 0.0;
  double vol_a_m = 0.0, vol_w_m = 0.0;
  foreach(reduction(+:pre_a_m) reduction(+:pre_w_m)
          reduction(+:vol_a_m) reduction(+:vol_w_m)) {
    pre_a_m = pre_a_m + (1.0-f[])*p[]*dv();
    pre_w_m = pre_w_m + (0.0+f[])*p[]*dv();
    vol_a_m = vol_a_m + (1.0-f[])*dv();
    vol_w_m = vol_w_m + (0.0+f[])*dv();
  }
  pre_a_m /= vol_a_m; pre_w_m /= vol_w_m;
  */

  double pre_a_m = 0.0;
  foreach_boundary(top,reduction(+:pre_a_m)) {
    double ff = 0.5*(f[0,1,0]+f[]);
    double pf = 0.5*(p[0,1,0]+p[]);
    pre_a_m = pre_a_m + (1.0-ff)*pf*sq(Delta);
  }
  pre_a_m = pre_a_m/sq(L0);

  double pre_w_m = 0.0;
  foreach_boundary(bottom,reduction(+:pre_w_m)) {
    double ff = 0.5*(f[0,-1,0]+f[]);
    double pf = 0.5*(p[0,-1,0]+p[]);
    pre_w_m = pre_w_m + (0.0+ff)*pf*sq(Delta);
  }
  pre_w_m = pre_w_m/sq(L0);

  double py_a_flux = 0.0, tx_a_flux = 0.0, ty_a_flux = 0.0;
  foreach_boundary (top,reduction(+:py_a_flux),reduction(+:tx_a_flux)
		    reduction(+:ty_a_flux)) {
    double pre_int = 0.5*(p[0,1,0]+p[]); // we interpolate to the face
    py_a_flux += sq(Delta)*( -(pre_int) );  
    tx_a_flux += sq(Delta)*mu2*(u.x[0,1,0]-u.x[])/Delta;  
    ty_a_flux += sq(Delta)*mu2*(u.y[0,1,0]-u.y[])/Delta;  
  }
  double py_w_flux = 0.0, ty_w_flux = 0.0;
  foreach_boundary (bottom,reduction(+:py_w_flux),reduction(+:ty_w_flux)) {
    double pre_int = 0.5*(p[0,-1,0]+p[]); // we interpolate to the face
    py_w_flux += sq(Delta)*( -(pre_int) );  
    ty_w_flux += sq(Delta)*mu1*(u.x[]-u.x[0,-1,0])/Delta;  
  }

  if (pid() == 0) {
  
    fflush(stderr);
    FILE * global_obs = fopen (fname, "a");
    fprintf (global_obs, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", 
                         t, 1.0*istep, area, eta_m, k_*sqrt(amp2), area_my, eta_my, k_*sqrt(amp2_my), pr_m,  
			 mf_px, mf_py, mf_pz, 
			 mp_px, mp_py, mp_pz, 
			 mf_vx, mf_vy, mf_vz, 
			 u_xi, u_yi, u_zi, 
			 ef_pn, ef_pp,
			 ef_fn, ef_fp,
			 ef_vn, ef_vp,
			 ep_gs, 
		         Sxx_gp1, Syy_gp1, Szz_gp1, 
			 Sxy_gp1, Sxz_gp1, Syz_gp1, 
			 Sxy_gp2, Sxz_gp2, Syz_gp2, 
		         Sxx_gu1, Syy_gu1, Szz_gu1, 
			 Sxy_gu1, Sxz_gu1, Syz_gu1, 
			 Sxy_gu2, Sxz_gu2, Syz_gu2, 
			 pre_a_m, pre_w_m, 
			 py_a_flux,tx_a_flux,ty_a_flux,py_w_flux,ty_w_flux);
    fclose(global_obs);

  }

}

event bulk_budgets (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_glo_my) {
//event bulk_budgets (i += 2) {

  fprintf(stderr, "I output the bulk energy and dissipation every t=%.10e\n", T0_/tout_glo_my);

  // --> Kinetic, Potential energy, dissipation of water and air (bulk)
  double keWat = 0., gpeWat = 0.;
  double keAir = 0., gpeAir = 0.;
  double keTot = 0., gpeTot = 0.;
  foreach(reduction(+:keWat),reduction(+:gpeWat) 
	  reduction(+:keAir),reduction(+:gpeAir)
	  reduction(+:keTot),reduction(+:gpeTot)) {
    double norm2 = 0.;
    foreach_dimension() {
      norm2 += sq(u.x[]);
    }
    keWat  += rho1*norm2*(0.0+f[])*dv();
    keAir  += rho2*norm2*(1.0-f[])*dv();
    keTot  += norm2*dv();
    gpeWat += rho1*g_*y*(0.0+f[])*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();
    gpeTot += g_*y*dv();
  }
  double rates[2];
  dissipation_rate(mu1,mu2,u,rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  
  // --> Mean velocity (all three components) for both phases
  double vol_a_m = 0., vol_w_m = 0.;
  coord vel_a_m = {0.,0.,0.}, vel_w_m = {0.,0.,0.};
  foreach(reduction(+:vol_a_m),reduction(+:vol_w_m), 
	  reduction(+:vel_a_m),reduction(+:vel_w_m)) {
    vol_a_m += (1.0-f[])*dv();
    vol_w_m += (0.0+f[])*dv();
    foreach_dimension() {
      vel_a_m.x += u.x[]*(1.0-f[])*dv(); 
      vel_w_m.x += u.x[]*(0.0+f[])*dv(); 
    }
  }

  // --> Mean momentum/velocity/force tension contribution/work by surface tension
  scalar phi_tmp[];
  curvature (f, phi_tmp);
  coord mom_m = {0.,0.,0.}, vel_m = {0.,0.,0}, stm_t = {0.,0.,0.};
  double ste_t = 0.0;
  foreach(reduction(+:mom_m),reduction(+:vel_m),reduction(+:stm_t),reduction(+:ste_t)) {
    foreach_dimension() {
      coord int_ft = {0.,0.,0}, int_et = {0.,0.,0.};
      if(interfacial(point,f) && phi_tmp[] != nodata) {
        int_ft.x = phi_tmp[]*(f[1,0,0]-f[-1,0,0])/(2.*Delta);
        int_et.x = phi_tmp[]*(f[1,0,0]-f[-1,0,0])/(2.*Delta)*u.x[];
      }
      mom_m.x += u.x[]*rho[]*dv();
      vel_m.x += u.x[]*dv();
      stm_t.x += f.sigma*int_ft.x*dv(); 
      ste_t += f.sigma*int_et.x*dv(); 
    }
  }

  // --> Pressure gradient integral, viscous term integral (air phase)
  coord grap_a_int = {0.,0.,0.}, divu_a_int = {0.,0.,0.};
  coord grap_w_int = {0.,0.,0.}, divu_w_int = {0.,0.,0.};
  foreach(reduction(+:grap_a_int),reduction(+:divu_a_int)
	  reduction(+:grap_w_int),reduction(+:divu_w_int)) {
    foreach_dimension() {
      double lap_u = (u.x[1,0,0]-2.0*u.x[]+u.x[-1,0,0] +
                      u.x[0,1,0]-2.0*u.x[]+u.x[0,-1,0] +
                      u.x[0,0,1]-2.0*u.x[]+u.x[0,0,-1])/sq(Delta);
      grap_a_int.x += -(p[1,0,0]-p[-1,0,0])/(2.*Delta)*(1.0-f[])*dv();
      divu_a_int.x += mu2*lap_u*(1.0-f[])*dv();
      grap_w_int.x += -(p[1,0,0]-p[-1,0,0])/(2.*Delta)*(0.0+f[])*dv();
      divu_w_int.x += mu1*lap_u*(0.0+f[])*dv();
    }
  }

  // --> quantities to be probed at the interface
  // n.b.: for the stress we can use just u since it contains derivatives
  scalar Sxx[]; scalar Syy[]; scalar Szz[];
  scalar Sxy[]; scalar Sxz[]; scalar Syz[];
  foreach() {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    Sxx[] = dudx+dudx; 
    Syy[] = dvdy+dvdy;
    Szz[] = dwdz+dwdz;
    Sxy[] = dudy+dvdx;
    Sxz[] = dudz+dwdx;
    Syz[] = dvdz+dwdy;
  }

  // --> Another version of the pressure term / viscous term integrals
  coord divc_a_int1 = {0.,0.,0.}, divc_w_int1 = {0.,0.,0.};
  coord grap_a_int1 = {0.,0.,0.}, grap_w_int1 = {0.,0.,0.};
  coord divu_a_int1 = {0.,0.,0.}, divu_w_int1 = {0.,0.,0.};
  foreach(reduction(+:divc_a_int1),reduction(+:divc_w_int1),
          reduction(+:grap_a_int1),reduction(+:grap_w_int1),
	  reduction(+:divu_a_int1),reduction(+:divu_w_int1)) {
    double div_con = 0;
    div_con = (u.x[1,0,0]*u.x[1,0,0]-u.x[-1,0,0]*u.x[-1,0,0])/(2.*Delta) +
	      (u.x[0,1,0]*u.y[0,1,0]-u.x[0,-1,0]*u.y[0,-1,0])/(2.*Delta) +
              (u.x[0,0,1]*u.z[0,0,1]-u.x[0,0,-1]*u.z[0,0,-1])/(2.*Delta);
    divc_a_int1.x += div_con*(1.0-f[])*dv();
    divc_w_int1.x += div_con*(0.0+f[])*dv();
    div_con = (u.y[1,0,0]*u.x[1,0,0]-u.y[-1,0,0]*u.x[-1,0,0])/(2.*Delta) +
	      (u.y[0,1,0]*u.y[0,1,0]-u.y[0,-1,0]*u.y[0,-1,0])/(2.*Delta) +
              (u.y[0,0,1]*u.z[0,0,1]-u.y[0,0,-1]*u.z[0,0,-1])/(2.*Delta);
    divc_a_int1.y += div_con*(1.0-f[])*dv();
    divc_w_int1.y += div_con*(0.0+f[])*dv();
    div_con = (u.z[1,0,0]*u.x[1,0,0]-u.z[-1,0,0]*u.x[-1,0,0])/(2.*Delta) +
	      (u.z[0,1,0]*u.y[0,1,0]-u.z[0,-1,0]*u.y[0,-1,0])/(2.*Delta) +
              (u.z[0,0,1]*u.z[0,0,1]-u.z[0,0,-1]*u.z[0,0,-1])/(2.*Delta);
    divc_a_int1.z += div_con*(1.0-f[])*dv();
    divc_w_int1.z += div_con*(0.0+f[])*dv();
    double div_pre = 0;
    div_pre = -(p[1,0,0]-p[-1,0,0])/(2.*Delta);
    grap_a_int1.x += div_pre*(1.0-f[])*dv();
    grap_w_int1.x += div_pre*(0.0+f[])*dv();
    div_pre = -(p[0,1,0]-p[0,-1,0])/(2.*Delta);
    grap_a_int1.y += div_pre*(1.0-f[])*dv();
    grap_w_int1.y += div_pre*(0.0+f[])*dv();
    div_pre = -(p[0,0,1]-p[0,0,-1])/(2.*Delta);
    grap_a_int1.z += div_pre*(1.0-f[])*dv();
    grap_w_int1.z += div_pre*(0.0+f[])*dv();
    double div_tau = 0;
    div_tau = (Sxx[1,0,0]-Sxx[-1,0,0]+Sxy[0,1,0]-Sxy[0,-1,0]+Sxz[0,0,1]-Sxz[0,0,-1])/(2.*Delta);
    divu_a_int1.x += mu2*div_tau*(1.0-f[])*dv();
    divu_w_int1.x += mu1*div_tau*(0.0+f[])*dv();
    div_tau = (Sxy[1,0,0]-Sxy[-1,0,0]+Syy[0,1,0]-Syy[0,-1,0]+Syz[0,0,1]-Syz[0,0,-1])/(2.*Delta);
    divu_a_int1.y += mu2*div_tau*(1.0-f[])*dv();
    divu_w_int1.y += mu1*div_tau*(0.0+f[])*dv();
    div_tau = (Sxz[1,0,0]-Sxz[-1,0,0]+Syz[0,1,0]-Syz[0,-1,0]+Szz[0,0,1]-Szz[0,0,-1])/(2.*Delta); 
    divu_a_int1.z += mu2*div_tau*(1.0-f[])*dv();
    divu_w_int1.z += mu1*div_tau*(0.0+f[])*dv();
  }

  // --> Compute the vorticity 2d
  scalar omega[];
  vorticity (u, omega);
  
  coord Fvor = {0.,0.,0.};
  foreach(reduction(+:Fvor)) {
    double lap_o = (omega[1,0,0]-2.0*omega[]+omega[-1,0,0])/sq(Delta) +
                   (omega[0,1,0]-2.0*omega[]+omega[0,-1,0])/sq(Delta) +
                   (omega[0,0,1]-2.0*omega[]+omega[0,0,-1])/sq(Delta);
    Fvor.x += -mu2*lap_o*y*(1.0-f[])*dv();
    Fvor.y += -mu2*lap_o*x*(1.0-f[])*dv();
  }

  // --> Compute the power vector
  vector t_vel[];
  foreach() {
    //
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx + dudx;
    double SDeformxy = dudy + dvdx;
    double SDeformxz = dudz + dwdx;
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy + dvdy;
    double SDeformyz = dvdz + dwdy;
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz + dwdz;
    t_vel.x[] = SDeformxx*u.x[] + SDeformxy*u.y[] + SDeformxz*u.z[]; 
    t_vel.y[] = SDeformyx*u.x[] + SDeformyy*u.y[] + SDeformyz*u.z[]; 
    t_vel.z[] = SDeformzx*u.x[] + SDeformzy*u.y[] + SDeformzz*u.z[]; 
  }

  // --> Pressure power integral, viscous power integral (per phase and for both kinetic and gravitational forces)
  double grap_pa_int = 0., grap_pw_int = 0.;
  double taup_pa_int = 0., taup_pw_int = 0.;
  foreach(reduction(+:grap_pa_int),reduction(+:grap_pw_int)
	  reduction(+:taup_pa_int),reduction(+:taup_pw_int)) {
    foreach_dimension() {
      double dpvel_dl = ( u.x[1,0,0]*p[1,0,0]-u.x[-1,0,0]*p[-1,0,0] )/(2.0*Delta);
      grap_pa_int += -dpvel_dl*(1.0-f[])*dv();
      grap_pw_int += -dpvel_dl*(0.0+f[])*dv();
      double dtvel_dl = ( t_vel.x[1,0,0]-t_vel.x[-1,0,0] )/(2.0*Delta);
      taup_pa_int += mu2*dtvel_dl*(1.0-f[])*dv();
      taup_pw_int += mu1*dtvel_dl*(0.0+f[])*dv();
    }
  }

  // --> tot_en
  scalar tot_en[]; 
  foreach() {
    double norm2 = 0.;
    foreach_dimension() {
      norm2 += sq(u.x[]);
    }
    tot_en[] = (norm2+g_*y)*dv();
  }

  // --> convective contribution
  double conv_ena = 0.0, conv_enw = 0.0; 
  foreach(reduction(+:conv_ena),reduction(+:conv_enw)) {
    foreach_dimension() {
      conv_ena += (u.x[1,0,0]*tot_en[1,0,0]-u.x[-1,0,0]*tot_en[-1,0,0])/(2.0*Delta)*(1.0-f[])*dv();
      conv_enw += (u.x[1,0,0]*tot_en[1,0,0]-u.x[-1,0,0]*tot_en[-1,0,0])/(2.0*Delta)*(0.0+f[])*dv();
    }
  }

  if (pid() == 0) { 

    char name[80];

    fflush(stderr);
    sprintf (name,"./budgets/vor_force.out");
    FILE * log_vor = fopen(name,"a");
    fprintf (log_vor, "%.10e %.10e %.10e %.10e\n", t, 1.0*i, 
		      Fvor.x, Fvor.y);
    fclose(log_vor);

    fflush(stderr);
    sprintf (name,"./budgets/en_bud_tot.out");
    FILE * log_en_bud = fopen(name,"a");
    fprintf (log_en_bud, "%.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, 
		         keTot, gpeTot, ste_t);
    fclose(log_en_bud);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_x_tot.out");
    FILE * log_mtx = fopen(name,"a");
    fprintf (log_mtx, "%.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, 
		      vel_m.x, mom_m.x, stm_t.x);
    fclose(log_mtx);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_y_tot.out");
    FILE * log_mty = fopen(name,"a");
    fprintf (log_mty, "%.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i,  
		      vel_m.y, mom_m.y, stm_t.y);
    fclose(log_mty);
    
    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_z_tot.out");
    FILE * log_mtz = fopen(name,"a");
    fprintf (log_mtz, "%.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, 
		      vel_m.z, mom_m.z, stm_t.z);
    fclose(log_mtz);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_x_air.out");
    FILE * log_max = fopen(name,"a");
    fprintf (log_max, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, vol_a_m, 
		      vel_a_m.x, grap_a_int.x, grap_a_int1.x, divu_a_int.x, divu_a_int1.x, divc_a_int1.x);
    fclose(log_max);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_y_air.out");
    FILE * log_may = fopen(name,"a");
    fprintf (log_may, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, vol_a_m, 
		      vel_a_m.y, grap_a_int.y, grap_a_int1.y, divu_a_int.y, divu_a_int1.y, divc_a_int1.y);
    fclose(log_may);
    
    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_z_air.out");
    FILE * log_maz = fopen(name,"a");
    fprintf (log_maz, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, vol_a_m, 
		      vel_a_m.z, grap_a_int.z, grap_a_int1.z, divu_a_int.z, divu_a_int1.z, divc_a_int1.z);
    fclose(log_maz);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_x_water.out");
    FILE * log_mwx = fopen(name,"a");
    fprintf (log_mwx, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, vol_w_m, 
		      vel_w_m.x, grap_w_int.x, grap_w_int1.x, divu_w_int.x, divu_w_int1.x, divc_w_int1.x);
    fclose(log_mwx);

    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_y_water.out");
    FILE * log_mwy = fopen(name,"a");
    fprintf (log_mwy, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, vol_w_m, 
		      vel_w_m.y, grap_w_int.y, grap_w_int1.y, divu_w_int.y, divu_w_int1.y, divc_w_int1.y);
    fclose(log_mwy);
    
    fflush(stderr);
    sprintf (name,"./budgets/mom_bud_z_water.out");
    FILE * log_mwz = fopen(name,"a");
    fprintf (log_mwz, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", t, 1.0*i, vol_w_m, 
		      vel_w_m.z, grap_w_int.z, grap_w_int1.z, divu_w_int.z, divu_w_int1.z, divc_w_int1.z);
    fclose(log_mwz);

    fflush(stderr);
    sprintf (name,"./budgets/ke_bud_air.out");
    FILE * fpair = fopen(name,"a");
    fprintf (fpair, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
             //t, 1.0*i, keAir/2., gpeAir + 0.125*L0*L0*L0, grap_pa_int, taup_pa_int, dissAir);
             t, 1.0*i, keAir/2., gpeAir, grap_pa_int, taup_pa_int, dissAir, conv_ena);
    fclose(fpair);

    fflush(stderr);
    sprintf (name,"./budgets/ke_bud_water.out");
    FILE * fpwater = fopen(name,"a");
    fprintf (fpwater, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
             //t, 1.0*i, keWat/2., gpeWat + 0.125*L0*L0*L0, grap_pw_int, taup_pw_int, dissWater);
             t, 1.0*i, keWat/2., gpeWat, grap_pw_int, taup_pw_int, dissWater, conv_enw); 
    fclose(fpwater);
    
  }

}

/**
   We want to compute some quantities at the interface. */

void output_int_qtn (char * fname, int istep, double time, scalar c, scalar p_a, double stp_eta) {

  // We first loop over all the interfacial points 
  // and we count them (per processor)
  
  int int_pt = 0;
  foreach(serial) { 
    if (interfacial (point, c)) {
      if (point.level == MAXLEVEL) {
        coord n = interface_normal(point, c), pp;
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
#if dimension > 2
        double zc = z + Delta*pp.z; // we keep here to avoid warning (otherwise: unused variable)
	Point point1 = locate (xc, yc, zc);
#else
	Point point1 = locate (xc, yc);
#endif
	if (point1.level > 0) {
	  POINT_VARIABLES;
	  int_pt++;
	}
      }
    }
  }

  //int tot_column = 21;
  int tot_column = 18;
  double t_mat[int_pt][tot_column];

  for (int j = 0; j < tot_column; j++) {
    for (int i = 0; i < int_pt; i++) {
      t_mat[i][j] = 0;
    }
  }
  fprintf(stderr, "First pass over interfacial points\n");

  scalar pos[];
#if dimension > 2
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#else
  coord G = {0.,1.}, Z = {0.,0.};
#endif
  position (c, pos, G, Z);

  scalar curv[];
  curvature (c, curv);

  // --> quantities to be probed at the interface
  // n.b.: for the stress we can use u since it contains derivatives
  scalar Sxx[]; scalar Syy[]; scalar Szz[];
  scalar Sxy[]; scalar Sxz[]; scalar Syz[];
  foreach() {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    Sxx[] = dudx+dudx; 
    Syy[] = dvdy+dvdy;
    Szz[] = dwdz+dwdz;
    Sxy[] = dudy+dvdx;
    Sxz[] = dudz+dwdx;
    Syz[] = dvdz+dwdy;
  }
  //boundary({Sxx,Syy,Szz,Sxy,Sxz,Syz}); // must be kept
  boundary({Sxx,Syy,Szz,Sxy,Sxz,Syz,p_a,u.x,u.y,u.z}); // must be kept

  int count = 0;
  foreach(serial) {
    if (interfacial (point, c)) {
      if (point.level == MAXLEVEL) {
        coord n = interface_normal(point, c), pp;
        double alpha1 = plane_alpha (c[], n);
        plane_area_center(n, alpha1, &pp);
	double eta = pos[];
	double cur = curv[];
	double norm_2 = sqrt(sq(n.x) + sq(n.y) + sq(n.z));
        double xc = x + Delta*pp.x;
        double yc = y + Delta*pp.y + stp_eta;
        double zc = z + Delta*pp.z;
#if dimension > 2
	Point point1 = locate (xc, yc, zc);
#else
	Point point1 = locate (xc, yc);
#endif
	if (point1.level > 0) {
	 
	  POINT_VARIABLES;

	  t_mat[count][0]  = xc;
#if dimension > 2
	  t_mat[count][1]  = zc;
#else
	  t_mat[count][1]  = 0;
#endif
	  t_mat[count][2]  = my_interpolation(point, p_a, xc, yc, zc);
	  t_mat[count][3]  = my_interpolation(point, Sxx, xc, yc, zc);
	  t_mat[count][4]  = my_interpolation(point, Syy, xc, yc, zc);
	  t_mat[count][5]  = my_interpolation(point, Szz, xc, yc, zc);
	  t_mat[count][6]  = my_interpolation(point, Sxy, xc, yc, zc);
	  t_mat[count][7]  = my_interpolation(point, Sxz, xc, yc, zc);
	  t_mat[count][8]  = my_interpolation(point, Syz, xc, yc, zc);
	  t_mat[count][9]  = my_interpolation(point, u.x, xc, yc, zc);
	  t_mat[count][10] = my_interpolation(point, u.y, xc, yc, zc);
#if dimension > 2
	  t_mat[count][11] = my_interpolation(point, u.z, xc, yc, zc);
#else
	  t_mat[count][11] = 0;
#endif

	  t_mat[count][12] = eta;           // already defined at the interface
	  t_mat[count][13] = cur;           // already defined at the interface
	  t_mat[count][14] = n.x/norm_2;    // already defined at the interface
	  t_mat[count][15] = n.y/norm_2;    // already defined at the interface
	  t_mat[count][16] = n.z/norm_2;    // already defined at the interface
	  t_mat[count][17] = Delta;         // no interpolation is meaningful
	  count++;

	}
      }
    }
  }
  fprintf(stderr, "Second pass over interfacial points\n");

  // We sort locally t_mat by the x coordinate (the first, i.e. 0-th, column) 

  qsort(t_mat, int_pt, tot_column*sizeof(double), compare);
  fprintf(stderr, "First sort\n");

#if _MPI

  // On multiple cores, we gather all the int_pt to the root pid 
  // and, then, we broadcast this information to all the processes
  
  int nproc;
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  int counts_it[nproc];
  if( pid() == 0 ) {
    MPI_Gather(&int_pt,1,MPI_INT,counts_it,1,MPI_INT,0,MPI_COMM_WORLD); // MPI_gather gathers by rank order
  }
  else {
    MPI_Gather(&int_pt,1,MPI_INT,NULL     ,1,MPI_INT,0,MPI_COMM_WORLD); // MPI_gather gathers by rank order
  }
  MPI_Bcast(counts_it,nproc,MPI_INT,0,MPI_COMM_WORLD);

  // Each processor knows the int_pt of the others. So we can compute 
  // the displacement, the total number of interfacial points and 
  // the number of elements owned by each processor

  int tot_int_pt = 0;
  int tot_el_p[nproc], disp_r[nproc], disp[nproc];
  for (int i = 0; i < nproc; i++) {
    disp_r[i]    = tot_int_pt;
    disp[i]      = disp_r[i]*tot_column;
    tot_int_pt  += counts_it[i];
    tot_el_p[i]  = counts_it[i]*tot_column;
  }

  // --> Gather to the root pid 
  // --> Sort by the first column

  double t_mat_tot[tot_int_pt][tot_column];

  if( pid() == 0 ) {
    int tot_el_p0[nproc], disp0[nproc];
    for (int i = 0; i < nproc; i++) {
      tot_el_p0[i] = tot_el_p[i];
      disp0[i]     = disp[i];
    }
    MPI_Gatherv(&t_mat,int_pt*tot_column,MPI_DOUBLE,t_mat_tot,tot_el_p0,disp0,MPI_DOUBLE,0,MPI_COMM_WORLD);
    fprintf(stderr, "GatherV 0\n");
    qsort(t_mat_tot, tot_int_pt, tot_column*sizeof(double), compare);
    fprintf(stderr, "Second sort 0\n");
  }
  else {
    int tot_el_p1[nproc], disp1[nproc];
    for (int i = 0; i < nproc; i++) {
      tot_el_p1[i] = tot_el_p[i];
      disp1[i]     = disp[i];
    }
    MPI_Gatherv(&t_mat,int_pt*tot_column,MPI_DOUBLE,NULL     ,tot_el_p1,disp1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
 
#else
 
  int tot_int_pt = int_pt;
  double t_mat_tot[tot_int_pt][tot_column];
  for (int j = 0; j < tot_column; j++) {
    for (int i = 0; i < int_pt; i++) {
      t_mat_tot[i][j] = t_mat[i][j];
    }
  }

#endif

  if (pid() == 0) {
    
    fflush(stderr);

    // Print in binary format
    FILE * fpver = fopen (fname,"w");
    for (int i = 0; i < tot_int_pt; i++) {
      for (int j = 0; j < tot_column; j++) {
	fwrite ( &t_mat_tot[i][j], sizeof(double), 1, fpver );
      }
    }
    fclose(fpver);
    fflush(fpver);

    // Log some info about the binary
    char name_1[100];
    sprintf (name_1, "./eta/global_int.out");
    FILE * global_int = fopen (name_1, "a");
    fprintf (global_int, "%.10e %09d %09d %09d\n", 
		          time, istep, tot_int_pt, tot_column);
    fclose(global_int);

  }
  fprintf(stderr, "Print\n");

}

//event tagging (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_tag_my) {
event tagging (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_tag_my) {
//event tagging (i += 2) {

  /*
  scalar f2[];
  scalar_clone (f2, f);
  foreach() {
    f2[] = f[];
  }
#if TREE
  f2.refine = f2.prolongation = fraction_refine;
  restriction({f2});
#endif

  if(do_remove) {
    remove_droplets (f2,min_size,threshold,false); // first remove droplets
    remove_droplets (f2,min_size,threshold,true);  // then remove bubbles
  }
  */

  fprintf(stderr, "I output the tagging every t=%.10e\n", T0_/tout_tag_my);

  if(do_tagging == 1) {

    char name_1[100], name_2[100], name_3[100];
    
    fprintf(stderr, "I do f1-tagging every t=%.10e\n", T0_/tout_cut_my);
    sprintf (name_1, "./tagging/num_bubbles_droplets_f1.out");
    sprintf (name_2, "./tagging/tagging_f1/droplets_f1_%09d.out", i);
    sprintf (name_3, "./tagging/tagging_f1/bubbles_f1_%09d.out", i);
    countDropsBubble (name_1, name_2, name_3, i, f);
 
    /*   
    fprintf(stderr, "I do f2-tagging every t=%.10e\n", T0_/tout_cut_my);
    sprintf (name_1, "./tagging/num_bubbles_droplets_f2.out");
    sprintf (name_2, "./tagging/tagging_f2/droplets_f2_%09d.out", i);
    sprintf (name_3, "./tagging/tagging_f2/bubbles_f2_%09d.out", i);
    countDropsBubble (name_1, name_2, name_3, i, f2); 
    */

  }
  //return 1;
  
}

event global_obs (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_glo_my) {
//event global_obs (i += 2) {

  scalar f2[];
  scalar_clone (f2, f);
  foreach() {
    f2[] = f[];
  }
#if TREE
  f2.refine = f2.prolongation = fraction_refine;
  restriction({f2});
#endif

  if(do_remove) {
    remove_droplets (f2,min_size,threshold,false); // first remove droplets
    remove_droplets (f2,min_size,threshold,true);  // then remove bubbles
  }

  fprintf(stderr, "I output the global observables file every t=%.10e\n", T0_/tout_glo_my);

  char name_0[100];

  // Note: 
  // --> we need a constant offset to ensure we can preserve 
  //     the same interface geometry
  // --> we need to a have a stp_eta = 0/1/2/4*Delta to ensure 
  //     the velocity gradients for the viscous contributions to the stress 
  //     are properly computed in the air side only

  double stp_eta_0 = 0.0*(L0/pow(2.0,MAXLEVEL)); // it corresponds to 0*Delta 
  fprintf(stderr, "stp_eta0 for gl. obs. = %.10e\n", stp_eta_0);
  sprintf (name_0, "./budgets/global_obs_ptot_0.out");
  output_global_obs_1 (name_0, i, f2, p  , stp_eta_0);

  double stp_eta_1 = 1.0*(L0/pow(2.0,MAXLEVEL)); // it corresponds to 1*Delta 
  fprintf(stderr, "stp_eta1 for gl. obs. = %.10e\n", stp_eta_1);
  sprintf (name_0, "./budgets/global_obs_ptot_1.out");
  output_global_obs_1 (name_0, i, f2, p  , stp_eta_1);

  double stp_eta_2 = 2.0*(L0/pow(2.0,MAXLEVEL)); // it corresponds to 2*Delta
  fprintf(stderr, "stp_eta2 for gl. obs. = %.10e\n", stp_eta_2);
  sprintf (name_0, "./budgets/global_obs_ptot_2.out");
  output_global_obs_1 (name_0, i, f2, p  , stp_eta_2);
  
  double stp_eta_3 = 4.0*(L0/pow(2.0,MAXLEVEL)); // it corresponds to 4*Delta
  fprintf(stderr, "stp_eta3 for gl. obs. = %.10e\n", stp_eta_3);
  sprintf (name_0, "./budgets/global_obs_ptot_3.out");
  output_global_obs_1 (name_0, i, f2, p  , stp_eta_3);

  /*
  // We shift f2 of my_stp_eta (we do at the end to avoid to pollute something)
  face vector uuf[]; coord ufv = {0, 1, 0};
  foreach_face()
    uuf.x[] = ufv.x;

  double dt2 = 0.25*CFL*L0/(1 << grid->maxdepth);
  double yp = 0, my_stp_eta = 4.0*(L0/pow(2.0,MAXLEVEL)); // it corresponds to 4*Delta
  int ps_i = 0;
  while (yp < my_stp_eta) {
    vof_advection2 (f2, dt2, uuf, i);
    yp += dt2*ufv.y;
    ps_i++;
    fprintf(stderr, "translation of the VoF along y%d\n", ps_i);
  }

  double stp_eta_0s = 0.0; // it corresponds to 0*Delta
  fprintf(stderr, "stp_eta0s for gl. obs. = %.10e\n", stp_eta_0s);
  sprintf (name_0, "./budgets/global_obs_ptot_0s.out");
  output_global_obs_1 (name_0, i, f2, p  , stp_eta_0s);
  */

}

event eta_loc (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_eta_my) {
//event eta_loc (i += 2) {

  if (do_eta_loc == 1) {

    //fprintf(stderr, "I output the eta_loc file every t=%.10e\n", T0_/tout_eta);
    fprintf(stderr, "I output the eta_loc file every t=%.10e\n", T0_/tout_eta_my);
   
    /*
    scalar f2[];
    scalar_clone (f2, f);
    foreach() {
      f2[] = f[];
    }
    f2.refine = f2.prolongation = fraction_refine;
    restriction({f2});
  
    if(do_remove) {
      //remove_droplets_all (f2,threshold,false); // first remove droplets
      //remove_droplets_all (f2,threshold,true);  // then remove bubbles
      remove_droplets (f2,min_size,threshold,false); // first remove droplets
      remove_droplets (f2,min_size,threshold,true);  // then remove bubbles
    }

    face vector uuf[]; coord ufv = {0, 1, 0};
    foreach_face()
      uuf.x[] = ufv.x;
    
    scalar f2[];
    scalar_clone (f2, f);
    foreach() {
      f2[] = f[];
    }
#if TREE
    f2.refine = f2.prolongation = fraction_refine;
    restriction({f2});
#endif

    double dt2 = 0.25*CFL*L0/(1 << grid->maxdepth);
    double yp = 0, my_stp_eta = 4.0*(L0/pow(2.0,MAXLEVEL)); // it corresponds to 4*Delta
    int ps_i = 0;
    while (yp < my_stp_eta) {
      vof_advection2 (f2, dt2, uuf, i);
      yp += dt2*ufv.y;
      ps_i++;
      fprintf(stderr, "translation of the VoF along y%d\n", ps_i);
    }
    */

    fflush(stderr);
    char eta_out[100];
    sprintf (eta_out, "./eta/eta_loc/eta_loc_t%09d.bin", i);
    //sprintf (eta_out, "./eta/eta_loc/eta_loc_t%.10e.bin", t);
    double stp_eta = 4.0*(L0/pow(2.0,MAXLEVEL)); // it corresponds to 4*Delta
    output_int_qtn (eta_out, i, t, f, p, stp_eta);
 
  }

}

/** 
   Dump every tout_res period */

event dumpstep (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_res_my) {
//event dumpstep (i += 2) {

  //fprintf(stderr, "I output the restarting files every t=%.10e\n", T0_/tout_res), fflush (stderr);
  fprintf(stderr, "I output the restarting files every t=%.10e\n", T0_/tout_res_my), fflush (stderr);
  
  if(counter < counter_max) {
    counter++;
  }
  else {
    counter = 1;
  }

  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  sprintf (dname, "dump_%d.bin", counter);
  dump (dname);
  //dump (dname, {f,u.x,u.y,u.z,g.x,g.y,g.z,rhov});

  /** 
     Add a symbolic link, log restarting info and size of the bin */

  if (pid () == 0) {

    char comm_2[80];
    sprintf (comm_2, "ln -sf dump_%d.bin restart.bin", counter);
    system(comm_2);

    fflush(stderr);
    char name_0[80];
    sprintf (name_0, "dump_%d.bin", counter);
    FILE * fp = fopen(name_0, "r");
    fseek(fp, 0L, SEEK_END);
    long int res = ftell(fp);
    fclose(fp);

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_restart.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %.10e\n", t, 1.0*i, 1.0*counter, 1.0*res);
    fclose(log_sim);

  }

}

event dump_backup (t = RELEASETIME; t <= T0_*end_sim; t += T0_/tout_rbk_my) {
//event dump_backup (i += 2) {

  fprintf(stderr, "I output the restarting files for backup every t=%.10e\n", T0_/tout_rbk_my), fflush (stderr);

  char dname[100];
  sprintf (dname, "./restart_bk/dump_%09d.bin", i);
  dump (dname);

  /**
     Log restarting info and size of the bin */

  if (pid () == 0) {

    fflush(stderr);
    FILE * fp = fopen(dname, "r");
    fseek(fp, 0L, SEEK_END);
    long int res = ftell(fp);
    fclose(fp);

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"./restart_bk/log_restart_bk.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %.10e\n", t, 1.0*i, t-RELEASETIME, 1.0*res);
    fclose(log_sim);

  }

}

/**
## End 

   We want to run up end_sim physical time. */

event end (t = end_sim*T0_) {

  dump ("end.bin");

  if ( pid() == 0 ) {

    char comm[80];
    sprintf (comm, "ln -sf end.bin restart.bin");
    system(comm);

  }

  return 1; // exit

}

/** 
   Adaptive function. We need a more strict criteria for water speed once the waves starts moving, i.e. t >= RELEASETIME. */ 
   
#if TREE
event adapt (i++) {

  if (t < RELEASETIME) {
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
  }

  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	u_water.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,uwemax,uwemax,uwemax}, MAXLEVEL, MINLEVEL);
  }

}
#endif
