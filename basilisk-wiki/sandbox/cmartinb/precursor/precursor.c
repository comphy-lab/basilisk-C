/** 
   This is the file to generate a precursor simulation for the two-phase TBL. */

#include "embed.h"
#include "navier-stokes/centered.h"   
#include "view.h"
#include "sandbox/perfs.h"
#include "sandbox/profile6.h"   // From Antoon
#include "sandbox/maxruntime.h"   
#include "sandbox/my_funcs.h"   

/**
   The input parameters are listed below. 
   Note that these input parameters are preinitialized to some values,
   which will be replaced with the ones passed by the command line. */

double Re_ast = 720.0;                 // Friction Reynolds number
double ak = 0.3;                       // initial steepness (relevant only for a Stokes wave)
double r_L0lam = 4.0;                  // ratio between the box size L0 and the (peak) wavelength
double rho_r = 1.225/1000.0;           // reference density (air)
double mu_r  = 2.2471881948940954e-06; // reference dynamic viscosity (air)
int MAXLEVEL = 9;                      // max level of the simulation
int MINLEVEL = 6;                      // min level of the simulation
double f_uemax = 0.3;                  // fraction of U_ast to define the refinement criteria for the velocity
double Tf_ = 3.14;                     // input physical time used to establish the dumping frequency
int do_en = 0;                         // dump or not different fields for ensemble simulations
int st_wave = 1;                       // initialize or not a Stokes wave
int MY_RAND = 2;                       // random number for a narrowbanded wave field (relevant if st_wave = 0)
int dump_now = 0;                      // dump a restarting file now 
int N_mod = 64;                        // number of mode for the initial spectrum
double theta_d = 0.;                   // misaligned angle (in degree)

/**
   We define some values: 
   --> the thermophyiscal properties in the airflow
   --> the friction velocity
   --> the water depth
   --> the (peak) wavelength, the (peak) wave number */

double mu2 = 1.0, rho2 = 1.0; // we change later
double U_ast = 1.0; // we change later
double h_  = 1.0; // we change later 
double lam = 1.0; // we change later
double k_  = 1.0; // we change later

/**
   Variables for the refinement criteria. */

double femax = 1.0e-4; // we change later 
double uemax = 0.0750; // we change later

/**
   For the restarting/ensemble step. */

int counter_max = 2;
int counter = 0;
int counter_max_ens = 20;
int counter_ens = 0;

/**
   We need to store some variables. */

face vector muv[], av[], alphav[];
scalar rhov[];

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
  return eta1 + ak*eta2 + sq(ak)*eta3 + h_;

}

/** 
   Add physical perturbations to promote faster transition to turbulence. */

double fy (double y) {
  double fy_f = (sq(1.0-sq(y)));
  return fy_f;
}

double dfy (double y) {
  double dfy_f = -4.0*y*(sq(1.0-sq(y)));
  return dfy_f;
}

double gxz (double x, double z) {
  double gxz_f = z*exp(-4.*(4.*sq(x)+sq(z)));
  return gxz_f;
}

double dgxz (double x, double z) {
  double dgxz_f = exp(-4.*(4.*sq(x)+sq(z)))*(1.-8.*sq(z));
  return dgxz_f;
}

/** 
   Set the wave field
   --> set_profile: 3rd order Stokes Wave;
   --> import_profile: broadbanded wave field generated in spectrum.h. */

void set_profile (scalar cs, face vector fs, vertex scalar phi) {
  foreach_vertex() {
    phi[] = y-WaveProfile(x,z);
  } 
  fractions (phi,cs,fs);
  fractions_cleanup (cs,fs); // remove inconsistencies
}

#include "sandbox/spectrum.h" // Used for new input method (the spectrum info)
void import_profile (int N_mod, scalar cs, face vector fs, vertex scalar phi) {
  power_input(N_mod);
  dkx_ = kx_[1] - kx_[0];
  dky_ = ky_[1] - ky_[0];
  fprintf (stderr, "dkx = %g, dky = %g\n", dkx_, dky_), fflush (stderr);
  foreach_vertex() {
    phi[] = y-(wave(x,z,N_mod)+h_);
  }
  fractions (phi,cs,fs);
  fractions_cleanup (cs,fs); // remove inconsistencies
}


/**
   Input parameters are passed in from the command line. */

int main(int argc, char *argv[]) {

  /**
     Provide the inputs */

  maxruntime (&argc, argv);

  if (argc > 1)
    Re_ast = atof(argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    r_L0lam = atof(argv[3]);
  if (argc > 4)
    rho_r = atof(argv[4]);
  if (argc > 5)
    mu_r = atof(argv[5]);
  if (argc > 6)
    MAXLEVEL = atoi(argv[6]);
  if (argc > 7)
    MINLEVEL = atoi(argv[7]);
  if (argc > 8)
    f_uemax = atof(argv[8]);
  if (argc > 9)
    Tf_ = atof(argv[9]);
  if (argc > 10)
    do_en = atoi(argv[10]);
  if (argc > 11)
    st_wave = atoi(argv[11]);
  if (argc > 12)
    MY_RAND = atoi(argv[12]);
  if (argc > 13)
    dump_now = atoi(argv[13]);
  if (argc > 14)
    N_mod = atoi(argv[14]);
  if (argc > 15)
    theta_d = atoi(argv[15]);
 
  if (argc < 16) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 16-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "maximum runtime = %8E seconds\n", _maxruntime), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " Re_ast = %8E\n ", Re_ast), fflush (stderr);
  fprintf(stderr, " a_0k = %8E\n ", ak), fflush (stderr);
  fprintf(stderr, " Ratio between the box size and the wavelength = %8E\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAXLEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MINLEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " Fraction of uemax = %8E\n ", f_uemax), fflush (stderr);
  fprintf(stderr, " Input physical time to output = %8E\n ", Tf_), fflush (stderr);
  fprintf(stderr, " do_en  = %d\n ", do_en), fflush (stderr);
  fprintf(stderr, " st_wave = %d\n ", st_wave), fflush (stderr);
  fprintf(stderr, " MY_RAND = %d\n ", MY_RAND), fflush (stderr);
  fprintf(stderr, " dump_now = %d\n ", dump_now), fflush (stderr);
  fprintf(stderr, " N_mode = %d\n ", N_mod), fflush (stderr);
  fprintf(stderr, " Misaligned angle (in degree) = %8E\n ", theta_d), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     The domain is a cubic box centered on the origin and of length
     $L_0=2*\pi$, periodic in the x- and z-directions. 
     Once $L_0$ is set, we automatically know 
     (1) the water depth, (2) the (peak) wavelength and (3) the (peak) wave number. */

  L0  = 2.0*pi;
  h_  = L0/(2.0*pi); // Water depth set following Wu et al. JFM2022
  lam = L0/r_L0lam;
  k_  = 2.0*pi/lam;

  origin (-L0/2., 0, -L0/2.);
  init_grid (1 << (MAXLEVEL - 2)); // we set a reasonable grid spacing just for initialization

  /**
     Set the boundary conditions */

  periodic (right);
  periodic (front);

  u.r[top] = neumann(0.);
  u.n[top] = dirichlet(0.);  
  u.t[top] = neumann(0.);
  p[top]   = neumann(0.);

  u.r[bottom] = neumann(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = neumann(0.);
  p[bottom]   = neumann(0.);

  uf.r[top] = neumann(0.);
  uf.n[top] = dirichlet(0.);  
  uf.t[top] = neumann(0.);
  pf[top]   = neumann(0.);

  uf.r[bottom] = neumann(0.);
  uf.n[bottom] = dirichlet(0.);
  uf.t[bottom] = neumann(0.);
  pf[bottom]   = neumann(0.);

  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);
  u.r[embed] = dirichlet(0.);
  p[embed]   = neumann(0.);

  uf.n[embed] = dirichlet(0.);
  uf.t[embed] = dirichlet(0.);
  uf.r[embed] = dirichlet(0.);
  pf[embed]   = neumann(0.);
 
  /**
     We set the thermophysical properties in the airflow 
     Note: this should match the ones employed for the two-phase simulations */

  rho2 = rho_r;
  mu2  = mu_r;
  
  /**
     We compute the U_ast value so that it matches the desired
     (1) Re_ast, (2) thermophysical properties and (3) size of the box. */

  U_ast = (mu2*Re_ast)/(rho2*(L0-h_));
  
  /**
     We reduce the tolerance of the solver */

  //TOLERANCE = 1e-6;
  /*
  if(st_wave != 1) {
    TOLERANCE = 1e-4;
    NITERMIN  = 2;
    //TOLERANCE = HUGE;
    //NITERMAX = 50;
  }
  */

  /**
     Give the address of some variables */

  a     = av;
  alpha = alphav;
  rho   = rhov;
  mu    = muv;

  /**
     Set the refinement criteria */

  uemax = f_uemax*U_ast;
  femax = 1.0e-4; // hard coded (works for most of the configurations)
  
  /**
     Run! */

  run();

}

/**
   User defined-function to append .ppm files. */

# define POPEN(name, mode) fopen (name ".ppm", mode)

/**
   Initial event. */

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " Re_ast = %8E\n ", rho2*U_ast*(L0-h_)/mu2), fflush (stderr);
  fprintf(stderr, " a_0k = %8E\n ", ak), fflush (stderr);
  fprintf(stderr, " Ratio between the box size and the wavelength = %8E\n ", r_L0lam), fflush (stderr);
  fprintf(stderr, " Reference density = %8E\n ", rho_r), fflush (stderr);
  fprintf(stderr, " Reference dynamic viscosity = %8E\n ", mu_r), fflush (stderr);
  fprintf(stderr, " MAXLEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MINLEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " Fraction of uemax  = %8E\n ", f_uemax), fflush (stderr);
  fprintf(stderr, " Input physical time to output = %8E\n ", Tf_), fflush (stderr);
  fprintf(stderr, " do_en  = %d\n ", do_en), fflush (stderr);
  fprintf(stderr, " st_wave = %d\n ", st_wave), fflush (stderr);
  fprintf(stderr, " MY_RAND = %d\n ", MY_RAND), fflush (stderr);
  fprintf(stderr, " dump_now = %d\n ", dump_now), fflush (stderr);
  fprintf(stderr, " N_mode = %d\n ", N_mod), fflush (stderr);
  fprintf(stderr, " Misaligned angle (in degree) = %8E\n ", theta_d), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, " Friction velocity = %8E\n ", U_ast), fflush (stderr);
  fprintf(stderr, " Re_ast based on the wavelength = %8E\n ", rho2*U_ast*lam/mu2), fflush (stderr);
  fprintf(stderr, " Refinement criteria for the velocity = %8E\n ", uemax), fflush (stderr);
  fprintf(stderr, " Refinement criteria for the cell/volume fraction = %8E\n ", femax), fflush (stderr);
  fprintf(stderr, " Water depth = %8E\n ", h_), fflush (stderr);
  fprintf(stderr, " (Peak) wavelength = %8E\n ", lam), fflush (stderr);
  fprintf(stderr, " (Peak) wave number = %8E\n ", k_), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     Create directories for better organization of the output files */

  fprintf(stderr, "Create directories (if needed) \n"), fflush (stderr);

  char comm[80];
  sprintf (comm, "mkdir -p profiles");
  system(comm);
  sprintf (comm, "mkdir -p restart_bk");
  system(comm);
  sprintf (comm, "mkdir -p restart_ens");
  system(comm);
  /*
  sprintf (comm, "mkdir -p phase");
  system(comm);
  */
  
  /**
     Restart from the latest dump files or initialize a new simulation */

  if (restore ("restart.bin")) {
    fprintf(stderr, "Simulation restarts from a dumped file\n"), fflush (stderr);
    for (scalar s in {u, g}) { 
      s.prolongation = refine_embed_linear;
    }
    ///*
    cs.prolongation = refine_injection;
    adapt_wavelet ({cs,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
    cs.prolongation = fraction_refine;
    //*/
    vertex scalar phi[]; // must be here after restore!
    if(st_wave == 1) {
      fprintf(stderr, "We use a 3rd-order Stokes wave!\n"), fflush (stderr);
      set_profile (cs,fs,phi);
    }
    else {
      fprintf(stderr, "We import a profile previously generated!\n"), fflush (stderr);
      allocate_arrays (N_mod);
      import_profile (N_mod,cs,fs,phi);
      free_arrays ();
    }
  }
  else {
    fprintf(stderr, "Initialization of all the variables\n"), fflush (stderr);
    refine (fabs(y-h_) < 0.2 && level < MAXLEVEL);
    vertex scalar phi[]; // must be here after the initial grid is set
    if(st_wave == 1) {
      fprintf(stderr, "We use a 3rd-order Stokes wave!\n"), fflush (stderr);
      set_profile (cs,fs,phi);
    }
    else {
      fprintf(stderr, "We import a profile previously generated!\n"), fflush (stderr);
      allocate_arrays (N_mod);
      import_profile (N_mod,cs,fs,phi);
      free_arrays ();
    }
    double y_ast = (L0-h_)/Re_ast;
    double Ubulk_ex = (mu2/(rho2*(L0-h_)))*(pow(0.5*Re_ast/0.09,(1.0/0.88))); // estimation based on Pope's relation 
    fprintf(stderr, "Ubulk_ex = %8E\n", Ubulk_ex), fflush (stderr);
    do {
      foreach() {
	if (phi[] > 0.05) {
	  u.x[] = cs[]*(U_ast/0.41)*(log(phi[]/y_ast));
          double x_n = 2.0*(x-0.0*L0)/L0;
          double y_n = 2.0*(y-0.0*L0)/L0-1.2;
          double z_n = 2.0*(z-0.0*L0)/L0;
          u.y[] = cs[]*(-1.0*gxz(z_n,x_n)*dfy(y_n)*Ubulk_ex*1.5);
          u.z[] = cs[]*(+1.0*fy(y_n)*dgxz(z_n,x_n)*Ubulk_ex*1.5);
	}
	else {
	  foreach_dimension() {
	    u.x[] = 0.0;
	  }
	}
      }
    }
#if TREE
    while (0);
#else
    while (0);
#endif
  }

  boundary({u.x,u.y,cs});
  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("cs");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.y", linear = true, n = {1,0,0}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  {
    static FILE * fp = POPEN ("3D_0", "a");
    save (fp = fp);
    fflush (fp);
  }

  clear();
  view (fov = 22.5, camera = "front", ty = -0.50,
        //theta = 0.0, phi = -pi/240.0, psi = 0.0, ty = -0.5,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("cs");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  {
    static FILE * fp = POPEN ("inter_0", "a");
    save (fp = fp);
    fflush (fp);
  }

  stats stat1 = statsf(cs); // log some stats
  double rms_cs = normf(cs).rms;

  if (pid() == 0) {

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"./log_stat_cs.out");
    FILE * log_c1 = fopen(name_1,"a");
    fprintf (log_c1, "%8E %8E %8E %8E %8E %8E\n", 
      	             t, 1.0*i, stat1.min, stat1.sum, stat1.max, rms_cs);
    fclose(log_c1);

  }

  if (dump_now == 1) {
    char dname[100];
    sprintf (dname, "restart.bin");
    dump (dname);
    return 1;
  }

}

/**
   We log the physical time, the number of time-steps and dt of the simulation. */

event log_simulation (i += 10) {
  
  int min_level = +100, max_level = -100;
  foreach(reduction(min:min_level) reduction(max:max_level)) {
    max_level = max(max_level,level);
    min_level = min(min_level,level);
  }

  if (pid() == 0) {

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_simulation.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %d %d\n", t, 1.0*i, dt, max_level, min_level);
    fclose(log_sim);

  }
 
  coord ua_m  = {0.,0.,0.};
  double vg_m = 0.;
  foreach(reduction(+:ua_m) reduction(+:vg_m)) {
    vg_m += (0.+cs[])*dv();
    foreach_dimension() {
      ua_m.x += u.x[]*(0.+cs[])*dv();
    }
  }

  coord Fp  = {0.,0.,0.};
  coord Fmu = {0.,0.,0.};
  embed_force_3d (p, u, mu, &Fp, &Fmu);

  double area = interface_area(cs);

  if (pid() == 0) {

    char name_1[80];
   
    fflush(stderr);
    sprintf (name_1,"log_forcing.out");
    FILE * log_sim1 = fopen(name_1,"a");
    fprintf (log_sim1, "%8E %8E %8E %8E %8E %8E\n", 
		        t, 1.0*i, vg_m,
			ua_m.x, ua_m.y, ua_m.z);
    fclose(log_sim1);

    fflush(stderr);
    sprintf (name_1,"log_force.out");
    FILE * log_sim2 = fopen(name_1,"a");
    fprintf (log_sim2, "%8E %8E %8E %8E %8E %8E %8E %8E %8E\n", 
		        t, 1.0*i, area,
		        Fp.x, Fp.y, Fp.z,
		        Fmu.x, Fmu.y, Fmu.z);
    fclose(log_sim2);

  }

}

void profile_output (int istep) {

  /**
     Compute the vorticity. */

  scalar omega[];
  vorticity (u, omega);

  char file[99];
  sprintf (file, "./profiles/prof_%09d.out", istep);
  scalar uxuy[],uxux[],uyuy[],uzuz[];
  foreach () {
    uxuy[] = u.x[]*u.y[];
    uxux[] = u.x[]*u.x[];
    uyuy[] = u.y[]*u.y[];
    uzuz[] = u.z[]*u.z[];
  }

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = y;
  }
  profiles ({u.x, u.y, u.z, cs, omega, uxuy, uxux, uyuy, uzuz, p}, phi, rf = 1.0, fname = file, n = 1 << MAXLEVEL);

}

event out_pro (t += 2.0*Tf_) {

  profile_output (i);
   
  if (pid() == 0) {
    char name_1[80];
    sprintf (name_1,"./profiles/log_pro.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
    fclose(log_sim);
  }

}

/**
   We update the density and the dynamic viscosity to match the desire Re_ast. */

event properties (i++) {
  foreach() {
    rhov[] = cm[]*rho2;
  }
  foreach_face() {
    alphav.x[] = fm.x[]/rho2;
    muv.x[]    = fm.x[]*mu2;
  }
}

/**
   We prescribe a uniform pressure gradient to drive the flow. */

event acceleration (i++) {

  double ampl    = sq(U_ast)/(L0-h_);
  double theta_r = (pi/180.)*theta_d; // convert into radiant
  
  coord th = {theta_r,0.,pi/2.-theta_r}; // th.x + th.z = pi/2.
  coord ev = {1.,0.,1.}; // forcing can be along x (streamwise) and/or z (spanwise), only
  foreach_face() {
    av.x[] += ampl*fm.x[]*fs.x[]*ev.x*cos(th.x);
  }

}

/** 
   Output video and field in uncompressed format 
   (we can compress later using convert .ppm to .mpg and open with mplayer) */

event movies (t += 5.0*Tf_) {
  
  boundary({u.x,u.y,u.z,cs});
  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("cs", color = "u.x");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.y", linear = true, n = {0,0,1}, alpha = 0.0);
  squares ("u.z", linear = true, n = {1,0,0}, alpha = 0.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
    fflush (fp);
  }

}

/** 
   Dump every tout_res period */

event dumpstep (t += 2.0*Tf_) {

  if(counter < counter_max) {
    counter++;
  }
  else {
    counter = 1;
  }

  char dname[100];
  sprintf (dname, "dump_%d.bin", counter);
  dump (dname);

  if (pid() == 0) {
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

  if(do_en == 1) {

    if(counter_ens < counter_max_ens) {
      counter_ens++;
    }
    else {
      return 1;
    }
    
    coord ua_m  = {0.,0.,0.};
    double vg_m = 0.;
    foreach(reduction(+:ua_m) reduction(+:vg_m)) {
      vg_m += (0.+cs[])*dv();
      foreach_dimension() {
        ua_m.x += u.x[]*(0.+cs[])*dv();
      }
    }

    if (pid() == 0) {
    
      fflush(stderr);
      char name_1[80];
      sprintf (name_1,"./restart_ens/log_ens_forcing.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%8E %8E %8E %8E %8E %8E %02d\n", 
          	         t, 1.0*i, vg_m,
          		 ua_m.x, ua_m.y, ua_m.z, counter_ens);
      fclose(log_sim);

    }
    
    sprintf (dname, "./restart_ens/dump_s%02d_ak%.2e_t%.10e.bin", counter_ens, ak, t);
    dump (dname);

  }

}

event dump_backup (t += 10.0*Tf_) {
 
  char dname[100];
  sprintf (dname, "./restart_bk/dump_%09d.bin", i);
  dump (dname);

}

/**
## End 

   We want to run up end_sim physical time. */

event end (t = 100000.0) {

  dump ("end.bin");

  if ( pid() == 0 ) {

    char comm[80];
    sprintf (comm, "ln -sf end.bin restart.bin");
    system(comm);

  }

  return 1; // exit

}

/** 
   Adaptive function */ 
   
#if TREE
event adapt (i++) {
  cs.prolongation = refine_injection;
  adapt_wavelet ({cs,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
  cs.prolongation = fraction_refine;
}
#endif