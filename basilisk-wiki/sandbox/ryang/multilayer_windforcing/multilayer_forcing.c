/**
# Field scale wave breaking (multilayer solver)
*/

/*   We test a phase correlated pressure forcing. */

// Modification of the face field to include pressure forcing
// Following the practice of Clement Robert for the surface tension term
// First we consider the most simple case where the pressure is pi/2 shifted
// and can be easily written as the _local_ gradient of eta
// the magnitude is defined through phat
// Later we'll replace p_surface with an external fft library
// NEED TO BE BEFORE hydro.h and hydro-tension.h
//

#define p_surface(eta, i) (- phat*(eta[i+1] - eta[i-1])/(2*Delta))
#define a_baro(eta, i) (G*gmetric(i)*(eta[i-1] - eta[i])/Delta)
#define a_baro_p(eta, i) (gmetric(i)*(p_surface (eta, i) - p_surface (eta, i - 1))/Delta)

double phat; // pressure magnitude

#include "grid/multigrid.h"
#include "hydro-tension.h" // For future extension of surface tension
#include "view.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "input.h" // Used for older input method (matrices of quantities)
#include "output_mpi.h" // Antoon's function for MPI compatible matrix output

/**
  Definition of input parameters. */ 

int NLAYER;     // number of layers
int LEVEL_data; // horizontal resolution
double TEND;    // end of the simulation
double nu;      // Viscosity for vertical diffusion
double L0;      // box size
double kp_;     // peak wavelength
double theta_H; // 
double g_;      // gravity
double Bo;      // Bond number
double eta_std; // root-mean-square of eta

/**
  Definition of some controlling parameters (we change later) */ 

double rho_w     = 1.0; // we keep fix
double sigma_val = 1.0; // surface tension
double h_        = 1.0; // depth of the water
double gpe_base  = 0.0; // gauge of potential energy
double T0_       = 1.0; // wave period

/**
  Define variables for the stress. */

double du0;
vector duv[];
scalar T;
#include "./spectrum_64.h" // Used for new input method (the spectrum info)

/** Function for writing fields at time t. */

int writefields (int i, const char *suffix) {
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50], filename5[50];
  vector u_temp;
  scalar w_temp, h_temp, T_temp;
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_%s_t%09d_l%02d", suffix, i, j);
    sprintf (filename2, "field/uy_%s_t%09d_l%02d", suffix, i, j);  
    sprintf (filename3, "field/uz_%s_t%09d_l%02d", suffix, i, j);  
    sprintf (filename4, "field/hl_%s_t%09d_l%02d", suffix, i, j); 
    sprintf (filename5, "field/T_%s_t%09d_l%02d", suffix, i, j); 
    if (j==0) {
      // The first layer is named u instead of u0
      sprintf (s, "u");
      u_temp = lookup_vector (s);
      sprintf (s, "w");
      w_temp = lookup_field (s);
      sprintf (s, "h");
      h_temp = lookup_field (s);
      sprintf (s, "T");
      T_temp = lookup_field (s);
    }
    else {
      sprintf (s, "u%d", j);
      u_temp = lookup_vector (s);
      sprintf (s, "w%d", j);
      w_temp = lookup_field (s);
      sprintf (s, "h%d", j);
      h_temp = lookup_field (s);
      sprintf (s, "T%d", j);
      T_temp = lookup_field (s);
    }
    FILE * fux = fopen (filename1, "w");
    output_matrix_mpi (u_temp.x, fux, N, linear = true);
    fclose (fux);
    FILE * fuy = fopen (filename2, "w");
    output_matrix_mpi (u_temp.y, fuy, N, linear = true);
    fclose (fuy);
    FILE * fuz = fopen (filename3, "w");
    output_matrix_mpi (w_temp, fuz, N, linear = true);
    fclose (fuz);
    FILE * fh = fopen (filename4, "w");
    output_matrix_mpi (h_temp, fh, N, linear = true);
    fclose (fh);
    FILE * fT = fopen (filename5, "w");
    output_matrix_mpi (T_temp, fT, N, linear = true);
    fclose (fT);    
  }
  return 0;

}

int main(int argc, char * argv[]) {

  /**
     Provide the inputs */

  if (argc > 1)
    NLAYER = atoi(argv[1]); // # of layers
  if (argc > 2)
    LEVEL_data = atoi(argv[2]); // Horizontal resolution
  if (argc > 3)
    TEND = atof(argv[3]); // End time
  if (argc > 4)
    nu = atof(argv[4]); // Viscosity for vertical diffusion
  if (argc > 5) 
    RANDOM = atoi(argv[5]); // An integer to seed the random number generator (don't use 0 or 1)
  if (argc > 6)
    L0 = atof(argv[6]); // Box size (not necessarily related to peak wave number)
  if (argc > 7) 
    kp_ = atof(argv[7]); // Peak wavelength 
  if (argc > 8)
    theta_H = atof(argv[8]); // Numerical parameter (blend explicit/implicit scheme to dump fast barotropic modes)
  if (argc > 9)
    g_ = atof(argv[9]); // gravitational acceleration
  if (argc > 10)
    Bo = atof(argv[10]); // Bond number

  if (argc < 11) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 11-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " NLAYER = %d\n ", NLAYER), fflush (stderr);
  fprintf(stderr, " Horizontal resolution = %d\n ", LEVEL_data), fflush (stderr);
  fprintf(stderr, " End time  = %8E\n ", TEND), fflush (stderr);
  fprintf(stderr, " Viscosity for vertical diffusion  = %8E\n ", nu), fflush (stderr);
  fprintf(stderr, " RANDOM integer = %d\n ", RANDOM), fflush (stderr);
  fprintf(stderr, " Box size = %8E\n ", L0), fflush (stderr);
  fprintf(stderr, " Peak wavelength = %8E\n ", kp_), fflush (stderr);
  fprintf(stderr, " theta_H = %8E\n ", theta_H), fflush (stderr);
  fprintf(stderr, " Gravity = %8E\n ", g_), fflush (stderr);
  fprintf(stderr, " Bond Number = %8E\n ", Bo), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     Create directories for better organization of the output */

  fprintf(stderr, "Create directories if needed\n"), fflush (stderr);

  char comm[80];
  sprintf (comm, "mkdir -p surface");
  system(comm);
  sprintf (comm, "mkdir -p field");
  system(comm);

  /**
     Set origin, grid, number of layers */

  origin (-L0/2., -L0/2.);
  N = 1 << LEVEL_data;
  nl = NLAYER;

  /**
     Define surface tension coefficient and compute the wave period */

  sigma_val = 0.0; //to remove Bond effect //rho_w*g_/(Bo*sq(kp_));
  T0_ = 2.*pi/sqrt(g_*kp_ + sigma_val*kp_*kp_*kp_/rho_w);
  fprintf(stderr, " Wave period = %8E\n ", T0_), fflush (stderr);

  /**
     Set boundary conditions */

  periodic (right);
  periodic (top);

  /**
     Set gravity, h and gpe */

  G = g_;
  h_ = 80.0; // set the water depth to be the peak wave length (the deep water assumption)
#if dimension == 2
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
#else
  gpe_base = -0.5*sq(h_)*L0*g_;
#endif

  /**
     Set CFL */

  CFL_H = HUGE;
  
  fprintf (stderr, "Read in parameters!\n");

  double us = 0; //Set the viscous stress at the interface if needed
  du0 = sq(us)/nu;
  fprintf (stderr, "du0 = %8E\n", du0);
  dut = duv;

  run();

}

/**
## Initialization
Read in the power spectrum and initialize the wave field. 
*/

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " NLAYER = %d\n ", NLAYER), fflush (stderr);
  fprintf(stderr, " Horizontal resolution = %d\n ", LEVEL_data), fflush (stderr);
  fprintf(stderr, " End time  = %8E\n ", TEND), fflush (stderr);
  fprintf(stderr, " Viscosity for vertical diffusion  = %8E\n ", nu), fflush (stderr);
  fprintf(stderr, " RANDOM integer = %d\n ", RANDOM), fflush (stderr);
  fprintf(stderr, " Box size = %8E\n ", L0), fflush (stderr);
  fprintf(stderr, " Peak wavelength = %8E\n ", kp_), fflush (stderr);
  fprintf(stderr, " theta_H = %8E\n ", theta_H), fflush (stderr);
  fprintf(stderr, " Gravity = %8E\n ", g_), fflush (stderr);
  //fprintf(stderr, " Bond number = %8E\n ", 1.0*g_/(sigma_val*sq(kp_))), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
         Set surface stress */

  foreach() {
    duv.x[] = du0;
  }
  // set tracer scalar
  T = new scalar[nl];
  tracers = list_append (tracers, T);

  if (!restore ("restart")) {
    power_input();
    dkx_ = kx_[1] - kx_[0];
    dky_ = ky_[1] - ky_[0];
    fprintf (stderr, "dkx = %g, dky = %g\n", dkx_, dky_);
    geometric_beta (1./3., true); // Varying layer thickness
    foreach() {
      zb[] = -h_;
      eta[] = wave(x, y);
      double H = wave(x, y) - zb[];
      foreach_layer() {
	h[] = H/nl;
      }
    }

    vertical_remapping (h, tracers);
    foreach() {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        u.x[] = u_x(x, y, z);
        u.y[] = u_y(x, y, z);
        w[] = u_z(x, y, z);
        T[] = 0.;
        z += h[]/2.;
      }
      for (int i=nl-10; i<nl; i++)
        T[0,0,i] = 1.;
    }
    fprintf (stderr,"Done initialization!\n");
    dump("initial");
    foreach() 
      sigma[] = sigma_val/rho_w;
  }
  else {
    geometric_beta (1./3., true); // when restarting, remember to specify the grid mapping method
    dtmax = 0.01;
    dt = dtnext (dtmax);
    char *suffix = "matrix";
    writefields (i, suffix);
    foreach() { 
      sigma[] = sigma_val/rho_w;
    }
  }

}

event face_fields (i++) {
  // Calculate the Std of eta
  double eta_sum = 0.;
  double eta_sum_2= 0.;
  double count   = 0.;
  foreach(reduction(+:eta_sum) reduction(+:eta_sum_2) reduction(+:count)){
    eta_sum   += eta[];
    eta_sum_2 += eta[]*eta[];
    count     += 1.;
  }
  eta_sum = eta_sum/count;
  eta_sum_2= eta_sum_2/count;
  eta_std = sqrt(eta_sum_2-eta_sum*eta_sum);

  phat = 0.1;// A constant wind pressure forcing
  //phat = (0.1613*0.1613-eta_std*eta_std)/0.1613*g_*2./pi*0.5; // A feedback forcing based on the target eta_std value 0.1613 or other values.

  double us    = 0;
  du0 = sq(us)/nu;
  foreach() {
    duv.x[] = du0;
  }
  dut = duv;

    dtmax = DT;
}

event face_fields (i++,first) {
  foreach_face (x) {
    double ax = a_baro_p (eta, 0);
    foreach_layer() {
      ha.x[] += hf.x[]*ax;
    }
  }
}

/** Add horizontal diffusion. */

event viscous_term (i++)
  horizontal_diffusion ((scalar *){u,w}, nu, dt);

/**
Outputs
*/

event energy_before_remap (i+=10, last) {

  if (i==10) {
    fprintf(stderr, "energy output before remap!\n");
    fflush(stderr);
  }
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer () {
      double norm2 = sq(w[]);
      foreach_dimension()
        norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  
  if (pid() == 0) {
    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"energy_before_remap.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %8E %8E\n", t, 1.0*i, ke/2.0, g_*gpe-gpe_base, eta_std);
    fclose(log_sim);
  }

}

/**
   Output 3-D field if needed */

event field_log (t = 0; t += 5.0; t <= TEND*T0_) {
  
  char *suffix = "matrix";
  writefields (t, suffix);
  //writefields (i, suffix);
  
  if (pid() == 0) {
    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_field.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
    fclose(log_sim);
  }

}

/**
   Output the surface layer. */

event surface_log (t = 0; t += 0.2; t <= TEND*T0_) {
//event surface_log (t = 0; t += 1.; t <= TEND*T0_) {
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50], filename5[50];
  sprintf (filename1, "surface/eta_matrix_%09d", (int) i);
  sprintf (filename2, "surface/ux_matrix_%09d", (int) i);
  sprintf (filename3, "surface/uy_matrix_%09d", (int) i);
  sprintf (filename4, "surface/uz_matrix_%09d", (int) i);
  sprintf (filename5, "surface/T_matrix_%09d", (int) i);
  FILE * feta = fopen (filename1, "w");
  // Might need to change to mpi function later
  output_matrix_mpi (eta, feta, N, linear = true);
  fclose (feta);
  sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi (u_temp.x, fux, N, linear = true);
  fclose (fux);
  FILE * fuy = fopen (filename3, "w");
  output_matrix_mpi (u_temp.y, fuy, N, linear = true);
  fclose (fuy);
  sprintf (s, "w%d", nl-1);
  scalar w_temp = lookup_field (s);
  FILE * fuz = fopen (filename4, "w");
  output_matrix_mpi (w_temp, fuz, N, linear = true);
  fclose (fuz); 
  sprintf (s, "T%d", nl-1);
  scalar T_temp = lookup_field (s);
  FILE * fT = fopen (filename5, "w");
  output_matrix_mpi (T_temp, fT, N, linear = true);
  fclose (fT);
  
  if (pid() == 0) {
    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_surface.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
    fclose(log_sim);
  }

}

event regulardump (t = 0; t +=T0_*5.0; t <= TEND*T0_) {
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
}

event endrun (t = TEND*T0_) {
  dump();
  delete ({T});
  return 1;
}