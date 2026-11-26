#include "embed.h"
#include "navier-stokes/centered.h"
#include "contact-fixed.h"
#define FILTERED 1
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

#define MINlevel 3                                              // minimum level

#define tsnap (1e0)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-3)                            // error tolerances in velocity


double  hf, tmax, Ldomain, Ca, mu_r, rho_r, theta0;
int MAXlevel;
double wind_speed;

// boundary conditions
u.t[bottom] = dirichlet(Ca);
u.n[bottom] = dirichlet(0.0);
p[bottom] = neumann(0.0);

u.t[right] = neumann(0.);
u.n[right] = dirichlet(-wind_speed);
p[right] = neumann(0.);

u.t[left] = neumann(0.);
u.n[left] = neumann(0.);
p[left] = dirichlet(0.0);

u.t[top] = neumann(0.);
u.n[top] = neumann(0.);
p[top] = dirichlet(0.0);

u.x[embed] = dirichlet(0.);
u.y[embed] = dirichlet(0.);


void restore_fs(){
  vertex scalar phi[], phi_vert[], phi_hor[];
  foreach_vertex(){
    phi_vert[] = -(x - 0.11 * L0);
    phi_hor[] = -(y - 0.11 * L0);
    phi[] = union(phi_hor[], phi_vert[]);
  }

  fractions (phi, cs, fs);
  boundary ({cs, fs});
  fractions_cleanup (cs, fs);
  boundary ({cs, fs});
  restriction ({cs, fs});
}

vector h[];
h.t[bottom] = contact_angle(theta0);

int main(int argc, char const *argv[]) {

  Ca         = argc > 1 ? atof(argv[1]) : 0.03;
  wind_speed = argc > 2 ? atof(argv[2]) : 1.;
  MAXlevel   = argc > 3 ? atoi(argv[3]) : 10;
  tmax       = argc > 4 ? atof(argv[4]) : 1e3;

  mu_r = 2e-2;
  rho_r = 1e-3;
  theta0 = 50.0*pi/180.0;

  // Ldomain = lr > 1 ? 32 : 32*lr;
  Ldomain = 2e1;
  hf = 0.5*Ldomain;

  fprintf(ferr, "Level %d tmax %g., hf %3.2f, Ca = %g\n", MAXlevel, tmax,  hf, Ca);

  L0=Ldomain;
  X0=-hf; Y0=0.;
  N = 1 << MAXlevel;
  char comm[80];
  sprintf (comm, "mkdir -p intermediate_%g facets_%g", Ca, Ca);
  system(comm);

  rho1 = 1e0; mu1 = 1e0;
  rho2 = rho1*rho_r; mu2 = mu1*mu_r;

  G.x = -1e0;

  f.sigma = 1.0;
  f.height = h;

  run();
}

event init(t = 0){
  char restart_f[80];
  sprintf(restart_f, "restart_%g", Ca);
  if(!restore (file = restart_f)) {
    vertex scalar phi[], phi_vert[], phi_hor[];
    foreach_vertex(){
      phi_vert[] = -(x - 0.11 * L0);
      phi_hor[] = -(y - 0.11 * L0);
      phi[] = union(phi_hor[], phi_vert[]);
    }

    fractions (phi, cs, fs);
    boundary ({cs, fs});
    fractions_cleanup (cs, fs);
    boundary ({cs, fs});
    restriction ({cs, fs});

    fraction(f, -x);
  }

  restore_fs();
}


scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
  adapt_wavelet((scalar *){f, u.x, u.y},
    (double[]){fErr, VelErr, VelErr},
    MAXlevel, MINlevel);
}

// Outputs
scalar pos[];
double x_cl_old = 0.;
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  char restart_file[80];
  sprintf(restart_file, "restart_%g", Ca);
  dump (file = restart_file);
  char nameOut[80];
  sprintf (nameOut, "intermediate_%g/snapshot-%5.4f", Ca, t);
  dump (file = nameOut);

  position(f, pos, {1,0 [0]});
  double max = statsf(pos).max;
  
  if (pid() == 0) {
    FILE * fp_inter_time = fopen ("cl_pos.dat", "a");

    fprintf(fp_inter_time, "%g %g\n", t, max);
    fclose(fp_inter_time);
  }

  fprintf(ferr, "dh = %g \n", max - x_cl_old );
  if (t > 0 && fabs(max - x_cl_old) < 1e-7){
    return 1;
  }
  if (max > 9.5){
    return 1;
  }
  x_cl_old = max;

}

event outputInterface(t += tsnap; t < tmax) {
  int N = 1 << MAXlevel;
  char facets_file[80];
  sprintf(facets_file, "facets_%g/facets-%d-t_%0.2f_%d", Ca, N, t, pid());
  FILE * fp_facets = fopen(facets_file, "w");

  output_facets (f,fp_facets);
  fclose(fp_facets);
  char command[200];
  sprintf(command, "LC_ALL=C  cat facets_%g/facets-%d-t_%0.2f_* > facets_%g/facets-%d-t_%0.2f", Ca, N, t, Ca, N, t);
  system(command);// allow to use linux command in the c code to concatenate our files

  char cleanup[80];
  sprintf(cleanup, "rm -f facets_%g/facets-%d-t_%0.2f_*", Ca, N, t);
  system(cleanup);
}