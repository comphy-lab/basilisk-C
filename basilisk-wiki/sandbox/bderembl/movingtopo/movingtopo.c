
/**
   Numerical clone of Dossmann et al 2016
   https://doi.org/10.1002/2016JC011990

   Inspired from 
   http://basilisk.fr/sandbox/popinet/movingcylinder.c

   Compile with    
   CC99='mpicc -std=c99' qcc -D_MPI=1 -D_NETCDF=1 -grid=multigrid3D -lm -lpnetcdf -O3 movingtopo.c -o movingtopo.e

   mpirun -np 4 movingtopo.e

   HPC:
   qcc -D_MPI=1  -D_NETCDF=1 -grid=multigrid3D -source movingtopo.c

*/

//#include "grid/multigrid.h"
//#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "../libs/boussinesq.h"
#include "../libs/extra.h"
#include "fractions.h"

#if _NETCDF
#include "../libs/pnetcdf_bas.h"
#endif


double mu0 = 0.;  //
double kappa0 = 0.;
double h0 = 0.06;
double lr = 0.05;
double xm = 0.;
double N2 = 1.;
double Ut = 1.;
double d_acc = 0.15;
double d_plate = 0.4;
double tau1;
double tau2;

int npy = 1;
int npz = 1;
int slip_bc = 1; // free slip by default

double tend = 5.;
double dtout = 0.1;

scalar topo[];
#if dimension > 2
#define MOUNTAIN (h0*exp(-(sq((x - xm)/lr))) - z)
#define STRATIFICATION  ((z*N2))
#else
#define MOUNTAIN (h0*exp(-(sq((x - xm)/lr))) - y)
#define STRATIFICATION  ((y*N2))
#endif

int main() {

  // declare user parameters
  params = array_new();
  add_param ("N", &N, "int");
  add_param ("npy", &npy, "int");
  add_param ("npz", &npz, "int");
  add_param ("slip_bc", &slip_bc, "int");
  add_param ("L0", &L0, "double");
  add_param ("h0", &h0, "double");
  add_param ("lr", &lr, "double");
  add_param ("Ut", &Ut, "double");
  add_param ("N2", &N2, "double");
  add_param ("xm", &xm, "double");
  add_param ("mu", &mu0, "double");
  add_param ("tend", &tend, "double");
  add_param ("dtout", &dtout, "double");
  add_param ("kappa", &kappa0, "double");

  // Read the configuration file params.in
  read_params("params.in");
  create_outdir();
  backup_config();

  // Topo CFL
  //dtmax = L0/N/Ut;

  dimensions(ny=npy);
#if dimension > 2
  dimensions(ny=npy, nz=npz);
#endif

  const face vector muc[] = {mu0, mu0, mu0};
  const face vector kappac[]  = {kappa0, kappa0, kappa0};

  if (mu0 > 0) mu = muc;
  if (kappa0 > 0) kappa = kappac;



  tau1 = 2*d_acc/Ut;
  tau2 = (L0 - d_plate/2 - 2*d_acc)/Ut;
//  tend = tau2 + 2*tau1;

  run();


  array_free (params);
}

event init (t = 0) {


  if (slip_bc == 0){
    u.t[top]    = dirichlet(0.);
    u.t[bottom] = dirichlet(0.);
    u.t[left]   = dirichlet(0.);
    u.t[right]  = dirichlet(0.);
#if dimension > 2
    u.t[front]  = dirichlet(0.);
    u.t[back]   = dirichlet(0.);

    u.r[top]    = dirichlet(0.);
    u.r[bottom] = dirichlet(0.);
    u.r[left]   = dirichlet(0.);
    u.r[right]  = dirichlet(0.);
    u.r[front]  = dirichlet(0.);
    u.r[back]   = dirichlet(0.);
#endif
  }

  coord vc = {Ut,0., 0.}; // the velocity of the topo
  xm = 0.;
  fraction (topo, MOUNTAIN);

  foreach() {
    b[] = STRATIFICATION;

    foreach_dimension()
      u.x[] = topo[]*vc.x + 0.05*Ut*noise();
  }

#if _NETCDF
  FILE * fp;
  if ((fp = fopen("restart.nc", "r"))) {
    read_nc({b,u}, "restart.nc");
    fclose(fp);
  }
#endif

#if _NETCDF
  char file_nc[90];
  sprintf (file_nc, "%svars.nc", dpath); 
  create_nc({b,u}, file_nc);
#endif
}

event moving_topo (i++) {

  double u_var;

  if (t<tau1){
    u_var = Ut*t/tau1;
    xm = Ut*t*t/(2*tau1);
  }
  else if((t>=tau1) & (t<tau1 + tau2)) {
    u_var = Ut;
    xm = d_acc + Ut*(t-tau1);
    }
  else if ((t>=tau1 + tau2) & (t<2*tau1 + tau2)) {
    double tprime = t - tau1 - tau2;
    u_var = Ut*(1 - tprime/tau1);
    xm = d_acc + Ut*tau2 + Ut*tprime - Ut*tprime*tprime/(2*tau1);
  }
  else {
    u_var = 0;
    xm = L0 - d_plate/2;
  }

  coord vc = {u_var,0., 0.}; // the velocity of the topo
  fraction (topo, MOUNTAIN);

  foreach()
    foreach_dimension()
      u.x[] = topo[]*vc.x + (1. - topo[])*u.x[];
  boundary ((scalar *){u});

}



event output (t = 0; t <= tend+1e-10;  t += dtout) {

#if _NETCDF
    write_nc();
#endif
}


event logfile (i++)
  fprintf (stderr, "%d %g %d %d %d %d\n", i, t,
	   mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
