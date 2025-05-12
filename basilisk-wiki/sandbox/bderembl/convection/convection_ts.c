/**
# Oceanic convection

This is a code for oceanic convection.

The ocean is forced at the surface by wind stress (tau0 in m2/s2), heat flux
(Qnet in W/m2), and fresh water flux (EmP in m/s). These fluxes impact the upper
layer velocity, temperature and salinity budgets.

The seawater (boussinesq) density is computed either with a linear (eos_type =
1) or non-linear (eos_type = 2) equation of state (UNESCO formula).

The linear equation of state is 

rho = rho0 (1 - alpha (T - T0) + beta (S - S0))

We prescribe an initial temperature and salinity profiles (dTdz0 and dSdz0)
respectively.

We use the full Coriolis force: f0_z is the vertical component of coriolis
parameter and f0_y the horizontal component.


## Compilation

for 2d code: -grid=multigrid
for 3d code: -grid=multigrid3D


qcc -grid=multigrid3D -D_NETCDF=1 -lm -O3 convection_ts.c -o cv.e

CC99='mpicc -std=c99' qcc -D_MPI=1  -D_NETCDF=1 -grid=multigrid3D -lm -lpnetcdf -O3 convection_ts.c -o cv.e

HPC:
qcc -D_MPI=1  -D_NETCDF=1 -grid=multigrid3D -source convection_ts.c

*/

double nu = 0.;
//double kappa = 1.e-5;
double Qnet = 0.; // heat flux in W/m2 (negative = cooling)
double EmP = 0.; // Evap minus Precip in m/s
double tau0 = 0.; // surface stress in m2/s2
double Tsurf = 20.; // surface temperature
double Ssurf = 0.; // surface salinity
double dTdz0 = 0.; // background vertical temperature gradient
double dSdz0 = 0.; // background vertical salinity gradient
double hm_i = 0.; // initial mixed layer thickness
double alphaT = 0.;  // thermal expansion computed around Tsurf
double betaS = 0.;   // haline contraction computed around Ssurf

// physical constants
double cp = 4000; // J/kg/K
double rho0 = 1000; //kg/m3
double gg = 9.8;  // m/s2

int NMOUT = 5000;
int N0; // bug in restore
double tend = 1.;
double DT_MAX = 100.;
double dtout = 1.;
double tout0 = 0;
char dpath[80]  = "./";
timer tel;
double tel_m = 1e10;
int npt = 0;
double nu0;
double T_noise = 1e-3;
double S_noise = 1e-3;
double tau_sponge = 0.;
double tol0 = 1e-7 [*];
double tol1 = 1e-4 [*];
int bottom_bc_TS = 1; // 1: neuman(dTdz), 0: neumann (0)
int eos_type = 1; // 1: linear, 2: non linear

// number of processor in z direction
int npz = 0;


#include "../libs/extra.h"
#include "navier-stokes/centered.h"
#include "../libs/coriolis3d.h"
#if _NETCDF
#include "../libs/pnetcdf_bas.h"
#endif


/**
   TS version of the boussinesq equations

   Qnet is the surface heat flux
   EmP is Evaporation minus precipitation. If EmP > 0 then dSdt > 0 (evaporation)
   https://mitgcm.readthedocs.io/en/latest/examples/global_oce_latlon/global_oce_latlon.html
*/

#include "tracer.h"
//#include "diffusion.h"
#include "../libs/sw_eos-80.h"

face vector av[];
scalar T[];
scalar S[];
scalar * tracers = {T,S};

event init (t = 0) {
  a = av;
}

event acceleration (i++) {

#if dimension > 2
  coord grav_dir = {0, 0, 1};
#else
  coord grav_dir = {0, 1};
#endif
  
/**
   Define the buoyancy variable.

   For now, we restrict to a linear equation of state because we need a boundary
   condition for b for the pressure solver. This is relatively easy to implement
   in the case of a linear equation of state.
*/
  scalar b[];
  if (bottom_bc_TS) {
#if dimension > 2
  b[back] = neumann(gg*(-dTdz0*alphaT -dSdz0*(-betaS)) );
#else
  b[bottom] = neumann(gg*(-dTdz0*alphaT -dSdz0*(-betaS)) );
#endif
  }
  if (eos_type == 1) {
    foreach() {
      b[] = gg*(alphaT*(T[] - Tsurf) - betaS*(S[] - Ssurf));
    }
  } else {
    foreach() {
      b[] = eos_b_s(T[],S[]);
    }
  }
  norm nb = normf(b);
    foreach() {
      b[] -= nb.avg;
    }

  foreach_face()
    av.x[] = grav_dir.x*face_value(b, 0);


}

/* event diffusion (i++) { */
/*   mgstats mgT; */
/*   mgstats mgS; */
/*   double kappa = 1e-2; */
/*   const face vector D[] = {kappa, kappa, kappa}; */
  
/*   mgT = diffusion (T, dt, D); */
/*   mgS = diffusion (S, dt, D); */
/* } */

/**
   
   Main code
*/


int main(int argc,char* argv[]) {
  tel = timer_start();

  params = array_new();
  add_param ("N", &N, "int");
  add_param ("npz", &npz, "int");
  add_param ("L0", &L0, "double");
  add_param ("nu", &nu, "double");
  add_param ("dTdz0", &dTdz0, "double");
  add_param ("dSdz0", &dSdz0, "double");
  add_param ("Tsurf", &Tsurf, "double");
  add_param ("Ssurf", &Ssurf, "double");
  add_param ("Qnet", &Qnet, "double");
  add_param ("EmP", &EmP, "double");
  add_param ("tau0", &tau0, "double");
  add_param ("hm_i", &hm_i, "double");
  add_param ("f0_x", &f0_x, "double");
  add_param ("f0_y", &f0_y, "double");
  add_param ("f0_z", &f0_z, "double");
  add_param ("alpha_H", &alpha_H, "double");
  add_param ("T_noise", &T_noise, "double");
  add_param ("S_noise", &S_noise, "double");
  add_param ("tau_sponge", &tau_sponge, "double");
  add_param ("DT_MAX", &DT_MAX, "double");
  add_param ("CFL", &CFL, "double");
  add_param ("tol0", &tol0, "double");
  add_param ("tol1", &tol1, "double");
  add_param ("tend", &tend, "double");
  add_param ("dtout", &dtout, "double");
  add_param ("tout0", &tout0, "double");
  add_param ("NMOUT", &NMOUT, "int");
  add_param ("bottom_bc_TS", &bottom_bc_TS, "int");
  add_param ("eos_type", &eos_type, "int");

  // origin
  X0 = 0;
#if dimension > 2
  Y0 = 0;
  Z0 = -L0; // will be adjsuted in init event in case non cubic domain
#else
  Y0 = -L0;
#endif

  N0 = N; // bug in restore

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in

  read_params(file_param);
  create_outdir();
  backup_config(file_param);
  
  DT = 1e-5;
  TOLERANCE = tol0;

//  TOLERANCE = 1e-4;
  periodic(right);
#if dimension > 2
  periodic(top);
  if (npz > 0){
#if _MPI
    int rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);
#endif
    dimensions(nz=npz);
  }
#endif

  u.x.gradient = minmod;
  u.y.gradient = minmod;
  u.z.gradient = minmod;
  run();
}

/**
## Initial conditions
Initial conditions correspond to a tanh profile for both buoyancy and velocity
*/



event init (t=0) {
#if _MPI
  Z0 = -L0*mpi_dims[2]/mpi_dims[0];
#endif


      foreach() {
#if dimension > 2
        T[] = z < -hm_i ? dTdz0*(z + hm_i) + Tsurf: Tsurf;
        S[] = z < -hm_i ? dSdz0*(z + hm_i) + Ssurf: Ssurf;
#else
        T[] = z < -hm_i ? dTdz0*(y + hm_i) + Tsurf: Tsurf;
        S[] = z < -hm_i ? dSdz0*(y + hm_i) + Ssurf: Ssurf;
#endif
        T[] += T_noise*noise();
        // no noise on S (S>0)
      }

      alphaT = eos_alpha_s(Tsurf,Ssurf);
      betaS  = eos_beta_s(Tsurf,Ssurf);

      foreach() {
        u.x[] = 0.;
        u.y[] = 0.;
#if dimension > 2
        u.z[] = 0.;
#endif
      }

  N = N0; // restore function changes N


#if dimension > 2
  if (bottom_bc_TS) {
    T[back] = neumann(-dTdz0);
    S[back] = neumann(-dSdz0);
  }
  u.n[front] = dirichlet(0);
  u.n[back] = dirichlet(0);
#else
  if (bottom_bc_TS) {
    T[bottom] = neumann(-dTdz0);
    S[bottom] = neumann(-dSdz0);
  }
#endif

  boundary (all);

#if _NETCDF
  char file_nc[80];
  sprintf (file_nc, "%svars.nc", dpath); 
  create_nc({T,S,u}, file_nc);
#endif
}


/**
   Surface forcing (temperature and salt flux)
*/

// trying this formulation instead of specification of neumann BC so that it is
// independent of nu
event surface_fluxes (i++) {
#if dimension > 2
  foreach_boundary(front) {
#else
  foreach_boundary(top) {
#endif
    T[]   += dt*Qnet/(rho0*cp*Delta);
    S[]   += dt*EmP*S[]/Delta;
  }
}

/**
   Surface forcing (momentum flux)
*/

event acceleration (i++) {

#if dimension > 2
  foreach_boundary(front) {
#else
  foreach_boundary(top) {
#endif
    u.x[] += dt*tau0/Delta;
  }

}

  
event bottom_sponge (i++) {
  
  if (tau_sponge > 0){

    foreach(){
      if (z < 0.9*Z0){
        T[] -= dt*(T[] - (dTdz0*(z + hm_i) + Tsurf))/tau_sponge;
        S[] -= dt*(S[] - (dSdz0*(z + hm_i) + Ssurf))/tau_sponge;
        foreach_dimension()
          u.x[] -= dt*u.x[]/tau_sponge;
      }
    }
  }
}
 
/**
   Output
*/

event output (t = tout0; t <= tend+1e-10;  t += dtout) {

  double ket = 0., vd = 0.;
  foreach(reduction(+:ket)) {
      ket += 0.5*(sq(u.x[]) + sq(u.y[]))*sq(Delta);
  }
  

  fprintf (stdout,"i = %d, dt = %g, t = %g, ke = %.2e\n",
           i, dt, t, ket);

  if (N < NMOUT){
#if _NETCDF
    write_nc();
#endif
  } else {
    char name[100];
    sprintf (name, "%soutput-%09d", dpath, i);
    fprintf(stdout,name);
    dump (name);
  }

}

event snapshot(i++){
  double elapsed = timer_elapsed(tel);

  if ( elapsed > tel_m) {
    char name[100];
    sprintf (name, "%ssnapshot-%09d", dpath, i);
    dump (name);
    return 1;
  }
}

/**
## Now trying to keep salinity under control....

if we have dSdz0 = 0 (uniform salinity), we expect that it will remain
uniform. But it does not. The primary reason is the 


This *back_to_defaults* event is added so that we can start with a very small timestep (*DT*) and tolerance for the Poisson problems (*TOLERANCE*), so that the iterative solver has some low-impact iterations before it finds the not-consistently-initialized pressure field.
*/
event back_to_defaults(i=1000){
  /* DT = DT2; */
//  DT = DT_MAX;
  TOLERANCE = tol1;
}
/* event back_to_defaults2(i=1000){ */
/*   /\* DT = DT2; *\/ */
/*   DT = 10; */
/* } */

event increase_dt(i++){
  /* DT = DT2; */
  stats s = statsf(S);

  if ( DT < DT_MAX)
    DT = 1.01*DT;
  fprintf (stdout,"i = %d, dt = %g, t = %g, min = %g, max = %g\n",
           i, dt, t, s.min, s.max);
}
