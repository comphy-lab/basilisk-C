/**

# Kelvin Helmholtz instability 

This routine is heavily inspired from http://basilisk.fr/sandbox/Antoonvh/kh.c
and the boussinesq module from
http://basilisk.fr/sandbox/acastillo/convection/convection_boussinesq.h


Compile with


qcc -grid=multigrid -lm -O3 kelvin_helmholtz.c -o kh.e

CC99='mpicc -std=c99' qcc -D_MPI=1 -grid=multigrid -lm -O3 kelvin_helmholtz.c -o kh.e

HPC:
qcc -D_MPI=1 -grid=multigrid -source kelvin_helmholtz.c -D_NETCDF=1
*/

#include "grid/multigrid.h"
#include "boussinesq.h"
#include "auxiliar_input.h"
#include "higher-order.h"
//#define interpolate_array interpolate_array_3
#define interpolate_linear interpolate_quadratic
#include "tracer-particles.h"
#if _NETCDF
#include "netcdf_bas.h"
#include "netcdf_particles.h"
#endif

char file_nc[100];
char file_nc_part[100];

Particles ptr;

int main() {
  tel = timer_start();

  read_params();
  create_outdir();
  backup_config();
  

  TOLERANCE = 1e-4;
  periodic(right);

  run();
}

/**
## Initial conditions
Initial conditions correspond to a tanh profile for both buoyancy and velocity
*/

event init (t=0) {
  if (!restore (file = "restart.bin")) {
    FILE * fp;
    if ((fp = fopen("b0.bas", "r"))) {
      input_matrixl ({b}, fp, ox=X0, oy=Y0);
      fclose(fp);
    } else {
      foreach()
        b[] = tanh(y) + N2b*y;
    }
    
    if ((fp = fopen("u0.bas", "r"))) {
      input_matrixl ((scalar *) {u}, fp, ox=X0, oy=Y0);
      fclose(fp);
    } else {
      foreach() {
        u.x[] = tanh(y) + 1e-1*noise();
        u.y[] = 0.;
    }
    }
  }
  N = N0; // restore function changes N

  boundary (all);

  if (npt > 0)
    ptr = init_tp_square (npt,0,0,0.8);

#if _NETCDF
  sprintf (file_nc, "%svars.nc", dpath); 
  create_nc({b,u}, file_nc, 1);

  if (npt > 0){
    sprintf (file_nc_part, "%sparts.nc", dpath); 
    create_nc_particles(file_nc_part);
  }
#endif
}

event output (t = 0; t <= tend+1e-10;  t += dtout) {

  double ket = 0., vd = 0.;
  foreach(reduction(+:ket) reduction(+:vd)) {
      ket += 0.5*(sq(u.x[]) + sq(ar*u.y[]))*sq(Delta)*ar;
      //      viscous dissipation // FIXME: aspect ratio
      vd += sq(Delta)*ar*(sq(u.x[1] - u.x[-1]) +
                          sq(u.x[0,1] - u.x[0,-1])/sq(ar) +
                          sq(u.y[1] - u.y[-1])*sq(ar) +
                          sq(u.y[0,1] - u.y[0,-1])
                          )/sq(2.*Delta);
  }
  vd *= 1./Re;
  

  fprintf (stdout,"i = %d, dt = %g, t = %g, ke = %.2e, vd = %.2e\n",
           i, dt, t, ket, vd);

  if (N0 < NMOUT){
#if _NETCDF
    write_nc({b,u}, file_nc, 1, t, linear = true);
#else
    char name[100];
    sprintf (name,"%silast.dat", dpath);
    FILE * fp = fopen (name, "w");
    fprintf(fp,"%d",i);
    fclose(fp);
  
    sprintf (name,"%sb%09d.bas", dpath, i);
    fp = fopen (name, "w");
    output_matrixl ({b}, fp, linear = true);
    fclose(fp);
  
    sprintf (name,"%su%09d.bas", dpath, i);
    fp = fopen (name, "w");
    output_matrixl ((scalar *) {u}, fp, linear = true);
    fclose(fp);
  
#if TREE
    if (tol_b || tol_u) {
      scalar l[];
      foreach()
        l[] = level;
      boundary({l});
      sprintf (name,"%slevel%09d.bas", dpath, i);
      fp = fopen (name, "w");
      output_matrixl ({l}, fp, linear = true);
      fclose(fp);
    }
#endif
#endif
  } else {
    char name[100];
    sprintf (name, "%soutput-%09d", dpath, i);
    fprintf(stdout,name);
    dump (name);
  }

  if (npt > 0) {

    particle_boundary(ptr);
#if _NETCDF
    write_nc_particles(ptr,file_nc_part);
#else

    scalar epsilon[];
    scalar ke[];
    scalar del2b[];
    foreach() {
      del2b[] = laplacian(b);
      ke[] = 0.5*(sq(u.x[]) + sq(ar*u.y[]));
      //      viscous dissipation // FIXME: aspect ratio
      epsilon[]= (sq(u.x[1] - u.x[-1]) +
                  sq(u.x[0,1] - u.x[0,-1])/sq(ar) +
                  sq(u.y[1] - u.y[-1])*sq(ar) +
                  sq(u.y[0,1] - u.y[0,-1])
                  )/sq(2.*Delta);
    }
    boundary ({del2b,ke,epsilon});

    char name[80];
    sprintf (name, "%spart-%09d-%04d.dat", dpath, i,pid());
    FILE * fpp = fopen (name, "w");
    if (!fpp) {
      fprintf (stdout,"could not open particle file");
      return 1;
    }

    foreach_particle() {
      fprintf(fpp, "%lu %g %g %g %g %g %g ", p().tag, x, y,
              interpolate_5 (b, x, y),
              interpolate_5 (del2b, x, y),
              interpolate_5 (epsilon, x, y),
              interpolate_5 (ke, x, y));
      fprintf(fpp, "%g %g \n", p().u.x, p().u.y);
    }
    fclose(fpp);
#endif
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
