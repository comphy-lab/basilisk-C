
/**
   This program writes and read a netcdf file
*/


// CC99='mpicc -std=c99' qcc -D_MPI=1 -O3 test_netcdf.c -o test_netcdf.e -lm -lnetcdf 


#include "grid/multigrid.h"
#include "run.h"
#include "../libs/netcdf_bas.h"

scalar omega[];
scalar psi[];

char fileout[100] = "test.nc";

int main() {
  N = 4;
  init_grid (N);
  foreach() {
    omega[] = x + y;
    psi[] = y;
  }


/**
   Create the netcdf file. we must specify which variable we want to store and the file name.
*/

  create_nc({omega, psi}, fileout);

/**
   Store a snapshot
*/

  write_nc();

  printf("Storing data in %s: \n", fileout);
  foreach() {
    printf ("x = %e, y = %e, psi = %e omega = %e\n", x,y, psi[], omega[]);
  }

/**
   Pretend we do a time step
*/

  t = 1;
  foreach() {
    omega[] += 1.;
    psi[] += 2.;
  }

/**
   Store new snapshot
*/

  write_nc();

/**
   Read the netcdf file. we can only read the first iteration
   TODO: make it more generic so that we can read any iteration.
*/

  foreach() {
    omega[] = 0.;
    psi[] = 0.;
  }
  printf("Reading data in %s: \n", fileout);

  read_nc({omega, psi}, "test.nc");

  foreach() {
    printf ("x = %e, y = %e, psi = %e omega = %e\n", x,y, psi[], omega[]);
  }
  
  run(); // this matters to cleanup the global list
  
}
