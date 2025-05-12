/* You need to run this code using MPI and see the out file for understanding the bug CC='mpicc -D_MPI=4' make foreach_bound_MPI.tst**/
#include "grid/tree.h"
#include "run.h"
#include "view.h"

#define dA sq(Delta)
double nusselt_top (scalar T)
{
  double nutop=0.;
  foreach_boundary (top,reduction(+:nutop))
    nutop += dA*(T[] - T[0,1])/Delta;
  return nutop;
}

scalar T[];

int main() {
  L0 = 1;
  init_grid (128);
  run();
}

T[top] = dirichlet(-0.5);

event init(t=0){
  foreach()
    T[] = 1.;
  boundary({T});
}

//It works
#if 0
event flux(t=0){
  foreach_boundary(top)
    printf("%.9g %.9g %.9g %.9g\n", x, y, T[0,1], T[]);
}
#endif


//It works
#if 0
event flux(t=0){
  char names[80];
  sprintf(names, "test.%d.dat", pid());
  FILE * fp = fopen (names, "w");
  foreach_boundary(top)
    fprintf(fp, "%.9g %.9g %.9g %.9g\n", x, y, T[0,1], T[]);
  fflush(fp);
  fclose(fp);
}
#endif


//I don't understand why these cases doesn't works 
//case 1 

#if 0
event flux(t=0){
  foreach_boundary(top)
    fprintf(stdout, "%.9g %g\n", T[0,1], T[]);
  fflush(stdout);
}
#endif

//case 2
#if 0
event flux(t=0){
  char names[80] = "test.dat";
  FILE * fp = fopen (names, "w");
  foreach_boundary(top)
    fprintf(fp, "%.9g %.9g %.9g %.9g\n", x, y, T[0,1], T[]);
  fflush(fp);
  fclose(fp);
}
#endif