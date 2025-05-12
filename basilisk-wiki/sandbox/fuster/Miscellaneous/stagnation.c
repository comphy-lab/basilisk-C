#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

vector upot[];
int Nx = 64;

double tend = 0.1;

double Reinv = 1.e-3;
u.t[top] = dirichlet(0.);
u.n[top] = dirichlet(-1.);

u.t[bottom] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);


p[right] = dirichlet(0.);

int main()
{

  // coordinates of lower-left corner
  origin (0., 0.);
  // number of grid points
  init_grid (Nx);
 
  L0 = 5;
  // maximum timestep
  //DT = 0.1;
  // CFL number
  //CFL = 0.8;
  //
  const face vector muc[] = {Reinv,Reinv};
  mu = muc;

  run();
}

event init ( t = 0) {


  foreach () {
    upot.x[] = x;
    upot.y[] = -y;
  }

}

event outputlog ( t += 0.0001;  t <= tend) {
  printf("%g %g \n", t, dt);
}

event output ( t = end ) {

  double xmax = L0/Nx;

  FILE * fp = fopen("output_left.dat", "w");
  FILE * fp2 = fopen("output_bottom.dat", "w");
  FILE * fpKE = fopen("output_KE.dat", "w");
  foreach () {
    if (x < xmax)
      fprintf(fp, "%g %g %g \n", y - Delta/2. , (upot.x[0] + upot.x[-1])/2., (upot.y[] + upot.y[-1])/2);

    if (y < xmax)
      fprintf(fp2, "%g %g %g %g %g \n", x, upot.x[], u.x[], upot.y[], u.y[] );

    fprintf(fpKE, "%g %g %g %g \n", x, y, 0.5*(sq(u.x[]) + sq(u.y[])), 0.5*(sq(upot.x[]) + sq(upot.y[])) );   
  }
  fclose(fp);
  fclose(fp2);
  fclose(fpKE);


}