#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "interpolate.h"
#include "constant.h"


int main() {
  int N = 10;
  int Nx = 1000;
  int Ny = 1000;
  double X=0.5;
  double Y=0.5;

  // Generate points in polar coords
  // Data is from 0 to 2 and full circle
  double rho[N];
  double theta[N];
  double fp[N*N];
  for (int i = 0; i < N; ++i) {
    rho[i] = 1.0*i/(N-1);
    for (int j=0; j <N;++j){
      theta[j] = 2.0 * PI * 1.0*j / (1.0*N);
      fp[i*N +j] = rho[i]*rho[i]; // f(x)=x^2
    }  
  }

  // Defining the new grid
  // Cartesian from -1 to 1 in both X and Y
  double xi[Nx];
  double yi[Ny];
  double Fi[Nx*Ny];
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      xi[i] = -1.0 + 2.0 * i / (Nx - 1); // from -1 to 1
      yi[j] = -1.0 + 2.0 * j / (Ny - 1); // from -1 to 1
    }
  }

  // interp
  double r = sqrt(X*X+Y*Y);
  double tt = cos(X/r);
  printf("r %f\n",r);
  printf("tt %f\n",tt);

  double FFi = interp_lin(rho, theta, N, N, r, tt, fp);
  printf("true value is %f\n", r*r);
  printf("interp value is %f\n", FFi);

    return 0;
}
