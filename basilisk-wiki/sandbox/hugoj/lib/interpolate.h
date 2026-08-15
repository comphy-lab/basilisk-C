/* Linear interpolation in C
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>


void find_ibounds(const double *arr, int length, double target, int* low, int* high){
  /* This function finds the bounds in 'arr' that closely match 'target'
   *
   * we assume that 'target' is inside the range of 'arr'
   *
  */


  int lo = 0;
  int hi = length - 1;

  while (hi - lo > 1) {
      int mid = (lo + hi) / 2;

      if (arr[mid] <= target)
          lo = mid;
      else
          hi = mid;
  }

  *low = lo;
  *high = hi;
}


double interp_lin(double x[], double y[], int Nx, int Ny, double xi, double yi, double F[]) {
  /* Interpolate linearly F on grid (x,y) at position xi,yi 
   *
   * Using eq 98 of https://pages.hmc.edu/ruye/MachineLearning/lectures/ch7/node7.html
   * */
  
  int j0, j1, i0, i1;
  double Fi;
  find_ibounds(x, Nx, xi, &i0, &i1);
  find_ibounds(y, Ny, yi, &j0, &j1);
  
  // assert(i0 >= 0);
  // assert(i1 < Nx);
  // assert(i0 < i1);
  // assert(j0 >= 0);
  // assert(j1 < Ny);
  // assert(j0 < j1);

  double dx = xi - x[i0];
  double dy = yi - y[j0];
  double deltaX = x[i1] - x[i0];
  double deltaY = y[j1] - y[j0];
  double f11 = F[i1*Ny+j1];
  double f10 = F[i1*Ny+j0];
  double f01 = F[i0*Ny+j1];
  double f00 = F[i0*Ny+j0];
  
  // printf("j0 %d, j1 %d, i0 %d, i1 %d\n", j0, j1, i0, i1);
  // printf("dx %f, dy %f\n", dx, dy);
  // printf("deltaY %f, deltaX %f\n", deltaY, deltaX);
  // printf("f11 %f, f10 %f, f01 %f, f00 %f\n", f11, f10, f01, f00);
  Fi = ( dx/deltaX * dy/deltaY * f11 +
      dy/deltaY * (1-dx/deltaX) * f01 +
      dx/deltaX * (1-dy/deltaY) * f10 +
      (1-dx/deltaX-dy/deltaY+ dx/deltaX*dy/deltaY) * f00 );
  
  return Fi;
}
