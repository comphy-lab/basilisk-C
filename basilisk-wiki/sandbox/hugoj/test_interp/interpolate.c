/* Linear interpolation in C
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "constant.h"

void cart2pol(double x, double y, double rho, double phi) {
    rho = sqrt(x * x + y * y);
    phi = atan2(y, x);
}

void pol2cart(double rho, double phi, double x, double y){
    x = rho * cos(phi);
    y = rho * sin(phi);
}


int find_index_closest(double arr[], int length, double target) {
  /*
   Find the index of the closest element in 'arr'.
   */
    int res = 0;
    int lo = 0, hi = length - 1;

    while (lo < hi) {

        int mid = lo + (hi - lo) / 2;
        // Update res if mid is closer to target
        if (fabs(arr[mid] - target) < fabs(arr[res] - target)) {
            res = mid;
        }

        if (arr[mid] == target) {
            return mid;
        }
        else if (arr[mid] < target) {
            lo = mid + 1;
        }
        else {
            hi = mid - 1;
        } 
        // printf("imid %d, ires %d, arr[res] %f\n",mid,  res, arr[res]);
    }
    return res;
}

void find_ibounds(double arr[], int length, double target, int* low, int* high){
  /* This function finds the bounds in 'arr' that closely match 'target'
   *
   * we assume that 'target' is inside the range of 'arr'
   *
  */

  int closest = find_index_closest(arr, length, target);

  if (closest == length) {
    *low = closest-1;
    *high = closest;
    return;
  }

  if (arr[closest]-target < 0) {
    *low = closest;
    *high = closest+1;
  }
  else {
    *low = closest-1;
    *high = closest;
  }
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
  
  double dx = xi - x[i0];
  double dy = yi - y[j0];
  double deltaX = x[i1] - x[i0];
  double deltaY = y[j1] - y[j0];
  double f11 = F[i1*Nx+j1];
  double f10 = F[i1*Nx+j0];
  double f01 = F[i0*Nx+j1];
  double f00 = F[i0*Nx+j0];
  
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
