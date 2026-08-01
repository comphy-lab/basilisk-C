/**

Test for my bilinear interpolation.

Target is x^2+y^2 at x=y=0.5 so target = 0.5

I test interpolating from a 2D grid to this target point.
I build a grid x=y=[-1,1] with NxN points
The grid has increasing number of points, which should result in increasing
accuracy.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "interpolate.h"
#include "constant.h"
 

int main() {
  int N;
  int aN[6] = {10,20,40,60,80,100};
  double aX[9]={0.5, -1.0, -1.0, 1.0, 1.0, -1.0 + 1e-14, 1.0-1e-14, 0.3, 0.3};
  double aY[9]={0.5, -1.0, 1.0, -1.0, 1.0, 0.3, 0.3, -1.0 +1e-14, 1.0-1e-14};
  double Fi;
  double exact, X, Y;
  



  for (int k=0; k<9; ++k){
    X = aX[k];
    Y = aY[k];
    printf("Testing at X=%g,Y=%g\n", X, Y);
    exact = X*X+Y*Y;
    printf("-> true result F(%f,%f)=%f\n",X, Y, X*X+Y*Y);
    for (int i=0; i<6; ++i) {
      N = aN[i];
      double xi[N];
      double yi[N];
      for (int i=0; i<N; ++i){
        xi[i] = -1.0 +2.0*i / (N-1);
        yi[i] = -1.0 +2.0*i / (N-1);
      } 
      double F[N*N]; // f(x)=x^2
      for (int i=0; i<N;++i){
        for (int j=0; j<N;++j){
          F[i*N+j] = xi[i]*xi[i] + yi[j]*yi[j];
        } 
      }
      // interp
      Fi = interp_lin(xi, yi, N, N, X, Y, F); 
      printf("  N=%d ",N);
      printf("  interpolated Fi = %f\n",Fi);

      Fi = interp_lin(xi, yi, N, N, 0.05, 0.05, F); 
      
    }
    printf("\n");
  }
  return 0;
}
