#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "interpolate.h"
#include "constant.h"
 

int main() {
  int N;
  int aN[6] = {10,20,40,60,80,100};
  double X=0.5;
  double Y=0.5;
  double Fi;
    
  printf("true result F(%f,%f)=%f\n",X, Y, X*X+Y*Y);

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
        F[i*N+j] = xi[i]*xi[i] + yi[i]*yi[i];
      } 
    }
    // interp
    Fi = interp_lin(xi, yi, N, N, X, Y, F); 
    printf("N=%d ",N);
    printf("interpolated Fi = %f\n",Fi);

    Fi = interp_lin(xi, yi, N, N, 0.05, 0.05, F); 

  }

  return 0;
}
