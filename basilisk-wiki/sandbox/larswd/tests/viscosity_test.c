/** 

# Demonstration of bug in vertical viscosity module

The vertical viscosity module does not handle dry cells. To illustrate, we initialize 
a uniform river flow in a rectangular channel

*/

#include "grid/multigrid.h"
#include "layered/hydro.h"

double W; 
double H0;
double U = 0.1 [1,-1];
#define CHECK_FOR_DRY_CELLS 0

int main(){
  // Channel dimensions
  L0   = 100;
  Y0   = L0/2.;
  W    = 20;
  H0   = 10;
  
  // Mesh parameters
  N    = 8;
  nl   = 2;
  
  // Turn on viscosity
  #if CHECK_FOR_DRY_CELLS
    double _nu = 1e-6;
  #else
    nu = 1e-6;
  #endif
  run();
}

//
event init(i=0){
  foreach(){
    // Channel bottom
    if (fabs(y) > W/2.){
      zb[] = 0.0;
    } else {
      zb[] = -H0;
    }
    
    // Uniform flow which is zero outside channel
    foreach_layer(){
      h[] = -zb[]/nl;
      u.x[] = U*h[]/H0;
    }
    
  }
  FILE *   fp = fopen("out.txt", "w");
  foreach(){
    foreach_layer(){
      fprintf(fp, "u.x[] = %g\n", u.x[]);
    }
  }
  fclose(fp);
}

event print_change(i=1){
  FILE *   fp = fopen("out.txt", "w");
  foreach(){
    foreach_layer(){
      fprintf(fp, "u.x[] = %g\n", u.x[]);
    }
  }
  fclose(fp);
}

