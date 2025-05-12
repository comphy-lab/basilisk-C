/**
# Dry cells lead to singularities with vertical diffusion

Currently, there are no checks checking if ```h[]``` is nonzero in the [vertical diffusion module](http://basilisk.fr/src/layered/diffusion.h). Hence, dry cells inevitably cause singularities if viscosity is set to be anything but zero. 
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"

/**
If viscous is set to 1, the program crashes with an Arithmetic exception. 
If viscous is set to 0, the program runs normally.
*/

#define viscous 1

int main(){
  N  = 8;
  nl = 1;
  #if viscous
    nu = 1;
  #else 
    nu = 0.0;
  #endif
  periodic(right);
  run();
}


/** 
We start by initializing 
a simple grid where the top half of the cells
are dry. 
*/

event init(i=0){
  foreach(){
    zb[] = -1.0;
    if ( y >= 0.5)
      zb[] = 0.0;
    foreach_layer(){
      h[]  = -zb[]/nl;
      u.x[] = 0.1;
    }
  }
  
  // Printing the initial h field values
  FILE * fp = fopen("out0.txt", "w");
  foreach(){
    foreach_layer(){
      fprintf(fp, "(x,y) = (%g,%g), h[] = %g\n", x,y,h[]);
    }
  }
  fclose(fp);
}
/**
Output before simulation without viscosity:

![](singularity_diffusion/out0.txt)
*/

/**
And after one timestep
*/
event next_timestep(i=1){
  FILE * fp = fopen("out1.txt", "w");
  foreach(){
    foreach_layer(){
      fprintf(fp, "(x,y) = (%g,%g), h[] = %g\n", x,y,h[]);
    }
  }
  fclose(fp);
}

/**
Output before simulation without viscosity:

![](singularity_diffusion/out1.txt)
*/


