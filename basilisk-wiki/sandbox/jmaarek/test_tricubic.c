/**
# Tricubic interpolation test case


This is a simple test case for tricubic interpolation, we use a 3D-sinusoidal
function as reference and interpolate it on a 4 times finer grid.

We also do a triquadratic interpolation and correct it using the interpolation
coefficient of the center value of the 3x3 stencil.
*/

#define BGHOSTS 2
#define BICUBIC 1
#include "grid/octree.h"
#include "alex_functions.h"


#define Pi 3.14159265358979323846
double ref_function(double x, double y, double z){
  // return (x-0.2)*(x-0.4)*(x-0.45)*(x+0.1);
  return cos(6.2*Pi*x+60)*cos(6.2*Pi*y+60)*cos(6.2*Pi*z+60)+Pi;
}

scalar s[];
s[left] = ref_function(x,y,z);
s[right] = ref_function(x,y,z);
s[top] = ref_function(x,y,z);
s[bottom] = ref_function(x,y,z);
s[front] = ref_function(x,y,z);
s[back] = ref_function(x,y,z);


int main(){
  for(int j = 0; j<=4; j++){
    int MAXLEVEL = 4+j;
    init_grid (1 << MAXLEVEL);
    origin(-L0/2.,-L0/2.);

    foreach(){
      s[] = ref_function(x,y,z);
    }
    boundary({s});

    double tn = 0.0;
    double error_tricubic = 0.0;
    double error_triquadratic = 0.0;

    foreach_vertex(reduction(+:tn) reduction(+:error_triquadratic) reduction(+:error_tricubic)){
      coord p_interp = {-0.5, -0.5, -0.5};

      if((point.i-1) > (1 << grid-> depth)){
        p_interp.x = 0.5;
      }

      if((point.j-1) > (1 << grid-> depth)){
        p_interp.y = 0.5;
      }

      if((point.k-1) > (1 << grid-> depth)){
        p_interp.z = 0.5;
      }

      int Stencil[2] = {-1,-1};
      //currently there is a bug on the +z boundary. Working on fixing this
      if ((p_interp.z == -0.5)){
        error_triquadratic += fabs(mytriquadratic(point , s, p_interp)-ref_function(x,y,z));
        error_tricubic += fabs(mytricubic(point , s, Stencil, p_interp)-ref_function(x,y,z));
      tn += 1;}


    }


/**
We output the L1-error.
*/
    fprintf(stdout, "%d %g %g\n", 1 << MAXLEVEL, error_tricubic/tn, error_triquadratic/tn);
  }
}
