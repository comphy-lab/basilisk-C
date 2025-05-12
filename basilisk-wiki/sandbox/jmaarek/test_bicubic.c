/**
# Bicubic interpolation test case


This is a simple test case for bicubic interpolation, we use a 2D-sinusoidal
function as reference and interpolate it on a 4 times finer grid.

We also do a biquadratic interpolation and correct it using the interpolation
coefficient of the center value of the 3x3 stencil.
*/

#define BGHOSTS 2
#define BICUBIC 1
#include "alex_functions.h"


#define Pi 3.14159265358979323846
double ref_function(double x, double y){
  // return (x-0.2)*(x-0.4)*(x-0.45)*(x+0.1);
  return cos(6.2*Pi*x+60)*cos(6.2*Pi*y+60);
}

scalar s[];
s[left] = ref_function(x,y);
s[right] = ref_function(x,y);
s[top] = ref_function(x,y);
s[bottom] = ref_function(x,y);


int main(){
  for(int j = 0; j<=4; j++){
    int MAXLEVEL = 4+j;
    init_grid (1 << MAXLEVEL);
    origin(-L0/2.,-L0/2.);

    foreach(){
      s[] = ref_function(x,y);
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

      /*if((point.k-1) > (1 << grid-> depth)){
        p_interp.z = 0.5;
      }*/

      int Stencil[2] = {-1,-1};

      /*if ((p_interp.z == -0.5)){
        error_triquadratic += fabs(mytriquadratic(point , s, p_interp)-ref_function(x,y));
        error_tricubic += fabs(mytricubic(point , s, Stencil, p_interp)-ref_function(x,y));
      tn += 1;}

      else{
        //error_triquadratic += fabs(mytriquadratic(point , s, p_interp)-ref_function(x,y));
        fprintf(stderr, "%g %g %g ", p_interp.x, p_interp.y, p_interp.z);
        double q = s[2,2,2];
        fprintf(stderr, "%g\n", q);

        tn += 1;
      }*/
	error_triquadratic += fabs(mybiquadratic(point , s, p_interp,0)-ref_function(x,y));
        error_tricubic += fabs(bicubic(point , s, Stencil, p_interp,0)-ref_function(x,y));
      tn += 1;


    }


/**
We output the L1-error.
*/
    fprintf(stdout, "%d %g %g\n", 1 << MAXLEVEL, error_tricubic/tn, error_triquadratic/tn);
  }
}