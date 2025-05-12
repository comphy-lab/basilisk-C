/**
# Weno Reconstruction
*/

#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "utils.h"

int main(){

 L0 = 4;
 origin(-2);
 init_grid(1<<5);

 scalar u[];
 vertex scalar u_weno[], u_simple[];
 double u_S1, u_S2, Beta1, Beta2, w1_tilda, w2_tilda, wsum_tilda, w1, w2;
 double epsilon = 1E-06, gamma1 = 1./3., gamma2 = 2./3.;

 foreach(){
   if(x<=0)
     u[] = 2. + x*x;
   else if (x>=0)
     u[] = -2. - x*x;
  }
 u[left] = dirichlet(6);
 u[right] = dirichlet(-6);
 boundary({u});

 foreach_vertex(){
  u_S1 =  3.*u[]/2. - 1.*u[-1]/2.;
  u_S2 =  1.*u[1]/2. + 1.*u[]/2.;
  
  Beta1 = 4.*sq(u[]-u[-1]);
  Beta2 = 4.*sq(u[1]-u[]);
  
  w1_tilda = gamma1/sq(epsilon + Beta1);
  w2_tilda = gamma2/sq(epsilon + Beta2);
  wsum_tilda = w1_tilda + w2_tilda;

  w1 = w1_tilda/wsum_tilda;
  w2 = w2_tilda/wsum_tilda;
  
  u_weno[] = w1*u_S1 + w2*u_S2;
  u_simple[] = gamma1*u_S1 + gamma2*u_S2;
 }

 FILE *fp;
 fp = fopen("Weno.dat","w");
 foreach_vertex()
    fprintf(fp,"%g %g %g \n",x,u_weno[],u_simple[]);

}

/**
## Results

~~~gnuplot Fluctuations in Interpolation
set output 'InterpolationFluctuations.png'
set xlabel 'X'
set ylabel 'U(X)'
set grid
plot 'Weno.dat' u 1:2 w p t 'WENO-SCHEME', 'Weno.dat' u 1:3 w p t 'Direct Polynomial Interpolation'
~~~

*/