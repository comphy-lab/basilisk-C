/**
# 1D Reconstruction of a Smooth function using WENO 3rd Order Schemes for spatial discretization 


We want to reconstruct the following cell averaged function to construct the point values at the cell faces

T(x,0) = 0, if x<=-1 ; T(x,0) = ( 1-x^2 )^5, if -1 <= x <= 1 ; T(x,0) = 0, if x>=1

*/

#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "utils.h"

void WENO3_Reconstruction(int depth){

 scalar X[];
 face vector u_weno[], u_simple[], Error[];
 double T_S1, T_S2, Beta1, Beta2;
 double epsilon = 1E-06, gamma1, gamma2, T_minus, w1_tilda, w2_tilda, wsum_tilda, w1, w2, lower, upper;

 foreach(){
   if(x<=-1)
     X[] = 0.;
   else if(x >=-1 && x<=1){
      upper = x+Delta/2.;
      lower  = x-Delta/2.;
      X[] = ((pow(upper,11)/11. - 5.*pow(upper,9)/9. + 10.*pow(upper,7)/7. - 2.*pow(upper,5) + 5.*pow(upper,3)/3. - upper) - (pow(lower,11)/11. - 5.*pow(lower,9)/9. + 10.*pow(lower,7)/7. - 2.*pow(lower,5) + 5.*pow(lower,3)/3. - lower))/Delta;     
     }
   else if(x>=1)
     X[] = 0.;
  }
 X[left] = dirichlet(0);
 X[right] = dirichlet(0);
 boundary({X});

 foreach_face(){
  
  T_S1 = -1.*X[-2]/2. + 3.*X[-1]/2.;
  Beta1 = sq(X[-1]-X[-2]);
  T_S2 = 1.*X[-1]/2. + 1.*X[]/2.;
  Beta2 = sq(X[]-X[-1]);
  gamma1 = 1./3.;
  gamma2 = 2./3.;
  w1_tilda = gamma1/sq(epsilon + Beta1);
  w2_tilda = gamma2/sq(epsilon + Beta2);
  wsum_tilda = w1_tilda + w2_tilda;
  w1 = w1_tilda/wsum_tilda;
  w2 = w2_tilda/wsum_tilda;
  T_minus = w1*T_S1 + w2*T_S2;

  u_weno.x[] = T_minus;
  if(x<-1)
    u_simple.x[] = 0.;
  else if( x>=-1 && x<=1)
    u_simple.x[] = pow((sq(x)-1),5);
  else if(x>1)
    u_simple.x[] = 0.;
 }  

 foreach_face()
     Error.x[] = u_weno.x[] - u_simple.x[];
  
 FILE *fp;
 fp = fopen("Weno.dat","w");
 foreach_face()
    fprintf(fp,"%g %g %g \n",x,u_weno.x[],Error.x[]);
 fclose(fp);
 fp = fopen("ErrorvsGrid.dat","a");
 fprintf(fp,"%g %g \n",log10(pow(2,depth)),log10(normf(Error.x).max));
 fclose(fp);
}


int main(){

 L0 = 4;
 origin(-2);
 int depth;
 system("rm -f ErrorvsGrid.dat");
 for(depth=8;depth<12;depth++){
   init_grid(1<<depth);
   WENO3_Reconstruction(depth);
 }
}


/**
## Results

~~~gnuplot Reconstructed Function
set output 'Interpolation.png'
set xlabel 'X'
set ylabel 'U(X)'
set grid
plot 'Weno.dat' u 1:2 w p t 'WENO-SCHEME'
~~~

~~~gnuplot Error in Interpolation
set output 'InterpolationError.png'
set xlabel 'X'
set ylabel 'Error(X)'
set grid
plot 'Weno.dat' u 1:3 w p t 'Error-WENO-SCHEME'
~~~

~~~gnuplot Error Scaling of the Weno 3rd order reconstruction
set output 'ErrorScaling_Weno3_Smooth.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(L1Norm-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'Error-Scaling-WENO3', f(x) t title_f(a,b)
~~~
*/