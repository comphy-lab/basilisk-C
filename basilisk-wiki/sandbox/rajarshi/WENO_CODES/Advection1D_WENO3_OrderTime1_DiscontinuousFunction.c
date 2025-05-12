/**
# 1D Advection of smooth function using WENO 3rd Order Schemes & First order time scheme


We want to solve the following equation.
$$
\partial_t\mathbf{T} + \partial_x\mathbf{(uT)} = 0
$$

u(x,t) = 2
\break

T(x,0) = 2 + x^2 , if x<=0 and T(x,0) = -2 -x^2, if x>0

*/

#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "utils.h"

void Weno_Reconstruction(scalar *Xl, vector *ul, vector *Wenol){

 double T_S1, T_S2, Beta1, Beta2, w1_tilda, w2_tilda, wsum_tilda, w1, w2, T_plus, T_minus;
 double epsilon,gamma1,gamma2;
 epsilon = 1E-06;
 
 scalar X = Xl[0];
 face vector u = ul[0];
 face vector Weno = Wenol[0];

 foreach_face(){
  
  T_S1 = 1.*X[-1]/2. + 1.*X[]/2.;
  Beta1 = sq(X[]-X[-1]);
  T_S2 = 2.*X[]/2. - 1.*X[-1]/2.;
  Beta2 = sq(X[1]-X[]);
  gamma1 = 2./3.;
  gamma2 = 1./3.;
  w1_tilda = gamma1/sq(epsilon + Beta1);
  w2_tilda = gamma2/sq(epsilon + Beta2);
  wsum_tilda = w1_tilda + w2_tilda;
  w1 = w1_tilda/wsum_tilda;
  w2 = w2_tilda/wsum_tilda;
  T_plus = w1*T_S1 + w2*T_S2;
  
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

/* The next step is to use T_plus & T_minus values compute the godunov flux (which is the Riemann problem) */

  if(u.x[] >= 0)
       Weno.x[] = u.x[]*T_minus;
  else if(u.x[] < 0)
       Weno.x[] = u.x[]*T_plus;
  }  
}

int main(){

 L0 = 4;
 origin(-2);
 init_grid(1<<8);
 int N = 256;
 /*
 char Prefix[] = "T-Data__Time-";
 char Suffix[] = ".dat";
 char FName[19];
 */
 scalar T[];
 face vector vel[], F_weno[];
 double t,dt;
 t = 0.;
 dt = L0/(N*40*20);
 
 foreach(){
   if(x<=-1)
     T[] = 2.;
   else if(x >=-1 && x<=0)
     T[] = 2. - sq(x+1);
   else if (x>=0  && x<=1)
     T[] = sq(x-1) - 2.;
   else if(x>=1)
     T[] = -2.;
   }

 foreach_face()
   vel.x[] = 2.;
  
 T[left] = neumann(0);
 T[right] = neumann(0);
 boundary({T,vel});

 FILE *Fp;
 Fp = fopen("WENO-Original.dat","w");
 foreach()
   fprintf(Fp,"%g %g\n",x,T[]);
 fclose(Fp);


 while(t<0.45){
   
   Weno_Reconstruction({T},{vel},{F_weno});
   foreach()
     T[] = T[] - dt*(F_weno.x[1] - F_weno.x[])/(Delta);
   T[left] = neumann(0);
   T[right] = neumann(0);
   boundary({T});
   
   t = t+dt;
 }

 t -= dt;
 scalar T_Analytical[];
 foreach(){
   if (x<=-1+2.*t)
     T_Analytical[] = 2;
   else if (x >= -1+2.*t && x <= 2.*t)
     T_Analytical[] = 2. - sq(x-2.*t+1);
   else if (x >= 2.*t && x <= 1+2.*t)
     T_Analytical[] = sq(x-2.*t-1) - 2;
   else if (x >= 1+2.*t)
     T_Analytical[] = -2.;
   }
 
 Fp = fopen("WENO.dat","w");
 foreach()
   fprintf(Fp,"%g %g %g %g \n",x,T[], T_Analytical[], T[]-T_Analytical[]);
 fclose(Fp);
}

/**
## Results

~~~gnuplot Advection 1D
set output 'Advecting_Tracer.png'
set xlabel 'X'
set ylabel 'T(X)'
set grid
plot 'WENO.dat' u 1:2 w p t 'Advected-Solution', 'WENO-Original.dat' u 1:2 w p t 'Initial-Condition'
~~~

~~~gnuplot Advection 1D
set output 'Observing_Diffusion.png'
set xlabel 'X'
set ylabel 'T(X)'
set grid
plot 'WENO.dat' u 1:2 w p t 'Numerical-Solution', 'WENO.dat' u 1:3 w p t 'Analytical-Solution'
~~~


~~~gnuplot Error 1D
set output '1DError.png'
set xlabel 'X'
set ylabel 'Error(X)'
set grid
plot 'WENO.dat' u 1:4 w p t 'Error'
~~~
*/
