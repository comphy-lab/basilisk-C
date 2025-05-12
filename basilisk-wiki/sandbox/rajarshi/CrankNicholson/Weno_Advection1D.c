/**
# 1D Advection Problem using WENO 3rd Order Schemes 


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

void Weno_Reconstruction(scalar *Xl, face vector *ul, face vector *Wenol){

 double T_S1, T_S2, Beta1, Beta2, w1_tilda, w2_tilda, wsum_tilda, w1, w2, T_plus, T_minus;
 double epsilon,gamma1,gamma2;
 
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
  T_plus = w1*u_S1 + w2*u_S2;
  
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
  T_minus = w1*u_S1 + w2*u_S2;

/* The next step is to use T_plus & T_minus values compute the godunov flux (which is the Riemann problem) */

  if(T_minus <= T_plus){
     if(u.x[]*T_minus <= u.x[]*T_plus)
       Weno.x[] = u.x[]*T_minus;
     else
       Weno.x[] = u.x[]*T_plus;
    }
  else if(T_minus > T_plus){
     if(u.x[]*T_minus <= u.x[]*T_plus)
       Weno.x[] = u.x[]*T_plus;
     else
       Weno.x[] = u.x[]*T_minus;
    }
  
  }  

}

int main(){

 L0 = 4;
 origin(-2);
 init_grid(1<<5);
 int N = 32;
 char Prefix[] = "T-Data__Time-";
 char Suffix[] = ".dat";
 char FName[19];
 scalar T[];
 face vector vel[], F_weno[];
 double t,dt;
 t = 0.;
 dt = L0/(N*4);
 
 foreach(){
   if(x<=0)
     T[] = 2. + x*x;
   else if (x>=0)
     T[] = -2. - x*x;
   }
 foreach_face()
   vel.x[] = 2.;
  
 T[left] = dirichlet(6);
 T[right] = dirichlet(-6);
 boundary({T,vel});
 FILE *Fp;
 while(t<0.5){
   Weno_Reconstruction({T},{vel},{F_weno});
   foreach()
     T[] = T[] - dt*(F_weno.x[1] -F_weno.x[])/(L0/N);
   t = t+dt;
   T[left] = dirichlet(6);
   T[right] = dirichlet(-6);
   boundary({T});
   sprintf(FName,"%s%g%s",Prefix,time,Suffix);
   Fp = fopen(FName,"w");
   foreach()
    fprintf(Fp,"%g %g \n",x,T[]);
   fclose(Fp);
 } 
}

