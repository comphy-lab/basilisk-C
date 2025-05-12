/**
# 1D Advection of a Discontinuous function using WENO 5th Order Schemes for spatial discretization and RK3 temporal discretization 

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

 double T_S1,T_S2,T_S3,Beta1,Beta2,Beta3,w1_tilda,w2_tilda,w3_tilda,wsum_tilda,w1,w2,w3,T_plus,T_minus;
 double epsilon,gamma1,gamma2,gamma3;
 epsilon = 1E-06;
 
 scalar X = Xl[0];
 face vector u = ul[0];
 face vector Weno = Wenol[0];
 face vector gradL[];
 
 foreach_face()
  gradL.x[] = (X[]-X[-1])/(Delta);

 foreach_face(){

  if(u.x[] > 0){

     T_S1 =   1.*(X[-2] - Delta*gradL.x[-2])/3. - 7.*X[-2]/6. + 11.*X[-1]/6.;
     Beta1 = 13.*sq(X[-2] - Delta*gradL.x[-2] - 2.*X[-2] + X[-1])/12. + sq(X[-2] - Delta*gradL.x[-2] - 4.*X[-2] + 3.*X[-1])/4.;
     T_S2 =  -1.*X[-2]/6. + 5.*X[-1]/6. + 1.*X[]/3.;
     Beta2 = 13.*sq(X[-2] - 2.*X[-1] + X[])/12. + sq(X[-2]-X[])/4.;
     T_S3 =   1.*X[-1]/3. + 5.*X[]/6. - 1.*X[1]/6.;
     Beta3 = 13.*sq(X[-1] - 2.*X[] + X[1])/12. + sq(3.*X[-1] - 4.*X[] + X[1])/4.;
     gamma1 = 1./10.;
     gamma2 = 3./5.;
     gamma3 = 3./10.;
     w1_tilda = gamma1/sq(epsilon + Beta1);
     w2_tilda = gamma2/sq(epsilon + Beta2);
     w3_tilda = gamma3/sq(epsilon + Beta3);
     wsum_tilda = w1_tilda + w2_tilda + w3_tilda;
     w1 = w1_tilda/wsum_tilda;
     w2 = w2_tilda/wsum_tilda;
     w3 = w3_tilda/wsum_tilda;
     T_minus = w1*T_S1 + w2*T_S2 + w3*T_S3; 
     Weno.x[] = u.x[]*T_minus;

   }
  
 else if(u.x[] <= 0 ){
 
     T_S1 =  -1.*X[-2]/6. + 5.*X[-1]/6. + 1.*X[]/3.;
     Beta1 = 13.*sq(X[-2] - 2.*X[-1] + X[])/12. + sq(X[-2] - 4.*X[-1] + 3.*X[])/4.;
     T_S2 =   1.*X[-1]/3. + 5.*X[]/6. - 1.*X[1]/6.;
     Beta2 = 13.*sq(X[-1] - 2.*X[] + X[1])/12. + sq(X[-1]-X[1])/4.;
     T_S3 =  11.*X[]/6. - 7.*X[1]/6. + 1.*X[2]/3.;
     Beta3 = 13.*sq(X[] - 2.*X[1] + X[2])/12. + sq(3.*X[] - 4.*X[1] + X[2])/4.;
     gamma1 = 3./10.;
     gamma2 = 3./5.;
     gamma3 = 1./10.;
     w1_tilda = gamma1/sq(epsilon + Beta1);
     w2_tilda = gamma2/sq(epsilon + Beta2);
     w3_tilda = gamma3/sq(epsilon + Beta3);
     wsum_tilda = w1_tilda + w2_tilda + w3_tilda;
     w1 = w1_tilda/wsum_tilda;
     w2 = w2_tilda/wsum_tilda;
     w3 = w3_tilda/wsum_tilda;
     T_plus = w1*T_S1 + w2*T_S2 + w3*T_S3;
     Weno.x[] = u.x[]*T_plus;

   }
  }
}


int main(){

 L0 = 4;
 origin(-2);
 init_grid(1<<8);
 int N = 256;

 scalar T[],T1[],T2[];
 face vector vel[], F_weno[];
 double t,dt,velocity;

 velocity = 2.;
 t = 0.;
 dt = L0/(N*20*10);
 
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
   vel.x[] = velocity;
  
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
   if (x<=-1+velocity*t)
     T_Analytical[] = 2.;
   else if(x >=-1+velocity*t && x<=velocity*t)
     T_Analytical[] = 2. - sq(x - velocity*t + 1);
   else if (x>=velocity*t  && x<=1+velocity*t)
     T_Analytical[] = sq(x - velocity*t - 1) - 2.;
   else if(x>=1+velocity*t)
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
