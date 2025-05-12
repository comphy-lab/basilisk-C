/**
# 1D Advection of a Smooth function using WENO 4th Order Schemes for spatial discretization and first order temporal discretization 


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

 double T_S1,T_S2,T_S3,Beta1,Beta2,Beta3,w1_tildaP,w2_tildaP,w3_tildaP,w1_tildaM,w2_tildaM,w3_tildaM,wsum_tildaP,wsum_tildaM, w1P,w1M,w2P,w2M,w3P,w3M,T_plus,T_minus;
 double gamma1p,gamma1m,gamma2p,gamma2m,gamma3p,gamma3m,SigmaP,SigmaM,TsplitP,TsplitM;
 double epsilon,gamma1,gamma2,gamma3;
 epsilon = 1E-06;
 
 scalar X = Xl[0];
 face vector u = ul[0];
 face vector Weno = Wenol[0];

 foreach_face(){
  
  if(u.x[] < 0){

    T_S1 =  -1.*X[-2]/6. + 5.*X[-1]/6. + 1.*X[]/3.;
    Beta1 = 13.*sq(X[-2]-2.*X[-1]+X[])/12. + sq(X[]-X[-2])/4.;
    T_S2 =  3.*X[-1]/2. - X[-2]/2.;
    Beta2 = sq(X[-1]-X[-2]);
    T_S3 =  3.*X[]/2. - X[1]/2.;
    Beta3 = sq(X[1]-X[]);

    gamma1 = 1.;
    gamma2 = -1./6.;
    gamma3 = 1./6.;
    gamma1p = 2.*gamma1;
    gamma1m = gamma1; 
    gamma2p = -1.*gamma2;
    gamma2m = -2.*gamma2;
    gamma3p = 2.*gamma3;
    gamma3m = gamma3;
    SigmaP = gamma1p + gamma2p + gamma3p;
    SigmaM = gamma1m + gamma2m + gamma3m;
    gamma1p /= SigmaP;
    gamma2p /= SigmaP;
    gamma3p /= SigmaP;
    gamma1m /= SigmaM;
    gamma2m /= SigmaM;
    gamma3m /= SigmaM;
  
    w1_tildaP = gamma1p/sq(epsilon + Beta1);
    w1_tildaM = gamma1m/sq(epsilon + Beta1);
    w2_tildaP = gamma2p/sq(epsilon + Beta2);
    w2_tildaM = gamma2m/sq(epsilon + Beta2);
    w3_tildaP = gamma3p/sq(epsilon + Beta3);
    w3_tildaM = gamma3m/sq(epsilon + Beta3);

    wsum_tildaP = w1_tildaP + w2_tildaP + w3_tildaP;
    wsum_tildaM = w1_tildaM + w2_tildaM + w3_tildaM;
  
    w1P = w1_tildaP/wsum_tildaP;
    w2P = w2_tildaP/wsum_tildaP;
    w3P = w3_tildaP/wsum_tildaP;
    w1M = w1_tildaM/wsum_tildaM;
    w2M = w2_tildaM/wsum_tildaM;
    w3M = w3_tildaM/wsum_tildaM;

    TsplitP = w1P*T_S1 + w2P*T_S2 + w3P*T_S3;
    TsplitM = w1M*T_S1 + w2M*T_S2 + w3M*T_S3;  
    T_plus = SigmaP*TsplitP - SigmaM*TsplitM;  
     
    Weno.x[] = u.x[]*T_plus;
   }
  
 else if(u.x[] >=0 ){
 
    T_S1 =  X[-1]/3. + 5.*X[]/6. -X[1]/6.;
    Beta1 = 13.*sq(X[-1]-2.*X[]+X[1])/12. + sq(X[1]-X[-1])/4.; 
    T_S2 =  3.*X[]/2. - X[1]/2.;
    Beta2 = sq(X[1]-X[]);
    T_S3 =  3.*X[-1]/2. - X[-2]/2.;
    Beta3 = sq(X[-1]-X[-2]);
  
    gamma1 = 1.;
    gamma2 = -1./6.;
    gamma3 = 1./6.;
    gamma1p = 2.*gamma1;
    gamma1m = gamma1; 
    gamma2p = -1.*gamma2;
    gamma2m = -2.*gamma2;
    gamma3p = 2.*gamma3;
    gamma3m = gamma3;
    SigmaP = gamma1p + gamma2p + gamma3p;
    SigmaM = gamma1m + gamma2m + gamma3m;
    gamma1p /= SigmaP;
    gamma2p /= SigmaP;
    gamma3p /= SigmaP;
    gamma1m /= SigmaM;
    gamma2m /= SigmaM;
    gamma3m /= SigmaM;

    w1_tildaP = gamma1p/sq(epsilon + Beta1);
    w1_tildaM = gamma1m/sq(epsilon + Beta1);
    w2_tildaP = gamma2p/sq(epsilon + Beta2);
    w2_tildaM = gamma2m/sq(epsilon + Beta2);
    w3_tildaP = gamma3p/sq(epsilon + Beta3);
    w3_tildaM = gamma3m/sq(epsilon + Beta3);

    wsum_tildaP = w1_tildaP + w2_tildaP + w3_tildaP;
    wsum_tildaM = w1_tildaM + w2_tildaM + w3_tildaM;
  
    w1P = w1_tildaP/wsum_tildaP;
    w2P = w2_tildaP/wsum_tildaP;
    w3P = w3_tildaP/wsum_tildaP;
    w1M = w1_tildaM/wsum_tildaM;
    w2M = w2_tildaM/wsum_tildaM;
    w3M = w3_tildaM/wsum_tildaM;

  
    TsplitP = w1P*T_S1 + w2P*T_S2 + w3P*T_S3;
    TsplitM = w1M*T_S1 + w2M*T_S2 + w3M*T_S3;
    T_minus = SigmaP*TsplitP - SigmaM*TsplitM;  
  
    Weno.x[] = u.x[]*T_minus;

   }
  }
}


int main(){

 L0 = 4;
 origin(-2);
 init_grid(1<<9);
 int N = 512;

 scalar T[];
 face vector vel[], F_weno[];
 double t,dt;
 t = 0.;
 dt = L0/(N*40*20);
 
 foreach(){
 
   if(x<=-1)
     T[] = 0.;
   else if(x >=-1 && x<=1)
     T[] = -pow((sq(x)-1),5);
   else if(x>=1)
     T[] = 0.;
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
     T_Analytical[] = 0.;
   else if (x >= -1+2.*t && x <= 1+2.*t)
     T_Analytical[] = -pow((sq(x-2.*t)-1),5);
   else if (x >= 1+2.*t)
     T_Analytical[] = 0.;
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
