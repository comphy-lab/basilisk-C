/**
#Advection of a 2D Periodic function on a uniform grid - WENO5 + RK4
*/

#include "grid/multigrid.h"
#include "utils.h"
#include "temporal.h"
#define dimension 2
#define BGHOSTS 2
  
void (* Flux_Computation)  (scalar *, vector *, vector *) = Weno_Flux_Order5;  
double (* Time_Marching) (scalar *, vector *, double, double, double) = RungeKutta4;

void Convergence(int depth){

 L0 = 4;
 origin(-2,-2);
 init_grid(1<<depth);

 scalar T[], T_Analytical[];
 vector vel[];

 double u,v;
 u = 3;
 v = 2;
 double t,dt,tfinal;
 t = 0.;
 dt = 1E-05;
 
 foreach()
   T[] = 4.*(cos(pi*x/2+pi*Delta/4)-cos(pi*x/2-pi*Delta/4))*(cos(pi*y/2+pi*Delta/4)-cos(pi*y/2-pi*Delta/4))/(pi*pi*Delta*Delta);
    
 foreach(){
    vel.x[] = u;
    vel.y[] = v;
   }
 
 periodic(left);
 periodic(top);

 boundary({T});
 boundary((scalar *){vel});
 
 FILE *Fp;
 Fp = fopen("Original-Profile.dat","w");
 foreach(){
   fprintf(Fp,"%g %g %g \n",x,y,T[]);
  }
 fclose(Fp);
 
 tfinal = Time_Marching({T},{vel},t,1E-02,dt);

 scalar Error[];
 foreach(){
   T_Analytical[] = 4.*(cos(pi*(x-u*tfinal)/2+pi*Delta/4)-cos(pi*(x-u*tfinal)/2-pi*Delta/4))*(cos(pi*(y-v*tfinal)/2+pi*Delta/4)-cos(pi*(y-v*tfinal)/2-pi*Delta/4))/(pi*pi*Delta*Delta);
   Error[] = T[] - T_Analytical[];
 }

 Fp = fopen("WENO.dat","w");
 foreach()
   fprintf(Fp,"%g %g %g\n",x,y,Error[]);
 fclose(Fp);

 Fp = fopen("Error-Scaling.dat","a");
 fprintf(Fp,"%g %g\n",log10(pow(2,depth)),log10(normf(Error).max));
 fclose(Fp);
}

int main(){
   system("rm -f Error-Scaling.dat");
   for(int depth=4; depth<=8; depth++)
      Convergence(depth);
}


/**
## Results
~~~gnuplot Error Convergence of WENO 5 scheme
set output 'ErrorConvergence.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(L2Norm-Error)'
set grid
f(x) = a*x + b
fit f(x) 'Error-Scaling.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'Error-Scaling.dat' u 1:2 w p t 'Error-Convergence-WENO5', f(x) t title_f(a,b)
~~~

~~~gnuplot Advection 2D - Error Distribution
set output '2DAdvection.png'
set xlabel 'x'
set ylabel 'y'
set zlabel 'E(x,y)'
splot 'WENO.dat'
~~~

*/
