/**
#Advection of a 1D Compact function on a uniform grid ( WENO 5 + RK4 )
*/

#include "grid/multigrid1D.h"
#include "utils.h"
#include "temporal.h"

#define dimension 1
#define velocity 4
#define BGHOSTS 2

double (* Time_Marching) (scalar *, vector *, double, double, double) = RungeKutta3;   /* options : Temporal1, RungeKutta3, RungeKutta4    */
void (* Flux_Computation)  (scalar *, vector *, vector *) = Weno_Flux_Order5;  /* options : Minmod_Flux_Order2, Weno_Flux_Order3, Weno_Flux_Order4, Weno_Flux_Order5 */


double compact_polynomial (double xu, double xl, double Delta){
  if(sq((xu+xl)/2.) <= 1)
     return ((1./Delta)*((xu - (7./3.)*pow(xu,3) + (21./5.)*pow(xu,5) - 5.*pow(xu,7) + (35./9.)*pow(xu,9) - (21./11.)*pow(xu,11) + (7./13.)*pow(xu,13) - (1./15.)*pow(xu,15)) - (xl - (7./3.)*pow(xl,3) + (21./5.)*pow(xl,5) - 5.*pow(xl,7) + (35./9.)*pow(xl,9) - (21./11.)*pow(xl,11) + (7./13.)*pow(xl,13) - (1./15.)*pow(xl,15))));
  else
     return(0);
}

void convergence(int depth){
   
 scalar s[];
 vector u[];
 face vector uf[], Flux_value[];
 
 foreach()
   s[] = compact_polynomial(x+Delta/2.,x-Delta/2.,Delta);

 foreach()
   u.x[] = velocity;
  
 s[left] = neumann(0);
 s[right] = neumann(0);
 boundary({s});
 boundary((scalar *){u});

 FILE *Fp;
 Fp = fopen("Initial_Condition.dat","w");
 foreach()
   fprintf(Fp,"%g %g\n",x,s[]);
 fclose(Fp);

 double dt, Tf, t_now; 
 dt = 1E-04;
 Tf = 0.1;
 t_now = Time_Marching({s},{u},0,Tf,dt);

 scalar S_Analytical[], Error[];
 foreach()
    S_Analytical[] = compact_polynomial(x-velocity*t_now+Delta/2.,x-velocity*t_now-Delta/2.,Delta);
 boundary({S_Analytical});
 
 Fp = fopen("Solution.dat","w");
 foreach()
    fprintf(Fp,"%g %g \n",x,s[]);
 fclose(Fp);

 Fp = fopen("Error-Distribution.dat","w");
 foreach(){
   Error[] = s[] - S_Analytical[];
   fprintf(Fp,"%g %g \n",x, Error[]);
 }
 fclose(Fp);

 Fp = fopen("ErrorvsGrid.dat","a");
 fprintf(Fp,"%g %g\n",log10(grid->tn),log10(normf(Error).max));
 fclose(Fp);
}


int main(){

 system("rm -f ErrorvsGrid.dat");
 L0 = 4;
 origin(-2);
 for(int depth = 7; depth <= 10; depth++){
    init_grid(1<<depth);
    convergence(depth);
   }
}
 
/**
## Results

~~~gnuplot Final-Solution
set output 'Final-Solution.png'
set xlabel 'X'
set ylabel 'T(X)'
set grid
plot 'Solution.dat' u 1:2 w p t 'Final Solution'
~~~

~~~gnuplot Error distribution 1D - Uniform Grid
set output '1DError.png'
set xlabel 'X'
set ylabel 'Error(X)'
set grid
plot 'Error-Distribution.dat' u 1:2 w p t 'Error-FullDomain'
~~~

~~~gnuplot Error Convergence of the Advection operator
set output 'Error-Convergence.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(L2Norm-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'Error-Convergence', f(x) t title_f(a,b)
~~~

*/