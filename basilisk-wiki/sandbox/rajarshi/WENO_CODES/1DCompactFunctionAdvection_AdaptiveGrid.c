/**
#ADVECTION OF A COMPACT TRACER ON AN ADAPTIVE GRID ( WENO-5 + RK-4 + PROLONGATION-5 )
*/

#include "grid/bitree.h"
#include "utils.h"
#include "temporal.h"

#define dimension 1
#define BGHOSTS 2
#define velocity 4

void (* Flux_Computation)  (scalar *, vector *, vector *) = Weno_Flux_Order5;  /* options : Minmod_Flux_Order2, Weno_Flux_Order3, Weno_Flux_Order4, Weno_Flux_Order5*/
double (* Time_Marching) (scalar *, vector *, double, double, double) = RungeKutta4;   /* options : Temporal1, RungeKutta3, RungeKutta4    */

double compact_polynomial (double xu, double xl, double Delta){
  if(sq((xu+xl)/2.) <= 1)
     return ((1./Delta)*((xu - (7./3.)*pow(xu,3) + (21./5.)*pow(xu,5) - 5.*pow(xu,7) + (35./9.)*pow(xu,9) - (21./11.)*pow(xu,11) + (7./13.)*pow(xu,13) - (1./15.)*pow(xu,15)) - (xl - (7./3.)*pow(xl,3) + (21./5.)*pow(xl,5) - 5.*pow(xl,7) + (35./9.)*pow(xl,9) - (21./11.)*pow(xl,11) + (7./13.)*pow(xl,13) - (1./15.)*pow(xl,15))));
  else
     return(0);
}

void convergence (double cmax) {

 L0 = 4;
 origin(-2);
 init_grid(1<<6);

 scalar s[];
 astats statistics; 
 s[left] = neumann(0);
 s[right] = neumann(0);
 s.prolongation = refine_weno;
 
 do { 
       foreach()
          s[] = compact_polynomial(x+Delta/2.,x-Delta/2.,Delta);
       boundary({s}); 
       statistics = adapt_wavelet ({s}, (double []){cmax}, maxlevel = 13, {s});
    } while (statistics.nf > 1);

 vector u[];
 foreach()
    u.x[] = velocity;
 u.n[left] = neumann(0);
 u.n[right] = neumann(0);
 for(scalar v in {u})
   v.prolongation = refine_weno;
 boundary((scalar *){u});
  
 FILE *Fp;
 double dt, Tf, t_now;
 dt = 5E-05;
 Tf = 0.01;
 t_now = 0.;

 while(t_now < Tf){
    t_now = Time_Marching({s}, {u}, t_now, t_now+dt, dt);
    statistics = adapt_wavelet ({s}, (double []){cmax}, maxlevel = 13);  
  }

 scalar s_analytical[], Error[];
 foreach()
   s_analytical[] =  compact_polynomial(x-velocity*t_now+Delta/2.,x-velocity*t_now-Delta/2.,Delta);
 boundary({s_analytical});
 Fp = fopen("Error-Distribution.dat","w");
 foreach(){
   Error[] = s[] - s_analytical[];
   fprintf(Fp,"%g %g %g %g \n",x, s[], s_analytical[], Error[]);
 }
 fclose(Fp);
 
 Fp = fopen("ErrorvsGrid.dat","a");
 fprintf(Fp,"%g %g\n",log10(grid->tn),log10(normf(Error).max));
 fclose(Fp);

 free_grid(); 
}

int main(){

 system("rm -f ErrorvsGrid.dat");
 for(double cmax = 0.01 ; cmax >= 0.0005; cmax/=2)
    convergence(cmax);

}

/**
## Results

~~~gnuplot Final-Solution-Comparison
set output 'Final-Solution-Comparison.png'
set xlabel 'X'
set ylabel 'T(X)'
set grid
plot 'Error-Distribution.dat' u 1:2 w p t 'Numerical-Solution', 'Error-Distribution.dat' u 1:3 w p t 'Analytical-Solution'
~~~

~~~gnuplot Error distribution 1D - Full Domain
set output '1DError.png'
set xlabel 'X'
set ylabel 'Error(X)'
set grid
plot 'Error-Distribution.dat' u 1:4 w p t 'Error-FullDomain'
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