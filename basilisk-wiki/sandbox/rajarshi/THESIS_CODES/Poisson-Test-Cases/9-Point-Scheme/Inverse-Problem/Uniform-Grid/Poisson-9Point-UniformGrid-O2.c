/**
#Poisson-Helmholtz solver on a Uniform Grid
*/

#include "grid/multigrid.h"
#include "utils.h"
#define dimension 2

#if HIGHER
   #include "../../../../Header_Files/poisson_O4.h"
   #define TOL 1E-08
   #define BGHOSTS 2
#else
   #include "poisson.h"
   #define TOL 1E-04
#endif

scalar A[],B[],A_an[],Error[];

double Poisson_LHS (double x, double y, double Delta){
  return ((cos(2.*pi*(x-Delta/2.))-cos(2.*pi*(x+Delta/.2)))*(cos(2.*pi*(y-Delta/2.))-cos(2.*pi*(y+Delta/2.)))/(sq(2.*pi*Delta)));
}

double Poisson_RHS (double x, double y, double Delta){
  return ( -2.*(cos(2.*pi*(x-Delta/2.))-cos(2.*pi*(x+Delta/.2)))*(cos(2.*pi*(y-Delta/2.))-cos(2.*pi*(y+Delta/2.)))/(sq(Delta)) );
}

void convergence (int depth){

  L0 = 1;
  origin(-0.5,-0.5);
  init_grid(1<<depth);
  int N = pow(2,depth);
  foreach_dimension()
     periodic (left);

  foreach(){
    B[] = Poisson_RHS(x,y,Delta);
    A_an[] = Poisson_LHS(x,y,Delta);
  }
  boundary({B,A_an});

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  poisson (A,B,tolerance=TOL,nrelax=4);
  update_perf();
  timer_print (perf.gt, iter, perf.tnc);

  foreach()
     Error[] = A[] - A_an[];
  stats s2 = statsf(Error);
  foreach()
    Error[] -= s2.sum/s2.volume;

  fprintf(stderr,"%g %g\n",pow(2,depth),normf(Error).max);
}

int main(){

  for(int depth = 5; depth <= 8; depth++)
     convergence(depth);
}

/**

~~~gnuplot Accuracy of the solution as a function of the level of refinement
reset
set term @PNG enhanced
set key font ",12"
set xtics font ", 12"
set ytics font ", 12"
set xlabel "x-units" font ",14"
set ylabel "y-units" font ",14"
set xlabel 'Spatial resolution'
set ylabel 'Max Error Norm'
set cbrange [1:1]
set logscale
set xtics 32,2,256
ftitle(a,b) = sprintf("order %4.2f", -b)
f1(x)=a1+b1*x
fit f1(x) 'log' u (log($1)):(log($2)) via a1,b1
f2(x)=a2+b2*x
fit f2(x) '../Poisson-9Point-UniformGrid-O4/log' u (log($1)):(log($2)) via a2,b2
set xrange [16:512]
plot exp (f1(log(x))) t ftitle(a1,b1) lc rgb "red", \
     'log' u 1:2 t '5P-Poisson-2' pt 5 ps 0.5 lc rgb "red", \
     exp (f2(log(x))) t ftitle(a2,b2) lc rgb "blue" lw 2, \
     '../Poisson-9Point-UniformGrid-O4/log' u 1:2 t '9P-Poisson-4' pt 7 ps 0.5 lc rgb "blue"
~~~

~~~gnuplot Error vs Computing Time
reset
set term @PNG enhanced
set key font ",12"
set xtics font ",12"
set ytics font ",12"
set xlabel 'Computing time (sec)' font ",14" 
set ylabel 'Error-Max' font ",14" offset 10
set cbrange [1:2]
set logscale
set grid ytics
set ytics format "%.0e"
! awk '($1 == "#") {print $5}' < out > time
! awk '($1 != "#") {print $1,$2}' < log > error
! awk '($1 == "#") {print $5}' < ../Poisson-9Point-UniformGrid-O4/out > time-O4
! awk '($1 != "#") {print $1,$2}' < ../Poisson-9Point-UniformGrid-O4/log > error-O4
set xrange [0.001:10]
plot '< paste time error' u 1:3 w lp t '5P-Poisson-2' pt 5 ps 0.5 lc rgb "red",            \
     '< paste time error' u 1:3:(sprintf("%g",$2)) with labels font ",12"  rotate by 45 offset char -1,-1 notitle, \
     '< paste time-O4 error-O4' u 1:3:(sprintf("%g",$2)) with labels font ",12"  rotate by 45 offset char -1,-1 notitle, \
     '< paste time-O4 error-O4' u 1:3 w lp t '9P-Poisson-4' pt 7 ps 0.5 lc rgb "blue"
~~~

~~~gnuplot Computing steps vs Resolution
reset
set term @PNG enhanced
set key font ",12"
set xtics font ",12"
set ytics font ",12"
set xlabel 'Grid Resolution' font ",14" 
set ylabel 'Steps.Points/sec' font ",14" offset 10
set logscale
set cbrange [1:2]
set xrange [16:512]
set xtics 32,2,256
set ytics format "%.0e"
! awk '($1 == "#") {print $9}' < out > time
! awk '($1 != "#") {print $1,$2}' < log > error
! awk '($1 == "#") {print $9}' < ../Poisson-9Point-UniformGrid-O4/out > time-O4
! awk '($1 != "#") {print $1,$2}' < ../Poisson-9Point-UniformGrid-O4/log > error-O4

plot '< paste time error' u 2:1 w lp t '5P-Poisson-2' pt 5 ps 0.5 lc rgb "red",            \
     '< paste time-O4 error-O4' u 2:1 w lp t '9P-Poisson-4' pt 7 ps 0.5 lc rgb "blue"
~~~
*/
