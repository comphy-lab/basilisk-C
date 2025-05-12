/**
#TEST CASE - ADVECTION OF A PERIODIC TRACER (BCG vs WENO vs O5)

In this preliminary test case, we demonstrate that the WENO and the O5 advection schemes in 1D, as well as the time marching 
RK4 scheme, which have been implemented in Basilisk are functional. We study the error convergence and the time performance of 
the 1D-advection solvers and compare it with the existing BCG tracer advection solver.

We use a periodic sinusoidal tracer function in a periodic domain and advect it using a constant velocity. We run the simulation over
one time-period and compare the advected solution to the initial condition to observe the error and performance of the solvers. 
*/


#define dimension 1
#include "grid/multigrid.h"
#include "utils.h"
#include "../../Header_Files/advection.h"
#define vel 0.1

scalar f[] , * tracers = {f};

static double Point_Value (double x) {
  return (sin(2.*pi*x));
}

#if WENO
static double volume_avg (double x, double Delta) {
  double q = (Delta/2.)*sqrt(3./5.);
  return ( (5./18.) * Point_Value (x-q) + (8./18.) * Point_Value (x) + (5./18.) * Point_Value (x+q));
}
#endif

void face_velocity(face vector u, double time){
  foreach_face()
     u.x[] = vel;
  boundary((scalar *){u});
} 

#if !WENO
event velocity (i++) {
  face_velocity (u, t);
}
#endif

int main() {

  L0=1.;
  origin(-0.5);
  periodic(left);
  DT=0.25;
  for (N=64; N<=512; N*=2){
     CFL=0.9; 
     run();
     CFL=0.5;
     run();
  }
}

event init ( i = 0 ) {     
  #if WENO
    foreach()
      f[] = volume_avg (x,Delta);
  #else
    foreach()
      f[] = Point_Value(x);
  #endif
  boundary({f});
}

event field (t = 10) {
  scalar e[];
  #if WENO
    foreach()
      e[] = f[] - volume_avg (x,Delta);
  #else
    foreach()
      e[] = f[] - Point_Value(x);
  #endif
  boundary ({e});
  norm n = normf (e);
  fprintf (stderr, "%g %d %g %g %g\n", CFL, N, n.avg, n.rms, n.max);
}

/**

~~~gnuplot Error-convergence with spatial resolution (BCG vs O5 vs WENO)
reset
set key bottom left
set key font ",8"
set xtics font ",8"
set ytics font ",8"
set xlabel "Grid Resolution" font ",8"
set ylabel "Error-Max" font ",8"
set cbrange [1:2]
! awk '($1 == "0.9") {print $2,$5}' < log > CFL-0.9
! awk '($1 == "0.5") {print $2,$5}' < log > CFL-0.5
! awk '($1 == "0.9") {print $2,$5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno/log > CFL-Weno-0.9
! awk '($1 == "0.5") {print $2,$5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno/log > CFL-Weno-0.5
! awk '($1 == "0.9") {print $2,$5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5/log > CFL-O5-0.9
! awk '($1 == "0.5") {print $2,$5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5/log > CFL-O5-0.5
ftitle(a,b) = sprintf("Order %4.2f",-b)
f1(x)=a1+b1*x
fit f1(x) 'CFL-0.5' u (log($1)):(log($2)) via a1,b1
f2(x)=a2+b2*x
fit f2(x) 'CFL-0.9' u (log($1)):(log($2)) via a2,b2
f3(x)=a3+b3*x
fit f3(x) 'CFL-Weno-0.5' u (log($1)):(log($2)) via a3,b3
f4(x)=a4+b4*x
fit f4(x) 'CFL-Weno-0.9' u (log($1)):(log($2)) via a4,b4
f5(x)=a5+b5*x
fit f5(x) 'CFL-O5-0.5' u (log($1)):(log($2)) via a5,b5
f6(x)=a6+b6*x
fit f6(x) 'CFL-O5-0.9' u (log($1)):(log($2)) via a6,b6
set logscale
set xrange [16:1024]
set xtics 32,2,512
plot 'CFL-0.5' u 1:2 t 'BCG(CFL=0.5)' ps 0.5 pt 5 lc rgb "red", exp(f1(log(x))) t ftitle(a1,b1) lc rgb "red", \
     'CFL-0.9' u 1:2 t 'BCG(CFL=0.9)' ps 0.5 pt 5 lc rgb "blue", exp(f2(log(x))) t ftitle(a2,b2) lc rgb "blue", \
     'CFL-Weno-0.5' u 1:2 t 'WENO(CFL=0.5)' ps 0.5 pt 7 lc rgb "green", exp(f3(log(x))) t ftitle(a3,b3) lc rgb "green", \
     'CFL-Weno-0.9' u 1:2 t 'WENO(CFL=0.9)' ps 0.5 pt 7 lc rgb "black", exp(f4(log(x))) t ftitle(a4,b4) lc rgb "black", \
     'CFL-O5-0.5' u 1:2 t 'O5(CFL=0.5)' ps 0.5 pt 9 lc rgb "orange", exp(f5(log(x))) t ftitle(a5,b5) lc rgb "orange", \
     'CFL-O5-0.9' u 1:2 t 'O5(CFL=0.9)' ps 0.5 pt 9 lc rgb "brown", exp(f6(log(x))) t ftitle(a6,b6) lc rgb "brown"
~~~

~~~gnuplot Computing time vs Error-Norm (BCG vs O5 vs WENO)
reset
set key font ",8"
set xtics font ",8"
set ytics font ",8"
set xlabel 'Computing time (sec)' font ",8" 
set ylabel 'Error-Max' font ",8"
set cbrange [1:2]
set logscale
set grid ytics
set ytics format "%.0e"
set xrange [0.001:1]
! awk '($1 == "#") {print $5}' < out > time
! awk '($1 != "#") {print $2,$5}' < log > error
! awk '($1 == "#") {print $5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno/out > time-weno
! awk '($1 != "#") {print $2,$5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno/log > error-weno
! awk '($1 == "#") {print $5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5/out > time-O5
! awk '($1 != "#") {print $2,$5}' < ../Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5/log > error-O5
plot '< paste time error' every 2::1::7 u 1:3 w lp t 'BCG(CFL=0.5)' pt 5 ps 0.5 lc rgb "blue", \
     '< paste time-weno error-weno' every 2::1::7 u 1:3:(sprintf("%g",$2)) with labels font ",12" offset char 2,0 notitle, \
     '< paste time-weno error-weno' every 2::1::7 u 1:3 w lp t 'WENO(CFL=0.5)' pt 7 ps 0.5 lc rgb "black", \
     '< paste time-O5 error-O5' every 2::1::7 u 1:3 w lp t 'O5(CFL=0.5)' pt 9 ps 0.5 lc rgb "orange"
~~~
*/
