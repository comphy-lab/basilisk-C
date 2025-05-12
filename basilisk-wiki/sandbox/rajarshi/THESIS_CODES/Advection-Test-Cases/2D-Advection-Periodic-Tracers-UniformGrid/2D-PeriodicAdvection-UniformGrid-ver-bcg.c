/**
#2D Advection of a Periodic tracer field


In this test case we wish to observe the convergence and performance of the WENO solver in 2D, where an additional gaussian-quadrature interpolation is
required to implement the higher-order scheme. For details see (Rajarshi-PhD-Thesis) section 2.8. We have a setup similar to the 
[1D-case](http://basilisk.fr/sandbox/rajarshi/THESIS_CODES/Advection-Test-Cases/1D-Advection-Periodic-Function/Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.c) here in 2D, and the convergence and performance is compared with the bcg solver.
*/

#define dimension 2
#include "grid/multigrid.h"
#include "utils.h"
#include "../../Header_Files/advection.h"
#define vel 0.1

scalar f[], * tracers = {f};

#if !WENO
static double Point_Value (double x, double y) {
  return (sin(2.*pi*x)*sin(2.*pi*y));
}
#endif

#if WENO
static double volume_avg (double x, double y, double Delta) {
  return ((cos(2.*pi*(x-Delta/2.))-cos(2.*pi*(x+Delta/2.)))*(cos(2.*pi*(y-Delta/2.))-cos(2.*pi*(y+Delta/2.)))/sq(2.*pi*Delta));
}
#endif

void face_velocity(face vector u, double time){
  foreach_face(x)
     u.x[] = vel;
  foreach_face(y)
     u.y[] = vel;
  boundary((scalar *){u});
}

#if !WENO
event velocity (i++) {
  face_velocity (u, t);
}
#endif

int main () {
  L0=1.;
  origin(-0.5,-0.5);
  foreach_dimension()
      periodic(left);

  for (N=32; N<=256; N*=2){ 
     DT=1/(N*vel);
     CFL=0.8; 
     run();
     CFL=0.1;
     run();
  }

}

event init ( i = 0 ) {     

  #if WENO
     foreach()
        f[] = volume_avg (x, y, Delta);
  #else
     foreach()
        f[] = Point_Value(x,y);
  #endif
  boundary ({f});

}

event field (t = 10) {
  scalar e[];
  #if WENO
    foreach()
      e[] = f[] - volume_avg (x, y, Delta);
  #else
    foreach()
      e[] = f[] - Point_Value(x, y);
  #endif
  boundary ({e});
  norm n = normf (e);
  fprintf (stderr, "%g %d %g %g %g\n", CFL, N, n.avg, n.rms, n.max);
}

/**
~~~gnuplot Error-norm vs spatial resolution
reset
set key bottom left
set key font ",8"
set xtics font ",8"
set ytics font ",8"
set xlabel "Grid Resolution" font ",10"
set ylabel "Error-Max" font ",10"
set cbrange [1:2]
! awk '($1 == "0.8") {print $2,$5}' < log > CFL-0.8
! awk '($1 == "0.1") {print $2,$5}' < log > CFL-0.1
! awk '($1 == "0.8") {print $2,$5}' < ../2D-PeriodicAdvection-UniformGrid-ver-weno/log > CFL-Weno-0.8
! awk '($1 == "0.1") {print $2,$5}' < ../2D-PeriodicAdvection-UniformGrid-ver-weno/log > CFL-Weno-0.1
ftitle(a,b) = sprintf("Order %4.2f",-b)
f1(x)=a1+b1*x
fit f1(x) 'CFL-0.1' u (log($1)):(log($2)) via a1,b1
f2(x)=a2+b2*x
fit f2(x) 'CFL-0.8' u (log($1)):(log($2)) via a2,b2
f3(x)=a3+b3*x
fit f3(x) 'CFL-Weno-0.1' u (log($1)):(log($2)) via a3,b3
f4(x)=a4+b4*x
fit f4(x) 'CFL-Weno-0.8' u (log($1)):(log($2)) via a4,b4
set logscale
set xrange [16:512]
set xtics 32,2,256
plot 'CFL-0.1' u 1:2 t 'BCG(CFL=0.1)' ps 0.5 pt 5 lc rgb "red", exp(f1(log(x))) t ftitle(a1,b1) lc rgb "red", \
     'CFL-0.8' u 1:2 t 'BCG(CFL=0.8)' ps 0.5 pt 5 lc rgb "blue", exp(f2(log(x))) t ftitle(a2,b2) lc rgb "blue", \
     'CFL-Weno-0.1' u 1:2 t 'WENO(CFL=0.1)' ps 0.5 pt 7 lc rgb "green", exp(f3(log(x))) t ftitle(a3,b3) lc rgb "green", \
     'CFL-Weno-0.8' u 1:2 t 'WENO(CFL=0.8)' ps 0.5 pt 7 lc rgb "black", exp(f4(log(x))) t ftitle(a4,b4) lc rgb "black"
~~~

~~~gnuplot Error-norm vs Computing time
reset
set key font ",8"
set xtics font ",8"
set ytics font ",8"
set xlabel 'Computing time (sec)' font ",10" 
set ylabel 'Error-Max' font ",10"
set cbrange [1:2]
set logscale
set grid ytics
set ytics format "%.0e"
! awk '($1 == "#") {print $5}' < out > time
! awk '($1 != "#") {print $2,$5}' < log > error
! awk '($1 == "#") {print $5}' < ../2D-PeriodicAdvection-UniformGrid-ver-weno/out > time-weno
! awk '($1 != "#") {print $2,$5}' < ../2D-PeriodicAdvection-UniformGrid-ver-weno/log > error-weno

plot '< paste time error' every 2::0::7 u 1:3 w lp t 'BCG(CFL=0.8)' pt 5 ps 0.5 lc rgb "red",            \
     '< paste time error' every 2::1::7 u 1:3 w lp t 'BCG(CFL=0.1)' pt 5 ps 0.5 lc rgb "blue", \
     '< paste time-weno error-weno' every 2::1::7 u 1:3:(sprintf("%g",$2)) with labels font ",12" offset char 2,0 notitle, \
     '< paste time-weno error-weno' every 2::0::7 u 1:3 w lp t 'WENO(CFL=0.8)' pt 7 ps 0.5 lc rgb "green",            \
     '< paste time-weno error-weno' every 2::1::7 u 1:3 w lp t 'WENO(CFL=0.1)' pt 7 ps 0.5 lc rgb "black"
~~~

~~~gnuplot Computing steps vs Resolution
reset
set key font ",8"
set xtics font ",8"
set ytics font ",8"
set xlabel 'Grid Resolution' font ",10" 
set ylabel 'Steps.Points/sec' font ",10"
set logscale
set cbrange [1:2]
set yrange [1e+05:1e+08]
set xrange [16:512]
set ytics 1e+05,10,1e+08
set xtics 32,2,256
set ytics format "%.0e"
! awk '($1 == "#") {print $9}' < out > time
! awk '($1 != "#") {print $2,$5}' < log > error
! awk '($1 == "#") {print $9}' < ../2D-PeriodicAdvection-UniformGrid-ver-weno/out > time-weno
! awk '($1 != "#") {print $2,$5}' < ../2D-PeriodicAdvection-UniformGrid-ver-weno/log > error-weno

plot '< paste time error' every 2::0::7 u 2:1 w lp t 'BCG(CFL=0.8)' pt 5 ps 0.5 lc rgb "red",            \
     '< paste time error' every 2::1::7 u 2:1 w lp t 'BCG(CFL=0.1)' pt 5 ps 0.5 lc rgb "blue",            \
     '< paste time-weno error-weno' every 2::0::7 u 2:1 w lp t 'WENO(CFL=0.8)' pt 7 ps 0.5 lc rgb "green",  \
     '< paste time-weno error-weno' every 2::1::7 u 2:1 w lp t 'WENO(CFL=0.1)' pt 7 ps 0.5 lc rgb "black"
~~~
*/
