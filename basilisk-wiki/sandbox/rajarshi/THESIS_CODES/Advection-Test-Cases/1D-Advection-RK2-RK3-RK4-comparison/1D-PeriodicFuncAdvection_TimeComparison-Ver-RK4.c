/**
#Higher-order advection - Comparing RK2 vs RK4 vs SSP-RK3 time-marching.

This test case is an extension of the [Periodic-Tracer-Advection](http://basilisk.fr/sandbox/rajarshi/THESIS_CODES/Advection-Test-Cases/1D-Advection-Periodic-Function/Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.c) test-case, where we have kept the advection scheme as the O5 scheme, and use different time marching schemes to observe their effects on the error-norms. We have implemented the classical RK2,RK4 and a strong-stability-preserving RK3 scheme in Basilisk, and will observe from
the error-norm plots as to where the error is dictated by the advection scheme discretizations and in which cases is it dictated by the temporal scheme discretizations.

The simulations have been run from a CFL number of 0.82 (equivalent time-step = 0.016) to 0.05 (equivalent time-step = 0.001). The equivalent time-steps have been plotted
on the graph as data-point-labels. From the graph, it becomes clear that the RK4 scheme shows a saturated error-value accruing from the temporal scheme even for high-CFL numbers, which require a shorter computation time, and hence we naturally choose RK4 schemes to be used along with WENO schemes (continuous tracer-fields), but we will
use SSP-RK3 scheme for tracer-fields with discontinuities as done in the test case on [Discontinuous-Tracer-Advection](http://basilisk.fr/sandbox/rajarshi/THESIS_CODES/Advection-Test-Cases/1D-Advection-Discontinuous-Tracers/discontinuousadvection1D-NoLimiter.c).  
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

static double volume_avg (double x, double Delta) {
  double q = (Delta/2.)*sqrt(3./5.);
  return ( (5./18.) * Point_Value (x-q) + (8./18.) * Point_Value (x) + (5./18.) * Point_Value (x+q) );
}

void face_velocity ( face vector u, double t ) {
  foreach_face ()
     u.x[] = vel;
  boundary ( (scalar *) {u} );
}

#if !WENO
event velocity (i++) {
  face_velocity (u, t);
}
#endif

int main () {
  
   L0 = 2.;
   origin (-1.);
   N = 512;
   foreach_dimension ()
      periodic (left);
   CFL = 0.8;

   for (DT = 0.016; DT >= 0.001; DT /= 2)
      run();
}

event init ( i = 0 ) {     
  foreach()
     f[] = volume_avg (x, Delta);
  boundary ({f});
}

event logfile (t={0,1}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %ld %.12f %g %g\n", t, grid->tn, s.sum, s.min, s.max);
}

event field (t = 20) {
  scalar e[];
  foreach()
    e[] = f[] - volume_avg (x, Delta);
  boundary ({e});
  norm n = normf (e);
  fprintf (stderr, "%g %g %g %g\n", DT, n.avg, n.rms, n.max);
}


/**

~~~gnuplot Computation time performance graph using RK2, SSP-RK3 and RK4 schemes. Grid resolution 512.
reset
set key font ",10"
set xtics font ",8"
set ytics font ",8"
set xlabel "Computing time (sec)" font ",10"
set ylabel "Error-Max" font ",10"
set logscale
set grid ytics
set ytics format "%.0e"
! awk '($1 == "#") {print $5}' < out > time
! awk '($1 != "#") {print $2,$4,$1}' < log > error
! awk '($1 == "#") {print $5}' < ../1D-PeriodicFuncAdvection_TimeComparison-Ver-RK2/out > time-RK2
! awk '($1 != "#") {print $2,$4,$1}' < ../1D-PeriodicFuncAdvection_TimeComparison-Ver-RK2/log > error-RK2
! awk '($1 == "#") {print $5}' < ../1D-PeriodicFuncAdvection_TimeComparison-Ver-SSPRK3/out > time-SSPRK3
! awk '($1 != "#") {print $2,$4,$1}' < ../1D-PeriodicFuncAdvection_TimeComparison-Ver-SSPRK3/log > error-SSPRK3
set cbrange [1:2]
set yrange [1e-10:1e-03]
set xrange [0.07:7]
set xtics 0.07,10,7
set key top right
plot '< paste time-RK2 error-RK2' u 1:3:(sprintf("%g",$4)) with labels font ",8" offset char 0.5,0.5 notitle, \
     '< paste time-RK2 error-RK2' u 1:3 w lp t 'WENO5-RK2' pt 7 ps 0.5 lw 2 lc rgb "blue", \
     '< paste time-RK2 error-SSPRK3' u 1:3:(sprintf("%g",$4)) with labels font ",8" offset char 3.2,0.7 notitle, \
     '< paste time-SSPRK3 error-SSPRK3' u 1:3 w lp t 'WENO5-SSP-RK3' pt 9 ps 0.5 lw 2 lc rgb "green", \
     '< paste time error' u 1:3:(sprintf("%g",$4)) with labels font ",8" offset char 0,-1 notitle, \
     '< paste time error' u 1:3 w lp t 'WENO5-RK4' pt 5 ps 0.5 lw 2 lc rgb "red"
~~~
*/ 
