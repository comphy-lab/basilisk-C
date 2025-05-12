/**
#TEST CASE - PROLONGATION OPERATOR ON CONTINUOUS FUNCTIONS

This test case compares the effect of different Prolongation-operators when applied on a 2D-continuous-function. We use the traditionally applied bi-linear
prolongation, and two different versions of the fifth-order prolongation function. We take a simple 2D-continuous function which is refined in a given
sub-domain. The prolongated values are compared with analytical values for all three prolongation operators. We display the error convergence and 
computational time required for each operator. We observe the increased computing cost of the Limiting-prolongation operator viz-a-viz the other two
operators on the error vs computing time plots.

A dis-continuous field will display an error convergence of order-three using the limiting version of the prolongation function, and the no-limiting version will give numerical oscillations near the region of dis-continuity as can be seen in the test case :- [PROLONGATION - DISCONTINUOUS FUNCTIONS](http://basilisk.fr/sandbox/rajarshi/THESIS_CODES/Adaptivity-Test-Cases/Limited-Version/Test_Case-Prolongation-Limiting-ver_bilinear.c)
*/

#include "grid/tree.h"
#define dimension 2
#include "utils.h"

#if HIGHER
   #define BGHOSTS 2
   #if LIMITED
       #include "../../Header_Files/Adapt_Limiters.h"
   #else
        #include "../../Header_Files/Adapt_No_Limiters.h"
   #endif
#endif

static double FuncValues (double x, double y) {
  if(sq(x)<=0.0625 && sq(y)<=0.0625) 
    return (pow(1-sq(x/0.25),7)*pow(1-sq(y/0.25),7));
  else
    return 0;
}

#if HIGHER
static double Volume_Avg (double x, double y, double Delta) {
  double sum=0.;
  double q = (Delta/2.)*sqrt(3./5.);
  sum += (5./18.)*( (5./18.)*FuncValues(x-q,y-q) + (8./18.)*FuncValues(x,y-q) + (5./18.)*FuncValues(x+q,y-q) );
  sum += (8./18.)*( (5./18.)*FuncValues(x-q,y) + (8./18.)*FuncValues(x,y) + (5./18.)*FuncValues(x+q,y) );
  sum += (5./18.)*( (5./18.)*FuncValues(x-q,y+q) + (8./18.)*FuncValues(x,y+q) + (5./18.)*FuncValues(x+q,y+q) );
  return (sum);
}
#endif

void convergence (int depth){

  L0 = 1;
  origin (-0.5,-0.5);
  init_grid(1<<depth);
  scalar s[],error[];

#if HIGHER  
  s.refine = refine_order5;
  s.prolongation = refine_order5;
  foreach()
     s[] = Volume_Avg (x,y,Delta);
#else
  foreach()
     s[] = FuncValues (x,y);
#endif

  foreach_dimension(){
     s[left] = neumann(0);
     s[right] = neumann(0);
   }
  boundary({s}); 

#if HIGHER
  FILE * fp;
  if(depth==7){
     fp = fopen("BeforeRefine.dat","w");
     foreach()
        fprintf(fp,"%g %g %g\n",x,y,s[]);
     fclose(fp);
   }
#endif

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  refine((level < depth+1 ) && (sq(x) + sq(y) <= 0.04) );
  update_perf();
  timer_print (perf.gt, iter, perf.tnc);

#if HIGHER
  if(depth==7){
     fp = fopen("AfterRefine.dat","w");
     foreach()
       fprintf(fp,"%g %g %g \n",x,y,s[]);
     fclose(fp);
   }
#endif

  foreach(){
     #if HIGHER
         error[] = s[] - Volume_Avg (x,y,Delta);
     #else
         error[] = s[] - FuncValues (x,y);
     #endif
  }
  fprintf(stderr,"%g %g \n",pow(2,depth),normf(error).max);
} 

int main(){

  for (int depth = 7; depth<=10; depth++)
     convergence(depth); 
}

/**
## Results

~~~gnuplot Test function Before Refinement
reset
set xtics font ",8" 0.25
set ytics font ",8" 0.25
set ztics font ",8" 0.25
set xlabel "x" font ",8"
set ylabel "y" font ",8"
set zlabel "f(x,y)" font ",8"
splot '../Test_Case-Prolongation-NoLimiting-ver_higherL/BeforeRefine.dat' u 1:2:3 w p notitle ps 0.2 lc rgb "brown"
~~~

~~~gnuplot Test function after refinement
reset
set xtics font ",8" 0.25
set ytics font ",8" 0.25
set ztics font ",8" 0.25
set xlabel "x" font ",8"
set ylabel "y" font ",8"
set zlabel "f(x,y)" font ",8"
splot '../Test_Case-Prolongation-NoLimiting-ver_higherL/AfterRefine.dat' u 1:2:3 w p notitle ps 0.2 lc rgb "brown"
~~~

~~~gnuplot Error Convergence of the Prolongation operator (bilinear vs higher-order-limiting vs higher-order-nolimiting)
reset
set key bottom left
set key font ",8"
set key spacing 1.25
set xtics font ",8"
set ytics font ",8"
set xlabel "Grid Resolution" font ",8"
set ylabel "Error-Max" font ",8"
set logscale
set cbrange [1:2]
set grid
set xtics 128,2,1024
ftitle(a,b) = sprintf("Order %4.2f",-b)
f1(x) = a1 + b1*x
fit f1(x) 'log' u (log($1)):(log($2)) via a1, b1
f2(x) = a2 + b2*x
fit f2(x) '../Test_Case-Prolongation-NoLimiting-ver_higherNL/log' u (log($1)):(log($2)) via a2, b2
f3(x) = a3 + b3*x
fit f3(x) '../Test_Case-Prolongation-NoLimiting-ver_higherL/log' u (log($1)):(log($2)) via a3, b3
set xrange [64:2048]
plot 'log' u 1:2 w p t 'Bilinear' pt 5 ps 0.75 lc rgb "red", exp(f1(log(x))) t ftitle(a1,b1) lc rgb "red", \
     '../Test_Case-Prolongation-NoLimiting-ver_higherNL/log' u 1:2 w p t 'Higher-NoLimiting' pt 7 ps 0.75 lc rgb "blue", exp(f2(log(x))) t ftitle(a2,b2) lc rgb "blue", \
     '../Test_Case-Prolongation-NoLimiting-ver_higherL/log' u 1:2 w p t 'Higher-Limiting' pt 9 ps 0.75 lc rgb "green", exp(f3(log(x))) t ftitle(a3,b3) lc rgb "green"  
~~~

~~~gnuplot Error vs Computing Time
reset
set key bottom left
set key font ",8"
set key spacing 1.25
set xtics font ",8"
set ytics font ",8"
set xlabel 'Computing time (sec)' font ",8" 
set ylabel 'Error-Max' font ",8"
set cbrange [1:2]
set logscale
set grid ytics
set yrange [1e-12:1]
set xrange [1e-2:10]
set ytics format "%.0e"
! awk '($1 == "#") {print $5}' < out > time
! awk '($1 != "#") {print $1,$2}' < log > error
! awk '($1 == "#") {print $5}' < ../Test_Case-Prolongation-NoLimiting-ver_higherNL/out > time-O5NL
! awk '($1 != "#") {print $1,$2}' < ../Test_Case-Prolongation-NoLimiting-ver_higherNL/log > error-O5NL
! awk '($1 == "#") {print $5}' < ../Test_Case-Prolongation-NoLimiting-ver_higherL/out > time-O5L
! awk '($1 != "#") {print $1,$2}' < ../Test_Case-Prolongation-NoLimiting-ver_higherL/log > error-O5L

plot '< paste time error' u 1:3 w lp t 'Bilinear' pt 5 ps 0.5 lc rgb "red",            \
     '< paste time error' u 1:3:(sprintf("%g",$2)) with labels font ",10"  rotate by 45 offset char 1,1 notitle, \
     '< paste time-O5NL error-O5NL' u 1:3 w lp t 'Higher-NoLimiting' pt 7 ps 0.5 lc rgb "blue", \
     '< paste time-O5NL error-O5NL' u 1:3:(sprintf("%g",$2)) with labels font ",10"  rotate by 45 offset char -1,-1 notitle, \
     '< paste time-O5L error-O5L' u 1:3 w lp t 'Higher-Limiting' pt 9 ps 0.5 lc rgb "green", \
     '< paste time-O5L error-O5L' u 1:3:(sprintf("%g",$2)) with labels font ",10"  rotate by 45 offset char 1,1 notitle
~~~
*/
