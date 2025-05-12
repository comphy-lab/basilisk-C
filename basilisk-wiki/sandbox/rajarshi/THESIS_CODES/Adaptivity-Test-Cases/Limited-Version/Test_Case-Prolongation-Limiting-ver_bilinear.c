/**
#TEST CASE - PROLONGATION FUNCTION ON DISCONTINUOUS FUNCTIONS

This test case compares the effect of different Prolongation-functions when applied to a discontinuous-function. We use the traditionally applied bi-linear
prolongation, and two different versions of the fifth-order prolongation function. We take a simple 1D-discontinuous function, and the grid is
refined by one level in the sub-domain : ( −0.25 < x < 0.25 ). The refinement region includes the field-discontinuity. The prolongated fine-cell values (marked by blue dots) along with the coarse cell values (marked by red dots) are plotted. We also demonstrate the third order convergence of the error-norm
which is caused due to the field discontinuity. A continuous field will display an error convergence of order-five using the same prolongation function
as can be seen in the test case :- [PROLONGATION - CONTINUOUS FUNCTIONS](http://basilisk.fr/sandbox/rajarshi/THESIS_CODES/Adaptivity-Test-Cases/Non-Limited-Version/Test_Case-Prolongation-NoLimiting-ver_bilinear.c)
*/

#define dimension 1
#include "grid/bitree.h"
#include "utils.h"

#if HIGHER
   #define BGHOSTS 2
   #if LIMITED
       #include "../../Header_Files/Adapt_Limiters.h"
   #else
        #include "../../Header_Files/Adapt_No_Limiters.h"
   #endif
#endif

static double FuncValues (double x) {
  if(x<0) 
    return ( 1 - sin(pi*x) );
  else
    return (-1 - sin(pi*x) );
}

#if HIGHER
#define BGHOSTS 2
static double Volume_Avg (double x, double Delta) {
  double q = (Delta/2.)*sqrt(3./5.);
  return ( (5./18.)*FuncValues(x-q) + (8./18.)*FuncValues(x) + (5./18.)*FuncValues(x+q) );
}
#endif

void convergence (int depth){

  L0 = 1;
  origin (-0.5);
  init_grid(1<<depth);
  scalar s[],error[];

#if HIGHER  
  s.refine = refine_order5;
  s.prolongation = refine_order5;
  foreach()
     s[] = Volume_Avg (x,Delta);
#else
  foreach()
     s[] = FuncValues (x);
#endif

  foreach_dimension(){
     s[left] = neumann(0);
     s[right] = neumann(0);
   }
  boundary({s}); 

  FILE * fp;
  if(depth==5){
     fp = fopen("BeforeRefine.dat","w");
     foreach()
        fprintf(fp,"%g %g\n",x,s[]);
     fclose(fp);
   }

  refine((level < depth+1 ) && (sq(x) + sq(y) <= 0.0625) );

  if(depth==5){
     fp = fopen("AfterRefine.dat","w");
     foreach()
       fprintf(fp,"%g %g \n",x,s[]);
     fclose(fp);
   }


  foreach(){
     #if HIGHER
         error[] = s[] - Volume_Avg (x,Delta);
     #else
         error[] = s[] - FuncValues (x);
     #endif
  }
  fprintf(stderr,"%g %g \n",pow(2,depth+1),normf(error).max);
} 

int main(){
  for (int depth = 4; depth<=9; depth++)
     convergence(depth); 
}

/**
## Results

~~~gnuplot BILINEAR-PROLONGATION
reset
set xlabel 'x' font ",8"
set ylabel 'f(x)' font ",8"
set grid
set xtics 0.5
plot 'AfterRefine.dat' u 1:2 w p notitle pt 5 ps 0.25 lc rgb "blue", \
     'BeforeRefine.dat' u 1:2 w p notitle pt 7 ps 0.35 lc rgb "red"
~~~

~~~gnuplot HIGHER-PROLONGATION (NON-LIMITED-VERSION)
reset
set xlabel 'x' font ",8"
set ylabel 'f(x)' font ",8"
set grid
set xtics 0.5
plot '../Test_Case-Prolongation-Limiting-ver_higherNL/AfterRefine.dat' u 1:2 w p notitle pt 5 ps 0.25 lc rgb "blue", \
     '../Test_Case-Prolongation-Limiting-ver_higherNL/BeforeRefine.dat' u 1:2 w p notitle pt 7 ps 0.35 lc rgb "red"
~~~

~~~gnuplot HIGHER-PROLONGATION (LIMITED-VERSION)
reset
set xlabel 'x' font ",8"
set ylabel 'f(x)' font ",8"
set grid
set xtics 0.5
plot '../Test_Case-Prolongation-Limiting-ver_higherL/AfterRefine.dat' u 1:2 w p notitle pt 5 ps 0.25 lc rgb "blue", \
     '../Test_Case-Prolongation-Limiting-ver_higherL/BeforeRefine.dat' u 1:2 w p notitle pt 7 ps 0.35 lc rgb "red"
~~~

~~~gnuplot ERROR CONVERGENCE OF THE LIMITED PROLONGATION FUNCTION
reset
set key font ",10"
set xtics font ",8"
set ytics font ",8"
set xlabel "Grid Resolution" font ",8"
set ylabel "Error-Max" font ",8"
set logscale
set cbrange [1:2]
set grid
set xtics 32,2,1024
ftitle(a,b) = sprintf("Order %4.2f",-b)
f2(x) = a2 + b2*x
fit f2(x) '../Test_Case-Prolongation-Limiting-ver_higherL/log' u (log($1)):(log($2)) via a2, b2
set xrange [16:2048]
plot '../Test_Case-Prolongation-Limiting-ver_higherL/log' u 1:2 w p t 'Higher-Limiting' pt 7 ps 0.5 lc rgb "red", exp(f2(log(x))) t ftitle(a2,b2) lc rgb "red"
~~~
*/
