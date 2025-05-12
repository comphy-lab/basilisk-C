/**
#Wave propagation - Implicit saint venant solver.
*/

#include "grid/multigrid1D.h"
#include "saint-venant-implicit.h"

double amp;

int main(){

   init_grid(512);
   periodic(right);
   G = 10.;
   DT = 1e-2;
   theta = 1.3;

   system("rm -f log.dat");
   amp = 0.1;
   run();
   

   amp = 0.4;
   run();

   amp = 1.;
   run();

  
   system("cat out-* > log");
}

event init (i=0) {
   foreach()
     h[] = 1. + amp*sin(2.*pi*x);
}

event logfile (i++) {
   char name[80];
   sprintf(name, "log-%g",amp);
   static FILE * fp = fopen(name,"w");
   fprintf(fp,"%g %g %g \n",t,dt,statsf(h).sum);
}

event output (t<=0.2; t+=0.1) {
   char name[80];
   sprintf(name, "out-%g",amp);
   static FILE * fp = fopen(name,"w");
   foreach()
      fprintf(fp, "%g %g \n",x,h[]);
   fprintf(fp,"\n");
   fflush(fp);
}

/**
For small wave amplitudes, non-linear effects are small and we get the
following wave evolution.

~~~gnuplot Solution for a wave amplitude of 0.1
set key bottom left
set xlabel 'x'
set ylabel 'z'
plot 'out-0.1' w l t 'implicit-O5', '../implicit_weno/out-0.1' w l t 'implicit-weno' 
~~~

~~~gnuplot Solution for a wave amplitude of 0.4
plot 'out-0.4' w l t 'implicit-O5', '../implicit_weno/out-0.4' w l t 'implicit-weno'
~~~

~~~gnuplot Solutions for a wave amplitude of 1
plot 'out-1' w l t 'implicit-O5', '../implicit_weno/out-1' w l t 'implicit-weno' 
~~~
*/

