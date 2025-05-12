/**
#Gaussian Quadrature Integrals

We test the order of the three point and two point scheme respectively
*/

#include "grid/multigrid.h"
#define dimension 2

#include "utils.h"

static double Volume_Avg_Analytical ( double x, double y, double Delta ) {

  return ( ( cos( 2.*pi*( x - Delta/2. ) ) - cos( 2.*pi*(x+Delta/2.) ) )*
	   ( cos( 2.*pi*(y-Delta/2.) ) - cos( 2.*pi*(y+Delta/2.) ) )
	   /sq(2.*pi*Delta) );

}


static double func ( double x, double y ) {

  return(sin(2.*pi*x)*sin(2.*pi*y));

}


#if THREEPOINT

static double VolumeAvg_GuassQuadrature ( double x, double y, double Delta ) {

  double sum = 0.;
  double q = Delta/2.*sqrt(3./5.);
  sum += 5.*(5.*func (x - q, y - q) +
	     8.*func (x, y - q) +
	     5.*func (x + q, y - q));
  sum += 8.*(5.*func (x - q, y) +
	     8.*func (x, y) +
	     5.*func (x + q, y));
  sum += 5.*(5.*func (x - q, y + q) +
	     8.*func (x, y + q) +
	     5.*func (x + q, y + q));
  return sum/sq(18.);

}

#else

static double VolumeAvg_GuassQuadrature ( double x, double y, double Delta ) {

  double sum = 0.;
  double q = Delta/2.*sqrt(1./3.);
  sum += func (x - q, y - q) + func (x + q, y - q);
  sum += func (x - q, y + q) + func (x + q, y + q);
  return sum/4;

}

#endif

clock_t start,end;

void convergence ( int depth ){

  start = clock();
  L0 = 1.;
  origin (-0.5,-0.5);
  init_grid (1<<depth);
  foreach_dimension()
    periodic (left);

  scalar s[], s_an[], error[];
  foreach(){
    s_an[]  = Volume_Avg_Analytical (x, y, Delta );
    s[]     = VolumeAvg_GuassQuadrature ( x, y, Delta );
    error[] = s[] - s_an[];
  } 
  boundary({s,s_an,error});
  end = clock();
  fprintf(stderr, "%g %g %g\n", pow(2,depth), statsf(error).max, (end-start)/(double)CLOCKS_PER_SEC );
 
}

int main () {

  for (int depth = 3; depth <= 7; depth++)
    convergence (depth); 

}

/**

   ~~~gnuplot Error convergence of the quadrature formulations
   reset
   set key font ",12"
   set key spacing 1.5
   set xtics font ",12"
   set ytics font ",12"
   set xlabel "Grid Resolution" font ",12"
   set ylabel "Error-Max" font ",12" offset 10
   ftitle(a,b) = sprintf("ORDER-%4.2f",-b)
   f2P(x)=a+b*x
   fit f2P(x) 'log' u (log($1)):(log($2)) via a,b
   f3P(x)=a2+b2*x
   fit f3P(x) '../Gaussian-three/log' u (log($1)):(log($2)) via a2,b2
   set logscale
   set xrange [4:256]
   set xtics 4,2,256
   set grid ytics
   set cbrange [1:2]
   set ytics format "%.0e"
   plot 'log' u 1:2 t 'Gaussian-two-point' pt 7 ps 0.7, \
        exp(f2P(log(x))) t ftitle(a,b) lc 1 lw 2, \
        '../Gaussian-three/log' u 1:2 t 'Gaussian-three-point' pt 5 ps 0.7, \
        exp(f3P(log(x))) t ftitle(a2,b2) lc 3 lw 2
   ~~~

   ~~~gnuplot Time performance of the quadrature formulation
   reset
   set key font ",12"
   set key spacing 1.5
   set xtics font ",12"
   set ytics font ",12"
   set xlabel "Computing Time (sec)" font ",12"
   set ylabel "Error-Max" font ",12" offset 10
   set logscale
   set grid ytics
   set ytics format "%.0e"
   set cbrange [1:2]
   set key top right
plot 'log' u 3:2:(sprintf("%d",$1)) with labels font ",12" offset char 0.8,0.8 notitle, \
     'log' u 3:2 w lp t 'Gaussian-two-point' pt 7 ps 0.7 lw 2 lc 1,            \
     '../Gaussian-three/log' u 3:2:(sprintf("%d",$1)) with labels font ",12" offset char 0.8,0.8 notitle, \
     '../Gaussian-three/log' u 3:2 w lp t 'Gaussian-three-point' pt 5 ps 0.7 lw 2 lc 3
~~~
*/
