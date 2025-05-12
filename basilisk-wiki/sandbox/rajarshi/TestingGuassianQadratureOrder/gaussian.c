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


void convergence ( int depth ){

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

  fprintf(stderr, "%g %g \n", pow(2,depth), statsf(error).max );
 
}

int main () {

  for (int depth = 3; depth <= 7; depth++)
    convergence (depth); 

}

/**

   ~~~gnuplot Error convergence of the quadrature formulations
   ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
   f2P(x)=a+b*x
   fit f2P(x) 'log' u (log($1)):(log($2)) via a,b
   f3P(x)=a2+b2*x
   fit f3P(x) '../gaussian-three/log' u (log($1)):(log($2)) via a2,b2
   set xlabel 'Resolution'
   set ylabel 'Error'
   set logscale
   set xrange [4:256]
   set xtics 4,2,256
   set grid ytics
   set cbrange [1:2]
   set ytics format "%.0e"
   set key outside
   plot 'log' u 1:2 t 'Gaussian-two-point', \
        exp(f2P(log(x))) t ftitle(a,b), \
        '../gaussian-three/log' u 1:2 t 'Gaussian-three-point', \
        exp(f3P(log(x))) t ftitle(a2,b2)
   ~~~

*/
