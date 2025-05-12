/**
#Advection of a discontinuous 1D tracer field.

In this test case we want to observe the qualitative performance of the BCG-advection scheme with different limiter options
viz. No-limiter / Minmod / Superbee and we compare the performance to the WENO-5 advection scheme which has third-order
limiters by construction of the scheme. 

We start with a tracer which has smooth but narrow combination of Gaussians, a square wave, a sharp triangular wave and a 
half ellipse. We advect it using a constant advection current over one time-period in a periodic domain, and observe the 
qualitative shapes of the advected tracer when compared to the initial tracer field to give us a picture of the dissipative
or oscillatory nature of different limited schemes. While the superbee limiter is the best choice among all the bcg-advection
schemes, the weno-5 scheme has lower dissipation compared to the superbee limiter.
*/

#include "grid/multigrid.h"
#define dimension 1

#include "utils.h"
#include "../../Header_Files/advection.h"
#define vel 0.1

scalar f[] , * tracers = {f};

double Function_F ( double x, double alpha, double a ) {
  if ( 1 - sq(alpha)*sq(x-a) >= 0 )
     return ( sqrt (1 - sq(alpha)*sq(x-a)) );
  else
     return (0);
}

double Function_G ( double x, double beta, double z ) {
  return ( exp( -1.*beta*sq(x-z) ) ); 
}

double Point_Value ( double x, double a, double z, double del, double alpha, double beta ) {
 
   if ( x >= -0.8 && x <= -0.6 )
       return ( ( Function_G(x,beta,z-del) + Function_G(x,beta,z+del) + 4.*Function_G(x,beta,z) )/6. );
   else if ( x >= -0.4 && x <= -0.2 )
       return (1.);
   else if ( x >= 0. && x <= 0.2 )
       return ( 1 - 10.*fabs(x-0.1) );
   else if ( x >= 0.4 && x <= 0.6 )
       return ( ( Function_F(x,alpha,a-del) + Function_F(x,alpha,a+del) + 4.*Function_F(x,alpha,a) )/6. );
   else
       return (0.);
  
}

double volume_avg ( double x, double Delta ) {
   
  double a = 0.5, z = -0.7, del = 0.005, alpha = 10.;
  double beta = log(2)/(36.*sq(del));
  double Quad = Delta*sqrt(3./5.);
  return ( (5./18.)*Point_Value (x-Quad,a,z,del,alpha,beta) + (8./18.)*Point_Value (x,a,z,del,alpha,beta) + (5./18.)*Point_Value (x+Quad,a,z,del,alpha,beta) );
   
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
   foreach_dimension ()
      periodic (left);

   DT = .1;
   CFL = 0.4;
   for (N = 512; N <= 512; N *= 2)
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
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);

   FILE * fp = fopen ("Data-512.dat","w");
   foreach ()
      fprintf (fp, "%g %g %g %g \n", x, volume_avg (x,Delta), f[], e[] );
   fclose (fp); 
}


/**

~~~gnuplot
reset
set xlabel 'x'
set ylabel 'Tracer ( x , t_{final} )'
set key font ",8"
set yrange [-0.2:1.25]
plot 'Data-512.dat' u 1:2 w l lt 2 lw 3 lc rgb "black" t 'Analytical', \
     '' u 1:3 w l lw 2 lt 1 lc rgb "blue" t 'bcg-NoLimiting', \
     '../discontinuousadvection1D-Minmod/Data-512.dat' u 1:3 w l lw 2 lt 1 lc rgb "green" t 'bcg-Minmod', \
     '../discontinuousadvection1D-Superbee/Data-512.dat' u 1:3 w l lw 2 lt 1 lc rgb "red" t 'bcg-Superbee'  
~~~

~~~gnuplot
reset
set xlabel 'x'
set ylabel 'Tracer ( x , t_{final} )'
set key font ",8"
plot 'Data-512.dat' u 1:2 w l lt 2 lw 3 lc rgb "black" t 'Analytical', \
     '../discontinuousadvection1D-weno/Data-512.dat' u 1:3 w l lw 2 lt 1 lc rgb "blue" t 'weno', \
     '../discontinuousadvection1D-Superbee/Data-512.dat' u 1:3 w l lw 2 lt 1 lc rgb "red" t 'bcg-Superbee'  
~~~
*/ 
