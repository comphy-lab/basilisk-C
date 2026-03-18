/** 
## Inclined free surface Poiseuille flow fluid with the layered solver
An example of 2D free surface flow  over a inclined plate  is presented here.
The configuration is periodic. 
*/
#define ML 1

# include "grid/multigrid1D.h"
# if !ML
#  include "saint-venant.h"
# else // ML
# include "layered/hydro.h"
//#  include "hydroMT.h" // need to be tested with  diffusionMT.h
#  include "layered/remap.h"
# endif // ML


const double HR = 1., NU = 0.1, T0 = 100;
double slope;
scalar uold[];
//FILE *file;

int main()
{
  periodic(right);
 
  slope = 0.359;  
  L0 = 1.;
  G  = 1;
  N  = 16; // no need to refine, the velocity varies only along the vertical direction
  nu = NU;
  nl = 64;

  // CFL_H = HUGE;
  // DT = 0.5;
   for (N = 8; N<= 32; N*= 2){
    for (nl = 4; nl <= 256; nl *= 2)
      run();
    fprintf (stderr,"\n\n");
  }
}

/**
We initialise the topography and the initial thickness of each layer *h*. */

event init (i = 0)
{
  foreach() {
    zb[] = 0.;
    //zb[] = -(x-X0)*slope; //;
#if !ML
    h[] = HR - zb[];
#else
    foreach_layer()
      h[] = (HR)/nl;
#endif
 	//eta[] = HR;  
  }

//Velocity initialization
  foreach (){
    uold[] = 0;  
    //foreach_layer(){
    //u.x[] = y/DT;
  //}
  }
}

#if !ML
event acc (i++, t<=T0){
  foreach () 
    for (vector u in ul) {
      u.x[] = u.x[] + G*sin(slope)*dt;
  }
}
#else
event acceleration (i++, t<=T0)
{
  foreach_face()
    foreach_layer()
      ha.x[] += G*sin(slope)*hf.x[]; 
}

#endif

#if 0
/** We check for convergence. */
event logfile (t += 0.1; t<=T0) {
//event logfile (i++;t<=T0) {
  double du = change (u.x,uold);
    if (i > 10 && du < 1e-11)
      return 1; /* stop */
}
#endif

/**
## Outputs

We save the velocity profile at regular intervals. */
#define uan( z ) ( (G*sin(slope)/(2*nu) )*(2-z)*z)

event profiles ( t = end )
{
  if ( N == 16 && nl == 128 ){
    foreach() {
      double z = zb[];
#if !ML
      for (int l = 0; l <= nl - 1; l++) {
        vector u = ul[l]; 
        fprintf (stdout,"%g %g %g %g %g %g\n", x, z, u.x[], uan( (z + h[]*layer[l]/2) ), dt, t);
        z += h[]*layer[l];
      }
#else
      foreach_layer(){
        fprintf (stdout,"%g %g %g %g %g %g\n", x, z, u.x[], uan( (z + h[]/2) ), dt, t);
        z += h[];

      }
#endif
      fprintf (stdout,"\n \n");   
    }
  }
}

/**
We also compute the error between the numerical solution and the analytical solution   */

event error (t = end)
{
  int i = 0;
  foreach() {
    if (i++ == N/2) {
      double z = zb[];
      double norm = 0, norm2 = 0, normax = 0; 
    #if ML
      foreach_layer() {
        double e = fabs(u.x[] - uan ( (z + h[]/2.) ) );
        norm += e*h[];
        norm2 += sq(e)*h[];
        normax = max( normax, e );
        z += h[];
      }
    #else
        for (int l = 0; l <= nl - 1; l++) {
          vector u = ul[l]; 
          double e = fabs(u.x[] - uan ( (z + h[]*layer[l]/2.) ));
          norm += e*h[]*layer[l];
          norm2 += sq(e)*h[]*layer[l];
          normax = max( normax,e );
          z += h[]*layer[l];
        }
    #endif
      norm = norm/z;
      norm2 = sqrt(norm2/z);
      fprintf (stderr, "%d %d %g %g %g %g %g\n", N, nl, norm, norm2, normax, dt, t);
    }
  }
}

/**
# Results
~~~gnuplot Velocity profiles for flow
set key bottom right
 set ylabel "z"
 set xlabel "u"
p 'out' u 3:2 w lp t'U computed', "" u 4:2 w l t 'U exact'
~~~

# Convergence
We track the value of the relative error on the velocity for various
number of layers $\text{nl}$. The error decrease with 2nd ordrer precision.

~~~gnuplot Variation of the relative error with number of layers, for different grid $N$.
set key bottom left
set xlabel 'nl'
set ylabel 'Max |e|'
set logscale
#set format x "%.0e"
set format y "%.2e"

set xrange [2:512]
set cbrange [1:2]
set xtics 2,2,512

set yrange [1e-6:1e-1]
#set cbrange [1:2]
#set xtics 1e-5,10,1e-1

set grid ytics

#set label 1 "N^{-1}" at 32,0.9*16**-1 font ',12' textcolor rgb 'red' offset 0,1
set label 1 "N^{-2}" at 8,0.01*16**-2 font ',12' textcolor rgb 'purple' offset 0,1

plot for [i=0:3] 'log' index i u 2:5 t "N=".columnhead(1) with lp lw 2 ps 1.5 lt i+2,\
         [4:8<<2] 0.02*x**-2 t '' w l lw 2 lc rgb 'purple'

         #[4:8<<2] 0.001*x**-1 t '' w l lw 2 lc rgb 'red', [4:8<<2] 0.001*x**-2 t '' w l lw 2 lc rgb 'purple'
~~~

*/ 





























