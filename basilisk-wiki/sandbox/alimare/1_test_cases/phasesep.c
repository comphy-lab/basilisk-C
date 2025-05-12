/**
# Phase separation in 1D

This tests the VdW EOS and spinodal decomposition */

#define MU 1e-1
#include "../vdw.h"

/**
We start with initial conditions
etc... 

![Animation of the density](phasesep/movie_normal.mp4)
*/

#define LEVEL 5

#define rho_c 1
#define R_g 1 
#define theta 0.95
#define p_c 1

double P0(double x)
{
  double rhop;
  rhop=x/rho_c;
  return p_c*rhop*theta*(8/(3-rhop) - 3*rhop/theta);
}
  
int main()
{
  mgu.nrelax = 20;
  TOLERANCE = 1.e-6;
  periodic (right); 
  // periodic(top);
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  DT = 1e-3;
  system ("rm -r gnuplot");
  system ("mkdir gnuplot");
  run();
}

event init (i = 0)
{
  lambda=0.01;
  foreach()
    {
      foreach_dimension()
        mu.x[] = MU;
      rho[] = rho_c + 0.4*sin(2*pi*x);
      q.x[] = q.y[] = 0.;
    }
   boundary(all);
}

event logfile (i += 200) {
  stats s = statsf (rho);
  stats s2 = statsf(u.y);
  stats s3 = statsf(u.x);
  fprintf (stderr, "%g %d %g %g %g %.8f %g %g %g %g\n", t, i, dt, s.min, s.max, s.sum,
    s2.min, s2.max, s3.min, s3.max);
}

int fileNumber = 0;

event gnuplot(t+=0.1,last){
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if(i==0){
    fprintf (fp,
     "unset key\n"
     "set xlabel 'x'\n"
     "set ylabel 'Density'\n"
     );
  }
  fprintf (fp,
     "set term png notransparent truecolor\n"
     "set yrange [0.4:1.8]\n"
     "set output 'gnuplot/plot-%06d.png'\n",fileNumber++);
     // "stats 'out' nooutput\n", 
     
  fprintf (fp,
     "set title 't = %03d'\n"
     "plot '-' u 1:2 t sprintf('t=%d')\n",
     fileNumber,fileNumber);
  foreach (serial) {
    fprintf (fp, "%g %g %g\n", x, rho[],q.x[]);
  }
  fprintf (fp, "e\n");
  fflush (fp);
}

event end(t= 40){
  system ("for f in gnuplot/plot-*.png; do"
    " convert $f ppm:- && rm -f $f; done | "
    "ppm2mp4 movie_normal.mp4");
  fprintf (stderr, "\n\nDone\n");
}
