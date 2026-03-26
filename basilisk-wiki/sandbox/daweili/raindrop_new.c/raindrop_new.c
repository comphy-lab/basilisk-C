/**
# Films do not break (on GPUs either)

![Inside distance function](raindrop_new/movie.mp4)(width="80%") 

~~~gnuplot
set ylabel 'Minimal distance'
set xlabel 'Time'
set logscale y
set grid
plot 'log' u 1:($2*2.) w l t ''
~~~
*/

#include "grid/gpu/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"

const double tend = 12;
const double We = 5.;
const double Re = 1000.;
const int maxlevel = 11;

double dmin0;
scalar Y[];

FILE * fp;

int main()
{
  periodic(right);

  dmin0 = 1.;
  L0 = 5.;
  X0 = - L0/2.;
  //Y0 = - L0/2. - 1.;
  Y0 = -L0 + 0.1 ;
  const scalar sigma[] = 1./We;
  d.sigmaf = sigma;

  fp = fopen("dmin.dat","w");
  rho1 = 1, rho2 = 0.1;
  mu1 = mu2 = 1./Re;

  N = 1 << maxlevel;
  run();
  
  fclose(fp);
}

event acceleration (i++)
{
  foreach_face (y)
    a.y[] -= 1.;
}

event init (t = 0) {
  foreach()
    d[] = - sq(y + 0.01*sin(2*pi*x/L0)) + sq(0.5);
}

event movies (i += 30; t <= tend)
{
  scalar df[];
  foreach()
    df[] = f[];
    //df[] = f[] > 0.5 ? d[] : nodata;
  output_ppm (df, spread = -1, file = "movie.mp4", n = 1024);
  //  output_ppm (df, spread = -1, file = "zoom.mp4", box = {{-1.25,-1.75},{1.25,2}});
}

scalar dmin[];

event logfile (i++)
{
  
  scalar dmin[];
  foreach ()
    dmin[] = y*(1.-f[]);
  double dmintmp = 0.1-statsf(dmin).max;
  

  fprintf (stderr, "%g %g\n", t, dmintmp );
  fprintf (fp, "%g %g\n", t, dmintmp );

}