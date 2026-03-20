/**
# Films do not break (on GPUs either)

![video](foam_gpu/movie.mp4)(width="50%") 


*/

#include "grid/gpu/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"

const double tend = 40;
const double We = 2.;
const double Re = 50;
const int maxlevel = 7;

double dmin0;
scalar Y[];

FILE * fp;

int main()
{
  periodic(right);

  
  L0 = 5.;
  X0 = - L0/2.;
  Y0 = -L0/2;
  const scalar sigma[] = 1./We;
  d.sigmaf = sigma;

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
  foreach() {
    d[] =  sin(2*pi*y/L0*5) + 0.01*sin(2*pi*x/L0);
  }
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

