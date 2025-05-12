/**
# Bouncing Saint-Venant bump with the Mc Cormack scheme in 2D */

#include "../vdw.h"
#include "view.h"

/**
We start with initial conditions
etc... as when using the standard 
[bump2D.c](/src/test/bump2D.c) .


![Density during simulation](bump-vdw/rho.mp4)
 */

#define LEVEL 8

#define rho_c 1
#define R_g 1 
#define theta 0.95
#define p_c 1

#define MU 1e-3

double P0(double x)
{
  double rhop;
  rhop=x/rho_c;
  return p_c*rhop*theta*(8/(3-rhop) - 3*rhop/theta);
}
  
int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  DT = 1e-3;
  run();
}

event init (i = 0)
{
  lambda=1.e-2;
  foreach()
    {
      mu[] = 1.e-5;
      rho[] = 0.1 + exp(-200.*(sq(x)+sq(y)));
      q.x[] = q.y[] = 0.;
    }
}

event adapt(i++){
  adapt_wavelet ({rho,q.x,q.y},
    (double[]){1.e-3,1.e-2,1.e-2},8, 5);
}


event logfile (i++;t<=0.4) {
  stats s = statsf (rho);
  fprintf (stderr, "%g %4.2g     %4.2g %4.2g %.8f\n", t,  dt, s.min,s.max,s.sum);
  if(s.min < 1.e-3){
    dump();
    exit(1);
  }
}


event outputfile (t+=5.e-3) {
  squares("rho");
  save("rho.mp4");
}

