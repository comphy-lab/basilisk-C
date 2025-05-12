/**
# Spinodal decomposition with a phase field method

 */
#define growthTimestep 0.1
#define MU 1e-5
#include "../vdw.h"
#include "view.h"
#include "../diverging_cool_warm.h"
/**



![Density during simulation, isoline is the critical density](spinodal/rho.mp4)

*/

#define LEVEL 8
#define LAMBDA 1.4e-1
#define rho_c 1.06344
#define R_g 1 
#define theta 0.9
#define p_c 1


// Calculate adimensionalized numbers.

double P0(double x)
{
  double rhop;
  rhop=x/rho_c;
  return p_c*rhop*theta*(8/(3-rhop) - 3*rhop/theta);
  // return 0.5*sq(x);
}

double mynoise(){
  return (1-2.*exp(noise())/exp(1));
}
  
int main()
{
  origin (-0.5, -0.5);
  periodic(right);
  periodic(top);
  init_grid (1 << LEVEL);
  DT = 1e-4;
  run();
}

event init (i = 0)
{
  mgu.nrelax = 10;
  lambda=LAMBDA;
  TOLERANCE = 1.e-6;
  foreach()
    {
      foreach_dimension()
        mu.x[] = MU;
      rho[] = rho_c *( 1.+0.2*noise());
      q.x[] = q.y[] = 0.;
    }
  boundary({mu,rho,q});
}


event end (t=1)
  dump();

event outputfile (t=1.e-2;t*=1.005) {
  face vector fs[];
  isoline("rho", rho_c);
  squares("rho",  map = mycoolwarm);
  save("rho.mp4");

  stats s = statsf (rho);
  stats s2 = statsf(u.x);
  fprintf (stderr, "%6.4g %4.2e    %6.3g %6.2g %4.2g %4.2g %ld\n", t,  dt,
    s.min,s.max,s2.min,s2.max, grid->tn);
}


event adapt(i++){
  adapt_wavelet ({rho,q.x,q.y},
    (double[]){5.e-3,3.e-3,3.e-3},LEVEL, 5);
}
