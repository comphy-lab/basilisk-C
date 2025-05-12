/**
In this case, we want to simulate the bouncing of an oil droplet in a stratified mixture of water and ethanol due to a Marangoni instability described in [this paper](https://doi.org/10.1017/jfm.2023.415).
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "../src/advdiff.h"

#define wtp(xnd) (0.01*xnd + 0.209513)
#define st(wte) (1. - 3.1434513*wte + 4.0225458*wte*wte - 1.8601032*wte*wte*wte)

int LEVEL = 8;

const double Oh = 0.88433, Ga = 0.023274;
const double tend = 120000.;
const double rhoe = 785./998.;
scalar sigmaf[];

int main()
{
  size (20);
  origin (-4.,0.);
  rho1 = 1.; // water
  rho2 = 966./998.; // drop
  mu1 = 1./96.; //water
  mu2 = 1.; //drop
  
  Diff_C1 = 2.38095*1.e-9/0.0966; //water
  Diff_C2 = 2.38095*1.e-12; //drop
  
  d.sigmaf = sigmaf;
    
  TOLERANCE = 1e-4 [*];

  /* double U_drop = - 2./((2. + 3.*mu1/mu2)*(2. + Diff_C2/Diff_C1))/mu2; */
  /* fprintf(ferr,"#%g \n",U_drop); */
  
    N = 1 << LEVEL;
    run();
}

/**
We initialize the signed distance *d* -ve d[] is drop and +ve d[] is water + ethanol and the surface tension field. */

event init (t = 0)
{
  if (!restore (file = "restart")) {
    foreach() {
      d[] = sqrt (sq(x) + sq(y)) - 1; // -ve d[] is drop and +ve d[] is water
      PhiC[] = wtp(x); // Misbehaves if PhiC is 0 in drop 
      sigmaf[] = st(PhiC[])/sq(Oh);
    }
  }
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume. */

double u_drop = 0.;

event logfile (i += 2)
{
  double xb = 0., vb = 0., sb = 0.;
  static double xb0 = 0., previous = 0.;
  if (t == 0.)
    previous = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vb += u.x[]*dv;
    xb += x*dv;
    sb += dv;
  }
  static double sb0 = 0.;
  if (i == 0) {
    sb0 = sb;
    fprintf (ferr, "\n#i dt t/t0 dsb xb vb ta u_drop dt perf.t perf.speed\n");
  }
  u_drop = t > previous ? (xb/sb - xb0)/(t - previous) : 0.;
  fprintf (ferr, "%d %g %g %g %g %g %g %g %g %g %g\n",
	   i,dt,t, (sb - sb0)/sb0, xb/sb, vb/sb,
	   (t + previous)/2., u_drop,
	   dt, perf.t, perf.speed);
  xb0 = xb/sb, previous = t;
}

/**
Surface-tension is updated on the basis of Phi field. The acceleration term is added with extra contribution from Boussinesq approximation due to stratification of water with ethanol.
*/

event acceleration (i++) {
  foreach(){
    sigmaf[] = st(PhiC[])/sq(Oh);
    f[] = clamp(f[], 0., 1.);
  }
  
  face vector av = a;
  foreach_face(x){
    av.x[] -=  Ga;
    av.x[] -=  Ga*(f[])*(rhoe - rho2)*(PhiC[])/(f[]*rho1+(1. - f[])*rho2);
  }
}

/**
Finally, we output the priliminary fields to a dump file.
 */
event movie (t += 10){
  char filename[60];
  sprintf(filename,"dump-%4.4f",t);
  dump(filename,{f,u,PhiC,sigmaf,d});

}

/**
The simulations until tend is reached.
*/
event end (t = tend)
{

}