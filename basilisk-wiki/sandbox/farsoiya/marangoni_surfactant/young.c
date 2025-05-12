/**
# Marangoni-induced translation due to a temperature gradient

This reproduces the test case in section 3.4 of [Al Saud et al.,
2018](#alsaud2018) which should be consulted for more details. 

![Final velocity field, interface, temperature and adaptive
 mesh](young/fields.png){ width="100%" }

*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
//vector chid[];
#include "marangoni.h"
#include "view.h"

int LEVEL = 8;

/**
See section 3.4 of [Al Saud et al., 2018](#alsaoud2018). Note that we
use a capillary number *Ca* 10 times larger than in Al Saud et al. to
make the computation faster, but the results are good also for
$\text{Ca} = 0.066$. */

const double R = 0.5 [1], NablaT = 0.133, Mu = 0.1, Rho = 0.2 [0];
const double Re = 0.066, Ca = 0.66;
const double Gamma_T = 0.1 ; //Re*sq(Mu)/(Rho*sq(R)*NablaT);
const double Gamma_0 = 0.1; //(Gamma_T*R*NablaT)/Ca;
const double t0 = Mu/(Gamma_T*NablaT);
const double Cdrop = 1., Cbulk = 1.;
double U_drop;

/**
We need a field to store the variable surface tension coefficient. */

// scalar sigmaf[];
int counter2 = 0;
int main()
{

  size (15.*R);
  origin (- L0/2.);
  rho1 = rho2 = Rho;
  mu1 = mu2 = Mu;
  // d.sigmaf = sigmaf;
  
  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 0.1e-4 [*];
  
  U_drop =  2./((2. + 3.*mu2/mu1)*(2. + Cdrop/Cbulk))*Gamma_T*R*NablaT/mu1;

  for (LEVEL = 7; LEVEL <= 9; LEVEL++) {
    N = 1 << LEVEL;
    counter2 = 0;
    run();
  }
}

/**
We initialize the signed distance *d* and the surface tension gradient. */
double temp_0 = 0;
scalar temp[];
event init (t = 0)
{
  fraction(f,sqrt (sq(x) + sq(y)) - R);

  foreach() {
    // d[] = sqrt (sq(x) + sq(y)) - R;
    temp[] = temp_0 + NablaT*x;

    sigmaf[] = Gamma_0 + Gamma_T*(temp_0 - temp[]);
  }
  boundary({sigmaf});
}

/**
A utility function to check the convergence of the multigrid
solvers. */

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stdout, "%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	     mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	     mg.nrelax);
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

double u_drop = 0.;

event logfile (i += 10)
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
    // fprintf (filename, "\nt dsb xb vb/U_drop ta u_drop/U_drop dt perf.t perf.speed\n");
  }
  char filename[200];
  sprintf(filename, "young-%d.dat",(1 << grid->maxdepth));
	static FILE * fileptr2 = fopen(filename,"w");

  u_drop = t > previous ? (xb/sb - xb0)/(t - previous) : 0.;
  fprintf (fileptr2, "%g %g\n", (t + previous)/2./t0, u_drop/U_drop); fflush(fileptr2);
  //  fprintf (stdout, "%g %g\n", (t + previous)/2./t0, u_drop/U_drop); fflush(stdout);

  xb0 = xb/sb, previous = t;
}

event graphics (t=0; t += 0.03*t0; t <= 3.*t0)
// event graphics (i=0;i+=1;i < 10)
{
   double U = - Gamma_T*R*NablaT/Mu;
   fprintf (stderr, "%d %g %g\n", N/16, u_drop/U, U_drop/U);
   if (LEVEL == 7) {
     view (fov = 30, near = 0.01, far = 1000,
 	  tx = 0.009, ty = -0.076, tz = -0.291,
 	  width = 1239, height = 575);
     draw_vof (c = "f", filled = - 1, fc = {1,1,1});
     draw_vof (c = "f", lw = 2);
     squares (color = "temp", spread = 0.8, linear = true);
     vectors (u = "u", scale = 1);
     cells ();
     save ("fields.png");
   }
     counter2++;
    printf("\n %d %g %g %g",i,dt,t,t0); fflush(stdout);

}


event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-2, 1e-5, 1e-5}, LEVEL);
  
}

/**
~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})	
import matplotlib.ticker as mticker


plt.figure()


ts128, us128 = np.loadtxt('young-128.dat',delimiter=' ',unpack=True)
ts256, us256 = np.loadtxt('young-256.dat',delimiter=' ',unpack=True)
ts512, us512 = np.loadtxt('young-512.dat',delimiter=' ',unpack=True)

sts512, sus512 = np.loadtxt('../stephane-512.dat',delimiter=' ',unpack=True)

anexp = np.ones(len(ts128));
plt.plot(sts512[1:-1:20],sus512[1:-1:20],'yo',label='Al. Saoud et. al. (2018) $1/32 R_0$');


plt.plot(ts128[1:-1:2],us128[1:-1:2],'bx',label='$1/8 R_0$');
plt.plot(ts256[1:-1:20],us256[1:-1:20],'gx',label='$1/16 R_0$');
plt.plot(ts512[1:-1:20],us512[1:-1:20],'kx',label='$1/32 R_0$');

plt.plot(ts128,anexp,'r',label='Young et. al. (1959)');



#plt.ylim(0,1.2)
#plt.xlim(0,3.5)

plt.legend();
plt.xlabel(r'$t\;\nabla T/\Gamma_T$')
plt.ylabel(r'$U(t)/U_{drop}$')
plt.tight_layout()

plt.savefig('youngetal.svg')


~~~

## References

~~~bib
@hal{alsaud2018, hal-01706565}
~~~
*/
