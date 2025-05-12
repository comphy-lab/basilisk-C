/**
# A So-called $p$-adaptive method

Rather than using grid refinement, the accuracy of the spatial approximations can locally be improved by using a higher order $p$olynominal method. As such, we *adaptively* switch between a second and fourth order accuracte formulation based on a wavelet-based error estimate. 
*/

#define BGHOSTS 2
#include "grid/multigrid.h"
#include "runge-kutta.h"
#include "run.h"
/**
## Advection Diffusion. 

The advection-diffusion equation is solved on a 2D equidistant grid, according to a advection velocity face field (`uf`) and diffusivity (`mudiff`)
*/
(const) face vector mudiff;
(const) face vector uf;
/**
The $p$ refimenement criterion is $\zeta$ and the error estimated is stored in field `w`. 
*/
scalar w[];
double zeta = 0.01;
/**
Functions are defined for 2nd and 4th order accurate representations of the netto flux for advection and diffusion. 
*/
double adv2 (Point point, scalar s) {
  double adv = 0;
  foreach_dimension() {
    adv += (uf.x[]*(s[] + s[-1]) -
	    uf.x[1]*(s[1] + s[0]))/2.;
  }
  return adv;
}

double adv4 (Point point, scalar s) {
  double adv = 0;
  foreach_dimension() {
    adv += (uf.x[]*(7.*(s[] + s[-1]) - (s[-2] + s[1])) -
	    uf.x[1]*(7.*(s[1] + s[]) - (s[-1] + s[2])))/12.;
  }
  return adv;
}

double diff2 (Point point, scalar s) {
  double diff = 0;
  foreach_dimension() {
    diff +=  (-mudiff.x[]*(s[] - s[-1]) +
	      mudiff.x[1]*(s[1] - s[]))/Delta; 
  }
  return diff;
}

double diff4 (Point point, scalar s) {
  double diff = 0;
  foreach_dimension() {
    diff +=   (-mudiff.x[]*(15.*(s[] - s[-1]) - (s[1] - s[-2])) +
	       mudiff.x[1]*(15.*(s[1] - s[]) - (s[2] - s[-1])))/(12.*Delta); 
  }
  return diff;
}
/**
The `update` function is formulated such that it can be used with the Runge-Kutta solver.
*/
static void adv_diff_update (scalar * s, double t, scalar *kl) {
  boundary (s);
  for (int j = 0; j < list_len(s); j++) {
    scalar m = s[j];
    scalar k = kl[j];
    wavelet (m, w);
    foreach() {
      if (fabs(w[]) < zeta)
	k[] = (adv2 (point, m) + diff2 (point, m))/Delta;
      else
	k[] = (adv4 (point, m) + diff4 (point, m))/Delta;
    }
  }
}

struct Adv_Diff{
  scalar s;       
  double dt;      // timestep
  face vector D;  // diffusivity, default 0
  face vector uf; // Advection velocity, default 0
};

double PECLET = 0.4; // Default Cell Peclet Number
trace
int adv_diff (struct Adv_Diff p) {
  scalar f = p.s;
  double dtmin = p.dt;
  if (p.D.x.i) mudiff = p.D;
  else mudiff = zerof;
  if (p.uf.x.i) uf = p.uf;
  else uf = zerof;
  /**
     The forward-in-time integrator is porne to numerical
     instabiliy. Therefore we introduce and ensure a limited
     cell-Peclet number ($Pe$,`PECLET`) and a maximum Courant number
     ($Co$, `CFL`);
     
     $$\mathrm{d}t_{\mathrm{Diff}}<Pe\frac{\Delta^2}{\kappa},$$
     
     $$\mathrm{d}t_{\mathrm{Adv}}<Co\frac{\Delta}{u},$$
     
     by cutting the requested time step `dt` into `it` equal
     parts. The default values of the `global double`s `PECLET`$=0.4$ and
     `CFL` $=0.6$ (see `src/common.h`).
  */
  if (!is_constant(mudiff.x) && !is_constant(uf.x)){
    foreach_face() {
      double DTmin = min(PECLET*sq(Delta)/mudiff.x[], CFL*Delta/fabs(uf.x[]));
      if (DTmin < dtmin)
	dtmin = DTmin;
    }
  } else { 
    foreach_face (reduction(min:dtmin)) {
      double kappa = max(mudiff.x[], 1E-30);
      double ufx = max(fabs(uf.x[]), 1E-30);
      double delt = (double)L0/(1 << grid->maxdepth);
      dtmin = min(PECLET*sq(delt)/kappa, CFL*delt/ufx);
      break;
    }
  }
  int it = 1;
  if (dtmin < p.dt)
    it = (int)((dt/dtmin) + .99);
  double dtD = dt/(double)(it);
  for (int m = 0; m < it; m++) 
    runge_kutta ({f}, t, dtD, adv_diff_update, 4);
  return it;
}
/**
## Setup 

We consider a scalar field $s$ (`s`). 
*/

scalar s[];

int main() {
  s.prolongation = refine_injection;
  L0 = 10;
  X0 = Y0 = -L0/2;
  /**
  On a $128 \times 128$ grid, we perform the computations with various runs using a decreasing refinement criterion. 
  */
  N = 128;
  for (zeta = 0.01; zeta > 0.001; zeta /= 2)
    run();
}
/**
The solution is initialized...
*/
event init (t = 0) {
  foreach()
    s[] = exp(-sq(x + 1));
  DT = 0.01;
  boundary ({s});
  wavelet (s, w);
}
/**
Time integration is performed along side with update to the error-estimate field.
*/
event advance (i++) {
  dt = dtnext(DT);
  printf("%g %d %d\n",t, i, adv_diff (s, dt, unityf, unityf));
}
/**
## Output

There is only visual output
*/
event movie (i += 2; t < 1) {
  output_ppm (s, file = "s.mp4", n = 256);
  output_ppm (w, file = "w.mp4", n = 256);
  scalar c[];
  foreach() 
    c[] = fabs(w[]) > zeta ? 1. : 0.;
  output_ppm (c, file = "c.mp4", n = 256, min = 0, max = 1);
}
/**
## Results

The 1D solutions look like this:

![](polynominal_adaptivity/s.mp4)

The cells that are update with 4th order methods are marked red in the movie below:

![](polynominal_adaptivity/c.mp4)

*/