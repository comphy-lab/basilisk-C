/** Riemann solver for St. Venant.  Uses wave propagation algorithm.
*/

#include "waveprop.h"

scalar h[];
vector hu[];
double grav = 1;

int mwaves = 3;
int limiters[3] = {4,4,4}, *mthlim = limiters;

scalar * scalars = {h,hu};

void wpa_rpsolver (int dir, int meqn, int mwaves, double *ql, double *qr, 
		   double *waves, double *speeds, double *amdq, double *apdq, 
		   double *flux)
{
  if (meqn != 3) {
    fprintf (stderr, "Error : meqn != 3\n");
    exit(0);
  }

  int mu = (dir == 0) ? 1 : 2;
  int mv = (dir == 0) ? 2 : 1;
  
  double  hl = ql[0],  hr  = qr[0];
  double hul = ql[mu], hur = qr[mu];
  double hvl = ql[mv], hvr = qr[mv];
  
  double delta[3];  
  delta[0] = qr[0]   - ql[0];
  delta[1] = qr[mu]  - ql[mu];
  delta[2] = qr[mv]  - ql[mv];
  
  if (hl <= 0 || hr <= 0) {
    fprintf (stderr, "hr, hl <= 0;   %g %g\n", hl, hr);
    exit (0);
  }

  double ubar, vbar, cbar, divsqrt;
  divsqrt = sqrt(hl) + sqrt(hr);
  
  ubar = (hul/sqrt(hl) + hur/sqrt(hr))/divsqrt;
  vbar = (hvl/sqrt(hl) + hvr/sqrt(hr))/divsqrt;
  cbar = sqrt(grav*(hr + hl)/2.0);
  
  double a1, a2, a3, cbar2;
  cbar2 = 2*cbar;
  a1 = (-delta[1] + (ubar + cbar) * delta[0])/cbar2;
  a2 = -vbar*delta[0] + delta[2];
  a3 = (delta[1] - (ubar - cbar) * delta[0])/cbar2;
  
  int mw = 0;
  waves[0  + mw*meqn] = a1;
  waves[mu + mw*meqn] = a1*(ubar - cbar);
  waves[mv + mw*meqn] = a1*vbar;
  speeds[mw] = ubar - cbar;
  
  mw = 1;
  waves[0  + mw*meqn] = 0;
  waves[mu + mw*meqn] = 0;
  waves[mv + mw*meqn] = a2;
  speeds[mw] = ubar;
  

  mw = 2;
  waves[0  + mw*meqn] = a3;
  waves[mu + mw*meqn] = a3*(ubar + cbar);
  waves[mv + mw*meqn] = a3*vbar;
  speeds[mw] = ubar + cbar;
  
  for(int m = 0; m < meqn; m++) {
    amdq[m] = 0;
    apdq[m] = 0;
  }
  
  for(int m = 0; m < meqn; m++)
    for(int mw = 0; mw < mwaves; mw++)
      if (speeds[mw] < 0) 
	amdq[m] += speeds[mw]*waves[m+mw*meqn];
      else
	apdq[m] += speeds[mw]*waves[m+mw*meqn];  
  
  flux[0] = hur;
  flux[mu] = hur*hur/hr + 0.5*grav*hr*hr;  
  flux[mv] = hur*hvr/hr;                
}
