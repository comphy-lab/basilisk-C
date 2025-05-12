/**
# Low-Storage RK-3 and RK-4 schemes

The default is 3-stage 3rd order method, if defined otherwise, a
5-stage 4th order method is used.

As suggested in,

van Heerwaarden, et al. *MicroHH 1.0: ...* Geoscientific Model
Development (2017).
*/
#include "run.h"
#ifndef RKORDER
#define RKORDER (3)
#endif
#if (RKORDER == 3)
// Williamson, J. H.: Low-Storage Runge-Kutta schemes, J.
// Comput.Phys., 35, 48–56, 1980.
#define STAGES (3)
double An[STAGES] = {0., -5./9., -153./128.};
double Bn[STAGES] = {1./3., 15./16., 8./15.};
#else
// Carpenter, M.  H.  and Kennedy, C.  A.: Fourth-order
// 2N-storageRunge-Kutta schemes, Tech. Rep. TM-109112, NASA
// LangleyResearch Center, 1994
#define STAGES (5)
double An[STAGES] = {0.,
		     -567301805773. /1357537059087.,
		     -2404267990393./2016746695238.,
		     -3550918686646./2091501179385.,
		     -1275806237668./842570457699.};
double Bn[STAGES] = {1432997174477./9575080441755. ,
		     5161836677717./13612068292357.,
		     1720146321549./2090206949498. ,
		     3134564353537./4481467310338. ,
		     2277821191437./14882151754819.};
#endif

scalar * dfl = NULL; 
void A_Time_Step (scalar * fl, double dt,
		  void (* Lu) (scalar * ul, scalar * dul)) {
  if (dfl == NULL)
    dfl = list_clone (fl);
  scalar * dfltmp = list_clone (dfl);
  for (int Stp = 0; Stp < STAGES; Stp++) {
    Lu (fl, dfltmp);
    scalar f, df, dftmp;
    foreach() {
      for (f, df, dftmp in fl, dfl, dfltmp) {
	df[]  = An[Stp]*df[] + dftmp[];
	f[]  += Bn[Stp]*df[]*dt;
      }
    }
  }
  delete (dfltmp); free (dfltmp); dfltmp = NULL;
}

event rm_dfl (t = end) {
  delete (dfl); free (dfl); dfl = NULL;
}
/**
## Test

* [4th order scheme](tlsrk.c)

## Usage

* [A co-located Navier--Stokes solver](nsrk.h)
 */
