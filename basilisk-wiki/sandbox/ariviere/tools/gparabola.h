/**
I slighly modified parabola.h of the main branch to add a fit for the gaussian curvature.
   **/
#include "utils.h"


// Define this to use a x^iy^j polynomial with i = 0...NP-1, j = 0...NP-1
// #define NP 3


static double parabola_fit_gcurvature (ParabolaFit * p,
				      double gkappamax)
{
  double kappa;
#if dimension == 2//just the mean curvature
  double dnm = 1. + sq(p->a[1]);
  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);
//  if (kmax)
//    *kmax = fabs (kappa);
#else /* 3D */
# ifdef NP
  double hxx = 2.*p->a[2*NP], hyy = 2.*p->a[2], hxy = p->a[NP + 1];
  double hx = p->a[NP], hy = p->a[1];
# else
  double hxx = 2.*p->a[0], hyy = 2.*p->a[1], hxy = p->a[2];
  double hx = p->a[3], hy = p->a[4];
# endif
  double dnm = 1. + sq(hx) + sq(hy);
  kappa = (hxx*hyy - hxy*hxy)/(dnm*dnm);
//  kappa = - (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)
//    /sqrt (dnm*dnm*dnm); //where does the - come from?
//  if (kmin) {
//      double km = - (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)
//    /sqrt (dnm*dnm*dnm);
//      double a = km*km/4. - kappa;
//      *kmin = fabs (km/2.);
//    if (a >= 0.)
//      *kmax += sqrt (a);
//  }
#endif /* 3D */
  if (fabs (kappa) > gkappamax) {
//    if (kmax)
//      *kmax = kappamax;
    return kappa > 0. ? gkappamax : - gkappamax;
  }
  return kappa;
}

#if AXI
static void parabola_fit_axi_gcurvature (const ParabolaFit * p,
					double r, double h,
					double * kappa)
{
  double nr = (p->m.x*p->a[1] + p->m.y)/sqrt (1. + sq(p->a[1]));
  /* limit the minimum radius to half the grid size */
  double kaxi = nr/max(r, h/2.);
  *kappa *= kaxi;//check that
//  if (kmax)
//    *kmax = max (*kmax, fabs (kaxi));
}
#endif /* 2D */
