/**
# Skalak law for the Eulerian elasticity framework.

In this file we implement the Skalak law to compute the elastic surface stresses
of biological capsules. The computation of the stresses rely on the Eulerian
quantities defined in [elasticity.h](elasticity.h). The force are transferred
to the fluid in the whole membrane region computed in [capsule.h](capsule.h).

This implementation is greatly inspired from
[log-conform.h](http://basilisk.fr/src/log-conform.h)
*/


#ifndef AREA_DILATATION_MODULUS
  #define AREA_DILATATION_MODULUS 1.
#endif

event acceleration (i++)
{
  /**
  We use Skalak's law to compute the shear stress in the membrane region
  */
  double p, q;
  foreach()
  {
    if (GAMMA)
    {
       // p = J[]*(G.x.x[] + G.y.y[]) - 1.; // For 3-dimensional simulations
       p = J[]*(G.x.x[] + G.y.y[]); // For 2-dimensional simulations
       q = J[]*(AREA_DILATATION_MODULUS*(sq(J[]) - 1.) - 1.);
       foreach_dimension()
       {
         T_s.x.x[] = ngcaps[]*(p*G.x.x[] + q*(1. - (sq(extended_n.x[]))));
         T_s.x.y[] = ngcaps[]*(p*G.x.y[] - q*(extended_n.x[]*extended_n.y[]));
       }
    else {
      foreach_dimension() {
        T_s.x.x[] = 0.;
        T_s.x.y[] = 0.;
      }
    }
  }
  boundary((scalar *){T_s});
}
