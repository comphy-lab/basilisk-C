/**
# Boussinesq buoyancy for isopycnal layers

This adds buoyancy to the [multilayer solver](hydro.h), assuming that
each layer has a constant density variation.

This is a special case of the more general [buoyancy module](dr.h)
which should be consulted for more details.

The density variations in each layer are defined by the user using the
`drho` array whose dimension must match the numbers of layers. */

extern double * drho;

event acceleration (i++)
{
  scalar q = new scalar[nl];
  
  foreach() {
    double ph = 0.;
    for (point.l = nl - 1; point.l >= 0; point.l--) {
      double dp = G*drho[point.l]*h[];
      ph += dp;
      q[] = ph;
    }
  }
  
  foreach_face()
    hpg (pg, q, 0,
	 ha.x[] += pg);

  delete ({q});
}
