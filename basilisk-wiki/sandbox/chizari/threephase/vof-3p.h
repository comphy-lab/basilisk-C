/**
The *vof-3p.h* is basically same as [vof.h](http://basilisk.fr/src/vof.h) here,
the only difference is the fact that we are using the
[fraction-3p.h](/sandbox/chizari/threephase/fractions-3p.h)
instead of [fraction.h](/src/fraction.h) since the recustruction
of the two-phase cells in the three-phase flow is modified.
Therefore, the include of the fraction.h shall be modified to
*fractions-3p.h* and the rest remains the same. */

attribute {
  scalar * tracers;
  bool inverse;
}

#include "fractions.h"

extern scalar * interfaces;
extern face vector uf;
extern double dt;

event defaults (i = 0)
{
#if TREE
  for (scalar c in interfaces)
    c.refine = c.prolongation = fraction_refine;
#endif
}

event stability (i++) {
  if (CFL > 0.5)
    CFL = 0.5;
}

foreach_dimension()
static void sweep_x (scalar c, scalar cc)
{
  vector n[];
  scalar alpha[], flux[];
  double cfl = 0.;
  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }
    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl) {
	double cl = c[-1], cc = c[], cr = c[1];
	if (t.inverse)
	  cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
	gf[] = 0.;
	static const double cmin = 0.5;
	if (cc >= cmin) {
	  if (cr >= cmin) {
	    if (cl >= cmin) {
	      if (t.gradient)
		gf[] = t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
	      else
		gf[] = (t[1]/cr - t[-1]/cl)/(2.*Delta);
	    }
	    else
	       gf[] = (t[1]/cr - t[]/cc)/Delta;
	  }
	  else if (cl >= cmin)
	    gf[] = (t[]/cc - t[-1]/cl)/Delta;
	}
      }
    }
    boundary (gfl);
  }

  reconstruction (c, n, alpha);
  
  foreach_face(x, reduction (max:cfl)) {

    double un = uf.x[]*dt/(Delta*fm.x[]), s = sign(un);
    int i = -(s + 1.)/2.;

    if (un*fm.x[]*s/cm[] > cfl)
      cfl = un*fm.x[]*s/cm[];

    double cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
      rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			  (coord){-0.5, -0.5, -0.5},
			  (coord){s*un - 0.5, 0.5, 0.5});
    
    flux[] = cf*uf.x[];
    
    scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = c[i];
      if (t.inverse)
	cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
	double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
	tflux[] = ff*cf1*uf.x[];
      }
      else
	tflux[] = 0.;
    }
  }
  delete (gfl); free (gfl);

#if TREE
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_append (fluxl, flux);
  for (int l = depth() - 1; l >= 0; l--)
    foreach_halo (prolongation, l) {
#if dimension == 1
      if (is_refined (neighbor(-1)))
	for (scalar fl in fluxl)
	  fl[] = fine(fl);
      if (is_refined (neighbor(1)))
	for (scalar fl in fluxl)
	  fl[1] = fine(fl,2);
#elif dimension == 2
      if (is_refined (neighbor(-1)))
	for (scalar fl in fluxl)
	  fl[] = (fine(fl,0,0) + fine(fl,0,1))/2.;
      if (is_refined (neighbor(1)))
	for (scalar fl in fluxl)
	  fl[1] = (fine(fl,2,0) + fine(fl,2,1))/2.;
#else // dimension == 3
      if (is_refined (neighbor(-1)))
	for (scalar fl in fluxl)
	  fl[] = (fine(fl,0,0,0) + fine(fl,0,1,0) +
		  fine(fl,0,0,1) + fine(fl,0,1,1))/4.;
      if (is_refined (neighbor(1)))
	for (scalar fl in fluxl)
	  fl[1] = (fine(fl,2,0,0) + fine(fl,2,1,0) +
		   fine(fl,2,0,1) + fine(fl,2,1,1))/4.;
#endif
    }
  free (fluxl);
#endif

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     cfl - 0.5), fflush (ferr);

  foreach() {
    c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
    scalar t, tflux;
    for (t, tflux in tracers, tfluxl)
      t[] += dt*(tflux[] - tflux[1])/(cm[]*Delta);
  }
  boundary ({c});
  boundary (tracers);

  delete (tfluxl); free (tfluxl);
}

void vof_advection (scalar * interfaces, int i)
{
  for (scalar c in interfaces) {
    scalar cc[];
    foreach()
      cc[] = (c[] > 0.5);

    void (* sweep[dimension]) (scalar, scalar);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    boundary ({c});
    for (d = 0; d < dimension; d++)
      sweep[(i + d) % dimension] (c, cc);
  }
}

event vof (i++)
  vof_advection (interfaces, i);
